// Build: cmake --build build --target column_management_test; run: bin/column_management_test
#include "cvrp.hpp"
#include "full_r1c_coeffs_getter.hpp"
#include "r1c_coeffs_getter.hpp"
#include "t_node_in.hpp"
#include "t_node_out.hpp"

namespace RouteOpt::Application::CVRP {
    bool tellIfInt(const std::vector<double> &x, const std::vector<SequenceInfo> &cols);

    namespace ImposeBranchingDetail {
        void addBranchCutToUnsolved(BbNode *node, BbNode *&node2,
                                   const std::pair<int, int> &edge, int num_buckets, bool symmetry);
        void addBranchCutToUnsolvedInEnu(BbNode *node, BbNode *&node2,
                                        const std::pair<int, int> &edge, int *pool);
    }
}

using namespace RouteOpt;
using namespace RouteOpt::Application::CVRP;

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void addColumn(BbNode &node, SequenceInfo column, double cost) {
    std::vector<int> indices;
    for (int customer: column.col_seq) indices.push_back(customer - 1);
    indices.push_back(3);
    std::vector<double> values(indices.size(), 1.0);
    int begin = 0;
    SAFE_SOLVER(node.refSolver().addVars(1, static_cast<int>(indices.size()), &begin,
        indices.data(), values.data(), &cost, nullptr, nullptr, nullptr, nullptr))
    node.refCols().push_back(std::move(column));
    SAFE_SOLVER(node.refSolver().updateModel())
}

void buildNode(Solver &environment, BbNode &node) {
    BbNode::buildModel(1, 4, &environment, &node);
    addColumn(node, {{1}, SINGLETON_COLUMN_IDX}, 100);
    addColumn(node, {{2}, SINGLETON_COLUMN_IDX}, 100);
    addColumn(node, {{3}, SINGLETON_COLUMN_IDX}, 100);
    addColumn(node, {{1, 2, 3}, 2}, 10);
    addColumn(node, {{2, 1, 3}, 2}, 200);
    SAFE_SOLVER(node.refSolver().reoptimize())
}

bool hasSingleton(const BbNode &node, int customer) {
    return std::any_of(node.getCols().begin() + 1, node.getCols().end(),
                       [customer](const SequenceInfo &column) {
                           return column.col_seq == std::vector<int>{customer};
                       });
}

void testCleanup(Solver &environment) {
    BbNode node;
    buildNode(environment, node);
    // A stabilization artificial variable with no row coefficients and positive reduced cost.
    int begin = 0;
    double cost = 300;
    SAFE_SOLVER(node.refSolver().addVars(1, 0, &begin, nullptr, nullptr, &cost,
        nullptr, nullptr, nullptr, nullptr))
    node.refCols().push_back({{0, 1}, STAB_ARTI_COLUMN_IDX});
    SAFE_SOLVER(node.refSolver().updateModel())
    node.cleanIndexColForNode(0.25);
    require(node.getCols().size() == 6, "Cleanup must delete the expensive ordinary route");
    require(hasSingleton(node, 1) && hasSingleton(node, 2) && hasSingleton(node, 3),
            "Cleanup deleted a singleton");
    require(node.getCols().back().forward_concatenate_pos == STAB_ARTI_COLUMN_IDX,
            "Cleanup deleted a stabilization artificial column");
    require(node.getCols().front().col_seq == std::vector<int>({1, 2, 3}),
            "Cleanup deleted artificial column zero");
}

void testFixing(Solver &environment) {
    BbNode node;
    buildNode(environment, node);
    SAFE_SOLVER(node.refSolver().getObjVal(&node.refValue()))
    // An optimal dual certificate giving every singleton reduced cost above the gap of 9.
    std::vector<double> duals{3, 3, 4, 0};
    node.cleanIndexColForNode(20, duals);
    require(!hasSingleton(node, 1) && !hasSingleton(node, 2) && !hasSingleton(node, 3),
            "Singleton tag must not prevent reduced-cost fixing");
    require(node.getCols().size() == 2, "Fixing should leave artificial zero and the cheap route");
}

void testBranching(Solver &environment, bool enumerated) {
    for (const auto edge: {std::pair{0, 1}, std::pair{1, 2}}) {
        BbNode node;
        buildNode(environment, node);
        node.refIfInEnumState() = enumerated;
        BbNode *child = nullptr;
        if (enumerated) {
            ImposeBranchingDetail::addBranchCutToUnsolvedInEnu(&node, child, edge, nullptr);
        } else {
            node.refAllForwardBuckets() = new Bucket *[4];
            for (int i = 0; i < 4; ++i) {
                node.refAllForwardBuckets()[i] = new Bucket[1];
                for (int j = 0; j < 4; ++j)
                    if (i != j) node.refAllForwardBuckets()[i][0].bucket_arcs.push_back(j);
            }
            ImposeBranchingDetail::addBranchCutToUnsolved(&node, child, edge, 1, true);
            Rank1Cuts::Rank1CutsDataShared data(4);
            Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter getter(data);
            child->preprocessEnumeration(nullptr, getter, std::vector<routeOptLong>(4), {}, 0);
        }
        std::unique_ptr<BbNode> true_child(child);
        if (edge.first == 0) {
            require(!hasSingleton(node, 1) && hasSingleton(node, 2),
                    "Forbidding a depot edge must delete its singleton");
            require(hasSingleton(*child, 1), "Requiring a depot edge must allow its singleton");
        } else {
            require(!hasSingleton(*child, 1) && !hasSingleton(*child, 2),
                    "Requiring a customer edge must delete incompatible singletons");
            require(hasSingleton(node, 1) && hasSingleton(node, 2),
                    "Forbidding a customer edge must allow singletons");
        }
    }
}

void testCheckpoint(Solver &environment) {
    BbNode node;
    buildNode(environment, node);
    const auto path = std::filesystem::temp_directory_path() /
                      ("routeopt-singleton-" + std::to_string(
                          std::chrono::steady_clock::now().time_since_epoch().count()));
    OutNodeNameSpace::outFile.open(path, std::ios::binary);
    OutNodeNameSpace::writeColSequence(node.refSolver(), node.getCols());
    OutNodeNameSpace::outFile.close();
    std::vector<double> cost;
    std::vector<SequenceInfo> columns;
    InNodeNameSpace::inFile.open(path, std::ios::binary);
    InNodeNameSpace::readColSequence(cost, columns);
    InNodeNameSpace::inFile.close();
    std::filesystem::remove(path);
    require(columns == node.getCols() && cost[1] == 100,
            "Checkpoint changed route content or cost");
    require(columns[1].forward_concatenate_pos == 0,
            "Checkpoint must keep the existing route-position format");
}

void testRouteInterpretation() {
    const SequenceInfo tagged{{1}, SINGLETON_COLUMN_IDX};
    const SequenceInfo ordinary{{1}, 0};
    require(tagged.getForwardConcatenatePos() == 0 && tagged == ordinary,
            "The tag must preserve singleton route semantics");
    require(tellIfInt({1.0}, {tagged}), "Tagged singleton rejected as an integer route");
    require(!tellIfInt({0.5}, {tagged}), "Fractional singleton accepted as integral");
    require(!tellIfInt({1.0}, {{{0, 1}, STAB_ARTI_COLUMN_IDX}}),
            "Stabilization artificial column accepted as an integer route");

    Rank1Cuts::Rank1CutsDataShared data(4);
    Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter getter(data);
    Rank1Cuts::R1c cut;
    cut.info_r1c = {{1, 2, 3}, 0};
    cut.arc_mem = {{1, 2}, {1, 3}, {2, 3}};
    Eigen::SparseMatrix<double> coefficients;
    getter.getR1CCoeffs(std::vector<SequenceInfo>{tagged, ordinary, {{1, 2}, 1}},
                       {cut}, nullptr, true, coefficients);
    require(coefficients.coeff(0, 0) == 0 && coefficients.coeff(0, 1) == 0 && coefficients.coeff(0, 2) == 1,
            "Singleton tag changed limited-memory Rank-1 coefficients");
}

void testFullMemoryCoefficients(Solver &environment) {
    BbNode node;
    BbNode::buildModel(1, 4, &environment, &node);
    addColumn(node, {{1, 2, 3}, 2}, 10);
    Rank1Cuts::Rank1CutsDataShared data(4);
    Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter getter(data);
    Rank1Cuts::R1c cut, singleton_cut;
    cut.info_r1c = {{1, 2, 3}, 0};
    singleton_cut.info_r1c = {{1}, 0};
    for (bool sparse: {false, true}) {
        if (sparse) addColumn(node, {{3}, SINGLETON_COLUMN_IDX}, 100);
        Eigen::SparseMatrix<double> actual, reference;
        getter.getR1CCoeffs(node.getCols(), {cut, singleton_cut}, &node.refSolver(), false, actual);
        getter.getR1CCoeffs(node.getCols(), {cut, singleton_cut}, nullptr, false, reference);
        require(actual.coeff(0, 0) == 1 && actual.coeff(1, 0) == 0,
                "Artificial column must have exactly the cut rhs coefficient");
        require((Eigen::MatrixXd(actual) - Eigen::MatrixXd(reference)).isZero(),
                "Solver and route-based full-memory coefficients differ");
    }
}

void testEnumerationRowDeletion(Solver &environment, bool with_branch) {
    BbNode node;
    BbNode::buildModel(1, 4, &environment, &node);
    addColumn(node, {{1, 2, 3}, 2}, 10);
    if (with_branch) {
        int begin = 0;
        int indices[]{0, 1};
        double values[]{1, 1}, rhs = 1;
        char sense = SOLVER_EQUAL;
        SAFE_SOLVER(node.refSolver().addConstraints(1, 2, &begin, indices, values, &sense, &rhs, nullptr))
        node.refBrCs().push_back({{1, 2}, 4, true});
    }
    // More cut rows than columns, as after enumeration removes non-elementary routes.
    for (int customer = 1; customer <= 3; ++customer) {
        int begin = 0;
        double rhs = 0;
        char sense = SOLVER_LESS_EQUAL;
        SAFE_SOLVER(node.refSolver().addConstraints(1, 0, &begin, nullptr, nullptr, &sense, &rhs, nullptr))
        Rank1Cuts::R1c cut;
        cut.info_r1c = {{customer}, 0};
        cut.idx_r1c = 3 + customer + static_cast<int>(with_branch);
        node.refR1Cs().push_back(cut);
    }
    std::reverse(node.refR1Cs().begin(), node.refR1Cs().end());
    SAFE_SOLVER(node.refSolver().updateModel())
    Rank1Cuts::Rank1CutsDataShared data(4);
    Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter getter(data);
    node.preprocessEnumeration(nullptr, getter, std::vector<routeOptLong>(4), {}, 0);
    int rows, columns;
    SAFE_SOLVER(node.refSolver().getNumRow(&rows))
    SAFE_SOLVER(node.refSolver().getNumCol(&columns))
    require(rows == 4 && columns == 2 && node.refR1Cs().empty(),
            "Enumeration did not delete precisely the obsolete cut rows");
    if (with_branch) require(node.refBrCs().front().idx_brc == INVALID_BRC_INDEX,
                             "Deleted branch row still has an LP index");
}

int main() {
    try {
        static_assert(SINGLETON_COLUMN_IDX == STAB_ARTI_COLUMN_IDX - 1);
        testRouteInterpretation();
        Solver environment{};
        SAFE_SOLVER(environment.loadEnv(nullptr))
        SAFE_SOLVER(environment.setEnvOutputFlag(0, true))
        SAFE_SOLVER(environment.setEnvThreads(1, true))
        testCleanup(environment);
        testFixing(environment);
        testBranching(environment, false);
        testBranching(environment, true);
        testCheckpoint(environment);
        testFullMemoryCoefficients(environment);
        testEnumerationRowDeletion(environment, false);
        testEnumerationRowDeletion(environment, true);
        environment.freeEnv();
        std::cout << "PASS: singleton tag, heuristic cleanup, fixing, branching, integrality, Rank-1 coefficients, checkpoint, enumeration cut buffers\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
