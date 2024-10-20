#include <utility>
#include "cvrp.hpp"
#include "branching.hpp"
#include "dynamics.hpp"
#include "machine_learning.hpp"
#include "train.hpp"
using namespace std;
using namespace std::chrono;
using namespace Eigen;


void BaseBranching::testCG(bool if_exact_CG, bool if_record_product_val,
                           bool if_record_improvement, bool if_force_complete) {
    if (current_branching_info.branch_pair.size() == 1) {
        verbose_call(cout << "testCG: current_branching_info.branch_pair.size() == 1, return!" << endl;)
        dynamic_call(node->getDynamicKNode() = Dynamics::opt_k;)
        return;
    }
    auto beg = high_resolution_clock::now();
    int BeforeNumRow = cvrp->getNumRow();
    auto org_val = node->getCurrentNodeVal();
    int cnt = 0;
    std::vector<std::pair<std::pair<int, int>, double> >
            Branch_Val_consider_ub(current_branching_info.branch_pair.size()),
            Branch_Val(current_branching_info.branch_pair.size());

    Brc bf;
    bf.idx_br_c = cvrp->getNumRow();
    node->getBrCs().emplace_back(bf);
    int col_start = cvrp->getNumCol();

    auto pseudo_down = &branching_history.heuristic_improvement_down;
    auto pseudo_up = &branching_history.heuristic_improvement_up;
    if (if_exact_CG) {
        pseudo_down = &branching_history.exact_improvement_down;
        pseudo_up = &branching_history.exact_improvement_up;
    }

    cvrp->preprocess4TestCG(node);
    vector<int> solver_ind;
    vector<double> solver_val;
    vector<int> const_for_branching(cvrp->getNumCol());
    iota(const_for_branching.begin(), const_for_branching.end(), cvrp->getNumCol());
    bool if_changed = false;

    ++cvrp->getNumRow();

    double max_product = 0;
    for (auto &edge: current_branching_info.branch_pair) {
        if (!if_force_complete) if (max_product > current_branching_info.branch_lp[edge]) continue;
        cout << MID_PHASE_SEPARATION;
        cout << "Evaluate on ( " << edge.first << " , " << edge.second << " )...\n";
        double temp_val;
        cvrp->getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
        if (!if_changed) {
            safe_solver(cvrp->addBranchConstraint(solver_ind,
                solver_val,
                SOLVER_EQUAL,
                0,
                nullptr,
                node->getSolver()))
            if_changed = true;
        } else {
            cvrp->changeBranchConstraint(solver_ind,
                                         solver_val,
                                         SOLVER_EQUAL,
                                         0,
                                         BeforeNumRow,
                                         node->getSolver());
        }
        safe_solver(node->getSolver().updateModel())
        node->getBrCs().back().edge = edge;
        node->getBrCs().back().br_dir = false;
        if (cvrp->getIfInEnuState()) {
            cvrp->addBranchConstraint2ColPoolInEnumByColMap(node, edge);
            cvrp->solveLPByInspection(node, true, !if_exact_CG, false);
        } else {
            if_test_cg = true;
            cvrp->solveLPInLabeling(node, true, if_exact_CG, false);
            if_test_cg = false;
        }
        node->getIfTerminated() = false;
        safe_solver(node->getSolver().getObjVal(&temp_val))
        auto dif1 = calculateDifference(temp_val, org_val);
        ml_call(MASTER_VALVE_ML != 0, MachineLearning::recordDiscrepancyLongInfo(edge, dif1, true))
        ml_call(MASTER_VALVE_ML == ML_GET_DATA_2, GetTrainingData::findCGScore4OneEdge(edge, dif1, true))
        auto dif1_consider_ub = calculateDifference(temp_val, org_val, true);
        cout << SMALL_PHASE_SEPARATION;
        int len = cvrp->getNumCol() - col_start;
        if (const_for_branching.size() < len) {
            int arr_beg = (int) const_for_branching.size();
            int arr_val = const_for_branching.back() + 1;
            const_for_branching.resize(len);
            iota(const_for_branching.begin() + arr_beg, const_for_branching.end(), arr_val);
        }
        cvrp->rmLPCols(node, vector<int>(const_for_branching.begin(), const_for_branching.begin() + len));
        safe_solver(cvrp->inverseLastBranchConstraint(SOLVER_EQUAL, 1, node->getSolver()))
        node->getBrCs().back().br_dir = true;
        safe_solver(node->getSolver().updateModel())
        safe_solver(node->getSolver().getNumCol(&cvrp->getNumCol()))
        if (cvrp->getIfInEnuState()) {
            cvrp->solveLPByInspection(node, true, !if_exact_CG, false);
        } else {
            if_test_cg = true;
            cvrp->solveLPInLabeling(node, true, if_exact_CG, false);
            if_test_cg = false;
        }
        node->getIfTerminated() = false;
        safe_solver(node->getSolver().getObjVal(&temp_val))
        auto dif2 = calculateDifference(temp_val, org_val);
        ml_call(MASTER_VALVE_ML != 0, MachineLearning::recordDiscrepancyLongInfo(edge, dif2, false))
        ml_call(MASTER_VALVE_ML == ML_GET_DATA_2, GetTrainingData::findCGScore4OneEdge(edge, dif2, false))
        auto dif2_consider_ub = calculateDifference(temp_val, org_val, true);
        auto product = dif1 * dif2;
        auto product_consider_ub = dif1_consider_ub * dif2_consider_ub;
        cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  real_pd= " << setw(6)
                << left << product << " if consider ub: " << dif1_consider_ub << " " << dif2_consider_ub << " "
                << product_consider_ub << endl;
        len = cvrp->getNumCol() - col_start;
        if (const_for_branching.size() < len) {
            int arr_beg = (int) const_for_branching.size();
            int arr_val = const_for_branching.back() + 1;
            const_for_branching.resize(len);
            iota(const_for_branching.begin() + arr_beg, const_for_branching.end(), arr_val);
        }
        cvrp->rmLPCols(node, vector<int>(const_for_branching.begin(), const_for_branching.begin() + len));
        safe_solver(node->getSolver().updateModel())
        safe_solver(node->getSolver().getNumCol(&cvrp->getNumCol()))
        Branch_Val[cnt] = {edge, product};
        Branch_Val_consider_ub[cnt++] = {edge, product_consider_ub};
        if (product_consider_ub > max_product) {
            max_product = product_consider_ub;
        }
        if (if_record_improvement) {
            (*pseudo_down)[edge].first += dif1;
            ++(*pseudo_down)[edge].second;
            (*pseudo_up)[edge].first += dif2;
            ++(*pseudo_up)[edge].second;
        }
    }
    Branch_Val.resize(cnt);
    Branch_Val_consider_ub.resize(cnt);
    safe_solver(node->getSolver().delConstraints(1, &BeforeNumRow))
    safe_solver(node->getSolver().updateModel())
    --cvrp->getNumRow();

    cvrp->postprocess4TestCG(node);
    safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))

    if (if_record_product_val) {
        sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
            return a.second > b.second;
        });
        current_branching_info.branch_pair_val = Branch_Val;
        auto it =
                max_element(Branch_Val_consider_ub.begin(), Branch_Val_consider_ub.end(),
                            [](const auto &a, const auto &b) {
                                return a.second < b.second;
                            });
        current_branching_info.branch_pair = {it->first};
        if (it->first != Branch_Val[0].first) cout << "The final branching decision consider upper bound!" << endl;
    } else {
        auto it =
                max_element(Branch_Val_consider_ub.begin(), Branch_Val_consider_ub.end(),
                            [](const auto &a, const auto &b) {
                                return a.second < b.second;
                            });
        current_branching_info.branch_pair = {it->first};
    }
    node->getBrCs().pop_back();
    dynamic_call(node->getDynamicKNode() = Dynamics::opt_k;)
    /**
     * is_terminated, value. these two cannot be deleted! since when exact cg could determine if this node is no promising to be solved
     */
    node->getIfTerminated() = false;
    node->getCurrentNodeVal() = org_val;
    auto end = high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    verbose_call(cout << "testCG spent= " << eps << "s" << endl)
    dynamic_call(Dynamics::solve_cg_time = eps;
        Dynamics::getAverageT4LPNHeuristic(eps / cnt / 2, false);
    )
}
