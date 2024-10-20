#ifndef CVRP_BBNODE_HPP
#define CVRP_BBNODE_HPP

#include <vector>
#include <unordered_map>
#include "macro.hpp"
#include "cuts.hpp"
#include <deque>
#include "Eigen/Sparse"
#include "solver.hpp"
#include "dynamics.hpp"

class CVRP;

class BbNode;

struct Bucket;

using sparseRowMatrixXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using sparseColMatrixXd = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using sparseRowMatrixXI = Eigen::SparseMatrix<int, Eigen::RowMajor>;
using denseRowMatrixXi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using denseColMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using RowVectorXT = Eigen::Matrix<size_t, 1, Eigen::Dynamic>;
using RowVectorXd = Eigen::Matrix<double, 1, Eigen::Dynamic>;

struct SequenceInfo {
    std::vector<int> col_seq; // exclude 0; forward sequence
    std::vector<double> main_res; // for each vertex, record the res used
    int forward_concatenate_pos{}; //record the last pos of forward partial label, == reached
    bool operator==(const SequenceInfo &other) const {
        return col_seq == other.col_seq &&
               main_res == other.main_res &&
               forward_concatenate_pos == other.forward_concatenate_pos;
    }
};

class BbNode {
public:
    bool is_terminated{};
    int size_enumeration_col_pool{};
    int valid_size{};
    int tree_level{}, index{};
    int num_edges{};
    int num_rows_in_bucket_graph{};
    int num_forward_bucket_arcs{}, num_forward_jump_arcs{};
    double last_gap{1};
    double value{};
    bool if_just_enter_enu{};
    double br_value_improved{};
    Solver solver{};
    std::vector<SequenceInfo> cols{};
    std::pair<std::vector<SequenceInfo>, std::vector<double> > only_frac_sol{}; //easier to extract sequence
    std::pair<std::vector<SequenceInfo>, std::vector<double> > all_lp_sol{}; //include integer sol
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int> >, PairHasher>
    edge_col_map{}; // col idx & visit times
    std::vector<bool> deleted_columns_in_enumeration_pool{};
    std::vector<int> edge_tail, edge_head;
    std::vector<double> edge_value{};
    Bucket **all_forward_buckets{};
    RowVectorXT index_columns_in_enumeration_column_pool;
    RowVectorXd cost_for_columns_in_enumeration_column_pool;
    std::deque<sparseRowMatrixXd> matrix_in_enumeration;
    sparseColMatrixXd basic_matrix;
    std::vector<Rcc> rccs;
    std::vector<R1c> r1cs;
    std::vector<Brc> brcs;

    void allocateMem(int num);

    explicit BbNode(int num, CVRP *cvrp);

    BbNode(BbNode *node, int idx, const Brc &bf, int num_buckets_per_vertex);

    BbNode(BbNode *node, int num_buckets_per_vertex, const bool *if_use_arc);

    BbNode(BbNode *node, int idx, const Brc &bf);

    double obtainParentNodeValue();

    virtual ~BbNode();

    Dynamics dynamics{};

    int num_backward_bucket_arcs{}, num_backward_jump_arcs{};
    Bucket **all_backward_buckets{};
    std::vector<std::vector<std::vector<int> > > topological_order_forward{};
    std::vector<std::vector<std::vector<int> > > topological_order_backward{};

    explicit BbNode(int num);

    int &getTreeLevel();

    double &getCurrentNodeVal();

    int &getNodeIdx();

    void printBrCInfo();

    double &getNodeBrValueImproved();

    bool &getIfJustEnterEnu();

    Solver &getSolver();

    int &getNumEdges();

    std::vector<double> &getEdgeValue();

    std::vector<int> &getEdgeTail();

    std::vector<int> &getEdgeHead();

    std::vector<SequenceInfo> &getCols();

    std::vector<Brc> &getBrCs();

    bool &getIfTerminated();

    int &getDynamicKNode();
};

#endif //CVRP_BBNODE_HPP
