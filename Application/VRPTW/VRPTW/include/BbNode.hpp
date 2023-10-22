
#ifndef CVRP_BBNODE_HPP
#define CVRP_BBNODE_HPP

#include <vector>
#include <unordered_map>
#include "MACRO.hpp"
#include "Cuts.hpp"
#include <deque>
#include <Eigen/Sparse>
#include "Solver.hpp"

class CVRP;

class BbNode;

struct Bucket;

using sparseRowMatrixXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using RowVectorXT = Eigen::Matrix<size_t, 1, Eigen::Dynamic>;
using RowVectorXd = Eigen::Matrix<double, 1, Eigen::Dynamic>;

struct BidirectionalLinkedList {
  BidirectionalLinkedList *p_node{nullptr}, *l_node{nullptr}, *r_node{nullptr};
  std::unordered_map<int, std::vector<int>> edge_to_cols{};

  void becomeParent(BbNode *node, CVRP *cvrp);
  void giveBirth(BidirectionalLinkedList *&lnode, BidirectionalLinkedList *&rnode);
  void deleteSelf();
  explicit BidirectionalLinkedList(BidirectionalLinkedList *pnode) { p_node = pnode; }
  ~BidirectionalLinkedList() = default;
};

class BbNode {
 public:
  bool is_integer{};
  bool is_terminated{};
  int size_enumeration_col_pool{};
  int valid_size{};
  int num_parent_cols{};
  int num_parent_cols_in_lp_solutions{};
  int tree_level{}, index{};
  int num_edges{};
  int num_rows_in_bucket_graph{};
  int num_forward_bucket_arcs{}, num_forward_jump_arcs{};
  double last_gap{1};
  double value{};
  Solver solver{};
  std::vector<bool> deleted_columns_in_enumeration_pool{};
  std::vector<int> edge_tail, edge_head;
  std::vector<double> edge_value{};
  std::vector<size_t> index_columns{};
  BidirectionalLinkedList *ptr{};
  Bucket **all_forward_buckets{};
  RowVectorXT index_columns_in_enumeration_column_pool;
  RowVectorXd cost_for_columns_in_enumeration_column_pool;
  std::deque<sparseRowMatrixXd> matrix_in_enumeration;
  std::vector<std::pair<size_t, double>> index_for_lp_solutions_in_column_pool;
  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher>
	  column_pool_mapping;
  std::vector<Rcc> rccs;
  std::vector<R1c> r1cs;
  std::vector<R1cMulti> r1cs_multi;
  std::vector<Brc> brcs;

  void allocateMem(int num);

  explicit BbNode(int num, int p_col, CVRP *cvrp);

  BbNode(BbNode *node, BidirectionalLinkedList *ptr, int p_col, int idx, const Brc &bf, int num_buckets_per_vertex);

  BbNode(BbNode *node, int num_buckets_per_vertex, int num_col, const bool *if_use_arc);

  BbNode(BbNode *node, int p_col, int idx, const Brc &bf);

  ~BbNode();

#ifdef USE_M_DYNAMICS
  double t_for_one_lp{};
  double geo_r_star{1};
  double c{1};
  std::unordered_map<std::pair<int, int>, std::tuple<double, double, int>, PairHasher> objective_change;
  double l_r_ratio{1};
  static void updateState(double new_value, double &old_value, int n);
  void calculateRStar(double lift, double &new_r_star, CVRP *cvrp);
#endif
  int num_backward_bucket_arcs{}, num_backward_jump_arcs{};
  Bucket **all_backward_buckets{};
#ifdef READ_ENUMERATION_TREES
  explicit BbNode(int num);
#endif
};

#endif //CVRP_BBNODE_HPP
