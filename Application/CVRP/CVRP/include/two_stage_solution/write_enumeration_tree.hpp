//
// Created by You, Zhengzhong on 4/25/24.
//

#ifndef INCLUDE_TWO_STAGE_SOLUTION_WRITE_ENUMERATION_TREE_HPP_
#define INCLUDE_TWO_STAGE_SOLUTION_WRITE_ENUMERATION_TREE_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"

class WriteEnumerationTree {
 public:
  static CVRP *cvrp;
  static BbNode *node;
  static std::unordered_map<yzzLong, std::tuple<std::vector<int>, double, size_t>> enumeration_col_idx;

  static void init(CVRP *cvrp);
  static void updateNode(BbNode *node);
  static void writeEnuTree();
  static void writeEnuCols();
  static void writeEnuCuts(const std::vector<std::pair<size_t, double>> &col_idx_cost);
  static void populateEnuCols(std::vector<std::pair<size_t, double>> &col_idx_cost,
							  const std::vector<SequenceInfo> &lp_col_info,
							  const std::vector<std::pair<size_t, double>> &enu_col_idx_cost);
};

#endif //INCLUDE_TWO_STAGE_SOLUTION_WRITE_ENUMERATION_TREE_HPP_
