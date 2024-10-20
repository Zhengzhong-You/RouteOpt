//
// Created by You, Zhengzhong on 4/25/24.
//

#ifndef INCLUDE_TWO_STAGE_SOLUTION_READ_ENUMERATION_TREE_HPP_
#define INCLUDE_TWO_STAGE_SOLUTION_READ_ENUMERATION_TREE_HPP_

#include "cvrp.hpp"

class ReadEnumerationTree {
 public:
  static CVRP *cvrp;
  static BbNode *node;
  static std::string tree_path;
  static std::string col_pool_path;

  static void init(CVRP *cvrp);
  static void getPath();

  static void restoreModel();
  static void readEnumTree(std::vector<std::pair<size_t, double>> &col_idx);
  static void readColumnPool(const std::vector<std::pair<size_t, double>> &col_idx);
  static void createNodeModel();
};

#endif //INCLUDE_TWO_STAGE_SOLUTION_READ_ENUMERATION_TREE_HPP_
