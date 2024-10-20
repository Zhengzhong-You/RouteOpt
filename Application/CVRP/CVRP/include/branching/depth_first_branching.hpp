//
// Created by You, Zhengzhong on 4/25/24.
//

#ifndef INCLUDE_DEPTH_FIRST_BRANCHING_DEPTH_FIRST_BRANCHING_HPP_
#define INCLUDE_DEPTH_FIRST_BRANCHING_DEPTH_FIRST_BRANCHING_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"
#include "branching.hpp"

class DepthFirstBranching : public BaseBranching {
 public:
  static StackTree bbt;
  static StackTree sub_bbt;
  static std::list<BbNode *> bbt_set;
  static std::list<BbNode *> sub_bbt_set;

  static void updateLowerBound();
  static void solve(StackTree &tree);
  static void solveOne(StackTree &tree);
  static void takeNodeOut(StackTree &tree);
  static void addNodeIn(StackTree &tree, BbNode *node);
};

#endif //INCLUDE_DEPTH_FIRST_BRANCHING_DEPTH_FIRST_BRANCHING_HPP_
