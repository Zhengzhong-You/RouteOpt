//
// Created by You, Zhengzhong on 5/1/24.
//

#ifndef INCLUDE_PRICING_DUAL_FASTER_DUAL_HPP_
#define INCLUDE_PRICING_DUAL_FASTER_DUAL_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"

class FasterDual {
 public:
  static CVRP *cvrp;
  static BbNode *node;
  static Solver sub_pricing_solver;
  static void init(CVRP *pr_cvrp);
  static void updateNode(BbNode *pr_node);
  static void reviseSubPricingModel();
  static void deleteCols4SubPricingModel(std::vector<int> &col_idx);
  static void changeModel4BetterDual();
  static void freeSubPricingModel();
  static void addCols2SubPricingModel(int ccnt_cnt, size_t nzcnt, std::vector<size_t> &solver_beg,
									  std::vector<int> &solver_ind, std::vector<double> &solver_val,
									  std::vector<double> &solver_obj);
};

#endif //INCLUDE_PRICING_DUAL_FASTER_DUAL_HPP_
