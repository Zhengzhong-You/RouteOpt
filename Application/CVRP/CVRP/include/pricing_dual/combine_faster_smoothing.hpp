//
// Created by You, Zhengzhong on 5/7/24.
//

#ifndef INCLUDE_PRICING_DUAL_COMBINE_FASTER_SMOOTHING_HPP_
#define INCLUDE_PRICING_DUAL_COMBINE_FASTER_SMOOTHING_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"
#include "faster_dual.hpp"
#include "dual_smoothing.hpp"

class CombineFasterSmoothing {
 public:
  static CVRP *cvrp;
  static BbNode *node;
  static double beta;
  static double gamma;
  static double lp2_value;
  static double lp2_norm;
  static std::vector<double> lp2_obj;
  static std::vector<double> pi_lp1;
  static std::vector<double> pi_lp2;

  static void init(CVRP *pr_cvrp);
  static void updateNode(BbNode *pr_node);
  static void getLP2Obj(const std::vector<double> &lp2_obj);
  static void getLP2Value(double &lp2_value);
  static void updatePiLP1();
  static void updatePiLP2();
  static void updateSmoothingState();
};

#endif //INCLUDE_PRICING_DUAL_COMBINE_FASTER_SMOOTHING_HPP_
