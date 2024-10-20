//
// Created by You, Zhengzhong on 6/27/24.
//

#ifndef INCLUDE_PRICING_DUAL_DUAL_TRIAL_HPP_
#define INCLUDE_PRICING_DUAL_DUAL_TRIAL_HPP_

#include "macro.hpp"
#include "solver_interface/solver.hpp"

class DualTrial {
 public:
  static int resolve_lp_method;
  static int num_iteration;
  static double rc_std;
};

#endif //INCLUDE_PRICING_DUAL_DUAL_TRIAL_HPP_
