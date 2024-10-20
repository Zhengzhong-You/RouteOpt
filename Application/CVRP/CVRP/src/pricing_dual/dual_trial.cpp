//
// Created by You, Zhengzhong on 6/27/24.
//


#include "dual_trial.hpp"

int DualTrial::resolve_lp_method = SOLVER_PRIMAL_SIMPLEX;
int DualTrial::num_iteration = 0;
double DualTrial::rc_std = 0.0;