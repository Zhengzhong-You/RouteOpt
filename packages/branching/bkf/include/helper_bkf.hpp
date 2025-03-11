/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HELPER_BKF_HPP
#define ROUTE_OPT_HELPER_BKF_HPP
namespace RouteOpt::Branching::BKF {
    double getLBTk(double omega, double alpha, double B, double k);

    double getUBTk(double omega, double alpha, double B, double k);

    double getRealTk(double omega, double alpha, double B, double k);

    void getRange(double omega,
                  double alpha,
                  double B,
                  double target,
                  double k1,
                  double est_m,
                  int &lb_k,
                  int &ub_k);

    void updateState(double new_value, double &old_value, int n);

    void updateStateWithWeights(double new_value, double &old_value, int n);

    void updateStateAverage(double new_value, double &old_value, int n);
}

#endif // ROUTE_OPT_HELPER_BKF_HPP
