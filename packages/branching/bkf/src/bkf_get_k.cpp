/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "bkf_macro.hpp"
#include "bkf_controller.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Branching::BKF {
    int BKFController::getBestK(const BKFDataShared &sharedData, double ub,
                                double lb) {
        int num = static_cast<int>(all_n);
        if (!if_init) {
            BKF_VERBOSE_EXEC(std::cout << "we directly test "<< num << std::endl;)
        QUIT:
            if_init = true;
            return num;
        }

        num = static_cast<int>(est_m);

        double beta;
        auto r_best = sharedData.getCurrentRBest();
        if (r_best < TOLERANCE) goto QUIT;
        beta = std::max((ub - lb - sharedData.getF()) / r_best, 0.);
        if (beta > 16) goto QUIT;


        double t_for_one_testing = tmp_if_b4 ? t_for_one_testing_b4 : t_for_one_testing_af;
        double c_for_node = tmp_if_b4 ? c_for_node_b4 : c_for_node_af;
        if (t_for_one_testing < TOLERANCE) goto QUIT;
        double omega = std::max(c_for_node / t_for_one_testing, 1.);
        double B = std::pow(2, beta);

        double alpha_beta = alpha * beta * std::log(std::sqrt(2));
        double k1 = std::sqrt(alpha_beta * (alpha_beta + 2 * alpha + 2 * omega - 2)) + alpha_beta + alpha - 1;
        auto k1_up = static_cast<int>(std::max(std::min(std::ceil(k1), static_cast<double>(est_m)), 1.));
        auto k1_down = static_cast<int>(std::max(std::min(std::floor(k1), static_cast<double>(est_m)), 1.));
        BKF_VERBOSE_EXEC(std::cout << "alpha= " << alpha << " beta= " << beta << " omega= " << omega;)

        if (k1_up == k1_down) {
            k1 = k1_up;
        } else {
            double k1_up_val = getLBTk(omega, alpha, B, k1_up);
            double k1_down_val = getLBTk(omega, alpha, B, k1_down);
            if (k1_up_val < k1_down_val) {
                k1 = k1_up;
            } else {
                k1 = k1_down;
            }
        }
        BKF_VERBOSE_EXEC(std::cout << " k1= " << k1 << std::endl)

        double ub_tk = getUBTk(omega, alpha, B, k1);
        int lb_k, ub_k;
        getRange(omega, alpha, B, ub_tk, k1, est_m, lb_k, ub_k);
        BKF_VERBOSE_EXEC(std::cout << "lb_k= " << lb_k << " ub_k= " << ub_k << std::endl)
        if (lb_k == ub_k) {
            num = lb_k;
        } else {
            double real_tk = getRealTk(omega, alpha, B, k1);
            getRange(omega, alpha, B, real_tk, k1, est_m, lb_k, ub_k);
            BKF_VERBOSE_EXEC(std::cout << "revise: lb_k= " << lb_k << " ub_k= " << ub_k << std::endl)
            if (lb_k == ub_k) {
                num = lb_k;
            } else {
                std::vector<double> tmp;
                for (int i = lb_k; i <= ub_k; ++i) {
                    tmp.emplace_back(getRealTk(omega, alpha, B, i));
                }
                auto it = std::min_element(tmp.begin(), tmp.end()) - tmp.begin() + lb_k;
                num = static_cast<int>(it);
            }
        }
        num = std::min(num, static_cast<int>(est_m));
        BKF_VERBOSE_EXEC(
            std::cout << "opt_k= " << num << std::endl;)
        goto QUIT;
    }
}
