/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <iostream>
#include <cmath>
#include <limits>
#include "route_opt_macro.hpp"
#include "bkf_macro.hpp"
#include "bkf_controller.hpp"
#include "helper_bkf.hpp"

namespace RouteOpt::Branching::BKF {
    //used before entering the branching phase
    void BKFDataShared::calculateRStar(double lift, int tree_level, bool dir, int node_idx,
                                       BKFController &controller) {
        if (lift < TOLERANCE) {
            PRINT_WARNING("lift is too small: "+std::to_string(lift));
        } else if (lift < -TOLERANCE) {
            THROW_RUNTIME_ERROR("lift is even negative: "+std::to_string(lift));
        }
        lift = std::max(lift, TOLERANCE);
        if (tree_level >= r_star_depth.size()) r_star_depth.resize(tree_level + 1);
        auto &r_star = r_star_depth[tree_level];

        auto &recordings = dir ? r_star.second.second : r_star.first.second;
        auto &increase = dir ? r_star.second.first : r_star.first.first;

        auto revised_lift = lift / (1 - controller.getAlpha() / (controller.getOptK(node_idx) + 1));
        updateState(revised_lift, increase, recordings);
        ++recordings;

        double new_r_star;
        if (r_star.second.second == 0 || r_star.first.second == 0) {
            new_r_star = std::numeric_limits<float>::max();
            for (auto &r: r_star_depth) {
                if (r.first.second == 0 || r.second.second == 0) continue;
                new_r_star = std::min(new_r_star, std::sqrt(r.first.first * r.second.first));
            }
            if (new_r_star == std::numeric_limits<float>::max()) new_r_star = revised_lift;
        } else {
            new_r_star = std::sqrt(r_star.first.first * r_star.second.first);
        }
        current_r_best = new_r_star * R_DISCOUNT;
        std::cout << "current_r_best= " << current_r_best << std::endl;
    }

    void BKFController::evaluateM1() {
        // if (if_adjust_m) return;
        auto max_num = MAX_ALPHA * all_n;
        est_m = std::min(est_m, max_num);
        alpha = est_m / all_n;
        // if_adjust_m = true;
    }
}
