/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <iostream>

#include "cvrp_macro.hpp"
#include "hgs_macro.hpp"
#include "hgs_controller.hpp"
#include "my_hgs.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Application::CVRP {
    void HGSController::runHGS() {
        if constexpr (app_type != APPLICATION_TYPE::CVRP) return;
        double time = std::max(MIN_RUN_TIME_HGS, B_HGS + A_HGS * std::pow((dim * DIM_DISCOUNT), POWER_FACTOR));
        if constexpr (IF_FIND_ALL_SOLUTIONS) {
            std::vector<std::pair<std::vector<std::vector<int> >, double> > sol;
            HGSWithAllSolution(f_name, true, SEED_HGS, time, sol);
            for (auto &route: sol) {
                for (auto &r: route.first) {
                    if (r.empty()) continue;
                    ip_opt_sol_ref.get().emplace_back(r);
                }
            }
        } else {
            std::vector<std::vector<int> > sol;
            double heur_UB;
            HGS(f_name, true, SEED_HGS, time, sol, heur_UB);
            if (ub_ref.get() > heur_UB + TOLERANCE) {
                ub_ref.get() = heur_UB;
                std::cout << "ub updated by HGS: " << ub_ref.get() << std::endl;
                ip_opt_sol_ref.get().reserve(sol.size());
                for (auto &route: sol) {
                    if (route.empty()) continue;
                    ip_opt_sol_ref.get().emplace_back(route);
                }
            }
        }
    }
}
