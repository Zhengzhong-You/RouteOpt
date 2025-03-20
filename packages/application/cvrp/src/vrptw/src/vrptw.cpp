/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "vrptw.hpp"
#include "vrptw_macro.hpp"

namespace RouteOpt::Application::CVRP {
    void VRPTW::checkSolutionFeasibility(const std::vector<double> &X,
                                         const std::vector<SequenceInfo> &cols,
                                         bool &feasible) {
        feasible = true;
        if (cols.size() == 1) {
            feasible = false;
            return;
        }
        for (auto &route: cols) {
            double local_cap = 0;
            for (auto &i: route.col_seq) {
                local_cap += getDemand()[i];
            }
            if (local_cap > getCap() + TOLERANCE) {
                std::cout << "cap for route ";
                for (auto &i: route.col_seq) {
                    std::cout << i << " ";
                }
                std::cout << " is " << local_cap << ", while cap= " << getCap() << std::endl;
                feasible = false;
                break;
            }
        }
    }
}
