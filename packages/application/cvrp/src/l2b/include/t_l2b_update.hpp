/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_UPDATE_HPP
#define ROUTE_OPT_T_L2B_UPDATE_HPP
#include <vector>
#include <iostream>
#include "route_opt_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::recordEdgeLongInfo(
        const std::unordered_map<BrCType, double, Hasher> &edge_map) {
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                auto &e = edge_long_info[{i, j}].aver_edge_lp;
                if (edge_map.find({i, j}) != edge_map.end()) {
                    e.first += edge_map.at({i, j});
                }
                ++e.second;
            }
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::recordDiscrepancyLongInfo(const BrCType &candidate, double cg_change,
                                                                           bool if_left) {
        auto &e =
                if_left
                    ? edge_long_info[candidate].aver_exact_lp_discrepancy_down
                    : edge_long_info[candidate].aver_exact_lp_discrepancy_up;
        auto lp = if_left ? edge_lp_change[candidate].first : edge_lp_change[candidate].second;
        double dis;
        if (lp < TOLERANCE) {
            PRINT_REMIND("lp change is near zero= " + std::to_string(lp));
            dis = 0;
        } else {
            dis = 1 - cg_change / lp;
        }
        e.first += dis;
        ++e.second;
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::cleanLastData() {
        edge_tmp_info.clear();
        edge_lp_change.clear();
    }
}

#endif // ROUTE_OPT_T_L2B_UPDATE_HPP
