/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "vrptw.hpp"
#include "vrptw_macro.hpp"

namespace RouteOpt::Application::CVRP {
    bool VRPTW::checkRouteTimeWindowFeasibility(const SequenceInfo &route) const {
        const auto &earliest_time = getEarliestTime();
        const auto &latest_time = getLatestTime();
        const auto &service_time = getServiceTime();
        const auto &cost_mat4_vertex = getDisMat();
        double now = 0;
        int prev = 0;
        for (auto node: route.col_seq) {
            now += service_time[prev] + cost_mat4_vertex[prev][node];
            now = std::max(now, earliest_time[node]);
            if (now > latest_time[node] + TOLERANCE) {
                return false;
            }
            prev = node;
        }
        now += service_time[prev] + cost_mat4_vertex[prev][0];
        return now <= latest_time[0] + TOLERANCE;
    }

    void VRPTW::checkSolutionFeasibility(const std::vector<SequenceInfo> &cols,
                                         bool &feasible) {
        feasible = true;
        for (const auto &route: cols) {
            if (!checkRouteDemandFeasibility(route)) {
                feasible = false;
                break;
            }
            if (!checkRouteTimeWindowFeasibility(route)) {
                feasible = false;
                break;
            }
        }
    }
}
