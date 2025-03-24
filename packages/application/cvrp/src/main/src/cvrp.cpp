/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include "cvrp.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRPSolver::getLowerBoundofMinimumNumberCars() {
        double sum_demand = std::accumulate(demand.data() + 1, demand.data() + dim, 0.0);
        int cap_k = std::ceil(sum_demand / cap);
        num_vehicle = cap_k;
        max_num_vehicle = dim - 1;
    }


    void CVRPSolver::initSolver() {
        pricing_controller.initLabelingMemory();
    }


    void CVRPSolver::printOptSol(std::ostream &os, int num_nodes, double lower_bound) {
        os << BIG_PHASE_SEPARATION;
        os << "\x1b[32m<Instance: " << ins_name << ">\n"
                << "<Solution>\n";

        for (const auto &route: ip_opt_sol) {
            if (route.empty())
                THROW_RUNTIME_ERROR("empty route found in optimal solution");
            for (int i = 0; i < route.size() - 1; i++) {
                os << route[i] << "-";
            }
            os << route.back() << "\n";
        }

        lower_bound = std::min(ceilTransformedNumberRelated(lower_bound), ub);

        double global_gap = 100.0 * (1.0 - lower_bound / ub);

        os << "<UB= " << ub << ">\n"
                << "<LB= " << lower_bound << ">\n"
                << "<Elapsed Time= " << roundTo(glob_timer.getTime()) << ">\n"
                << "<Nodes Explored= " << num_nodes << ">\n"
                << "<Global Gap= " << roundTo(global_gap, 3) << "%>\x1b[0m\n";
        os << BIG_PHASE_SEPARATION;

        if constexpr (IF_WRITE_NODE_OUT) {
            TwoStageController::updateUB(ins_name, ub);
            if (!tree_path.empty()) {
                auto node_name = std::string(NODE_FOLDER) + "/" + tree_path;
                TwoStageController::deleteInFile(node_name);
            }
        }
    }
}
