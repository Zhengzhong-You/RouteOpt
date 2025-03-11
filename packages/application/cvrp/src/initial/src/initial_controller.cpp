/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <iostream>
#include "cvrp_macro.hpp"
#include "vrptw_macro.hpp"
#include "solver.hpp"
#include "initial_controller.hpp"
#include "initial_macro.hpp"

namespace RouteOpt::Application::CVRP {
    double transformCost(double x) {
        if constexpr (app_type == APPLICATION_TYPE::VRPTW) {
            return (std::floor(10 * (x + VRPTW_DISTANCE_TOLERANCE))) / 10;
        }
        return std::floor(x + 0.5);
    }

    void InitialController::initialProcessing() {
        cost_mat4_vertex_ref.get().resize(dim, std::vector<double>(dim, 0));

        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                cost_mat4_vertex_ref.get()[i][j] = transformCost(
                    std::sqrt(static_cast<float>((info_vertex_ref.get()[i][1] - info_vertex_ref.get()[j][1]) * (
                                                     info_vertex_ref.get()[i][1] - info_vertex_ref.get()[j][1]) +
                                                 (info_vertex_ref.get()[i][2] - info_vertex_ref.get()[j][2]) * (
                                                     info_vertex_ref.get()[i][2] - info_vertex_ref.get()[j][2]))));
            }
        }

        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                cost_mat4_vertex_ref.get()[j][i] = cost_mat4_vertex_ref.get()[i][j];
            }
        }


        demand_ref.get().resize(dim);
        for (int i = 0; i < dim; ++i) {
            demand_ref.get()[i] = info_vertex_ref.get()[i][3];
        }

        getLowerBoundofMinimumNumberCarsRef();
        setSolverEnv();
    }

    void InitialController::setSolverEnv() const {
        auto &solver = solver_ref.get();
        SAFE_SOLVER(solver.loadEnv(nullptr))
        SAFE_SOLVER(solver.setEnvThreads(NUM_THREADS_LP, true))
        SAFE_SOLVER(solver.setEnvOutputFlag(0, true))
        SAFE_SOLVER(solver.setEnvInfUnbdInfo(1, true))
        SAFE_SOLVER(solver.setEnvMIPGap(MIP_GAP_TOLERANCE, true))
        SAFE_SOLVER(solver.setEnvFeasibilityTol(FeasibilityTol, true))
    }
}
