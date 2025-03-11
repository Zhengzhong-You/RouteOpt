/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */


/*
 * @file solver_grb.hpp
 * @brief Gurobi solver interface for RouteOpt.
 * Please refer to the Gurobi C API documentation for more details.
 */

#ifndef ROUTE_OPT_SOLVER_GRB_HPP
#define ROUTE_OPT_SOLVER_GRB_HPP
#include "solver_macro.hpp"
#include "gurobi_c.h"

namespace RouteOpt {
#if SOLVER_TYPE == 0
    using SOLVERmodel = GRBmodel;
    using SOLVERenv = GRBenv;
    constexpr auto SOLVER_EQUAL = GRB_EQUAL;
    constexpr auto SOLVER_GREATER_EQUAL = GRB_GREATER_EQUAL;
    constexpr auto SOLVER_LESS_EQUAL = GRB_LESS_EQUAL;
    constexpr auto SOLVER_BINARY = GRB_BINARY;
    constexpr auto SOLVER_INTEGER = GRB_INTEGER;
    constexpr auto SOLVER_CONTINUOUS = GRB_CONTINUOUS;
    constexpr auto SOLVER_OPTIMAL = 2;
    constexpr auto SOLVER_SUBOPTIMAL = 13;
    constexpr auto SOLVER_NUMERIC = 12;
    constexpr auto SOLVER_MIP_INFEASIBLE = 3;
    constexpr auto SOLVER_INF_OR_UNBD = 4;
    constexpr auto SOLVER_UNBOUNDED = 5;
    constexpr auto SOLVER_PRIMAL_SIMPLEX = 0;
    constexpr auto SOLVER_DUAL_SIMPLEX = 1;
    constexpr auto SOLVER_BARRIER = 2;
    constexpr auto SOLVER_CROSSOVER_DOWN = 0;
    constexpr auto SOLVER_CROSSOVER_DEFAULT = -1;
    constexpr auto SOLVER_MAX_SENSE = -1;
    constexpr auto SOLVER_TIME_LIMIT = 9;
    constexpr auto SOLVER_INFINITY = GRB_INFINITY;
#endif
}

#endif // ROUTE_OPT_SOLVER_GRB_HPP
