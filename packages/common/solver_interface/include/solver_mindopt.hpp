/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file solver_mindopt.hpp
 * @brief Header file for the MindOpt solver interface.
 *
 * This header defines constants and types specific to the Mindo solver, including variable types,
 * constraint types, and solver status codes. (Still in development)
 */

#ifndef ROUTE_OPT_SOLVER_MINDOPT_HPP
#define ROUTE_OPT_SOLVER_MINDOPT_HPP
#include "solver_macro.hpp"
#if SOLVER_TYPE==2
#include "Milp\LinearModel\LinearModel.hpp"
#endif

namespace RouteOpt {
#if SOLVER_TYPE == 2
    using VARTYPE = Mindo::VariableType;
    constexpr auto SOLVER_EQUAL = 0;
    constexpr auto SOLVER_GREATER_EQUAL = 2;
    constexpr auto SOLVER_LESS_EQUAL = 1;
    constexpr auto SOLVER_BINARY = 2;
    constexpr auto SOLVER_INTEGER = 1;
    constexpr auto SOLVER_CONTINUOUS = 0;
    constexpr auto SOLVER_MIP_INFEASIBLE = 1;
    constexpr auto SOLVER_UNBOUNDED = 2;
    constexpr auto SOLVER_INF_OR_UNBD = 3;
    constexpr auto SOLVER_PRIMAL_SIMPLEX = 0;
    constexpr auto SOLVER_DUAL_SIMPLEX = 1;
    constexpr auto SOLVER_MAX_SENSE = -1;
    constexpr auto SOLVER_TIME_LIMIT = 9;
#endif
}
#endif // ROUTE_OPT_SOLVER_MINDOPT_HPP
