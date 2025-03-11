/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file solver_cplex.hpp
 * @brief CPLEX solver interface for RouteOpt.
 * Please refer to the CPLEX C API documentation for more details.
 */

#ifndef ROUTE_OPT_SOLVER_CPLEX_HPP
#define ROUTE_OPT_SOLVER_CPLEX_HPP
#include "solver_macro.hpp"
#if SOLVER_TYPE == 1
#include <ilcplex/cplex.h>
#include <ilcplex/cplexx.h>
#include <cstdlib>
#endif

namespace RouteOpt {
#if SOLVER_TYPE == 1
    using SOLVERmodel = CPXLPptr;
    using SOLVERenv = CPXENVptr;

    constexpr auto SOLVER_EQUAL = 'E';
    constexpr auto SOLVER_GREATER_EQUAL = 'G';
    constexpr auto SOLVER_LESS_EQUAL = 'L';
    constexpr auto SOLVER_BINARY = 'B';
    constexpr auto SOLVER_INTEGER = 'I';
    constexpr auto SOLVER_CONTINUOUS = 'C';
    constexpr auto SOLVER_OPTIMAL = CPX_STAT_OPTIMAL;
    constexpr auto SOLVER_SUBOPTIMAL = CPX_STAT_FEASIBLE;
    constexpr auto SOLVER_NUMERIC = CPX_STAT_NUM_BEST;
    constexpr auto SOLVER_MIP_INFEASIBLE = CPXMIP_INFEASIBLE;
    constexpr auto SOLVER_INF_OR_UNBD = CPX_STAT_UNBOUNDED;
    constexpr auto SOLVER_UNBOUNDED = CPX_STAT_UNBOUNDED;
    constexpr auto SOLVER_PRIMAL_SIMPLEX = CPX_ALG_PRIMAL;
    constexpr auto SOLVER_DUAL_SIMPLEX = CPX_ALG_DUAL;
    constexpr auto SOLVER_BARRIER = CPX_ALG_BARRIER;
    constexpr auto SOLVER_CROSSOVER_DOWN = CPX_NONBASIC_SOLN;
    constexpr auto SOLVER_CROSSOVER_DEFAULT = 0;
    constexpr auto SOLVER_MAX_SENSE = CPX_MAX;
    constexpr auto SOLVER_TIME_LIMIT = CPX_PARAM_TILIM;
    constexpr auto SOLVER_INFINITY = CPX_INFBOUND;

#endif
}
#endif // ROUTE_OPT_SOLVER_CPLEX_HPP
