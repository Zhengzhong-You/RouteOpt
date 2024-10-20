//
// Created by You, Zhengzhong on 8/2/24.
//

#ifndef INCLUDE_SOLVER_INTERFACE_MINDOPT_HPP_
#define INCLUDE_SOLVER_INTERFACE_MINDOPT_HPP_
#include "technique_control.hpp"

#if SOLVER_TYPE == 2
#include "Milp\LinearModel\LinearModel.hpp"
#define SOLVER_EQUAL 0
#define SOLVER_GREATER_EQUAL 2
#define SOLVER_LESS_EQUAL 1
#define SOLVER_BINARY 2
#define SOLVER_INTEGER 1
#define SOLVER_CONTINUOUS 0
#define VARTYPE Mindo::VariableType
#define SOLVER_MIP_INFEASIBLE 1
#define SOLVER_UNBOUNDED 2
#define SOLVER_INF_OR_UNBD 3
#define SOLVER_PRIMAL_SIMPLEX 0
#define SOLVER_DUAL_SIMPLEX 1
#define SOLVER_MAX_SENSE (-1)
#define SOLVER_TIME_LIMIT 9
#endif
#endif //VRPTW_INCLUDE_SOLVER_INTERFACE_MINDOPT_HPP_
