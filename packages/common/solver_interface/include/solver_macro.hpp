/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file solver_macro.hpp
 * @brief Macros and utility functions for RouteOpt's solver interface.
 *
 * This header defines macros and utility functions used in the RouteOpt solver interface,
 * including error handling and solver-specific configurations.
 */

#ifndef ROUTE_OPT_SOLVER_MACRO_HPP
#define ROUTE_OPT_SOLVER_MACRO_HPP
#include <iostream>
#include <cstdlib>

namespace RouteOpt {
#define SOLVER_TYPE 0

#define SAFE_SOLVER(call) {                                      \
int solver_err = (call);                                         \
if (solver_err != 0) {                                           \
std::cerr << "\x1b[91m" << __FILE__ << ":" << __LINE__          \
<< ": error in " << #call                                       \
<< "\nERROR CODE = " << solver_err << "\x1b[0m" << std::endl;   \
std::exit(EXIT_FAILURE);                                         \
}                                                                \
}
}

#endif // ROUTE_OPT_SOLVER_MACRO_HPP
