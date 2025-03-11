/*
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file bkf_macro.hpp
 * @brief Macros and constants for the BKF (Best K Formula) algorithm.
 *
 * This header defines constants and macros used in the BKF algorithm, including tolerances,
 * maximum iterations, discount factors, and verbosity settings.
 */

#ifndef ROUTE_OPT_BKF_MACRO_HPP
#define ROUTE_OPT_BKF_MACRO_HPP

namespace RouteOpt::Branching::BKF {
    // BOOST_TOLERANCE_BIT: Bit-level tolerance used in boosting operations,
    // e.g., for adjusting precision in calculations.
    constexpr int BOOST_TOLERANCE_BIT = 6;

    // MAX_ITERATION: Maximum number of iterations allowed for solving equations in BKF computations.
    constexpr int MAX_ITERATION = 1000;

    // R_DISCOUNT: Discount factor applied to a parameter "R" during BKF calculations.
    constexpr double R_DISCOUNT = 0.72;

    // MAX_ALPHA: Maximum alpha value allowed in BKF computations.
    constexpr double MAX_ALPHA = 0.8;

    // Define BKF_VERBOSE to enable verbose execution in BKF debugging.
#define BKF_VERBOSE

#ifdef BKF_VERBOSE
    // BKF_VERBOSE_EXEC: Executes the provided code block when BKF_VERBOSE is defined.
    // This macro is useful for inserting debug or verbose logging code.
#define BKF_VERBOSE_EXEC(...) __VA_ARGS__;
#else
    // When BKF_VERBOSE is not defined, BKF_VERBOSE_EXEC does nothing.
    BKF_VERBOSE_EXEC();
#endif
} // namespace RouteOpt::Branching::BKF

#endif // ROUTE_OPT_BKF_MACRO_HPP
