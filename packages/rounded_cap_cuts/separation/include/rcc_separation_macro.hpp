/*
* Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rcc_separation_macro.hpp
 * @brief Macro definitions for RCC (Rounded Cap Cuts) separation in RouteOpt.
 *
 * This header defines constants that are used in the separation process for RCCs.
 * These constants include the maximum number of cuts allowed and the maximum violation tolerance.
 */

#ifndef ROUTE_OPT_RCC_SEPARATION_MACRO_HPP
#define ROUTE_OPT_RCC_SEPARATION_MACRO_HPP

namespace RouteOpt::RCCs::Separation {
    // Maximum number of cuts allowed in the RCC separation process.
    constexpr int MAX_NUM_OF_CUTS = 100;

    // Maximum violation (tolerance) allowed when evaluating cuts.
    constexpr double MAX_VIO = 1e-4;
} // namespace RouteOpt::RCCs::Separation

#endif // ROUTE_OPT_RCC_SEPARATION_MACRO_HPP
