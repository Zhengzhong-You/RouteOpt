/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HGS_MACRO_HPP
#define ROUTE_OPT_HGS_MACRO_HPP
namespace RouteOpt::Application::CVRP {
    constexpr bool IF_FIND_ALL_SOLUTIONS = false;
    constexpr int SEED_HGS = 1;
    constexpr double MIN_RUN_TIME_HGS = 1200.;
    constexpr double DIM_DISCOUNT = 0.05;
    constexpr double POWER_FACTOR = 2;
    constexpr double A_HGS = 2.;
    constexpr double B_HGS = 10.; //ax+b
}
#endif // ROUTE_OPT_HGS_MACRO_HPP
