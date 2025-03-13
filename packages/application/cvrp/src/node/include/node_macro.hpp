/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_NODE_MACRO_HPP
#define ROUTE_OPT_NODE_MACRO_HPP

namespace RouteOpt::Application::CVRP {
    constexpr double OBJ_ARTIFICIAL = 1e6;
    constexpr double LeftThresholdRCFixing4EnumerationPool = 0.6;
    constexpr double NodeLPDensityEstimation = 0.1;
}
#endif // ROUTE_OPT_NODE_MACRO_HPP
