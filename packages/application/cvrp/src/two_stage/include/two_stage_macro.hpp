/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_TWO_STAGE_MACRO_HPP
#define ROUTE_OPT_TWO_STAGE_MACRO_HPP
#include <string_view>

namespace RouteOpt::Application::CVRP {
    constexpr std::string_view NODE_FOLDER = "nodes";
    constexpr std::string_view NODE_FILE_SUFFIX = ".node";
    constexpr std::string_view UB_FOLDER = "ubs";
    constexpr std::string_view UB_FILE_SUFFIX = ".ub";
    constexpr double READ_MEMORY_RATIO = 2.;

    enum class COL_POOL_TYPE {
        TINY = 200000, // keep solving
        SMALL = 400000, // write small memory out
        MID = 800000, // write middle memory out
        LARGE = 1000000 // write large memory out
    };

    enum class OUT_NODE_MEMORY_USE {
        SMALL = 8,
        MID = 16,
        LARGE = 32
    };
}


#endif // ROUTE_OPT_TWO_STAGE_MACRO_HPP
