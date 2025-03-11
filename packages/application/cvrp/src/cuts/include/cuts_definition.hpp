/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_CUTS_DEFINITION_HPP
#define ROUTE_OPT_CUTS_DEFINITION_HPP
#include "rank1_macro.hpp"
#include "rcc_macro.hpp"

namespace RouteOpt::Application::CVRP {
    using R1c = Rank1Cuts::R1c;
    using Rcc = RCCs::Rcc;

    constexpr int INVALID_BRC_INDEX = -1;
    constexpr int INVALID_ROW_INDEX = -1;

    struct Brc {
        std::pair<int, int> edge;
        int idx_brc{INVALID_BRC_INDEX};
        bool br_dir;
    };
}

#endif // ROUTE_OPT_CUTS_DEFINITION_HPP
