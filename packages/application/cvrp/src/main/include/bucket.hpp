/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_BUCKET_HPP
#define ROUTE_OPT_BUCKET_HPP
#include <vector>
#include "label.hpp"


namespace RouteOpt::Application::CVRP {
    /**
     * Bucket structure for bucket graph representation
     * Contains regular arcs and jump arcs for a vertex in a resource bucket
     */
    struct Bucket {
        std::vector<int> bucket_arcs{};
        std::vector<std::pair<res_int, int> > jump_arcs{};//res and j
        int i{};
    };
}

#endif // ROUTE_OPT_BUCKET_HPP
