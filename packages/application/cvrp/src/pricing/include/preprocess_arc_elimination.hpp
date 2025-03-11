/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_PREPROCESS_ARC_ELIMINATION_HPP
#define ROUTE_OPT_PREPROCESS_ARC_ELIMINATION_HPP
#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    template<bool dir>
    void CVRP_Pricing::populateTellWhichBin4ArcElimination() {
        /**
         * one time calculation except for regenerate the bucket graph
         */
        if constexpr (dir) {
            if (!tell_which_bin4_arc_elimination_in_forward_sense.empty()) return;
        } else {
            if (!tell_which_bin4_arc_elimination_in_backward_sense.empty()) return;
        }
        size_t size = dim * dim * num_buckets_per_vertex;
        int dim_sq = dim * dim;
        if constexpr (dir) tell_which_bin4_arc_elimination_in_forward_sense.reserve(size);
        else tell_which_bin4_arc_elimination_in_backward_sense.reserve(size);
        std::unordered_map<int, int> map_bj_B;
        std::vector<std::pair<int, int> > vec_bj_B(num_buckets_per_vertex);
        map_bj_B.reserve(num_buckets_per_vertex + 1);
        Resource tmp_Resource, res_tuple;
        for (int b = 0; b <= num_buckets_per_vertex; ++b) map_bj_B[b] = dir ? -1 : num_buckets_per_vertex;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                int base = i * dim + j;
                for (auto &key_value: map_bj_B) key_value.second = dir ? -1 : num_buckets_per_vertex;
                for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0;
                     dir ? ++b : --b) {
                    res_tuple = dir ? Resource{} : resource;
                    auto res = dir ? b * step_size : std::min((b + 1) * step_size, resource.resources[0]);
                    res_tuple.resources[0] = res;
                    if (dir
                            ? !increaseMainResourceConsumption(res_tuple, tmp_Resource, i, j)
                            : !decreaseMainResourceConsumption(res_tuple, tmp_Resource, i, j))
                        break;
                    auto bj = static_cast<int>(tmp_Resource.resources[0] / step_size);
                    if (dir ? map_bj_B[bj] < b : map_bj_B[bj] > b) map_bj_B[bj] = b;
                }
                for (int bj = dir ? num_buckets_per_vertex - 1 : 0; dir ? bj >= 0 : bj < num_buckets_per_vertex;
                     dir ? --bj : ++bj) {
                    if (map_bj_B[bj] != (dir ? -1 : num_buckets_per_vertex)) {
                        for (int bj2 = dir ? bj + 1 : bj - 1; dir ? bj2 < num_buckets_per_vertex : bj2 >= 0;
                             dir ? ++bj2 : --bj2)
                            map_bj_B[bj2] = map_bj_B[bj];
                        break;
                    }
                }
                for (int bj = 0; bj < num_buckets_per_vertex; ++bj) {
                    (dir
                         ? tell_which_bin4_arc_elimination_in_forward_sense[base + bj * dim_sq]
                         : tell_which_bin4_arc_elimination_in_backward_sense[base + bj * dim_sq]) = map_bj_B[bj];
                }
            }
        }
    }
}

#endif // ROUTE_OPT_PREPROCESS_ARC_ELIMINATION_HPP
