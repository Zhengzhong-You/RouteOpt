/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRP_Pricing::updatePtr(Bucket **all_forward_buckets, Bucket **all_backward_buckets,
                                 std::vector<std::vector<std::vector<int> > > *topological_order_forward_ptr,
                                 std::vector<std::vector<std::vector<int> > > *topological_order_backward_ptr) {
        this->all_forward_buckets = all_forward_buckets;
        this->all_backward_buckets = all_backward_buckets;
        this->topological_order_forward_ptr = topological_order_forward_ptr;
        this->topological_order_backward_ptr = topological_order_backward_ptr;
    }


    void CVRP_Pricing::initializeOnceB4WholePricing() {
        ratio_dominance_checks_non_dominant = {};
    }


    void CVRP_Pricing::getNG() {
        ng_mem4_vertex.resize(dim, 0);
        initial_ng_size = std::min(InitialNGSize, dim - 1);
        std::vector<std::pair<int, double> > cost(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 1; j < dim; ++j) {
                cost[j].first = j;
                cost[j].second = cost_mat4_vertex_ref.get()[i][j];
            }
            cost[0].first = 0;
            cost[0].second = std::numeric_limits<float>::max();
            std::stable_sort(cost.begin(), cost.end(),
                             [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
                                 return a.second < b.second;
                             });
            routeOptLong &vst = ng_mem4_vertex[i];
            for (int k = 0; k < initial_ng_size; ++k) {
                vst.set(cost[k].first);
            }
        }
    }

    void CVRP_Pricing::assignMemory() {
        double aver = getAverageRouteLength();
        mem4_pricing = static_cast<size_t>(aver * MAX_ROUTE_PRICING);
        col_pool4_pricing = new int[mem4_pricing];

        label_assign = LABEL_ASSIGN / 2;
        reallocateLabel();
    }
}
