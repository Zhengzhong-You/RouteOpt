/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HEURISTIC_LABELING_HPP
#define ROUTE_OPT_HEURISTIC_LABELING_HPP
namespace RouteOpt::Application::CVRP {
    template<bool if_symmetry, PRICING_LEVEL pricing_level>
    int CVRP_Pricing::generateColumnsByHeuristic() {
        if (pricing_level == PRICING_LEVEL::HEAVY) {
            num_col_generated_ub = MaxNumRoutesInHeavierHeur;
        } else if (pricing_level == PRICING_LEVEL::LIGHT) {
            num_col_generated_ub = MaxNumRoutesInLighterHeur;
        } else {
            THROW_RUNTIME_ERROR("Pricing level does not match");
        }

        runLabeling<true, false, false, if_symmetry, pricing_level>(std::numeric_limits<float>::max());
        if (if_short_memory) {
            THROW_RUNTIME_ERROR("Lack of memory even for heuristic labeling");
        }

        if (!if_symmetry) {
            runLabeling<false, false, false, if_symmetry, pricing_level>(std::numeric_limits<float>::max());
            if (if_short_memory) {
                THROW_RUNTIME_ERROR("Lack of memory even for heuristic labeling");
            }
        }

        int ccnt = concatenateCols_prior_forward<if_symmetry>();

        return ccnt;
    }
}

#endif // ROUTE_OPT_HEURISTIC_LABELING_HPP
