/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_EXACT_HPP
#define ROUTE_OPT_EXACT_HPP
#include "pricing_macro.hpp"
#include "../../main/include/cvrp_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<bool if_symmetry>
    int CVRP_Pricing::generateColumnsByExact(double time_limit) {
        int ccnt = 0;
    RE_TRY:
        runLabeling<true, false, false, if_symmetry, PRICING_LEVEL::EXACT>(time_limit);

        if (if_short_memory) {
            reallocateLabel();
            initializeLabels<if_symmetry>();
            goto RE_TRY;
        }
        if (!if_exact_labeling_finished) {
            goto QUIT;
        }

        if (!if_symmetry) {
        RE_TRY2:
            runLabeling<false, false, false, if_symmetry, PRICING_LEVEL::EXACT>(time_limit);

            if (if_short_memory) {
                reallocateLabel();
                initializeLabels<if_symmetry>();
                goto RE_TRY2;
            }
            if (!if_exact_labeling_finished) {
                goto QUIT;
            }
        }

        ccnt = concatenateCols_prior_forward<if_symmetry>();
        updateDominanceStatics<if_symmetry>();

    QUIT:
        if_exact_cg_finished = if_exact_labeling_finished;
        return ccnt;
    }


    template<bool if_symmetry>
    void CVRP_Pricing::updateDominanceStatics() {
        NumExistedLabels = 0;
        NumExistedLabel_back = 0;

        for (int i = 1; i < dim; ++i) {
            for (int b = 0; b < num_buckets_per_vertex; ++b) {
                NumExistedLabels += static_cast<double>(label_array_in_forward_sense[i][b].size());
                if constexpr (!if_symmetry) {
                    NumExistedLabel_back += static_cast<double>(label_array_in_backward_sense[i][b].size());
                }
            }
        }

        ratio_dominance_checks_non_dominant.first +=
                if_symmetry
                    ? num_dominance_checks / NumExistedLabels
                    : num_dominance_checks
                      / (NumExistedLabels + NumExistedLabel_back);
        ++ratio_dominance_checks_non_dominant.second;
        // std::cout << "\x1b[94mratio_dominance_checks_non_dominant= "
        //    << ratio_dominance_checks_non_dominant.first / ratio_dominance_checks_non_dominant.second
        //    << "\x1b[0m" << std::endl;
        //
        // std::cout << "\x1b[94midx_glo= " << idx_glo << "\x1b[0m" << std::endl;
    }

    template<bool if_symmetry>
    void CVRP_Pricing::adjustResourceMeetPointInPricing() {
        if (!if_exact_labeling_cg || !if_exact_labeling_finished) return;
        if constexpr (!if_symmetry) {
            double dif = std::abs(NumExistedLabels - NumExistedLabel_back);
            double over = dif / std::min(NumExistedLabels, NumExistedLabel_back);
            if (over > NumberOfOverLabelsInMeetPoint) {
                if (NumExistedLabels > NumExistedLabel_back) {
                    if (last_tag == 1) {
                        meet_point_factor /= MeetPointFactorDecayFactor;
                    }
                    last_tag = -1;
                    meet_point_resource_in_bi_dir *= (1 - meet_point_factor);
                } else {
                    if (last_tag == -1) {
                        meet_point_factor /= MeetPointFactorDecayFactor;
                    }
                    last_tag = 1;
                    meet_point_resource_in_bi_dir *= (1 + meet_point_factor);
                }
            }
        }
    }
}
#endif // ROUTE_OPT_EXACT_HPP
