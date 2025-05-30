/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRPSolver::callEnumeration(BbNode *node) {
        if constexpr (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2) {
            return;
        }
        auto if_suc = pricing_controller.enumerateMIP<!IF_SYMMETRY_PROHIBIT>(node->getRCCs(), node->getR1Cs(),
                                                                             node->getBrCs(), optimal_dual_vector,
                                                                             ub,
                                                                             node->calculateOptimalGap(ub),
                                                                             node->refNumForwardBucketArcs(),
                                                                             node->refNumBackwardBucketArcs(),
                                                                             node->refIfInEnumState(),
                                                                             node->
                                                                             refIndexColumnsInEnumerationColumnPool(),
                                                                             node->
                                                                             refCostForColumnsInEnumerationColumnPool(),
                                                                             node->refValidSize(),
                                                                             app_type == APPLICATION_TYPE::VRPTW
                                                                                 ? demand
                                                                                 : std::vector<double>{},
                                                                             app_type == APPLICATION_TYPE::VRPTW
                                                                                 ? cap
                                                                                 : 0);

        if (!if_suc) return;

        node->preprocessEnumeration(pricing_controller.getColumnPoolPtr(),
                                    rank1_coefficient_getter, pricing_controller.getNGMem(),
                                    app_type == APPLICATION_TYPE::VRPTW
                                        ? demand
                                        : std::vector<double>{},
                                    app_type == APPLICATION_TYPE::VRPTW ? cap : 0);

        glob_timer.report();
    }
}
