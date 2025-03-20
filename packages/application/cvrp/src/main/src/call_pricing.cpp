/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "cvrp_macro.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRPSolver::callPricingAtBeg(BbNode *node) {
        setEnv(node);
        if constexpr (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2) {
            GetTrainingData<BbNode, std::pair<int, int>,
                PairHasher>::checkIfStopGeneratingData(node->getTreeSize(), node->refIfTerminate());
            if (node->getIfTerminate()) return;
        }
        callPricing(node, std::numeric_limits<float>::max());
        if (node->getIfRootNode() && app_type == APPLICATION_TYPE::VRPTW)
            augmentNGRound(
                node, pricing_controller.refNG());
    }

    void CVRPSolver::callPricing(BbNode *node, double labeling_time_limit) {
        if (node->getIfInEnumState()) {
        ENU:
            callInspection(node);
            BbNode::regenerateEnumMat(node, nullptr, false, optimal_dual_vector);
        } else {
            callLabeling(node, labeling_time_limit);
            if (!node->getIfInEnumState()
                && ml_type != ML_TYPE::ML_GET_DATA_1 && ml_type != ML_TYPE::ML_GET_DATA_2
            ) {
                if (pricing_controller.getIfCompleteCG()) {
                    std::vector<int> cstr_index;
                    node->findNonActiveCuts(optimal_dual_vector, cstr_index);
                    SAFE_SOLVER(node->refSolver().reoptimize(SOLVER_DUAL_SIMPLEX));
                }
            } else {
                goto ENU;
            }
        }

        if (glob_timer.getTime() > TIME_LIMIT) {
            PRINT_REMIND("time limit reached!");
            node->refIfTerminate() = true;
        }

        if constexpr (ml_type != ML_TYPE::ML_NO_USE) {
            if (!node->getIfTerminate()) l2b_controller.recordEdgeLongInfo(BbNode::obtainSolEdgeMap(node));
        }
    }

    void CVRPSolver::callInspection(BbNode *node) {
        constexpr bool if_update_column_pool = true;
        constexpr bool if_allow_delete_col = true;
        solveLPByInspection(node, if_update_column_pool, if_allow_delete_col);
    };

    void CVRPSolver::callLabeling(BbNode *node, double labeling_time_limit) {
        constexpr bool if_open_heur = true;
        constexpr bool if_open_exact = true;
        constexpr bool if_update_node_val = true;
        constexpr bool if_possible_terminate_early = false;
        constexpr bool if_fix_row = true;
        constexpr bool if_allow_delete_col = true;

        bool if_consider_regenerate_bucket_graph = node->getIfRootNode();
        bool if_fix_meet_point = (!node->getR1Cs().empty()) &&
                                 (rank1_separation_controller.getIfOnceUseNoSymmetryMem());

        solveLPInLabeling(node, if_open_heur, if_open_exact, if_update_node_val,
                          if_consider_regenerate_bucket_graph, if_possible_terminate_early,
                          if_fix_row, if_fix_meet_point, if_allow_delete_col, labeling_time_limit);


        if (node->getIfTerminate() || !pricing_controller.getIfCompleteCG()) return;

        if (!node->getIfRootNode() && (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2)) return;
        pricing_controller.eliminateArcs<!IF_SYMMETRY_PROHIBIT>(node->getRCCs(), node->getR1Cs(),
                                                                node->getBrCs(), optimal_dual_vector, ub,
                                                                node->calculateOptimalGap(ub),
                                                                node->refLastGap(),
                                                                node->refNumForwardBucketArcs(),
                                                                node->refNumBackwardBucketArcs(),
                                                                node->refNumForwardJumpArcs(),
                                                                node->refNumBackwardJumpArcs());

        if constexpr (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2) return;
        if (!node->getRCCs().empty() || !node->getR1Cs().empty()) callEnumeration(node);
    };
}
