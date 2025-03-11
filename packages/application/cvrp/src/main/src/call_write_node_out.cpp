/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "two_stage_controller.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRPSolver::callWriteNodeOut(BbNode *node,
                                      const Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                                      const Branching::BKF::BKFDataShared &bkf_data_shared) const {
        if (node->getIfInEnumState() && node->getColPoolSize() < static_cast<int>(COL_POOL_TYPE::TINY)) return;
        auto node_name = tree_path.empty() ? ins_name : tree_path;
        TwoStageController::writeNodeOut<!IF_SYMMETRY_PROHIBIT>(node_name,
                                                                node->refSolver(),
                                                                node->getCols(),
                                                                node->getRCCs(),
                                                                node->getR1Cs(),
                                                                node->getBrCs(),
                                                                node->getValue(),
                                                                node->getIdx(),
                                                                node->getIfInEnumState(),
                                                                dim,
                                                                pricing_controller.getNumBucketPerVertex(),
                                                                pricing_controller.getStepSize(),
                                                                pricing_controller.getExistLabelsInForward(),
                                                                node->getAllForwardBuckets(),
                                                                pricing_controller.getExistLabelsInBackward(),
                                                                node->getAllBackwardBuckets(),
                                                                node->getDeletedColumnsInEnumerationColumnPool(),
                                                                node->getIndexColPool(),
                                                                node->getCostColPool(),
                                                                pricing_controller.getColumnPoolPtr(),
                                                                pricing_controller.getMaxEnumerationSuccessGap(),
                                                                pricing_controller.getSuccessEnumerationGap(),
                                                                pricing_controller.getMinEnumerationFailGap(),
                                                                pricing_controller.getMaxGap2TryEnumeration(),
                                                                history,
                                                                bkf_data_shared);

        node->refIfTerminate() = true;
    }
}
