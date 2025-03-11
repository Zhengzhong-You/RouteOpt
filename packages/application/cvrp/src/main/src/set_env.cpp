/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "two_stage_controller.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRPSolver::setEnv(BbNode *node) {
        std::cout << "idx= " << node->getIdx() << " value= " << node->getValue() << " enu= " << node->
                getIfInEnumState() << std::endl;
        //print all brcs
        for (auto &brc: node->getBrCs()) {
            std::cout << "brc= " << brc.edge.first << "-" << brc.edge.second << " : " << (brc.br_dir ? "+" : "-")
                    << std::endl;
        }

        if constexpr (IF_WRITE_NODE_OUT) TwoStageController::updateUB(ins_name, ub);

        pricing_controller.updatePtr(node->getAllForwardBuckets(), node->getAllBackwardBuckets(),
                                     node->ptrTopologicalOrderForward(), node->ptrTopologicalOrderBackward());
        pricing_controller.initializeBucketGraphForNode<!IF_SYMMETRY_PROHIBIT>(
            node->refAllForwardBuckets(), node->refAllBackwardBuckets(),
            node->refNumForwardBucketArcs(), node->refNumBackwardBucketArcs());
        add_column_controller.updatePtr(node->refCols(), &node->refSolver(), node->getRCCs(), node->getR1Cs(),
                                        node->getBrCs(), node->getIndexColPool(), node->getCostColPool(),
                                        node->getMatrixColPool(),
                                        pricing_controller.getColumnPoolPtr());
    }
}
