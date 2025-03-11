/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "cvrp_macro.hpp"

namespace RouteOpt::Application::CVRP {
    std::pair<int, int> CVRPSolver::callMLCandidateSelection(BbNode *node,
                                                             Branching::BranchingHistory<std::pair<int, int>,
                                                                 PairHasher> &history,
                                                             Branching::BranchingDataShared<std::pair<int, int>,
                                                                 PairHasher> &data_shared,
                                                             Branching::CandidateSelector::BranchingTesting<BbNode,
                                                                 std::pair<int, int>, PairHasher> &tester) {
        if constexpr (ml_type == ML_TYPE::ML_NO_USE)
            THROW_RUNTIME_ERROR("machine learning module is not enabled, but called");
        int num_col;
        SAFE_SOLVER(node->refSolver().getNumCol(&num_col))
        std::vector<double> x(num_col);
        SAFE_SOLVER(node->refSolver().getX(0, num_col, x.data()))
        std::vector<int> route_length(num_col);
        for (int i = 0; i < num_col; ++i) {
            route_length[i] = static_cast<int>(node->getCols()[i].col_seq.size());
        }

        data_shared.refCandidateMap() = BbNode::obtainSolEdgeMap(node);

        switch (ml_type) {
            case ML_TYPE::ML_GET_DATA_1:
                l2b_train_controller.generateModelPhase1(node, history, data_shared, tester, node->getTreeSize(),
                                                         x, route_length);
                break;
            case ML_TYPE::ML_GET_DATA_2:
                l2b_train_controller.generateModelPhase2(node, l2b_predict_controller, history, data_shared, tester,
                                                         node->getTreeSize(), node->calculateOptimalGap(ub), x,
                                                         route_length);
                break;
            case ML_TYPE::ML_USE_MODEL:
                l2b_predict_controller.useMLInGeneralFramework(node, history, data_shared, tester,
                                                               node->getTreeSize(), node->calculateOptimalGap(ub), x,
                                                               route_length);
                break;
            default: ;
        }
        l2b_controller.cleanLastData();
        return data_shared.refBranchPair().front();
    }
}
