/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "two_stage_controller.hpp"

namespace RouteOpt::Application::CVRP {
    namespace InNodeNameSpace {
        //build LP;
        void buildLP(int dim,
                     int min_num_vehicle,
                     bool if_enu,
                     const std::vector<double> &cost,
                     Solver &solver,
                     CVRP_AddColumnController &add_column_controller,
                     CVRP_Pricing &pricing,
                     const std::vector<SequenceInfo> &cols,
                     const std::vector<Rcc> &rccs,
                     const std::vector<R1c> &r1cs,
                     const std::vector<Brc> &brcs) {
            //set up the LP by adding constraints;
            int num_row = static_cast<int>(rccs.size() + r1cs.size());
            for (auto &brc: brcs) {
                if (brc.idx_brc != INVALID_BRC_INDEX) num_row++;
            }
            num_row += dim;
            std::vector<char> sense(num_row, SOLVER_EQUAL);
            std::vector<double> rhs(num_row, 1.0);
            //set-partitioning constraints: pass
            //vehicle constraints
            sense[dim - 1] = SOLVER_GREATER_EQUAL;
            rhs[dim - 1] = min_num_vehicle;
            //rccs
            for (auto &rcc: rccs) {
                auto idx = rcc.idx_rcc;
                sense[idx] = rcc.form_rcc == static_cast<int>(RCCs::RCCForm::RCC_FORM_3)
                                 ? SOLVER_GREATER_EQUAL
                                 : SOLVER_LESS_EQUAL;
                rhs[idx] = rcc.rhs;
            }
            //r1cs
            for (auto &r1c: r1cs) {
                auto idx = r1c.idx_r1c;
                sense[idx] = SOLVER_LESS_EQUAL;
                rhs[idx] = r1c.rhs;
            }
            //brcs
            for (auto &brc: brcs) {
                auto idx = brc.idx_brc;
                if (idx == INVALID_BRC_INDEX) continue;
                sense[idx] = SOLVER_EQUAL;
                rhs[idx] = brc.br_dir ? 1.0 : 0.0;
            }
            SAFE_SOLVER(solver.newModel(MODEL_NAME, 0, nullptr, nullptr, nullptr, nullptr, nullptr))
            SAFE_SOLVER(solver.addConstraints(
                    num_row,
                    0,
                    nullptr,
                    nullptr,
                    nullptr,
                    sense.data(),
                    rhs.data(),
                    nullptr)
            )
            // new_cols_ref
            auto ccnt = static_cast<int>(cols.size());
            pricing.refNewCols() = cols; //need write this!
            std::vector<SequenceInfo> existing_cols;
            add_column_controller.updatePtr(
                existing_cols,
                &solver,
                rccs,
                r1cs,
                brcs,
                {},
                {},
                {},
                nullptr
            );
            add_column_controller.addColumns(ccnt, {}, false, if_enu);
            SAFE_SOLVER(solver.updateModel())

            //change the coefficients of the first column to rhs
            std::vector<int> cind(num_row);
            std::iota(cind.begin(), cind.end(), 0);
            std::vector<int> vind(num_row, 0);

            double v = cost[0];

            SAFE_SOLVER(solver.changeObj(0, 1, &v))
            SAFE_SOLVER(solver.changeCoeffs(num_row, cind.data(), vind.data(), rhs.data()))
            SAFE_SOLVER(solver.optimize())
        }

        //build Enu Matrix (if applicable);
        void buildEnuMatrix(
            BbNode *node,
            const int *col_pool,
            const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter) {
            node->refValidSize() = static_cast<int>(node->getIndexColPool().size());
            node->generateVertex2IndexColsAndEdge2IndexCols(col_pool, rank1_coefficient_getter, false);
        }
    }

    void CVRPSolver::callReadNodeIn(BbNode *node,
                                    Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                                    Branching::BKF::BKFDataShared &bkf_data_shared) {
        if (tree_path.empty()) {
            pricing_controller.initLabelingMemory();
            BbNode::buildModel(num_vehicle, dim, &solver, node);
            return;
        } {
            BbNode::setDim(dim);
            node->refSolver().getEnv(&solver);
        }
        std::vector<double> cost;
        auto node_name = std::string(NODE_FOLDER) + "/" + tree_path;
        TwoStageController::readNodeIn<!IF_SYMMETRY_PROHIBIT>(node_name,
                                                              cost,
                                                              node->refCols(),
                                                              node->refRCCs(),
                                                              node->refR1Cs(),
                                                              node->refBrCs(),
                                                              node->refValue(),
                                                              node->refIdx(),
                                                              node->refIfInEnumState(),
                                                              dim,
                                                              pricing_controller.refNumBucketPerVertex(),
                                                              pricing_controller.refStepSize(),
                                                              pricing_controller.refExistLabelsInForward(),
                                                              node->refAllForwardBuckets(),
                                                              pricing_controller.refExistLabelsInBackward(),
                                                              node->refAllBackwardBuckets(),
                                                              node->refDeletedColumnsInEnumerationColumnPool(),
                                                              node->refIndexColPool(),
                                                              node->refCostColPool(),
                                                              pricing_controller.refColumnPoolPtr(),
                                                              pricing_controller.refMem4Pricing(),
                                                              pricing_controller.refPoolBeg4Pricing(),
                                                              pricing_controller.refMaxEnumerationSuccessGap(),
                                                              pricing_controller.refSuccessEnumerationGap(),
                                                              pricing_controller.refMinEnumerationFailGap(),
                                                              pricing_controller.refMaxGap2TryEnumeration(),
                                                              history,
                                                              bkf_data_shared);
        InNodeNameSpace::buildLP(dim, num_vehicle,
                                 node->getIfInEnumState(),
                                 cost,
                                 node->refSolver(),
                                 add_column_controller,
                                 pricing_controller,
                                 node->getCols(),
                                 node->getRCCs(),
                                 node->getR1Cs(),
                                 node->getBrCs());

        if (node->getIfInEnumState()) {
            InNodeNameSpace::buildEnuMatrix(
                node,
                pricing_controller.getColumnPoolPtr(),
                rank1_coefficient_getter);
        } else {
            pricing_controller.initLabelingMemory();
        }
    }
}
