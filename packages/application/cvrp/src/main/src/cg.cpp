/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <iostream>
#include <iomanip>
#include "cvrp.hpp"

namespace RouteOpt::Application::CVRP {
    void printInfoLabeling(int iter, int num_col, int num_row,
                           double &mt, double &spt, double lp, double ub, bool if_force) {
        if (iter == 0 || iter % PRINT_LABELING_STEP_SIZE != 0 && !if_force) {
            return;
        }
        std::ios init(nullptr);
        init.copyfmt(std::cout);
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "it= " << std::setw(4) << std::left << iter
                << "ncol= " << std::setw(6) << std::left << num_col
                << "ncstr= " << std::setw(5) << std::left << num_row
                << "mt= " << std::setw(6) << std::left << mt
                << "spt= " << std::setw(6) << std::left << spt
                << "lpval= " << std::setw(10) << std::left << lp
                << "ub= " << std::setw(9) << std::left << ub
                << std::endl;
        std::cout.copyfmt(init);
        mt = 0;
        spt = 0;
    }

    void printTypeLabeling(PRICING_LEVEL cg_mode) {
        std::string message;
        switch (cg_mode) {
            case PRICING_LEVEL::EXACT:
                message = "EXACT";

                break;
            case PRICING_LEVEL::LIGHT:
                message = "LIGHT HEUR";

                break;
            case PRICING_LEVEL::HEAVY:
                message = "HEAVY HEUR";

                break;
            default:
                THROW_RUNTIME_ERROR("Unknown pricing level!");
        }
        printHeadLines("Pricing Level: " + message);
    }

    bool tellIfInt(const std::vector<double> &X) {
        bool is_integer = true;
        for (double i: X) {
            if (i > TOLERANCE && std::abs(i - 1) > TOLERANCE) {
                is_integer = false;
                break;
            }
        }
        return is_integer;
    }

    void updateIPOptSol(const std::vector<double> &X, const std::vector<SequenceInfo> &cols,
                        std::vector<std::vector<int> > &ip_opt_sol) {
        ip_opt_sol.clear();
        for (int i = 0; i < cols.size(); ++i) {
            if (cols[i].col_seq.empty()) continue;
            if (X[i] > 0.5) {
                ip_opt_sol.emplace_back(cols[i].col_seq);
            }
        }
    }

    void CVRPSolver::updateIntegerSolution(double val, const std::vector<double> &X,
                                           const std::vector<SequenceInfo> &cols, bool &if_integer, bool &if_feasible) {
        if (tellIfInt(X)) if_integer = true;
        else {
            if_integer = false;
            if_feasible = false;
            return;
        }
        if (val + TOLERANCE < ub) {
            std::vector<double> sol;
            std::vector<SequenceInfo> sols;
            for (int i = 0; i < cols.size(); ++i) {
                if (X[i] > 0.5) {
                    sol.emplace_back(X[i]);
                    sols.emplace_back(cols[i]);
                }
            }
            checkSolutionFeasibility(sol, sols, if_feasible);
            if (if_feasible) {
                ub = val;
                updateIPOptSol(X, cols, ip_opt_sol);
                std::cout << "\x1b[36mUpdated UB: " << ub << "\x1b[0m" << std::endl;
            }
        } else {
            if_feasible = true;
        }
    }

    void CVRPSolver::runColumnGenerationType(BbNode *node, PRICING_LEVEL cg_mode, double time_limit,
                                             bool if_possible_terminate_early,
                                             bool if_fix_row, bool if_fix_meet_point, bool if_allow_delete_col) {
        if (node->getIfTerminate()) return;

        int num_row;
        SAFE_SOLVER(node->refSolver().getNumRow(&num_row))
        std::vector<double> pi4_labeling(num_row);
        bool if_integer, if_feasible;
        int ccnt = 0;
        int iter = 0;
        double prior_value = node->getValue();
        int lp_method = SOLVER_BARRIER;

        double mt = 0, spt = 0;
        int num_col;
        std::vector<double> X;

        printTypeLabeling(cg_mode);


        constexpr bool if_symmetry = !if_symmetry_prohibit(app_type);
        pricing_controller.initializeOnceB4WholePricing();

        while (true) {
        START_OVER:
            if (node->getIfTerminate()) goto BREAK;

            mt += TimeSetter::measure([&]() {
                node->optimizeLPForOneIteration(prior_value, if_allow_delete_col, lp_method);
            });

            lp_method = SOLVER_PRIMAL_SIMPLEX;

            SAFE_SOLVER(node->refSolver().getNumCol(&num_col))
            X.resize(num_col);
            SAFE_SOLVER(node->refSolver().getX(0, num_col, X.data()))


            updateIntegerSolution(prior_value, X, node->getCols(), if_integer, if_feasible);

            if (if_integer && !if_feasible && !if_fix_row) {
                addFeasibilityCuts(num_row, X, node->getCols(), node->refRCCs(),
                                   node->refSolver());
                SAFE_SOLVER(node->refSolver().getNumRow(&num_row))
                pi4_labeling.resize(num_row);
                lp_method = SOLVER_DUAL_SIMPLEX;
                goto START_OVER;
            }


            printInfoLabeling(iter, num_col, num_row, mt, spt, prior_value, ub, false);

            if (if_possible_terminate_early) {
                if (cg_mode == PRICING_LEVEL::LIGHT && num_col > LP_COL_FINAL_LIMIT *
                    HEURISTIC_LIGHT_TESTING_MAX_COLUMN_RATIO) {
                    std::cout << "terminate early: " << num_col << " columns is too large for lighter heuristic!" <<
                            std::endl;
                    goto BREAK;
                }
                if (cg_mode == PRICING_LEVEL::HEAVY && num_col > LP_COL_FINAL_LIMIT *
                    HEURISTIC_HEAVY_TESTING_MAX_COLUMN_RATIO) {
                    std::cout << "terminate early: " << num_col << " columns is too large for heavier heuristic!" <<
                            std::endl;
                    goto BREAK;
                }
            }


            SAFE_SOLVER(node->refSolver().getDual(0, num_row, pi4_labeling.data()))

            pricing_controller.priceConstraints(node->getRCCs(),
                                                node->getR1Cs(),
                                                node->getBrCs(),
                                                pi4_labeling);

            ++iter;

            spt += TimeSetter::measure([&]() {
                switch (cg_mode) {
                    case PRICING_LEVEL::LIGHT:
                        ccnt = pricing_controller.generateColumnsByHeuristic<if_symmetry, PRICING_LEVEL::LIGHT>();
                        break;
                    case PRICING_LEVEL::HEAVY:
                        ccnt = pricing_controller.generateColumnsByHeuristic<if_symmetry, PRICING_LEVEL::HEAVY>();
                        break;
                    case PRICING_LEVEL::EXACT:
                        ccnt = pricing_controller.generateColumnsByExact<if_symmetry>(time_limit);
                        if (!if_fix_meet_point) pricing_controller.adjustResourceMeetPointInPricing<if_symmetry>();
                        pricing_controller.setTerminateMarker(prior_value, ub, node->refIfTerminate());
                        break;
                }
            });
            if (cg_mode == PRICING_LEVEL::EXACT && !pricing_controller.getIfCompleteCG())goto BREAK;

            if (!node->getIfTerminate())add_column_controller.addColumns(ccnt, pi4_labeling, true);

            if (ccnt == 0) {
                if (cg_mode == PRICING_LEVEL::EXACT) optimal_dual_vector = pi4_labeling;
                goto BREAK;
            }
            continue;

        BREAK:
            printInfoLabeling(iter, num_col, num_row, mt, spt, prior_value, ub, true);
            break;
        }
    }

    void CVRPSolver::solveLPInLabeling(BbNode *node, bool if_open_heur, bool if_open_exact,
                                       bool if_update_node_val,
                                       bool if_consider_regenerate_bucket_graph, bool if_possible_terminate_early,
                                       bool if_fix_row, bool if_fix_meet_point, bool if_allow_delete_col,
                                       double time_limit) {
        std::cout << SMALL_PHASE_SEPARATION;
        if (if_open_heur) {
            runColumnGenerationType(node, PRICING_LEVEL::LIGHT, std::numeric_limits<float>::max(),
                                    if_possible_terminate_early, if_fix_row,
                                    if_fix_meet_point, if_allow_delete_col);

            runColumnGenerationType(node, PRICING_LEVEL::HEAVY, std::numeric_limits<float>::max(),
                                    if_possible_terminate_early, if_fix_row,
                                    if_fix_meet_point, if_allow_delete_col);
        }

        if (if_open_exact) {
            runColumnGenerationType(node, PRICING_LEVEL::EXACT, time_limit, if_possible_terminate_early,
                                    if_fix_row,
                                    if_fix_meet_point, if_allow_delete_col);
            if (!pricing_controller.getIfCompleteCG()) goto QUIT;
            int num_col;
            SAFE_SOLVER(node->refSolver().getNumCol(&num_col))
            if (num_col > LP_COL_FINAL_LIMIT && if_allow_delete_col) node->cleanIndexColForNode();
        }

        if (if_update_node_val) {
            if (node->getIfTerminate()) {
                node->refValue() = ub;
            } else {
                if (pricing_controller.getIfCompleteCG()) SAFE_SOLVER(node->refSolver().getObjVal(& node->refValue()))
            }
        }

    QUIT:
        glob_timer.report();
        if (node->getIfTerminate()) return;
        if (if_consider_regenerate_bucket_graph && pricing_controller.getIfCompleteCG()) {
            pricing_controller.considerRegenerateBucketGraph<!IF_SYMMETRY_PROHIBIT>(
                node->refAllForwardBuckets(),
                node->refAllBackwardBuckets(),
                node->refNumForwardBucketArcs(),
                node->refNumForwardJumpArcs(),
                node->refNumBackwardBucketArcs(),
                node->refNumBackwardJumpArcs());
        }
    }


    int BbNode::generateColumnsByInspection(const std::vector<double> &pi, std::vector<int> &col_added,
                                            RowVectorXd &rc) {
        col_added.clear();
        auto size_pool = static_cast<int>(index_columns_in_enumeration_column_pool.size());
        if (size_pool == 0) return 0;
        auto &mat = matrix_in_enumeration;


        rc = cost_for_columns_in_enumeration_column_pool;

        int num = 0;
        for (auto &it: mat) {
            RowVectorXd dual(it.rows());
            for (int i = 0; i < it.rows(); ++i)dual(i) = pi[num++];
            rc -= dual * it;
        }

        std::vector<std::pair<int, double> > candidates;
        candidates.reserve(size_pool);

        for (int i = 0; i < size_pool; ++i) {
            if (!deleted_columns_in_enumeration_pool[i] && rc(i) < RC_TOLERANCE) {
                candidates.emplace_back(i, rc(i));
            }
        }

        if (static_cast<int>(candidates.size()) > MAX_NUM_ROUTES_Exact) {
            std::nth_element(
                candidates.begin(),
                candidates.begin() + MAX_NUM_ROUTES_Exact,
                candidates.end(),
                [](auto &lhs, auto &rhs) {
                    return lhs.second < rhs.second;
                }
            );
            candidates.resize(MAX_NUM_ROUTES_Exact);
        }

        col_added.resize(candidates.size());
        std::transform(candidates.begin(), candidates.end(), col_added.begin(), [](auto &p) { return p.first; });

        return static_cast<int>(col_added.size());
    }


    void CVRPSolver::solveLPByInspection(BbNode *node, bool if_update_column_pool, bool if_allow_delete_col) {
        if (node->getIfTerminate()) return;

        int num_row;
        SAFE_SOLVER(node->refSolver().getNumRow(&num_row))
        std::vector<double> pi4_labeling(num_row);
        //
        bool if_integer, if_feasible;
        int ccnt = 0;
        int iter = 0;
        double prior_value = node->getValue();
        int lp_method = SOLVER_BARRIER;
        //
        double mt = 0, spt = 0;
        int num_col;
        std::vector<double> X;
        RowVectorXd rc;
        auto &deleted_columns_in_enumeration_pool = node->refDeletedColumnsInEnumerationColumnPool();

        std::vector<int> col_added;

        while (true) {
        START_OVER:
            if (node->getIfTerminate()) goto BREAK;

            mt += TimeSetter::measure([&]() {
                node->optimizeLPForOneIteration(prior_value, lp_method);
            });


            if (lp_method == SOLVER_BARRIER) lp_method = SOLVER_PRIMAL_SIMPLEX;

            SAFE_SOLVER(node->refSolver().getNumCol(&num_col))
            X.resize(num_col);
            SAFE_SOLVER(node->refSolver().getX(0, num_col, X.data()))


            updateIntegerSolution(prior_value, X, node->getCols(), if_integer, if_feasible);


            printInfoLabeling(iter, num_col, num_row, mt, spt, prior_value, ub, false);

            SAFE_SOLVER(node->refSolver().getDual(0, num_row, pi4_labeling.data()))

            ++iter;


            spt += TimeSetter::measure([&]() {
                ccnt = node->generateColumnsByInspection(pi4_labeling, col_added, rc);
                if (if_update_column_pool) {
                    for (const auto &i: col_added) {
                        deleted_columns_in_enumeration_pool[i] = true;
                    }
                }
            });


            if (ccnt == 0) {
                if (if_update_column_pool) {
                    SAFE_SOLVER(node->refSolver().getObjVal(&node->refValue()))
                    auto opt_gap = node->calculateOptimalGap(ub);
                    if (opt_gap < RC_TOLERANCE) {
                        node->refValue() = ub;
                        node->refIfTerminate() = true;
                        goto BREAK;
                    }
                    auto size_pool = static_cast<int>(node->getIndexColPool().size());
                    for (int i = 0; i < size_pool; ++i) {
                        if (rc(i) > opt_gap) {
                            deleted_columns_in_enumeration_pool[i] = true;
                        }
                    }
                    optimal_dual_vector = pi4_labeling;
                    if (if_allow_delete_col)node->cleanIndexColForNode(ub, optimal_dual_vector);
                }
                goto BREAK;
            }

            if (!node->getIfTerminate())add_column_controller.addColumnsByInspection(col_added);
            continue;

        BREAK:
            printInfoLabeling(iter, num_col, num_row, mt, spt, prior_value, ub, true);
            break;
        }
        glob_timer.report();
    }
}
