/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include "cvrp.hpp"
#include "vrptw.hpp"
#include "rank1_data_shared.hpp"

namespace RouteOpt::Application::CVRP {
    namespace CuttingDetail {
        bool if_pure_rcc_tail = false;
        bool if_node_memory = true;
        constexpr bool if_print_tail_reason = true;

        inline void trySetInitialIfNodeMemory(const std::vector<R1c> &r1cs, int dim) {
            if (!if_node_memory) return;
            for (auto &r1c: r1cs) {
                for (auto &m: r1c.arc_mem) {
                    if (m.first.empty()) continue;
                    if (m.first.size() != dim) {
                        if_node_memory = false;
                        return;
                    }
                }
            }
        }

        inline void rollbackEasyWay(BbNode *node, const Solver &rollback_solver,
                                    const std::vector<SequenceInfo> &rollback_cols,
                                    const std::vector<Rcc> &rollback_rccs,
                                    const std::vector<R1c> &rollback_r1cs,
                                    const std::vector<Brc> &rollback_brcs) {
            std::cout << BIG_PHASE_SEPARATION;
            if (node->getIfInEnumState())
                THROW_RUNTIME_ERROR("rollback in enum state!");
            PRINT_REMIND("roll back to the previous cutting iteration!");

            if (rollback_solver.model) {
                node->refSolver().freeModel();
                node->refSolver().model = rollback_solver.copyModel();
                SAFE_SOLVER(node->refSolver().reoptimize(SOLVER_BARRIER))
                node->refCols() = rollback_cols;
                node->refRCCs() = rollback_rccs;
                node->refR1Cs() = rollback_r1cs;
                node->refBrCs() = rollback_brcs;
            } else {
                THROW_RUNTIME_ERROR("no rollback model can be used!");
            }
            std::cout << BIG_PHASE_SEPARATION;
        }

        inline bool increaseCondition1(double gap_improved, int &tail_off_counter) {
            if (gap_improved < NOMINAL_IMPROVE_THRESHOLD) {
                tail_off_counter = CutsTolerance;
                if constexpr (if_print_tail_reason)
                    PRINT_REMIND("1: gap improved= " + std::to_string(gap_improved) + " threshold= " +
                    std::to_string(NOMINAL_IMPROVE_THRESHOLD));
                return true;
            }
            return false;
        }

        inline bool increaseCondition2(double ub, double lb, double gap_improved, int &tail_off_counter) {
            double gap = 1 - lb / ub;

            if (gap > GAP_TOO_CLOSE_TO_ONE)
                return false;

            if ((gap_improved / (ub - lb)) < RELATIVE_IMPROVE_THRESHOLD) {
                ++tail_off_counter;
                if constexpr (if_print_tail_reason)
                    PRINT_REMIND("2: gap improved= " + std::to_string(gap_improved) + " ratio= " + std::to_string(
                    gap_improved / (ub - lb)) + " threshold= " + std::to_string(RELATIVE_IMPROVE_THRESHOLD));

                return true;
            }
            return false;
        }

        inline bool increaseCondition3(double lb, double gap_improved, int &tail_off_counter) {
            if ((gap_improved / lb) < LB_RELATIVE_IMPROVE_THRESHOLD) {
                ++tail_off_counter;
                if constexpr (if_print_tail_reason)
                    PRINT_REMIND("3: gap improved= " + std::to_string(gap_improved) + " ratio= " + std::to_string(
                    gap_improved / lb) + " threshold= " + std::to_string(LB_RELATIVE_IMPROVE_THRESHOLD));
                return true;
            }
            return false;
        }

        // Try to increase the tail-off counter if any condition is satisfied.
        inline bool tryIncreaseTailOffCounter(double ub, double lb, double gap_improved, int &tail_off_counter) {
            if (increaseCondition1(gap_improved, tail_off_counter) ||
                increaseCondition2(ub, lb, gap_improved, tail_off_counter) ||
                increaseCondition3(lb, gap_improved, tail_off_counter)) {
                return true;
            };
            return false;
        }

        inline void setIfTailOff(double eps, double ub, double now_val, double past_val, double br_value_improved,
                                 bool if_only_check_counter,
                                 double &eps_max,
                                 bool &if_tail_off,
                                 int &tail_off_counter) {
            double gap_improved = now_val - past_val;
            if (gap_improved < -TOLERANCE * now_val) {
                THROW_RUNTIME_ERROR("lp value decreased!");
            } else if (gap_improved < TOLERANCE * now_val) {
                if_tail_off = true;
                goto QUIT;
            }


            if (!if_only_check_counter && !equalFloat(br_value_improved, 0.)) {
                if (equalFloat(eps_max, 0.)) eps_max = std::numeric_limits<float>::max();
                if (eps > CUTTING_BRANCHING_RATIO * gap_improved / br_value_improved * eps_max) {
                    if_tail_off = true;
                    if constexpr (if_print_tail_reason)
                        PRINT_REMIND("eps= " + std::to_string(eps) + " threshold= " + std::to_string(
                        CUTTING_BRANCHING_RATIO * gap_improved / br_value_improved * eps_max));

                    goto QUIT;
                }
                if (gap_improved < CUTTING_BRANCHING_RATIO * br_value_improved) {
                    if_tail_off = true;
                    if constexpr (if_print_tail_reason)
                        PRINT_REMIND("gap improved= " + std::to_string(gap_improved) + " threshold= " + std::to_string(
                        CUTTING_BRANCHING_RATIO* br_value_improved));
                    goto QUIT;
                }
                if (eps > TailOffHardTime) {
                    if_tail_off = true;
                    if constexpr (if_print_tail_reason)
                        PRINT_REMIND("eps= " + std::to_string(eps) + " tail off hard time= " + std::to_string(
                        TailOffHardTime));

                    goto QUIT;
                }
            }
            if (tryIncreaseTailOffCounter(ub, now_val, gap_improved, tail_off_counter)) {
                if (tail_off_counter >= CutsTolerance) {
                    if_tail_off = true;
                    if constexpr (if_print_tail_reason)
                        PRINT_REMIND("tail off counter= " + std::to_string(tail_off_counter) + " threshold= " +
                        std::to_string(CutsTolerance));
                    goto QUIT;
                }
            }
        QUIT:
            if (eps_max == std::numeric_limits<float>::max()) eps_max = 0.;
            eps_max = std::max(eps_max, eps);
            if (if_tail_off) {
                if (!if_pure_rcc_tail) {
                    if_pure_rcc_tail = true;
                    if_tail_off = false;
                }
                tail_off_counter = 0;
            }
        }

        inline void configurePricingOptions(bool if_in_enu_state,
                                            bool if_root_node,
                                            bool if_non_robust_cuts,
                                            double &time_limit) {
            if (if_in_enu_state) {
                time_limit = std::numeric_limits<float>::max();
            } else {
                if (if_non_robust_cuts) {
                    if (if_root_node) {
                        if (if_node_memory) {
                            time_limit = CutGenTimeThresholdInPricingInitial;
                        } else {
                            time_limit = HardTimeThresholdInPricing;
                        }
                    } else {
                        time_limit = HardTimeThresholdInPricing;
                    }
                } else {
                    time_limit = std::numeric_limits<float>::max();
                }
            }
        }

        inline void configureCutsOptions(bool if_in_enu_state,
                                         bool if_root_node,
                                         bool &if_keep_rcc,
                                         bool &if_elementary,
                                         Rank1Cuts::Separation::MemoryType &limited_memory_type,
                                         Rank1Cuts::PRICING_HARD_LEVEL &pricing_hard_level,
                                         bool &if_search_mem,
                                         bool &if_select_cuts,
                                         bool &if_collect_sol) {
            if_keep_rcc = false;
            if (if_in_enu_state) {
                if_elementary = true;
                limited_memory_type = Rank1Cuts::Separation::MemoryType::NO_MEMORY;
                if_search_mem = false;
                if_select_cuts = false;
                if_collect_sol = false;
                pricing_hard_level = Rank1Cuts::PRICING_HARD_LEVEL::EASY;
            } else {
                if_elementary = false;
                limited_memory_type = if_node_memory
                                          ? Rank1Cuts::Separation::MemoryType::NODE_MEMORY
                                          : Rank1Cuts::Separation::MemoryType::ARC_MEMORY;
                if_search_mem = true;
                if_select_cuts = true;
                pricing_hard_level = if_node_memory
                                         ? Rank1Cuts::PRICING_HARD_LEVEL::EASY
                                         : Rank1Cuts::PRICING_HARD_LEVEL::EXTREMELY_HARD;
                if (if_root_node && limited_memory_type == Rank1Cuts::Separation::MemoryType::NODE_MEMORY) {
                    if_collect_sol = true;
                } else {
                    if_collect_sol = false;
                }
            }
        }


        inline void getSols(const std::vector<double> &x, const std::vector<SequenceInfo> &cols,
                            std::vector<double> &sol_x, std::vector<SequenceInfo> &sols) {
            if (x.size() != cols.size())
                THROW_RUNTIME_ERROR("col size is not equal to sol size");
            sol_x.resize(x.size());
            sols.resize(x.size());
            int sol_cnt = 0;
            for (int i = 0; i < x.size(); ++i) {
                if (x[i] > SOL_X_TOLERANCE) {
                    sol_x[sol_cnt] = x[i];
                    sols[sol_cnt] = cols[i];
                    ++sol_cnt;
                }
            }
            sol_x.resize(sol_cnt);
            sols.resize(sol_cnt);
        }

        inline void addCstr(const sparseRowMatrixXd &mat, std::vector<char> &sense,
                            std::vector<double> &rhs, Solver &solver) {
            size_t nz = mat.nonZeros();
            std::vector<size_t> solver_beg(mat.rows() + 1);
            std::vector<int> solver_ind(nz + sense.size());
            std::vector<double> solver_val(nz + sense.size());

            nz = 0;
            solver_beg[0] = 0;
            for (int k = 0; k < mat.outerSize(); ++k) {
                sparseRowMatrixXd::InnerIterator it(mat, k);
                if (it) {
                    if (equalFloat(rhs[k], 0)) {
                        if (it.col() != 0) {
                            solver_ind[nz] = static_cast<int>(it.col());
                            solver_val[nz] = it.value();
                            ++nz;
                        }
                    } else {
                        solver_ind[nz] = 0;
                        solver_val[nz] = rhs[k];
                        ++nz;
                        if (it.col() != 0) {
                            solver_ind[nz] = static_cast<int>(it.col());
                            solver_val[nz] = it.value();
                            ++nz;
                        }
                    }
                }

                ++it;

                for (; it; ++it) {
                    solver_ind[nz] = static_cast<int>(it.col());
                    solver_val[nz] = it.value();
                    ++nz;
                }
                solver_beg[k + 1] = nz;
            }

            solver_ind.resize(nz);
            solver_val.resize(nz);


            SAFE_SOLVER(solver.XaddConstraints(
                mat.rows(),
                nz,
                solver_beg.data(),
                solver_ind.data(),
                solver_val.data(),
                sense.data(),
                rhs.data(),
                nullptr
            ))
            SAFE_SOLVER(solver.updateModel())
        }

        inline void callRCC(
            int dim,
            double cap,
            const std::vector<double> &demand,
            bool if_keep_rcc,
            bool if_elementary,
            int &num_row,
            const std::vector<double> &sol_x,
            const std::vector<SequenceInfo> &sols,
            const std::vector<SequenceInfo> &cols,
            std::vector<Rcc> &existing_rccs,
            Solver &solver
        ) {
            std::vector<Rcc> new_rccs;
            RCCs::Separation::RCCSeparationController::generateRCCs(dim, cap, demand,
                                                                    if_keep_rcc,
                                                                    if_elementary,
                                                                    sol_x, sols,
                                                                    existing_rccs,
                                                                    new_rccs);

            if (new_rccs.empty()) return;
            printHeadLines("Separate RCCs");

            sparseRowMatrixXd mat;
            RCCs::CoefficientGetter::RCCCoefficientController::getCoefficientRCC(
                cols, new_rccs, if_elementary, mat);


            std::vector<char> sense(mat.rows(), SOLVER_LESS_EQUAL);
            if (if_elementary) std::fill(sense.begin(), sense.end(), SOLVER_GREATER_EQUAL);
            std::vector<double> rhs(mat.rows());

            for (int i = 0; i < new_rccs.size(); ++i) {
                auto &rcc = new_rccs[i];
                rcc.idx_rcc = num_row;
                ++num_row;
                if (if_keep_rcc) rcc.if_keep = true;
                rhs[i] = rcc.rhs;
            }

            existing_rccs.resize(existing_rccs.size() + new_rccs.size());
            std::transform(new_rccs.begin(), new_rccs.end(), existing_rccs.end() - static_cast<int>(new_rccs.size()),
                           [](const Rcc &rcc) { return rcc; });

            addCstr(mat, sense, rhs, solver);
        }


        inline void callRank1(Rank1Cuts::Separation::MemoryType limited_memory_type,
                              Rank1Cuts::PRICING_HARD_LEVEL pricing_hard_level,
                              bool if_collect_sol,
                              bool if_mem,
                              bool if_select_cuts,
                              const std::vector<double> &x,
                              const std::vector<SequenceInfo> &sols,
                              const std::vector<R1c> &old_cuts,
                              const std::vector<SequenceInfo> &cols,
                              Solver &solver,
                              Rank1Cuts::Separation::Rank1SeparationController &rank1_separation_controller,
                              Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter,
                              int &num_row,
                              std::vector<R1c> &existing_r1cs
        ) {
            std::vector<Rank1Cuts::Separation::RouteInfo> sol(sols.size());
            for (int i = 0; i < sols.size(); ++i) {
                auto &s = sol[i];
                const auto &sol_ = sols[i];
                s.col_seq = sol_.col_seq;
                s.forward_concatenate_pos = sol_.forward_concatenate_pos;
                s.frac_x = x[i];
            }
            rank1_separation_controller.updateInfo(limited_memory_type, pricing_hard_level, if_collect_sol, sol,
                                                   old_cuts);

            std::vector<R1c> cuts;
            rank1_separation_controller.separateRank1Cuts(cuts, if_mem, if_select_cuts);
            if (cuts.empty()) return;
            printHeadLines("Separate Rank1Cuts");

            sparseRowMatrixXd mat;
            rank1_coefficient_getter.getR1CCoeffs(cols, cuts, &solver, if_mem, mat);

            std::vector<char> sense(mat.rows(), SOLVER_LESS_EQUAL);
            std::vector<double> rhs(mat.rows());

            for (int i = 0; i < cuts.size(); ++i) {
                auto &cut = cuts[i];
                cut.idx_r1c = num_row;
                ++num_row;
                rhs[i] = cut.rhs;
            }

            existing_r1cs.resize(existing_r1cs.size() + cuts.size());
            std::transform(cuts.begin(), cuts.end(), existing_r1cs.end() - static_cast<int>(cuts.size()),
                           [](const R1c &r1c) { return r1c; });

            addCstr(mat, sense, rhs, solver);
        }
    }

    void CVRPSolver::callCutting(BbNode *node) {
        if (!node->getIfRootNode() && (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type ==
                                       ML_TYPE::ML_GET_DATA_2)) {
            return;
        }

        CuttingDetail::trySetInitialIfNodeMemory(node->getR1Cs(), dim);
        std::vector<double> x;
        std::vector<double> sol_x;
        std::vector<SequenceInfo> sols;
        int num_row;
        int old_row;
        bool if_keep_rcc, if_elementary, if_search_mem, if_select_cuts, if_collect_sol,
                if_tail_off = false;
        bool old_enu_state;
        double time_limit;
        Rank1Cuts::Separation::MemoryType limited_memory_type;
        Rank1Cuts::PRICING_HARD_LEVEL pricing_hard_level;
        auto &cols = node->getCols();
        double old_val;
        double eps_max = 0., eps = 0.;
        int tail_off_counter = 0;
        std::vector<double> gap_improved_vec;
        Solver rollback_solver{};
        rollback_solver.getEnv(&solver);
        std::vector<SequenceInfo> rollback_cols;
        std::vector<R1c> rollback_r1cs;
        std::vector<Rcc> rollback_rccs;
        std::vector<Brc> rollback_brcs;
        CuttingDetail::if_pure_rcc_tail = !node->getIfRootNode();
    CUTTING:

        if (node->getIfTerminate()) goto QUIT;

        old_val = node->getValue();

        x.resize(cols.size());
        SAFE_SOLVER(node->refSolver().getX(0, cols.size(), x.data()))
        CuttingDetail::getSols(x, cols, sol_x, sols);


        SAFE_SOLVER(node->refSolver().getNumRow(&num_row))
        old_row = num_row;

        CuttingDetail::configureCutsOptions(
            node->getIfInEnumState(),
            node->getIfRootNode(),
            if_keep_rcc,
            if_elementary,
            limited_memory_type,
            pricing_hard_level,
            if_search_mem,
            if_select_cuts,
            if_collect_sol
        );

        if (limited_memory_type == Rank1Cuts::Separation::MemoryType::ARC_MEMORY) {
            if (rollback_solver.model) rollback_solver.freeModel();
            rollback_solver.model = node->refSolver().copyModel();
            rollback_cols = cols;
            rollback_rccs = node->getRCCs();
            rollback_r1cs = node->getR1Cs();
            rollback_brcs = node->getBrCs();
        }

        CuttingDetail::callRCC(
            dim,
            cap,
            demand,
            if_keep_rcc,
            if_elementary,
            num_row,
            sol_x,
            sols,
            cols,
            node->refRCCs(),
            node->refSolver()
        );


        if (!CuttingDetail::if_pure_rcc_tail) goto PRICING;


    RANK1: {
            CuttingDetail::callRank1(
                limited_memory_type,
                pricing_hard_level,
                if_collect_sol,
                if_search_mem,
                if_select_cuts,
                sol_x,
                sols,
                node->getR1Cs(),
                cols,
                node->refSolver(),
                rank1_separation_controller,
                rank1_coefficient_getter,
                num_row,
                node->refR1Cs()
            );
        }
    PRICING:
        if (old_row == num_row) {
            goto SET_TAIL_OFF;
        } {
            bool if_care_lb_improvement = !node->getIfInEnumState();
            auto if_continue = node->validateCuts(old_row, if_care_lb_improvement);
            if (!if_continue) {
                goto QUIT;
            }
        }

        node->mergeR1Cs();

        if (node->getIfInEnumState())node->changeEnumMatByCuts(rank1_coefficient_getter);

        CuttingDetail::configurePricingOptions(
            node->getIfInEnumState(),
            node->getIfRootNode(),
            !node->getR1Cs().empty(),
            time_limit
        );

    ENTER:
        old_enu_state = node->getIfInEnumState();
        eps = TimeSetter::measure([&]() {
            callPricing(node, time_limit);
        });
        if (node->getIfTerminate()) goto QUIT;
        if (!node->getIfInEnumState()) {
            if (!pricing_controller.getIfCompleteCG()) {
                if (CuttingDetail::if_node_memory) {
                    rank1_separation_controller.convert2ArcMemory(node->refR1Cs());
                    CuttingDetail::if_node_memory = false;
                    time_limit = std::numeric_limits<float>::max();
                    tail_off_counter = 0;
                    old_val = -std::numeric_limits<float>::max(); //reset
                    goto ENTER;
                }
                CuttingDetail::rollbackEasyWay(node, rollback_solver, rollback_cols, rollback_rccs, rollback_r1cs,
                                               rollback_brcs);
                goto QUIT;
            }
        }
        if (node->getIfInEnumState() != old_enu_state) tail_off_counter = 0;

    SET_TAIL_OFF:
        CuttingDetail::setIfTailOff(eps, ub, node->getValue(), old_val, node->getBrValueImproved(),
                                    node->getIfInEnumState(), eps_max,
                                    if_tail_off, tail_off_counter);

        if (node->tellIf4MIP()) {
            terminateByMIP(node);
        }
        if (node->getIfTerminate()) {
            goto QUIT;
        }
        if (!if_tail_off) goto CUTTING;
    QUIT:
        if (rollback_solver.model) rollback_solver.freeModel();
        if (!node->getIfTerminate() && !node->getIfInEnumState()) {
            node->cleanIndexColForNode();
        }
    }

    void VRPTW::addFeasibilityCuts(int &num_row,
                                   const std::vector<double> &x,
                                   const std::vector<SequenceInfo> &cols,
                                   std::vector<Rcc> &existing_rccs,
                                   Solver &solver) {
        std::vector<double> sol_x;
        std::vector<SequenceInfo> sols;
        for (int i = 0; i < x.size(); ++i) {
            if (x[i] > SOL_X_TOLERANCE) {
                sol_x.emplace_back(x[i]);
                sols.emplace_back(cols[i]);
            }
        }
        CuttingDetail::callRCC(getDim(), getCap(), getDemand(), true, false,
                               num_row, sol_x, sols, cols, existing_rccs, solver);
    }
}
