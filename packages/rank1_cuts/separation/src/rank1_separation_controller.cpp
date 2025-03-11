/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "rank1_separation_controller.hpp"
#include "rank1_data_shared.hpp"
#include "solver.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    Rank1SeparationController::Rank1SeparationController(
        const Rank1CutsDataShared &rank1CutsDataShared,
        int max_row_rank1,
        int max_num_r1c3_per_round,
        int max_num_r1c_per_round,
        Solver &solver,
        const std::vector<std::vector<double> > &
        cost_mat4_vertex) : sharedData(cost_mat4_vertex),
                            cutGen(rank1CutsDataShared, sharedData),
                            memGen(rank1CutsDataShared, sharedData),
                            cutSelector(rank1CutsDataShared, sharedData),
                            rank1CutsDataShared_ref(rank1CutsDataShared) {
        // Set the data using setter methods
        cutGen.setMaxRowRank1(max_row_rank1);
        cutGen.setMaxNumR1c3PerRound(max_num_r1c3_per_round);
        cutGen.setMaxNumR1cPerRound(max_num_r1c_per_round);
        cutGen.setMaxHeuristicInitialSeedSetSizeRowRank1c(max_row_rank1 + 1);
        memGen.setSolver(solver);

        // Call additional methods to std::generate or std::set up further data if required
        cutGen.preprocess();
    }

    void Rank1SeparationController::setRhs(std::vector<R1c> &cuts) {
        for (auto &cut: cuts) {
            cut.rhs = rank1CutsDataShared_ref.get().getRhs(static_cast<int>(cut.info_r1c.first.size()),
                                                           cut.info_r1c.second);
        }
    }

    void Rank1SeparationController::separateRank1Cuts(
        std::vector<R1c> &cuts, bool if_mem,
        bool if_select_cuts) {
        if (if_mem) memGen.fillMemory();
        cutGen.generateRank1Cuts();
        if (if_mem) memGen.findMemory4Cuts();
        if (if_select_cuts) cutSelector.selectR1CsByVioNMemory();
        setRhs(sharedData.refCuts());
        cuts = sharedData.refCuts();
        cleanData();
    }

    void Rank1SeparationController::updateInfo(MemoryType limited_memory_type,
                                               PRICING_HARD_LEVEL pricing_hard_level,
                                               bool if_collect_sol,
                                               const std::vector<RouteInfo> &sol,
                                               const std::vector<R1c> &old_cuts) {
        // Update using setter methods
        cleanData();
        if (limited_memory_type == MemoryType::ARC_MEMORY) {
            if_once_use_no_symmetry_mem = true;
        }
        sharedData.setLimitedMemoryType(limited_memory_type);
        memGen.setPricingHardLevel(pricing_hard_level);
        auto sort_sol = sol;
        std::sort(sort_sol.begin(), sort_sol.end(),
                  [](const RouteInfo &a, const RouteInfo &b) {
                      return a.frac_x > b.frac_x;
                  });


        if (limited_memory_type == MemoryType::ARC_MEMORY) {
            //since hard pricing, then we carefully treat memory
            for (auto it = sort_sol.begin(); it != sort_sol.end();) {
                if (it->frac_x < SOL_X_RANK1_TOLERANCE || it->frac_x > 1 - SOL_X_RANK1_TOLERANCE) {
                    it = sort_sol.erase(it);
                } else {
                    ++it;
                }
            }
        }

        if (if_collect_sol && limited_memory_type == MemoryType::NODE_MEMORY) {
            sol_vec.emplace_back(sort_sol);
            if (sol_vec.size() > SOL_KEEPER_LIMIT) {
                PRINT_WARNING(
                    "The solution is retained but not yet utilized.\n"
                    "Are you sure it still needs to be kept?\n"
                    "The limit is " + std::to_string(SOL_KEEPER_LIMIT) +
                    ", and the earliest solution has been removed.");
                sol_vec.erase(sol_vec.begin());
            }
        } else sol_vec.clear();

        sharedData.setSol(sort_sol);
        sharedData.setOldCuts(old_cuts);
        sharedData.refCuts() = {};
    }

    void Rank1SeparationController::convert2ArcMemory(std::vector<R1c> &existing_cuts) {
        PRINT_REMIND("convert to arc memory since the pricing is too hard!");
        for (auto &r1c: existing_cuts) r1c.arc_mem.clear();
        for (auto &sol: sol_vec) {
            cleanData();
            if_once_use_no_symmetry_mem = true;
            sharedData.setLimitedMemoryType(MemoryType::ARC_MEMORY);
            memGen.setPricingHardLevel(PRICING_HARD_LEVEL::EXTREMELY_HARD);
            sharedData.setSol(sol);
            sharedData.refCuts() = existing_cuts;
            memGen.findMemory4Cuts();
            const auto &cuts = sharedData.refCuts();
            if (cuts.size() != existing_cuts.size()) {
                THROW_RUNTIME_ERROR("#cuts is changed!");
            }
            for (int i = 0; i < cuts.size(); ++i) {
                if (existing_cuts[i].info_r1c != cuts[i].info_r1c) {
                    THROW_RUNTIME_ERROR("existing cuts are changed!");
                }
                existing_cuts[i].arc_mem = cuts[i].arc_mem;
            }
        }
        sol_vec.clear();
        cleanData();
    }
}
