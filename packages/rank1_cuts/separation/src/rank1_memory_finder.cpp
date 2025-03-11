/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include <iostream>
#include "rank1_macro.hpp"
#include "rank1_memory_finder.hpp"
#include "helper_rank1_node_based_memory.hpp"
#include "helper_rank1_arc_based_memory.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    void MemGenerator::findMemory4Cuts() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        sharedData.getLimitedMemoryType() == MemoryType::NODE_MEMORY
            ? constructMemoryVertexBased(rank1CutsDataShared, sharedData, pricing_hard_level, solver,
                                         rank1_multi_mem_plan_map)
            : constructMemoryArcBased(rank1CutsDataShared, sharedData);
    }

    void MemGenerator::fillMemory() const {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        if (sharedData.getLimitedMemoryType() == MemoryType::NO_MEMORY) return;
        if (!sharedData.refCuts().empty())
            THROW_RUNTIME_ERROR("cuts should be empty at first!");
        std::vector<std::unordered_map<int, int> > local_v_r_map(rank1CutsDataShared.getDim());
        const auto &sol = sharedData.getSol();
        const auto &old_cuts = sharedData.getOldCuts();
        int r_idx = 0;
        for (const auto &i: sol) {
            for (const auto &j: i.col_seq) {
                ++local_v_r_map[j][r_idx];
            }
            ++r_idx;
        }
        int num_add_mem = 0;
        std::vector<double> vis_times(sol.size());
        int idx = 0;

        auto &cuts = sharedData.refCuts();
        auto &cut_record = sharedData.refCutRecord();

        for (auto &r1c: old_cuts) {
            memset(&vis_times[0], 0, sizeof(double) * sol.size());
            const auto &coeff = rank1CutsDataShared.getMultiplier(static_cast<int>(r1c.info_r1c.first.size()),
                                                              r1c.info_r1c.second);
            auto deno = rank1CutsDataShared.getDenominator(static_cast<int>(r1c.info_r1c.first.size()),
                                                           r1c.info_r1c.second);
            auto rhs = rank1CutsDataShared.getRhs(static_cast<int>(r1c.info_r1c.first.size()),
                                                  r1c.info_r1c.second);
            int cnt = 0;
            for (const auto &v: r1c.info_r1c.first) {
                for (auto &[fst, snd]: local_v_r_map[v]) {
                    vis_times[fst] += snd * coeff[cnt];
                }
                ++cnt;
            }

            std::transform(vis_times.begin(), vis_times.end(), sol.begin(), vis_times.begin(),
                           [deno](auto &a, auto &b) {
                               return static_cast<int>(a / deno + RANK1_TOLERANCE) * b.frac_x;
                           });

            if (auto vio = std::accumulate(vis_times.begin(), vis_times.end(), -static_cast<double>(rhs));
                vio > RANK1_TOLERANCE
            ) {
                cuts.emplace_back(r1c);
                cut_record[r1c.info_r1c.first].insert(r1c.info_r1c.second);
                ++num_add_mem;
            }
            ++idx;
        }
        RANK1_VERBOSE_EXEC(std::cout << "num_add_mem= " << num_add_mem << std::endl;)
    }
} // namespace RouteOpt::Rank1Cuts::Separation
