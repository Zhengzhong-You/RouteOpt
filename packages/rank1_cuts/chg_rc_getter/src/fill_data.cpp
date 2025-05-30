/*
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <unordered_set>

#include "rank1_rc_controller.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Rank1Cuts::RCGetter {
    void Rank1RCController::getRank1DualsInCG(
        const std::vector<R1c> &cuts,
        const std::vector<double> &pi_vector) {
        auto dim = data_shared_ref.get().getDim();


        std::vector<std::pair<int, double> > cut_dual;
        for (int i = 0; i < cuts.size(); ++i) {
            auto &r1c = cuts[i];
            double dual = std::abs(pi_vector[r1c.idx_r1c]);
            if (dual < DUAL_TOLERANCE) {
                continue;
            }
            cut_dual.emplace_back(i, dual);
        }

        std::sort(cut_dual.begin(), cut_dual.end(),
                  [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
                      return a.second > b.second;
                  });

        cg_v_cut_map.assign(dim, R1CUseStates(cut_dual.size()));
        if (cg_v_v_use_states.empty()) cg_v_v_use_states.resize(dim, std::vector<r1cIndex>(dim));
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                cg_v_v_use_states[i][j].reset();
            }
        }

        rank1_dual.resize(cut_dual.size());
        cg_r1c_denominator.resize(cut_dual.size());
        revised_rank1_dual.resize(cut_dual.size());

        int num = 0;
        for (auto &pr: cut_dual) {
            int i = pr.first;
            auto &r1c = cuts[i];
            auto cut_size = static_cast<int>(r1c.info_r1c.first.size());
            const auto &multi = data_shared_ref.get().getMultiplier(cut_size,
                                                                    r1c.info_r1c.second);
            auto denominator = data_shared_ref.get().getDenominator(cut_size,
                                                                    r1c.info_r1c.second);
            rank1_dual[num] = pi_vector[r1c.idx_r1c];
            cg_r1c_denominator[num] = denominator;
            revised_rank1_dual[num] = rank1_dual[num]; //reserved for future use

            std::unordered_set<int> mem_set;

            for (auto &[fst, snd]: r1c.arc_mem) {
                cg_v_v_use_states[fst][snd].set(num);
                cg_v_v_use_states[snd][fst].set(num); // keep symmetry
                mem_set.emplace(fst);
                mem_set.emplace(snd);
            }


            for (int j = 0; j < r1c.info_r1c.first.size(); ++j) {
                int n = r1c.info_r1c.first[j];
                mem_set.emplace(n);
                auto &tmp_n = cg_v_cut_map[n];
                int add = multi[j];
                tmp_n.add_vec.first[num] += add;
                tmp_n.add_vec.second.emplace_back(num);
            }

            for (auto &n: mem_set) {
                cg_v_cut_map[n].v_union_mem.emplace_back(num);
            }

            ++num;
        }
    }
}
