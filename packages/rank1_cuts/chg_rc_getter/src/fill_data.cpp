/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "rank1_rc_controller.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Rank1Cuts::RCGetter {
    void Rank1RCController::getRank1DualsInCG(
        const std::vector<R1c> &cuts,
        const std::vector<double> &pi_vector) {
        auto dim = data_shared_ref.get().getDim();

        if (cg_v_cut_map.empty()) cg_v_cut_map.resize(dim);
        if (cg_v_v_use_states.empty()) cg_v_v_use_states.resize(dim, std::vector<std::vector<int> >(dim));

        std::fill(cg_v_cut_map.begin(), cg_v_cut_map.end(), R1CUseStates());
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

        for (int i = 0; i < dim; ++i) {
            for (int j = 1; j < dim; ++j) {
                cg_v_v_use_states[i][j].assign(cut_dual.size(), RANK1_INVALID);
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
            revised_rank1_dual[num] = rank1_dual[num]; //reserved for future sue

            for (int j = 0; j < r1c.info_r1c.first.size(); ++j) {
                int n = r1c.info_r1c.first[j];
                auto &tmp_n = cg_v_cut_map[n];
                int add = multi[j];
                tmp_n.sparse.emplace_back(num, add);
                tmp_n.sparse_map.set(num);
                tmp_n.v_union_mem.emplace_back(num);
                tmp_n.union_map.set(num);
                for (int k = 0; k < dim; ++k) {
                    cg_v_v_use_states[k][n][num] = add;
                }
            }
            for (auto &m: r1c.arc_mem) {
                cg_v_cut_map[m.second].v_union_mem.emplace_back(num);
                cg_v_cut_map[m.second].union_map.set(num);
                for (auto &k: m.first) {
                    cg_v_v_use_states[k][m.second][num] = 0;
                }
            }
            ++num;
        }
    }
}
