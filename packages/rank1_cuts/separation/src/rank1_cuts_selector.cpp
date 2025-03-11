/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "rank1_cuts_selector.hpp"
#include "rank1_macro.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    void CutSelector::selectR1CsByVioNMemory() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        std::unordered_map<int, std::vector<std::pair<int, int> > > map_v_cut; // key: v; val: idx & multiplier
        const auto &cuts = sharedData.refCuts();
        const auto &old_cuts = sharedData.getOldCuts();
        const auto &sol = sharedData.getSol();
        map_v_cut.reserve(cuts.size());
        std::vector<int> rhs(cuts.size());
        for (int i = 0; i < cuts.size(); ++i) {
            auto &c = cuts[i].info_r1c;
            auto &multipliers = rank1CutsDataShared.getMultiplier(static_cast<int>(c.first.size()), c.second);
            for (int j = 0; j < c.first.size(); ++j) {
                map_v_cut[c.first[j]].emplace_back(i, multipliers[j]);
            }
            rhs[i] = rank1CutsDataShared.getRhs(static_cast<int>(c.first.size()), c.second);
        }
        std::vector<Eigen::RowVectorXd>
                cuts_coeffs(cuts.size(), Eigen::RowVectorXd::Zero(static_cast<int>(sol.size())));
        for (int i = 0; i < sol.size(); ++i) {
            for (const auto &j: sol[i].col_seq) {
                for (auto &k: map_v_cut[j]) {
                    cuts_coeffs[k.first][i] += k.second;
                }
            }
        }

        for (int i = 0; i < cuts.size(); ++i) {
            auto &c = cuts[i].info_r1c;
            const auto denominator = rank1CutsDataShared.getDenominator(static_cast<int>(c.first.size()), c.second);
            for (auto &j: cuts_coeffs[i]) j = static_cast<int>(j / denominator + RANK1_TOLERANCE);
        }

        std::vector<std::pair<double, int> > cuts_vio(cuts.size());
        for (int i = 0; i < cuts.size(); ++i)cuts_vio[i] = {0, i};
        Eigen::RowVectorXd frac_routes_vio(sol.size());
        for (int i = 0; i < sol.size(); ++i)frac_routes_vio[i] = sol[i].frac_x;
        for (int i = 0; i < cuts.size(); ++i)cuts_vio[i].first = cuts_coeffs[i].dot(frac_routes_vio) - rhs[i];

        std::sort(cuts_vio.begin(), cuts_vio.end(), [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
            return a.first > b.first;
        });

        int left_cuts = MAX_NUM_R1CS_IN_PRICING - static_cast<int>(old_cuts.size()) -
                        MIN_CUTS_ADDED_PER_ROUND;
        if (left_cuts < 0) {
            sharedData.refCuts() = {};
            return;
        }

        cuts_vio.resize(std::min(left_cuts, static_cast<int>(cuts_vio.size())));

        std::vector<int> vertex_related_r1c(rank1CutsDataShared.getDim(), MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX - 1);

        for (auto &r1c: old_cuts) {
            for (auto &it: r1c.info_r1c.first) {
                --vertex_related_r1c[it];
            }
            for (auto &it: r1c.arc_mem) {
                --vertex_related_r1c[it.second];
            }
        }

        std::vector<R1c> select_cut(cuts_vio.size());
        int num = 0;
        for (auto &i: cuts_vio) {
            int idx = i.second;
            auto &c = cuts[idx];
            bool if_keep = true;
            auto tmp = vertex_related_r1c;
            for (auto j: c.info_r1c.first) {
                if (tmp[j] <= 0) {
                    if_keep = false;
                    break;
                }
                --tmp[j];
            }
            for (auto &j: c.arc_mem) {
                if (tmp[j.second] <= 0) {
                    if_keep = false;
                    break;
                }
                --tmp[j.second];
            }
            if (!if_keep) continue;
            vertex_related_r1c = tmp;
            select_cut[num++] = c;
        }
        select_cut.resize(num);
        sharedData.refCuts() = select_cut;
    }
}
