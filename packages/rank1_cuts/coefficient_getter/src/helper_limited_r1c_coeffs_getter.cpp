/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "rank1_coefficient_controller.hpp"


namespace RouteOpt::Rank1Cuts::CoefficientGetter {
    void Rank1CoefficientGetter::getLimitedR1CPre(
        const std::vector<R1c> &cuts) {
        const auto &rank1CutsDataShared = data_shared_ref.get();
        lp_r1c_denominator.resize(cuts.size());
        auto dim = rank1CutsDataShared.getDim();
        if (lp_v_cut_map.empty()) lp_v_cut_map.resize(dim);
        if (lp_v_v_use_states.empty())lp_v_v_use_states.resize(dim, std::vector<r1cIndex>(dim));
        std::fill(lp_v_cut_map.begin(), lp_v_cut_map.end(),
                  std::make_pair(std::vector<int>(cuts.size(), 0), std::vector<int>()));
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                lp_v_v_use_states[i][j].reset();
            }
        }

        int num = 0;
        for (const auto &r1c: cuts) {
            // std::vector<int> multi;
            // int denominator, rhs;
            int size = static_cast<int>(r1c.info_r1c.first.size());
            auto denominator = rank1CutsDataShared.getDenominator(size, r1c.info_r1c.second);
            const auto &multi = rank1CutsDataShared.getMultiplier(size, r1c.info_r1c.second);
            lp_r1c_denominator[num] = denominator;

            for (auto &[fst, snd]: r1c.arc_mem) {
                lp_v_v_use_states[fst][snd].set(num);
                lp_v_v_use_states[snd][fst].set(num); // keep symmetry
            }

            for (int j = 0; j < r1c.info_r1c.first.size(); ++j) {
                int n = r1c.info_r1c.first[j];
                auto &tmp_n = lp_v_cut_map[n];
                int add = multi[j];
                tmp_n.first[num] += add;
                tmp_n.second.emplace_back(num);
            }
            ++num;
        }
    }

    void Rank1CoefficientGetter::getCoefficientExtendR1C(std::vector<int> &states,
                                                         std::vector<int> &sparse_rep,
                                                         std::unordered_map<int, int> &cnt,
                                                         int &valid_sparse_num,
                                                         int from,
                                                         int to
    ) {
        const auto &cut_mem = lp_v_v_use_states[from][to];
        const auto &add_vec = lp_v_cut_map[to].first;
        r1cIndex record = 0;
        for (int i = 0; i < valid_sparse_num;) {
            int idx = sparse_rep[i];
            if (cut_mem.test(idx)) {
                auto state = states[idx] + add_vec[idx];
                record.set(idx);
                if (state < lp_r1c_denominator[idx]) {
                    ++i;
                    states[idx] = state;
                } else if (state == lp_r1c_denominator[idx]) {
                    ++cnt[idx];
                    states[idx] = 0;
                    sparse_rep[i] = sparse_rep[--valid_sparse_num];
                } else {
                    ++cnt[idx];
                    ++i;
                    states[idx] = state - lp_r1c_denominator[idx];
                }
            } else {
                states[idx] = 0;
                sparse_rep[i] = sparse_rep[--valid_sparse_num];
            }
        }

        for (auto &pr: lp_v_cut_map[to].second) {
            if (!record.test(pr)) {
                sparse_rep[valid_sparse_num++] = pr;
                states[pr] = add_vec[pr];
            }
        }
    }
}
