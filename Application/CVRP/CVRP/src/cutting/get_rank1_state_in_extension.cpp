//
// Created by You, Zhengzhong on 11/25/23.
//

#include "cvrp.hpp"
#include "cutting/get_rank1_matrix.hpp"
using namespace std;

void CVRP::getLimitedR1CPre(BbNode *node,
                            const std::vector<int> &idx) {
    /**
     * get std::vector<std::vector<int>> lp_v_union_mem{};
            std::vector<std::vector<R1CUseStates>> lp_v_v_use_states{};
            std::vector<int> lp_r1c_denominator{};
     * to satisfy the requirement of getLimitedR1CCoeffs
     */
    lp_r1c_denominator.resize(idx.size());
    fill(lp_v_cut_map.begin(), lp_v_cut_map.end(), make_pair(R1CINDEX(), vector<int>()));

    for (int i = 0; i < dim; ++i) {
        for (int j = 1; j < dim; ++j) {
            lp_v_v_use_states[i][j].assign(idx.size(), RANK1_INVALID);
        }
    }

    int num = 0;
    for (auto i: idx) {
        auto &r1c = node->r1cs[i];
        vector<int> multi;
        int denominator, rhs;
        Rank1CutsSeparator::getMapPlanInfo(multi, denominator, rhs, static_cast<int>(r1c.info_r1c.first.size()),
                                      r1c.info_r1c.second);
        lp_r1c_denominator[num] = denominator;
        for (int j = 0; j < r1c.info_r1c.first.size(); ++j) {
            int n = r1c.info_r1c.first[j];
            auto &tmp_n = lp_v_cut_map[n];
            int add = multi[j];
            tmp_n.first.set(num);
            tmp_n.second.emplace_back(num);
            for (int k = 0; k < dim; ++k) {
                lp_v_v_use_states[k][n][num] = add;
            }
        }
        for (auto &m: r1c.arc_mem) {
            for (auto &k: m.first) {
                lp_v_v_use_states[k][m.second][num] = 0;
            }
        }
        ++num;
    }
}

void CVRP::getCoefficientExtendR1C(std::vector<int> &states,
                                   std::vector<int> &sparse_rep,
                                   std::unordered_map<int, int> &cnt,
                                   int &valid_sparse_num,
                                   int from,
                                   int to
) {
    auto &v_cut_map = lp_v_v_use_states[from][to];
    auto record = lp_v_cut_map[to].first;
    for (int i = 0; i < valid_sparse_num;) {
        int idx = sparse_rep[i];
        int add = v_cut_map[idx];
        bool is_add = true;
        if (add > 0) {
            int state = states[idx] + add;
            record.reset(idx);
            if (state == lp_r1c_denominator[idx]) {
                ++cnt[idx];
                states[idx] = 0;
                is_add = false;
            } else if (state < lp_r1c_denominator[idx]) {
                is_add = true;
                states[idx] = state;
            } else {
                ++cnt[idx];
                is_add = true;
                states[idx] = state - lp_r1c_denominator[idx];
            }
        } else if (add == RANK1_INVALID) {
            is_add = false;
            states[idx] = 0;
        }
        if (is_add) ++i;
        else {
            sparse_rep[i] = sparse_rep[--valid_sparse_num];
        }
    }
    for (auto &pr: lp_v_cut_map[to].second) {
        if (record.test(pr)) {
            sparse_rep[valid_sparse_num++] = pr;
            states[pr] = v_cut_map[pr];
        }
    }
}
