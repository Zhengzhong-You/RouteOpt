/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_DOMINANCE_HPP
#define ROUTE_OPT_DOMINANCE_HPP
#include "pricing_macro.hpp"
#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    /**
     * Check reduced cost dominance using R1C cuts
     */
    template<bool dir>
    bool CVRP_Pricing::doRCTermDominance(Label *ki, Label *kj) {
        double gap = kj->rc - ki->rc;
        if (!rank1_rc_controller_ref.get().doR1CDominance(gap, ki->r1c, kj->r1c)) return false;
        return true;
    }

    /**
     * Core dominance test between two labels
     * Checks resource consumption and visited nodes
     */
    template<bool dir, PRICING_LEVEL pricing_level>
    bool CVRP_Pricing::dominanceCore(Label *ki, Label *kj) {
        if constexpr (pricing_level == PRICING_LEVEL::EXACT) {
            //exact
            if (dir
                    ? tellResTupleRelations<'p'>(ki->res, kj->res)
                    : tellResTupleRelations<'p'>(kj->res, ki->res))
                return
                        false;
            if (((ki->pi & kj->pi) ^ (ki->pi)).any()) return false;
            if (!doRCTermDominance<dir>(ki, kj)) return false;
        } else if (pricing_level == PRICING_LEVEL::HEAVY) {
            if (dir
                    ? tellResTupleRelations<'p'>(ki->res, kj->res)
                    : tellResTupleRelations<'p'>(kj->res, ki->res))
                return
                        false;
        } else {
        }
        return true;
    }

    /**
     * Perform dominance check when adding new label
     * Removes dominated labels and checks if new label is dominated
     */
    template<bool dir, PRICING_LEVEL pricing_level>
    void CVRP_Pricing::doDominance(Label *ki, int j, int bj, bool &if_suc) {
        auto &labelList_j = dir ? label_array_in_forward_sense[j][bj] : label_array_in_backward_sense[j][bj];
        auto new_label = all_label + idx_glo;
        auto &tmp_rc = new_label->rc;
        double len;

        double tmp_rc_add = tmp_rc + RC_TOLERANCE
                , tmp_rc_sub = tmp_rc - RC_TOLERANCE;
        if_suc = true;
        if constexpr (CHECK_PRICING_LABELS) len = 0;
        auto it = labelList_j.begin();
        for (; it != labelList_j.end(); ++it) {
            ++num_dominance_checks;
            if constexpr (CHECK_PRICING_LABELS) ++len;
            auto kj = *it;
            if (kj->rc < tmp_rc_add) {
                if (dominanceCore<dir, pricing_level>(kj, new_label)) {
                HERE1:
                    labelList_j.splice(labelList_j.begin(), labelList_j, it);
                    if constexpr (CHECK_PRICING_LABELS) {
                        inner_bin_len.first += len;
                        inner_bin_len.second++;
                    }
                    if_suc = false;
                    return;
                }
            } else if (kj->rc > tmp_rc_sub) {
                if (dominanceCore<dir, pricing_level>(new_label, kj)) {
                HERE2:
                    kj->is_extended = true;
                    it = labelList_j.erase(it);
                    break;
                }
            } else {
                if (dominanceCore<dir, pricing_level>(kj, new_label)) {
                    goto HERE1;
                } else if (dominanceCore<dir, pricing_level>(new_label, kj)) {
                    goto HERE2;
                }
            }
        }

        for (; it != labelList_j.end();) {
            auto kj = *it;
            if (kj->rc < tmp_rc_add) {
                ++it;
                continue;
            }
            if (dominanceCore<dir, pricing_level>(new_label, kj)) {
                kj->is_extended = true;
                it = labelList_j.erase(it);
            } else ++it;
        }

        labelList_j.push_front(new_label);

        new_label->p_label = ki;
        new_label->is_extended = false;
        auto &bucket = dir
                           ? if_exist_extra_labels_in_forward_sense[j][bj]
                           : if_exist_extra_labels_in_backward_sense[j][bj];
        bucket.first[bucket.second++] = new_label;
        if (bucket.second == bucket.first.size()) {
            bucket.first.resize(bucket.first.size() * 2);
        }
    }

    /**
     * Check if a label is dominated by existing labels in previous buckets
     * Used to eliminate labels before extension
     */
    template<bool dir, PRICING_LEVEL pricing_level>
    void CVRP_Pricing::checkIfDominated(Label *&ki, int i, int b,
                                        bool &if_suc) {
        if_suc = true;
        double tmp_ki_rc_sub = ki->rc - RC_TOLERANCE;
        int len;
        if constexpr (CHECK_PRICING_LABELS) len = 0;
        for (int b4_b = (dir ? b - 1 : b + 1); dir ? b4_b >= 0 : b4_b < num_buckets_per_vertex; dir ? --b4_b : ++b4_b) {
            auto &b4_label_list = dir ? label_array_in_forward_sense[i][b4_b] : label_array_in_backward_sense[i][b4_b];
            if ((dir ? rc2_till_this_bin_in_forward_sense[i][b4_b] : rc2_till_this_bin_in_backward_sense[i][b4_b])
                > tmp_ki_rc_sub)
                break;
            ++num_dominance_checks;
            for (auto &p: b4_label_list) {
                if constexpr (CHECK_PRICING_LABELS) ++len;
                if (p->rc > tmp_ki_rc_sub) break;
                if (dominanceCore<dir, pricing_level>(p, ki)) {
                    if_suc = false;
                    if constexpr (CHECK_PRICING_LABELS) {
                        outer_bin_len.first += len;
                        outer_bin_len.second++;
                    }
                    return;
                }
            }
        }
        if constexpr (CHECK_PRICING_LABELS) {
            outer_bin_but_keep_len.first += len;
            outer_bin_but_keep_len.second++;
        }
    }
}

#endif // ROUTE_OPT_DOMINANCE_HPP
