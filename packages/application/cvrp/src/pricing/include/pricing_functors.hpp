/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_PRICING_FUNCTORS_HPP
#define ROUTE_OPT_PRICING_FUNCTORS_HPP
#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    inline bool CVRP_Pricing::increaseMainResourceConsumption(const Resource &nowResource,
                                                              Resource &newResource, int start, int end) {
        newResource = nowResource + resource_across_arcs_in_forward_sense[start][end];
        if (tellResTupleRelations<'p'>(newResource, ub4_vertex[end])) return false;
        newResource.takeLarger(lb4_vertex[end]);
        if (tellResTupleRelations<'c'>(newResource, resource_across_arcs_in_forward_sense[end][0])) return false;
        return true;
    }

    inline bool CVRP_Pricing::decreaseMainResourceConsumption(const Resource &nowResource,
                                                              Resource &newResource, int start, int end) {
        newResource = nowResource - resource_across_arcs_in_backward_sense[start][end];
        if (tellResTupleRelations<'p'>(lb4_vertex[end], newResource)) return false;
        newResource.takeSmaller(ub4_vertex[end]);
        /**
         * not type c
         */
        if (tellResTupleRelations<'p'>({}, newResource - resource_across_arcs_in_backward_sense[end][0])) return false;
        return true;
    }

    template<char type>
    bool CVRP_Pricing::tellResTupleRelations(const Resource &res1, const Resource &res2) const {
        if constexpr (type == 'c') {
            /**
             * concatenate: all is no larger than max_res, return false
             */
            if ((res1 + res2) <= resource) return false;
            return true;
        } else if (type == 'p') {
            /**
             * compare: all is no larger, return false
             */
            if (res1 <= res2) return false;
            return true;
        } else if (type == 'e') {
            /**
             * equal: both are equal, return true
             */
            if (res1 == res2) return true;
            return false;
        }
        THROW_RUNTIME_ERROR("tellResTupleRelations: unknown type");
    }

    inline void CVRP_Pricing::addPathByRC(double path_rc, Label *ki, Label *kj, int num) {
        if (path_rc < rc_std) {
            negative_rc_label_tuple.emplace_back(ki, kj, path_rc);
            if (negative_rc_label_tuple.size() >= num)
                rc_std = std::get<2>(negative_rc_label_tuple[negative_rc_label_tuple.size() - num]);
        }
    }

    template<bool dir>
    void CVRP_Pricing::sortLabelsInBinByRC(int i, int b) {
        auto &label_bin = (dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b]);
        label_bin.sort([](Label *a, Label *b) { return a->rc < b->rc; });
        auto
                &rc2_till_this_bin = (dir
                                          ? rc2_till_this_bin_in_forward_sense[i][b]
                                          : rc2_till_this_bin_in_backward_sense[i][b]);
        auto &rc2_bin = (dir ? rc2_bin_in_forward_sense[i][b] : rc2_bin_in_backward_sense[i][b]);
        rc2_bin = label_bin.empty() ? std::numeric_limits<float>::max() : label_bin.front()->rc;
        if constexpr (dir) {
            rc2_till_this_bin = (b
                                     ? std::min(rc2_till_this_bin_in_forward_sense[i][b - 1], rc2_bin)
                                     : rc2_bin);
        } else {
            rc2_till_this_bin = (b < num_buckets_per_vertex - 1
                                     ? std::min(rc2_till_this_bin_in_backward_sense[i][b + 1],
                                                rc2_bin)
                                     : rc2_bin);
        }
    }

    template<bool dir, bool if_clear_concatenate>
    void CVRP_Pricing::cleanAllPointers() {
        auto &label_array = dir ? label_array_in_forward_sense : label_array_in_backward_sense;
        auto &if_exist_extra_labels = dir
                                          ? if_exist_extra_labels_in_forward_sense
                                          : if_exist_extra_labels_in_backward_sense;

        for (int i = 0; i < dim; ++i) {
            for (int b = 0; b < num_buckets_per_vertex; ++b) {
                label_array[i][b].clear();
                if_exist_extra_labels[i][b].second = 0;
            }
        }

        if (if_clear_concatenate) {
            dir ? concatenate_labels_in_forward_cg.clear() : concatenate_labels_in_backward_cg.clear();
        }
    }

    template<bool dir, bool if_symmetry, bool if_reset_label_point, bool if_clear_all, bool if_clear_concatenate>
    void CVRP_Pricing::initializeLabels() {
        if constexpr (if_reset_label_point) {
            idx_glo = if_symmetry ? dim : 2 * dim - 1;
            rc_std = RC_TOLERANCE;
            num_dominance_checks = 0;
            negative_rc_label_tuple.clear();
        }


        if constexpr (if_clear_all) {
            cleanAllPointers<dir, if_clear_concatenate>();
        }


        if constexpr (dir) {
            can_leave_depot_forward.reset();
            for (auto j: all_forward_buckets[0][0].bucket_arcs) can_leave_depot_forward.set(j);
            for (int i = 1; i < dim; ++i) {
                if (!can_leave_depot_forward.test(i)) continue;
                auto &new_label = all_label[i];
                new_label.rc = chg_cost_mat4_vertex[0][i];
                new_label.is_extended = false;
                rank1_rc_controller_ref.get().updateR1CStates(new_label.rc, new_label.r1c, all_label->r1c, 0, i);
                auto bin = static_cast<int>(new_label.res.resources[0] / step_size);
                auto &bucket = label_array_in_forward_sense[i][bin];
                bucket.push_front(all_label + i);
                auto &bucket2 = if_exist_extra_labels_in_forward_sense[i][bin];
                bucket2.first[bucket2.second++] = all_label + i;
                auto rc_return = new_label.rc + chg_cost_mat4_vertex[i][0];
                if (adjust_brc_dual4_single_route.find(i) != adjust_brc_dual4_single_route.end()) {
                    rc_return += adjust_brc_dual4_single_route[i];
                }
                addPathByRC(rc_return, all_label + i, nullptr, num_col_generated_ub);
            }
        } else {
            can_leave_depot_backward.reset();
            for (auto j: all_backward_buckets[0][0].bucket_arcs) can_leave_depot_backward.set(j);
            int max_num = 2 * dim - 1;
            for (int i = dim; i < max_num; ++i) {
                int point = i - dim + 1;
                if (!can_leave_depot_backward.test(point)) continue;
                auto &new_label = all_label[i];
                new_label.rc = chg_cost_mat4_vertex[0][point];
                new_label.is_extended = false;
                rank1_rc_controller_ref.get().updateR1CStates(new_label.rc, new_label.r1c, all_label->r1c, 0, point);
                auto bin = static_cast<int>(new_label.res.resources[0] / step_size);
                auto &bucket = label_array_in_backward_sense[point][bin];
                bucket.push_front(all_label + i);
                auto &bucket2 = if_exist_extra_labels_in_backward_sense[point][bin];
                bucket2.first[bucket2.second++] = all_label + i;
                auto rc_return = new_label.rc + chg_cost_mat4_vertex[point][0];
                if (adjust_brc_dual4_single_route.find(point) != adjust_brc_dual4_single_route.end()) {
                    rc_return += adjust_brc_dual4_single_route[point];
                }
                addPathByRC(rc_return, nullptr, all_label + i, num_col_generated_ub);
            }
        }
    }

    template<bool dir>
    void CVRP_Pricing::populateRC2TillThisBinNRC2Bin() {
        for (int i = 1; i < dim; ++i) {
            for (int b = 0; b < num_buckets_per_vertex; ++b) {
                sortLabelsInBinByRC<dir>(i, b);
            }
        }
    }

    template<bool dir, bool if_symmetry, bool if_std_optgap>
    void CVRP_Pricing::concatenateOneLabelWithOtherLabels(Label *ki, int j, int arr_bj, double tmp_rc,
                                                          const Resource &tmp_res,
                                                          int &if_state) {
        double path_rc;
        double &which_rc = if_std_optgap ? opt_gap : rc_std;
        auto ptr_rc_till_this_bin = &rc2_till_this_bin_in_forward_sense[j][arr_bj];
        auto ptr_rc_bin = &rc2_bin_in_forward_sense[j][arr_bj];
        if constexpr ((dir && !if_symmetry) || (!dir && if_symmetry)) {
            ptr_rc_till_this_bin = &rc2_till_this_bin_in_backward_sense[j][arr_bj];
            ptr_rc_bin = &rc2_bin_in_backward_sense[j][arr_bj];
        }

        if (*ptr_rc_till_this_bin + tmp_rc > which_rc) {
            if_state = static_cast<int>(CONCATENATE_STATE::TOTALLY_FAIL); //no need to continue
            return;
        }
        if_state = static_cast<int>(CONCATENATE_STATE::FAIL);
        if (*ptr_rc_bin + tmp_rc < which_rc) {
            //most_negative_rc_in_this_bin
            auto &label_arr =
                    ((!dir && !if_symmetry) || (dir && if_symmetry))
                        ? label_array_in_forward_sense[j][arr_bj]
                        : label_array_in_backward_sense[j][arr_bj];
            for (auto &kj: label_arr) {
                path_rc = kj->rc + tmp_rc;
                if (path_rc > which_rc) break;
                if constexpr (if_symmetry) {
                    if (tellResTupleRelations<'c'>(tmp_res, kj->res)) continue;
                } else {
                    if (dir
                            ? tellResTupleRelations<'p'>(tmp_res, kj->res)
                            : tellResTupleRelations<'p'>(kj->res, tmp_res))
                        continue;
                }

                if ((ki->pi & kj->pi).any()) continue;

                if (!rank1_rc_controller_ref.get().concatenateR1CStates(path_rc, which_rc, ki->r1c, kj->r1c,
                                                                        ki->end_vertex, j))
                    continue;

                if constexpr (!if_std_optgap) {
                    addPathByRC(path_rc, ki, kj, num_col_generated_ub);
                } else {
                    if (path_rc < opt_gap) {
                        if_state = arr_bj; //what bin to concatenate
                        break;
                    }
                }
            }
        }
    }

    template<bool dir, bool if_last_half, bool if_complete, bool if_symmetry, bool if_std_optgap, bool if_res_updated,
        PRICING_LEVEL pricing_level>
    void CVRP_Pricing::updateLabel(const Resource &res, Label *ki, int i, int j, int &bj,
                                   bool &if_suc) {
        if_suc = false;
        if (ki->pi[j]) return;
        auto new_label = all_label + idx_glo;
        auto &tmp_res = new_label->res;
        auto &tmp_rc = new_label->rc;
        tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j]; //real rc


        if constexpr (!if_complete) {
            if constexpr (!if_res_updated) {
                if constexpr (dir) {
                    if (!increaseMainResourceConsumption(res, tmp_res, i, j)) return;
                    if constexpr (!if_last_half) {
                        if (tmp_res.resources[0] > meet_point_resource_in_bi_dir) {
                            concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_res);
                            return;
                        }
                    }
                } else {
                    if (!decreaseMainResourceConsumption(res, tmp_res, i, j)) return;
                    if constexpr (!if_last_half) {
                        if (tmp_res.resources[0] < meet_point_resource_in_bi_dir) {
                            concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_res);
                            return;
                        }
                    }
                }
            }
            if constexpr (if_last_half) {
                int concate_bj;
                int if_state;
                if constexpr (if_symmetry) {
                    concate_bj = static_cast<int>((resource.resources[0] - tmp_res.resources[0]) / step_size);
                } else {
                    concate_bj = static_cast<int>(tmp_res.resources[0] / step_size);
                }
                constexpr bool if_dif = dir ^ if_symmetry;
                for (; if_dif ? concate_bj < num_buckets_per_vertex : concate_bj >= 0;
                       if_dif ? ++concate_bj : --concate_bj) {
                    concatenateOneLabelWithOtherLabels<dir, if_symmetry, true>(ki,
                                                                               j,
                                                                               concate_bj,
                                                                               tmp_rc,
                                                                               tmp_res,
                                                                               if_state);
                    if (if_state == static_cast<int>(CONCATENATE_STATE::TOTALLY_FAIL)) break;
                    if (if_state >= 0) goto outside;
                }
            outside:
                if (if_state < 0) {
                    // this must be zero!
                    return;
                }
            }
        } else {
            if constexpr (dir) {
                if (!increaseMainResourceConsumption(res, tmp_res, i, j)) return;
            } else {
                if (!decreaseMainResourceConsumption(res, tmp_res, i, j)) return;
            }
        }

        bj = static_cast<int>(tmp_res.resources[0] / step_size);
        auto &tmp_PI = new_label->pi;
        tmp_PI = (ki->pi) & (ng_mem4_vertex[j]);
        tmp_PI.set(j);
        new_label->end_vertex = j;

        rank1_rc_controller_ref.get().updateR1CStates(tmp_rc, new_label->r1c, ki->r1c, i, j);
        if_suc = true;
    }


    template<typename T, bool dir, bool if_last_half, bool if_complete, bool if_symmetry, PRICING_LEVEL pricing_level>
    PRICING_STATE CVRP_Pricing::extendKernel4Exact(Label *ki,
                                                   int i,
                                                   Resource res,
                                                   const std::vector<T> &arc) {
        auto state = PRICING_STATE::NORMAL;
        int bj;
        bool if_suc;

        for (auto &pr: arc) {
            int j;
            if constexpr (std::is_same<T, int>::value) {
                j = pr;
            } else {
                j = pr.second;
                res.resources[0] = pr.first;
            }

            updateLabel<dir, if_last_half, if_complete, if_symmetry, false, false, pricing_level>(
                res, ki, i, j, bj,
                if_suc);

            if (!if_suc) continue;

            doDominance<dir, pricing_level>(ki, j, bj, if_suc);

            if (!if_suc) continue;

            if constexpr (!if_last_half) {
                if constexpr (if_symmetry) {
                    if constexpr (dir) {
                        if (can_leave_depot_forward.test(j)) {
                            addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
                                        all_label + idx_glo,
                                        nullptr,
                                        num_col_generated_ub);
                        }
                    } else {
                        if (can_leave_depot_backward.test(j)) {
                            addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
                                        nullptr,
                                        all_label + idx_glo,
                                        num_col_generated_ub);
                        }
                    }
                } else {
                    if constexpr (dir) {
                        if (can_leave_depot_backward.test(j))
                            addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
                                        all_label + idx_glo,
                                        nullptr,
                                        num_col_generated_ub);
                    } else {
                        if (can_leave_depot_forward.test(j))
                            addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
                                        nullptr,
                                        all_label + idx_glo,
                                        num_col_generated_ub);
                    }
                }
            }

            ++idx_glo;
            if (idx_glo == label_assign) {
                if_short_memory = true;
                state = PRICING_STATE::OUT_OF_MEMORY;
                goto QUIT;
            }
        }
    QUIT:

        return state;
    }
}


#endif // ROUTE_OPT_PRICING_FUNCTORS_HPP
