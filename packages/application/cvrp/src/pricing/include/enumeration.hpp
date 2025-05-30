/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_ENUMERATION_HPP
#define ROUTE_OPT_ENUMERATION_HPP
#include <valarray>

#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    namespace EnumerationDetail {
        template<bool if_symmetry>
        void adjustEnumerationStdBucketArcs(
            std::pair<double, double> &max_bucket_arc_suc_enumeration,
            std::pair<double, double> &min_bucket_arc_fail_enumeration,
            int num_forward_bucket_arcs,
            int num_backward_bucket_arcs,
            bool if_suc) {
            if (if_suc) {
                max_bucket_arc_suc_enumeration.first =
                        std::max(max_bucket_arc_suc_enumeration.first, static_cast<double>(num_forward_bucket_arcs));
                if constexpr (!if_symmetry) {
                    max_bucket_arc_suc_enumeration.second =
                            std::max(max_bucket_arc_suc_enumeration.second,
                                     static_cast<double>(num_backward_bucket_arcs));
                }
            } else {
                min_bucket_arc_fail_enumeration.first =
                        std::min(min_bucket_arc_fail_enumeration.first, static_cast<double>(num_forward_bucket_arcs));
                if constexpr (!if_symmetry) {
                    min_bucket_arc_fail_enumeration.second =
                            std::min(min_bucket_arc_fail_enumeration.second,
                                     static_cast<double>(num_backward_bucket_arcs));
                }
                min_bucket_arc_fail_enumeration.first =
                        std::max(max_bucket_arc_suc_enumeration.first + 1, min_bucket_arc_fail_enumeration.first);
                if constexpr (!if_symmetry) {
                    min_bucket_arc_fail_enumeration.second =
                            std::max(max_bucket_arc_suc_enumeration.second + 1,
                                     min_bucket_arc_fail_enumeration.second);
                }
            }
        }

        inline double getGapStdTryEnumeration(const std::pair<double, int> &success_enumeration_gap,
                                              double max_gap2_try_enumeration
                                              , double max_enumeration_success_gap
                                              , double min_enumeration_fail_gap) {
            if (success_enumeration_gap.second == 0) {
                return max_gap2_try_enumeration;
            }
            double mean = success_enumeration_gap.first / success_enumeration_gap.second;
            mean = std::max(mean / EnumerationFailFactor, mean + min_enumeration_exploration_added);
            mean = (mean + success_enumeration_gap.first) / (success_enumeration_gap.second + 1);
            mean = std::max(max_enumeration_success_gap, mean);
            mean = std::min(min_enumeration_fail_gap - TOLERANCE, mean);

            return mean;
        }

        inline void adjustEnumerationStdGap(double ub, double opt_gap, bool if_suc,
                                            double &gap_tolerance4_arc_elimination_n_enumeration,
                                            double &max_gap2_try_enumeration,
                                            std::pair<double, int> &success_enumeration_gap,
                                            double &max_enumeration_success_gap,
                                            double &min_enumeration_fail_gap) {
            auto gap = opt_gap / ub;
            if (if_suc) {
                if (gap > gap_tolerance4_arc_elimination_n_enumeration) {
                    gap_tolerance4_arc_elimination_n_enumeration = gap;
                }

                max_enumeration_success_gap = std::max(max_enumeration_success_gap, gap);

                double mean = 0;
                if (success_enumeration_gap.second) {
                    mean = success_enumeration_gap.first / success_enumeration_gap.second;
                }
                if (gap > mean) {
                    success_enumeration_gap.first = gap;
                    success_enumeration_gap.second = 1;
                } else {
                    success_enumeration_gap.first = getGapStdTryEnumeration(success_enumeration_gap,
                                                                            max_gap2_try_enumeration,
                                                                            max_enumeration_success_gap,
                                                                            min_enumeration_fail_gap) * (
                                                        success_enumeration_gap.second + 1);
                    ++success_enumeration_gap.second;
                }
                max_gap2_try_enumeration = getGapStdTryEnumeration(success_enumeration_gap,
                                                                   max_gap2_try_enumeration,
                                                                   max_enumeration_success_gap,
                                                                   min_enumeration_fail_gap);
            } else {
                min_enumeration_fail_gap = std::min(min_enumeration_fail_gap, gap);
                min_enumeration_fail_gap = std::max(max_enumeration_success_gap + TOLERANCE, min_enumeration_fail_gap);
                if (success_enumeration_gap.second) {
                    double mean = success_enumeration_gap.first / success_enumeration_gap.second;
                    mean *= EnumerationFailFactor; //reduce the mean
                    success_enumeration_gap.first += mean;
                    ++success_enumeration_gap.second; //increase the denominator, keep average the same
                } else {
                    max_gap2_try_enumeration = std::min(gap, max_gap2_try_enumeration * EnumerationFailFactor);
                }
            }
            if (success_enumeration_gap.second)
                std::cout << "average gap of successful enumeration= " << success_enumeration_gap.first /
                        success_enumeration_gap.
                        second
                        << std::endl;

            std::cout << "next try value > "
                    << ub * (1 - getGapStdTryEnumeration(success_enumeration_gap,
                                                         max_gap2_try_enumeration,
                                                         max_enumeration_success_gap,
                                                         min_enumeration_fail_gap)) << std::endl;
        }


        inline double getThresholdN(double m, int a0, double q) {
            if (q <= 1 + TOLERANCE) return m / a0 + TOLERANCE;
            return std::log(1 + m * (q - 1) / a0) / std::log(q) + TOLERANCE;
        }

        inline bool testIfNLinearPass(double max_label, double left_bin, int left_labels) {
            return left_labels > max_label * left_bin;
        }

        template<bool dir>
        void setLabelLimitMiddleCheck(const std::vector<int> &aq,
                                      int num_label,
                                      int max_label_in_enumeration,
                                      double stop_b,
                                      ENUMERATION_STATE &status) {
            if (num_label < NumCheckLabelInEnumeration) return;
            double q = 0.0;
            bool hasPrev = false;
            int prev = 0;
            int max_idx = -1;
            int max_val = -1;

            for (int i = (dir ? 0 : static_cast<int>(aq.size() - 1)); dir ? (i < aq.size()) : (i >= 0);
                 dir ? (++i) : (--i)) {
                int v = aq[i];
                if (v == 0) continue;
                if (hasPrev) {
                    q += v / static_cast<double>(prev);
                }
                prev = v;
                hasPrev = true;
                if (v > max_val) {
                    max_val = v;
                    max_idx = i;
                }
            }


            if (q == 0.0) return;

            if (!testIfNLinearPass(max_val, std::abs(stop_b - max_idx),
                                   max_label_in_enumeration - num_label)) {
                status = dir ? ENUMERATION_STATE::LABEL_FOR_LIMIT : ENUMERATION_STATE::LABEL_BACK_LIMIT;
                return;
            }

            double n = getThresholdN(MaxNumLabelsCheckBinFactorInEnumeration * max_label_in_enumeration, num_label,
                                     q) * (
                           dir ? (max_idx) : (aq.size() - max_idx - 1)) * ToleranceFactorInEnumerationLabelChecking;
            if constexpr (dir) {
                if (n < stop_b) {
                    status = ENUMERATION_STATE::LABEL_FOR_LIMIT;
                }
            } else {
                if (static_cast<double>(aq.size()) - n > stop_b) {
                    status = ENUMERATION_STATE::LABEL_BACK_LIMIT;
                }
            }
        }
    }

    template<bool if_symmetry>
    bool CVRP_Pricing::determineIfEnumeration(double ub, double opt_gap, int num_forward_bucket_arcs,
                                              int num_backward_bucket_arcs) {
        if_force_enumeration_suc = false;
        if (!if_arc_elimination_succeed) return false;

        auto gap = opt_gap / ub;
        if (gap > max_gap2_try_enumeration_enumeration || gap > EnumerationDetail::getGapStdTryEnumeration(
                success_enumeration_gap, max_gap2_try_enumeration_enumeration, max_enumeration_success_gap,
                min_enumeration_fail_gap)) {
            return false;
        }
        if (gap < max_enumeration_success_gap) {
            if_force_enumeration_suc = true;
            return true;
        }

        if (num_forward_bucket_arcs < max_bucket_arc_suc_enumeration.first) {
            if constexpr (!if_symmetry) {
                if (num_backward_bucket_arcs < max_bucket_arc_suc_enumeration.second) {
                    if_force_enumeration_suc = true;
                    return true;
                }
            } else {
                if_force_enumeration_suc = true;
                return true;
            }
        }

        if (num_forward_bucket_arcs > min_bucket_arc_fail_enumeration.first) {
            return false;
        }

        if constexpr (!if_symmetry) {
            if (num_backward_bucket_arcs > min_bucket_arc_fail_enumeration.second)
                return false;
        }
        return true;
    }

    template<bool if_symmetry>
    bool CVRP_Pricing::enumerateMIP(const std::vector<Rcc> &rccs,
                                    const std::vector<R1c> &r1cs,
                                    const std::vector<Brc> &brcs,
                                    const std::vector<double> &optimal_dual_vector,
                                    double ub,
                                    double opt_gap,
                                    int num_forward_bucket_arcs,
                                    int num_backward_bucket_arcs,
                                    bool &if_in_enu_state,
                                    RowVectorXT &index_columns_in_enumeration_column_pool,
                                    RowVectorXd &cost_for_columns_in_enumeration_column_pool,
                                    int &valid_size,
                                    const std::vector<double> &optional_demand_testifier,
                                    double optional_cap_testifier) {
        if (!determineIfEnumeration<if_symmetry>(ub, opt_gap, num_forward_bucket_arcs, num_backward_bucket_arcs))
            return false;

        std::cout << SMALL_PHASE_SEPARATION;
        printHeadLines("Try Enumeration");


        priceConstraints(rccs,
                         r1cs,
                         brcs,
                         optimal_dual_vector);

        this->opt_gap = opt_gap;

        ENUMERATION_STATE status;
        auto if_succeed = enumerateRoutes<if_symmetry>(if_in_enu_state,
                                                       status,
                                                       index_columns_in_enumeration_column_pool,
                                                       cost_for_columns_in_enumeration_column_pool,
                                                       valid_size,
                                                       optional_demand_testifier,
                                                       optional_cap_testifier);

        EnumerationDetail::adjustEnumerationStdBucketArcs<if_symmetry>(max_bucket_arc_suc_enumeration,
                                                                       min_bucket_arc_fail_enumeration,
                                                                       num_forward_bucket_arcs,
                                                                       num_backward_bucket_arcs,
                                                                       if_succeed);

        EnumerationDetail::adjustEnumerationStdGap(ub, opt_gap, if_succeed,
                                                   gap_tolerance4_arc_elimination_n_enumeration,
                                                   max_gap2_try_enumeration_enumeration,
                                                   success_enumeration_gap,
                                                   max_enumeration_success_gap,
                                                   min_enumeration_fail_gap);


        glob_timer.report();

        if (if_succeed) return true;
        return false;
    }

    template<bool if_symmetry>
    bool CVRP_Pricing::enumerateRoutes(
        bool &if_in_enu_state,
        ENUMERATION_STATE &status,
        RowVectorXT &index_columns_in_enumeration_column_pool,
        RowVectorXd &cost_for_columns_in_enumeration_column_pool,
        int &valid_size,
        const std::vector<double> &optional_demand_testifier,
        double &optional_cap_testifier) {
        if_short_memory = false;
        if (if_force_enumeration_suc) {
            max_label_in_enumeration = std::numeric_limits<int>::max();
            max_route_in_enumeration = std::numeric_limits<int>::max();
        } else {
            max_label_in_enumeration = MaxNumLabelInEnumeration;
            max_route_in_enumeration = MaxNumRouteInEnumeration;
        }
        int num_routes_now = 0;
        auto &cost_m = cost_for_columns_in_enumeration_column_pool;
        auto &ptr = index_columns_in_enumeration_column_pool;

        int index;
        std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > Tags;
        Tags.reserve(MaxNumRouteInEnumeration);

        auto copy_Forward_bucket = new std::vector<Label *> *[dim];
        for (int i = 0; i < dim; ++i) {
            copy_Forward_bucket[i] = new std::vector<Label *>[num_buckets_per_vertex];
            for (int b = 0; b < num_buckets_per_vertex; ++b) {
                copy_Forward_bucket[i][b].assign(label_array_in_forward_sense[i][b].begin(),
                                                 label_array_in_forward_sense[i][b].end());
            }
        }

        std::vector<Label *> **copy_Backward_bucket{};
        if constexpr (!if_symmetry) {
            copy_Backward_bucket = new std::vector<Label *> *[dim];
            for (int i = 0; i < dim; ++i) {
                copy_Backward_bucket[i] = new std::vector<Label *>[num_buckets_per_vertex];
                for (int b = 0; b < num_buckets_per_vertex; ++b) {
                    copy_Backward_bucket[i][b].assign(label_array_in_backward_sense[i][b].begin(),
                                                      label_array_in_backward_sense[i][b].end());
                }
            }
        }


        auto eps = TimeSetter::measure([&]() {
            if constexpr (!if_symmetry) {
                status = enumerateHalfwardRoutes<true, false>(Tags, copy_Backward_bucket, num_routes_now);
            } else {
                status = enumerateHalfwardRoutes<true, true>(Tags, copy_Forward_bucket, num_routes_now);
            }
        });
        printTimeMessage("labeling 1", eps);

        if (status != ENUMERATION_STATE::NORMAL) {
            return false;
        }

        if constexpr (!if_symmetry) {
            eps = TimeSetter::measure([&]() {
                status = enumerateHalfwardRoutes<false, false>(
                    Tags,
                    copy_Forward_bucket,
                    num_routes_now);
            });
            printTimeMessage("labeling 2", eps);

            if (status != ENUMERATION_STATE::NORMAL) {
                return false;
            }
        }

        eps = TimeSetter::measure([&]() {
            status = concatenateRoutesPriorForwardInEnumeration<if_symmetry>(Tags, num_routes_now);
        });
        printTimeMessage("concatenate", eps);

        if (status != ENUMERATION_STATE::NORMAL) {
            return false;
        }


        cost_m.resize(num_routes_now);
        ptr.resize(num_routes_now);


        auto num = size_t(num_routes_now * getAverageRouteLength());
        reallocatePricingPool(num);
        auto PricingWarning = static_cast<size_t>(PricingWarningThreshold * static_cast<double>(mem4_pricing));
        pool_beg4_pricing = 0;
        index = 0;
        Label *ki, *li, *p;
        std::vector<int> seq(dim + 1);
        for (auto &tag: Tags) {
            ki = std::get<0>(tag.second);
            li = std::get<1>(tag.second);

            cost_m[index] = std::get<2>(tag.second);
            ptr[index++] = pool_beg4_pricing;

            int cnt = 0;
            p = ki;
            while (p) {
                seq[cnt++] = p->end_vertex;
                p = p->p_label;
            }
            for (int k = 0; k < cnt; ++k) {
                col_pool4_pricing[pool_beg4_pricing++] = seq[cnt - 1 - k];
            }
            if (li) {
                p = li;
                while (p) {
                    col_pool4_pricing[pool_beg4_pricing++] = p->end_vertex;
                    p = p->p_label;
                }
            } else {
                col_pool4_pricing[pool_beg4_pricing++] = 0;
            }
            if (!optional_demand_testifier.empty()) {
                double local_cap = 0;
                for (size_t i = ptr[index - 1] + 1; i < pool_beg4_pricing; ++i) {
                    local_cap += optional_demand_testifier[col_pool4_pricing[i]];
                }
                if (local_cap > optional_cap_testifier + TOLERANCE) {
                    --index;
                    pool_beg4_pricing = ptr[index];
                    continue;
                }
            }
            resizePoolWarning(PricingWarning);
        }

        num_routes_now = index;
        cost_m.conservativeResize(num_routes_now);
        ptr.conservativeResize(num_routes_now);
        valid_size = num_routes_now;
        if_in_enu_state = true;
        std::cout << "#routes= " << num_routes_now << std::endl;
        return true;
    }

    template<bool dir, bool if_symmetry>
    ENUMERATION_STATE CVRP_Pricing::enumerateHalfwardRoutes(
        std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &Tags,
        std::vector<Label *> **copy_bucket,
        int &num_routes_now) {
        ENUMERATION_STATE status = ENUMERATION_STATE::NORMAL;
        (dir ? num_forward_labels_in_enu : num_backward_labels_in_enu) = 0;

        initializeLabels<dir, if_symmetry, false, true, true>();


        step_size_label_check_in_enumeration = resource.resources[0] / LabelsCheckBinInEnumeration;
        label_per_bin_in_enumeration.assign(
            static_cast<int>(resource.resources[0] / (step_size_label_check_in_enumeration - 1)) + 1,
            0);
        double stop_b = meet_point_resource_in_bi_dir / step_size_label_check_in_enumeration;

        for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
            for (auto &comp: dir ? topological_order_forward_ptr->at(b) : topological_order_backward_ptr->at(b)) {
                int index = 0;
            STILL_EXIST:
                for (; index < comp.size(); ++index) {
                    int i = comp[index];
                    auto &valid_num = (dir
                                           ? if_exist_extra_labels_in_forward_sense[i][b].second
                                           : if_exist_extra_labels_in_backward_sense[i][b].second);
                    if (!valid_num) continue;
                    auto &label_array = (dir
                                             ? if_exist_extra_labels_in_forward_sense[i][b].first
                                             : if_exist_extra_labels_in_backward_sense[i][b].first);
                    for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
                        auto &ki = label_array[vec_index];
                        if (ki->is_extended) continue;
                        ki->is_extended = true;
                        extendKernel4Enumeration<dir, if_symmetry>(
                            i, b, ki, copy_bucket, Tags, num_routes_now,
                            status);
                        if (status != ENUMERATION_STATE::NORMAL) goto outside;
                    }
                    EnumerationDetail::setLabelLimitMiddleCheck<dir>(label_per_bin_in_enumeration,
                                                                     dir
                                                                         ? num_forward_labels_in_enu
                                                                         : num_backward_labels_in_enu,
                                                                     max_label_in_enumeration,
                                                                     stop_b,
                                                                     status);
                    if (status != ENUMERATION_STATE::NORMAL) goto outside;
                    valid_num = 0;
                }
                for (index = 0; index < comp.size(); ++index) {
                    int i = comp[index];
                    if (dir
                            ? if_exist_extra_labels_in_forward_sense[i][b].second
                            : if_exist_extra_labels_in_backward_sense[i][b].second) {
                        goto STILL_EXIST;
                    }
                }
            }
        }
    outside:
        for (int i = 0; i < dim; ++i) {
            delete[]copy_bucket[i];
        }
        delete[] copy_bucket;
        std::cout << "#labels= " << static_cast<int>(dir ? num_forward_labels_in_enu : num_backward_labels_in_enu) <<
                " #routes= " <<
                num_routes_now << std::endl;
        return status;
    }

    template<bool dir, bool if_symmetry>
    void CVRP_Pricing::extendKernel4Enumeration(int i, int b, Label *ki,
                                                std::vector<Label *> **copy_bucket,
                                                std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> >
                                                &
                                                Tags,
                                                int &num_routes_now, ENUMERATION_STATE &status) {
        status = ENUMERATION_STATE::NORMAL;
        constexpr bool if_dif = dir ^ if_symmetry;
        for (int j: (dir
                         ? all_forward_buckets[i][b].bucket_arcs
                         : all_backward_buckets[i][b].bucket_arcs)) {
            if (ki->pi[j]) continue;
            auto new_label = all_label + idx_glo;
            auto &tmp_Resource = new_label->res;
            if (dir
                    ? !increaseMainResourceConsumption(ki->res, tmp_Resource, i, j)
                    : !decreaseMainResourceConsumption(ki->res, tmp_Resource, i, j))
                continue;
            auto &tmp_rc = new_label->rc;
            tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j]; //real rc
            auto if_keep = false;
            auto arr_bj = (if_symmetry
                               ? static_cast<int>((resource.resources[0] - tmp_Resource.resources[0]) / step_size)
                               : static_cast<int>(tmp_Resource.resources[0] / step_size));
            if constexpr (dir) {
                if (tmp_Resource.resources[0] > meet_point_resource_in_bi_dir) {
                    if constexpr (if_symmetry) {
                        if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap)
                            concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_Resource);
                    } else {
                        if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap)
                            concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_Resource);
                    }
                    continue;
                }
            } else {
                if (tmp_Resource.resources[0] < meet_point_resource_in_bi_dir) {
                    if constexpr (if_symmetry) {
                        if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap)
                            concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_Resource);
                    } else {
                        if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap)
                            concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_Resource);
                    }
                    continue;
                }
            }

            for (; if_dif ? arr_bj < num_buckets_per_vertex : arr_bj >= 0;
                   if_dif ? ++arr_bj : --arr_bj) {
                if (tmp_rc
                    + (if_dif
                           ? rc2_till_this_bin_in_backward_sense[j][arr_bj]
                           : rc2_till_this_bin_in_forward_sense[j][arr_bj])
                    > opt_gap)
                    break;
                if (tmp_rc + (if_dif ? rc2_bin_in_backward_sense[j][arr_bj] : rc2_bin_in_forward_sense[j][arr_bj])
                    > opt_gap)
                    continue;
                for (auto &kkj: copy_bucket[j][arr_bj]) {
                    double path_rc = tmp_rc + kkj->rc;

                    if (path_rc > opt_gap) break;

                    if constexpr (if_symmetry) {
                        if (tellResTupleRelations<'c'>(tmp_Resource, kkj->res)) continue;
                    } else {
                        if constexpr (dir) {
                            if (tellResTupleRelations<'p'>(tmp_Resource, kkj->res)) continue;
                        } else {
                            if (tellResTupleRelations<'p'>(kkj->res, tmp_Resource)) continue;
                        }
                    }

                    if ((ki->pi & kkj->pi).any()) continue;


                    if (!rank1_rc_controller_ref.get().concatenateR1CStates(path_rc, opt_gap, ki->r1c, kkj->r1c,
                                                                            i, j))
                        continue;


                    if_keep = true;
                    goto OUT;
                }
            }

        OUT:
            if (!if_keep) continue;

            int bj;
            updateEnumerationLabel(ki, i, j, bj);

            doDominanceEnumerationLabel<dir>(j, bj, if_keep);


            if (!if_keep) continue;


            rank1_rc_controller_ref.get().updateR1CStates(new_label->rc, new_label->r1c, ki->r1c, i, j);

            new_label->p_label = ki;
            new_label->is_extended = false;
            auto &bucket =
            (dir
                 ? if_exist_extra_labels_in_forward_sense[j][bj]
                 : if_exist_extra_labels_in_backward_sense[j][bj]);
            bucket.first[bucket.second++] = new_label;
            if (bucket.second == bucket.first.size()) {
                bucket.first.resize(bucket.first.size() * 2);
            }

            if constexpr (dir) {
                if (tmp_rc + chg_cost_mat4_vertex[j][0] < opt_gap) {
                    auto path_cost = new_label->cost + cost_mat4_vertex_ref.get()[j][0];
                    auto &tmp_PI = new_label->pi;
                    if (Tags.find(tmp_PI) == Tags.end()) {
                        Tags[tmp_PI] = {all_label + idx_glo, nullptr, path_cost};
                        ++num_routes_now;
                    } else if (std::get<2>(Tags[tmp_PI]) > path_cost) {
                        Tags[tmp_PI] = {all_label + idx_glo, nullptr, path_cost};
                    }
                }
            }

            if ((dir ? ++num_forward_labels_in_enu : ++num_backward_labels_in_enu) > max_label_in_enumeration) {
                status = dir ? ENUMERATION_STATE::LABEL_FOR_LIMIT : ENUMERATION_STATE::LABEL_BACK_LIMIT;
                goto QUIT;
            }
            ++label_per_bin_in_enumeration[static_cast<int>(
                tmp_Resource.resources[0] / step_size_label_check_in_enumeration)];
            ++idx_glo; //can be put here, because once go outside, the function will end
            if (idx_glo == label_assign) {
                status = ENUMERATION_STATE::OUT_OF_MEMORY;
                if_short_memory = true;
                goto QUIT;
            }
        }
    QUIT:;
    }


    inline void CVRP_Pricing::updateEnumerationLabel(Label *ki, int i, int j, int &bj) {
        auto new_label = all_label + idx_glo;
        auto &tmp_Resource = new_label->res;
        bj = static_cast<int>(tmp_Resource.resources[0] / step_size);
        auto &tmp_PI = new_label->pi;
        auto &tmp_Cost = new_label->cost;
        tmp_PI = ki->pi;
        tmp_PI.set(j);
        tmp_Cost = ki->cost + cost_mat4_vertex_ref.get()[i][j];
        new_label->end_vertex = j;
    }

    template<bool dir>
    void CVRP_Pricing::doDominanceEnumerationLabel(int j, int bj, bool &if_suc) {
        auto new_label = all_label + idx_glo;
        auto &labelList_j = dir ? label_array_in_forward_sense[j][bj] : label_array_in_backward_sense[j][bj];

        if_suc = true;
        auto it = labelList_j.begin();

        for (; it != labelList_j.end(); ++it) {
            auto kj = *it;
            if (kj->res.resources[MarkerSameNodeSetResourceIdx] != new_label->res.resources[
                    MarkerSameNodeSetResourceIdx])
                continue;
            if ((kj->pi ^ new_label->pi).none()) {
                auto state = dominanceCoreInEnumeration<dir>(new_label, kj);
                if (state == dominanceCoreInEnumeration_STATE::KI_DOMINATE_KJ) {
                    kj->is_extended = true;
                    *it = new_label;
                    dir ? --num_forward_labels_in_enu : --num_backward_labels_in_enu;
                    --label_per_bin_in_enumeration[kj->res.resources[0] / step_size_label_check_in_enumeration];
                    break;
                }
                if (state == dominanceCoreInEnumeration_STATE::KJ_DOMINATE_KI) {
                    if_suc = false;
                    return;
                }
            }
        }

        if (it == labelList_j.end()) {
            labelList_j.push_back(new_label); //here is push back
        }
    }

    template<bool dir>
    dominanceCoreInEnumeration_STATE CVRP_Pricing::dominanceCoreInEnumeration(Label *ki, Label *kj) {
        if (ki->cost < kj->cost) {
            if (!(dir ? tellResTupleRelations<'p'>(ki->res, kj->res) : tellResTupleRelations<'p'>(kj->res, ki->res)))
                return dominanceCoreInEnumeration_STATE::KI_DOMINATE_KJ;
        } else {
            if (!(dir ? tellResTupleRelations<'p'>(kj->res, ki->res) : tellResTupleRelations<'p'>(ki->res, kj->res)))
                return dominanceCoreInEnumeration_STATE::KJ_DOMINATE_KI;
        } //do not consider equal in this case, unlike rc
        return dominanceCoreInEnumeration_STATE::NO_DOMINANCE;
    }

    template<bool if_symmetry>
    ENUMERATION_STATE CVRP_Pricing::concatenateRoutesPriorForwardInEnumeration(
        std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &Tags,
        int &num_routes_now) {
        ENUMERATION_STATE status = ENUMERATION_STATE::NORMAL;

        if constexpr (!if_symmetry) populateRC2TillThisBinNRC2Bin<false>();
        else populateRC2TillThisBinNRC2Bin<true>();

        for (auto &label_list: concatenate_labels_in_forward_cg) {
            int i = label_list.first.first;
            int j = label_list.first.second;
            auto &label_vec = label_list.second;
            for (auto &pr: label_vec) {
                auto &ki = pr.first;
                auto &tmp_Resource = pr.second;
                double tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
                if constexpr (!if_symmetry) {
                    auto arr_bj = static_cast<int>((tmp_Resource.resources[0]) / step_size);
                    if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap)
                        continue;
                    for (; arr_bj < num_buckets_per_vertex; ++arr_bj) {
                        if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap)
                            break;

                        if (rc2_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap ||
                            (std::find(all_backward_buckets[j][arr_bj].bucket_arcs.begin(),
                                       all_backward_buckets[j][arr_bj].bucket_arcs.end(), i)
                             == all_backward_buckets[j][arr_bj].bucket_arcs.end()))
                            continue;

                        auto &label_arr = label_array_in_backward_sense[j][arr_bj];
                        checkGroupInner<if_symmetry>(label_arr, tmp_rc, tmp_Resource, ki, i, j, Tags,
                                                     num_routes_now,
                                                     status);
                        if (status != ENUMERATION_STATE::NORMAL) goto QUIT;
                    }
                } else {
                    auto arr_bj = static_cast<int>((resource.resources[0] - tmp_Resource.resources[0]) / step_size);
                    if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
                        continue;

                    for (; arr_bj >= 0; --arr_bj) {
                        if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
                            break;

                        if (rc2_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
                            continue;

                        auto &label_arr = label_array_in_forward_sense[j][arr_bj];
                        checkGroupInner<if_symmetry>(label_arr, tmp_rc, tmp_Resource, ki, i, j, Tags,
                                                     num_routes_now,
                                                     status);
                        if (status != ENUMERATION_STATE::NORMAL) goto QUIT;
                    }
                }
            }
        }
    QUIT:
        return status;
    }

    template<bool if_symmetry>
    void CVRP_Pricing::checkGroupInner(
        ListLabel &label_arr,
        double &tmp_rc, Resource &tmp_Resource,
        Label *ki, int i, int j,
        std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &Tags,
        int &num_routes_now, ENUMERATION_STATE &status
    ) {
        double path_cost;
        routeOptLong tmp_PI;
        for (auto &kj: label_arr) {
            auto path_rc = kj->rc + tmp_rc;
            if (path_rc > opt_gap) break;

            if constexpr (!if_symmetry) {
                if (tellResTupleRelations<'p'>(tmp_Resource, kj->res))continue;
            } else {
                if (tellResTupleRelations<'c'>(tmp_Resource, kj->res))continue;
            }

            if ((ki->pi & kj->pi).any()) continue;

            if (!rank1_rc_controller_ref.get().concatenateR1CStates(path_rc, opt_gap, ki->r1c, kj->r1c,
                                                                    i, j))
                goto HERE;

            path_cost = ki->cost + cost_mat4_vertex_ref.get()[i][j] + kj->cost;
            tmp_PI = ki->pi | kj->pi;
            if (Tags.find(tmp_PI) == Tags.end()) {
                Tags[tmp_PI] = std::make_tuple(ki, kj, path_cost);
                ++num_routes_now;
                if (num_routes_now > max_route_in_enumeration) {
                    status = ENUMERATION_STATE::ROUTE_LIMIT;
                    goto QUIT;
                }
            } else if (std::get<2>(Tags[tmp_PI]) > path_cost) {
                Tags[tmp_PI] = std::make_tuple(ki, kj, path_cost);
            }
        HERE:;
        }
    QUIT:;
    }
}

#endif // ROUTE_OPT_ENUMERATION_HPP
