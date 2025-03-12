/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_ARC_ELIMINATION_HPP
#define ROUTE_OPT_ARC_ELIMINATION_HPP
#include <enumeration.hpp>

namespace RouteOpt::Application::CVRP {
    template<bool dir, bool if_symmetry>
    void CVRP_Pricing::concatenatePhaseInArcElimination() {
        int bj;
        bool if_suc;


        auto beg = std::chrono::high_resolution_clock::now();

        for (auto &label_list: dir ? concatenate_labels_in_forward_cg : concatenate_labels_in_backward_cg) {
            int i = label_list.first.first;
            int j = label_list.first.second;
            auto &label_vec = label_list.second;
            for (auto &pr: label_vec) {
                auto &ki = pr.first;
                auto &tmp_Resource = all_label[idx_glo].res;
                tmp_Resource = pr.second; //do not change here! since res will not be updated later
                updateLabel<dir, true, false, if_symmetry, true, true, PRICING_LEVEL::EXACT>(
                    tmp_Resource, ki, i, j, bj, if_suc);
                if (!if_suc) continue;

                doDominance<dir, PRICING_LEVEL::EXACT>(ki, j, bj, if_suc);
                if (!if_suc) continue;


                ++idx_glo;
                if (idx_glo == label_assign) {
                    if_short_memory = true;
                    if_arc_elimination_succeed = false;
                    goto QUIT;
                }
            }
            if (
                std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - beg).count() >
                HardTimeThresholdInArcEliminationMidConcatenate) {
                if_arc_elimination_succeed = false;
            }
        }
    QUIT:;
    }

    template<bool if_symmetry>
    void CVRP_Pricing::runLabelingForArcElimination() {
        if_arc_elimination_succeed = true;
        if_short_memory = false;

        if constexpr (!if_symmetry) {
            concatenatePhaseInArcElimination<true, false>();
            if (!if_arc_elimination_succeed) goto QUIT;
            runLabeling<true, true, false, false, PRICING_LEVEL::EXACT>(arc_elimination_time);
            if (!if_exact_labeling_finished) goto QUIT;
            concatenatePhaseInArcElimination<false, false>();
            if (!if_arc_elimination_succeed) goto QUIT;
            runLabeling<false, true, false, false, PRICING_LEVEL::EXACT>(arc_elimination_time);
            if (!if_exact_labeling_finished) goto QUIT;
        } else {
            concatenatePhaseInArcElimination<true, true>();
            if (!if_arc_elimination_succeed) goto QUIT;
            runLabeling<true, true, false, true, PRICING_LEVEL::EXACT>(arc_elimination_time);
            if (!if_exact_labeling_finished) goto QUIT;
        }
    QUIT:
        if (!if_exact_labeling_finished || if_short_memory) if_arc_elimination_succeed = false;
        if (if_arc_elimination_succeed) {
            gap_improved_4_arc_elimination_n_enumeration += GapBarIncreased4ArcEliminationNEnumeration;
            arc_elimination_time += ArcEliminationTimeIncreased;
        }
    }

    inline bool CVRP_Pricing::determineIfArcElimination(double ub, double opt_gap, double &last_gap) {
        if_arc_elimination_succeed = false;
        if (!if_exact_labeling_finished || if_stop_arc_elimination) {
            return false;
        }

        auto now_gap = opt_gap / ub;
        if (!equalFloat(old_ub, ub)) {
            goto QUIT;
        }

        if (last_gap / now_gap > gap_improved_4_arc_elimination_n_enumeration
            || now_gap < gap_tolerance4_arc_elimination_n_enumeration) {
            goto QUIT;
        }
        return false;
    QUIT:
        old_ub = ub;
        last_gap = now_gap;
        return true;
    }

    template<bool if_symmetry>
    void CVRP_Pricing::eliminateArcs(
        const std::vector<Rcc> &rccs,
        const std::vector<R1c> &r1cs,
        const std::vector<Brc> &brcs,
        const std::vector<double> &optimal_dual_vector,
        double ub,
        double opt_gap,
        double &last_gap,
        int &num_forward_bucket_arcs,
        int &num_backward_bucket_arcs,
        int &num_forward_jump_arcs,
        int &num_backward_jump_arcs
    ) {
        if (!determineIfArcElimination(ub, opt_gap, last_gap)) return;

        std::cout << SMALL_PHASE_SEPARATION;
        printHeadLines("Run Arc Elimination");

        populateTellWhichBin4ArcElimination<true>();

        if constexpr (!if_symmetry) {
            populateTellWhichBin4ArcElimination<false>();
        }

        priceConstraints(rccs,
                         r1cs,
                         brcs,
                         optimal_dual_vector);

        this->opt_gap = opt_gap;


        auto eps = TimeSetter::measure([&]() { runLabelingForArcElimination<if_symmetry>(); });
        printTimeMessage("arc elimination", eps);

        if (!if_arc_elimination_succeed) goto QUIT;

        eps = TimeSetter::measure([&]() {
            eliminateBucketArcs<if_symmetry>(num_forward_bucket_arcs, num_backward_bucket_arcs);
        });
        printTimeMessage("eliminate bucket arcs", eps);


        eps = TimeSetter::measure([&]() {
            obtainJumpArcs<if_symmetry>(num_forward_jump_arcs, num_backward_jump_arcs);
        });
        printTimeMessage("obtain jump arcs", eps);
        if_arc_elimination_succeed = true;
    QUIT:
        std::cout << SMALL_PHASE_SEPARATION;
        glob_timer.report();
    }

    template<typename T, bool dir, bool if_symmetry>
    void CVRP_Pricing::concatenateTestKernelInArcElimination(int i,
                                                             int b,
                                                             const std::vector<T> &arc,
                                                             int dim_sq,
                                                             bool *stateBetween2Buckets,
                                                             int *latest_bucket) {
        Resource tmp_Resource;
        double tmp_rc;
        int if_state, arr_bj;
        for (auto &pr: arc) {
            int j;
            if constexpr (std::is_same<T, int>::value) {
                j = pr;
            } else {
                j = pr.second;
                Resource res = dir ? Resource{} : resource;
                res.resources[0] = pr.first;
                if (dir
                        ? !increaseMainResourceConsumption(res, tmp_Resource, i, j)
                        : !decreaseMainResourceConsumption(res, tmp_Resource, i, j))
                    continue;
            }
            int map1 = i * dim + j;
            int &latest = latest_bucket[map1];
            int old_latest = latest;
            auto &label_array = dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b];
            constexpr bool if_dif = dir ^ if_symmetry;
            for (auto ki: label_array) {
                if constexpr (std::is_same<T, int>::value) {
                    if (dir
                            ? !increaseMainResourceConsumption(ki->res, tmp_Resource, i, j)
                            : !decreaseMainResourceConsumption(ki->res, tmp_Resource, i, j))
                        continue;
                }
                tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
                if constexpr (dir) {
                    arr_bj =
                            if_symmetry
                                ? std::min(
                                    static_cast<int>((resource.resources[0] - tmp_Resource.resources[0]) / step_size),
                                    latest - 1)
                                : std::max(static_cast<int>((tmp_Resource.resources[0]) / step_size), latest + 1);
                } else {
                    arr_bj =
                            if_symmetry
                                ? std::max(
                                    static_cast<int>((resource.resources[0] - tmp_Resource.resources[0]) / step_size),
                                    latest + 1)
                                : std::min(static_cast<int>((tmp_Resource.resources[0]) / step_size), latest - 1);
                }

                for (int bj = if_dif ? num_buckets_per_vertex - 1 : 0; if_dif ? bj >= arr_bj : bj <= arr_bj;
                     if_dif ? --bj : ++bj) {
                    concatenateOneLabelWithOtherLabels<dir, if_symmetry, true>(ki,
                                                                               j,
                                                                               bj,
                                                                               tmp_rc,
                                                                               tmp_Resource,
                                                                               if_state);
                    if (if_state >= 0) {
                        latest = if_state;
                        goto outside;
                    }
                }
            outside:;
            }
            if (latest != old_latest) {
                int bi, chg_latest;
                if constexpr (if_symmetry) {
                    chg_latest = static_cast<int>((resource.resources[0] - latest * step_size) / step_size);
                    if (chg_latest < 0) chg_latest = 0;
                    else if (chg_latest >= num_buckets_per_vertex) chg_latest = num_buckets_per_vertex - 1;
                } else chg_latest = latest;
                bi = (dir
                          ? tell_which_bin4_arc_elimination_in_forward_sense[map1 + chg_latest * dim_sq]
                          : tell_which_bin4_arc_elimination_in_backward_sense[map1 + chg_latest * dim_sq]);
                int map2 = map1 + b * dim_sq;
                for (int k = b; dir ? k <= bi : k >= bi; dir ? ++k : --k, dir ? map2 += dim_sq : map2 -= dim_sq)
                    stateBetween2Buckets[map2] = true;
            }
        }
    }

    template<bool dir, bool if_symmetry>
    void CVRP_Pricing::eliminateBucketArcs(
        int dim_sq,
        bool *stateBetween2Buckets,
        int *latest_bucket,
        int &num_forward_bucket_arcs,
        int &num_backward_bucket_arcs) {
        memset(stateBetween2Buckets, 0, dim_sq * num_buckets_per_vertex * sizeof(bool));
        constexpr bool if_dif = dir ^ if_symmetry;
        if constexpr (if_dif) {
            memset(latest_bucket, -1, dim_sq * sizeof(int));
        } else {
            std::fill_n(latest_bucket, dim_sq, num_buckets_per_vertex);
        }
        int num_bucket_arcs = 0;
        for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
            for (int i = 1; i < dim; ++i) {
                concatenateTestKernelInArcElimination<std::pair<res_int, int>, dir, if_symmetry>(i,
                    b,
                    (dir
                         ? all_forward_buckets[i][b].jump_arcs
                         : all_backward_buckets[i][b].jump_arcs),
                    dim_sq,
                    stateBetween2Buckets,
                    latest_bucket);
                concatenateTestKernelInArcElimination<int, dir, if_symmetry>(i,
                                                                             b,
                                                                             (dir
                                                                                  ? all_forward_buckets[i][b].
                                                                                  bucket_arcs
                                                                                  : all_backward_buckets[i][b].
                                                                                  bucket_arcs),
                                                                             dim_sq,
                                                                             stateBetween2Buckets,
                                                                             latest_bucket);
            }
        }

        std::vector<int> tmp_vec;
        tmp_vec.reserve(dim);
        int map1 = 0;
        for (int b = 0; b < num_buckets_per_vertex; ++b, map1 += dim_sq) {
            int map2 = map1 + dim;
            for (int i = 1; i < dim; ++i, map2 += dim) {
                tmp_vec.clear();
                for (int j: dir
                                ? all_forward_buckets[i][b].bucket_arcs
                                : all_backward_buckets[i][b].bucket_arcs) {
                    if (stateBetween2Buckets[map2 + j]) {
                        tmp_vec.emplace_back(j);
                        ++num_bucket_arcs;
                    }
                }
                (dir ? all_forward_buckets[i][b].bucket_arcs : all_backward_buckets[i][b].bucket_arcs) =
                        tmp_vec;
            }
        }

        eliminateBuketArc4Depot<dir, if_symmetry>();
        double ratio;
        if (dir) {
            ratio = static_cast<double>(num_bucket_arcs) / static_cast<double>(max_num_forward_graph_arc);
        } else {
            ratio = static_cast<double>(num_bucket_arcs) / static_cast<double>(max_num_backward_graph_arc);
        }
        std::cout << "previous= "
                << static_cast<double>(num_bucket_arcs) / (dir ? num_forward_bucket_arcs : num_backward_bucket_arcs)
                * 100
                << "% maximum= "
                << ratio * 100 << "%"
                << std::endl;
        (dir ? num_forward_bucket_arcs : num_backward_bucket_arcs) = num_bucket_arcs;
        getTopologicalOrder<if_symmetry>();
    }

    template<bool dir, bool if_symmetry>
    void CVRP_Pricing::eliminateBuketArc4Depot() {
        bool constexpr if_dif = dir ^ if_symmetry;
        auto &allBuckets = dir ? all_forward_buckets : all_backward_buckets;
        auto &RC2TillThisBin = if_dif ? rc2_till_this_bin_in_backward_sense : rc2_till_this_bin_in_forward_sense;
        auto &RC2Bin = if_dif ? rc2_bin_in_backward_sense : rc2_bin_in_forward_sense;
        auto &labelArray = if_dif ? label_array_in_backward_sense : label_array_in_forward_sense;

        std::unordered_set<int> depot_set;
        for (auto &arc: allBuckets[0][0].bucket_arcs) depot_set.emplace(arc);

        for (int i = 1; i < dim; ++i) {
            if (depot_set.find(i) == depot_set.end()) continue;
            bool if_delete = true;

            for (int b = (if_dif ? 0 : (num_buckets_per_vertex - 1)); if_dif ? (b < num_buckets_per_vertex) : (b >= 0);
                 if_dif ? ++b : --b) {
                if (RC2TillThisBin[i][b] + chg_cost_mat4_vertex[i][0] > opt_gap) break;
                if (RC2Bin[i][b] + chg_cost_mat4_vertex[i][0] > opt_gap) continue;

                if (!labelArray[i][b].empty()) {
                    auto &labels = labelArray[i][b];
                    for (auto &label: labels) {
                        if (label->rc + chg_cost_mat4_vertex[i][0] > opt_gap) break;
                        auto tmp_res = label->res; //cannot use &
                        if (if_dif
                                ? (decreaseMainResourceConsumption(tmp_res, tmp_res, i, 0))
                                : (
                                    increaseMainResourceConsumption(tmp_res, tmp_res, i, 0))) {
                            if_delete = false;
                            break;
                        }
                    }
                }
            }

            if (if_delete) {
                depot_set.erase(i);
            }
        }

        if (depot_set.size() != allBuckets[0][0].bucket_arcs.size()) {
            allBuckets[0][0].bucket_arcs.clear();
            allBuckets[0][0].bucket_arcs.assign(depot_set.begin(), depot_set.end());
            std::sort(allBuckets[0][0].bucket_arcs.begin(), allBuckets[0][0].bucket_arcs.end());
        }
    }

    template<bool if_symmetry>
    void CVRP_Pricing::eliminateBucketArcs(int &num_forward_bucket_arcs, int &num_backward_bucket_arcs) {
        int dim_sq = dim * dim;
        auto stateBetween2Buckets = new bool[dim_sq * num_buckets_per_vertex];
        auto latest_bucket = new int[dim_sq];

        if constexpr (!if_symmetry) {
            eliminateBucketArcs<true, false>(dim_sq, stateBetween2Buckets, latest_bucket, num_forward_bucket_arcs,
                                             num_backward_bucket_arcs);
            eliminateBucketArcs<false, false>(dim_sq, stateBetween2Buckets, latest_bucket, num_forward_bucket_arcs,
                                              num_backward_bucket_arcs);
        } else {
            eliminateBucketArcs<true, true>(dim_sq, stateBetween2Buckets, latest_bucket, num_forward_bucket_arcs,
                                            num_backward_bucket_arcs);
        }

        delete[]stateBetween2Buckets;
        delete[]latest_bucket;
    }

    template<bool if_symmetry>
    void CVRP_Pricing::obtainJumpArcs(int &num_forward_jump_arcs, int &num_backward_jump_arcs) {
        auto tmp = new std::bitset<2> *[dim];
        for (int i = 0; i < dim; ++i) tmp[i] = new std::bitset<2>[num_buckets_per_vertex];
        obtainJumpArcs<true>(tmp, num_forward_jump_arcs, num_backward_jump_arcs);
        if constexpr (!if_symmetry) {
            obtainJumpArcs<false>(tmp, num_forward_jump_arcs, num_backward_jump_arcs);
        }
        for (int i = 0; i < dim; ++i) {
            delete[]tmp[i];
        }
        delete[]tmp;
    }

    template<bool dir>
    void CVRP_Pricing::obtainJumpArcs(std::bitset<2> **bitMap,
                                      int &num_forward_jump_arcs,
                                      int &num_backward_jump_arcs) {
        int num_jump_arcs = 0;
        bool if_used;

        for (int i = 1; i < dim; ++i) {
            for (int b = 0; b < num_buckets_per_vertex; ++b) {
                if constexpr (dir) {
                    all_forward_buckets[i][b].jump_arcs.clear();
                } else {
                    all_backward_buckets[i][b].jump_arcs.clear();
                }
                for (int j = 1; j < dim; ++j) bitMap[j][b] = 2;
                bitMap[i][b] = 1;
                for (int j: (dir
                                 ? all_forward_buckets[i][b].bucket_arcs
                                 : all_backward_buckets[i][b].bucket_arcs))
                    bitMap[j][b] = 0;
            }
            for (int b = (dir ? 0 : num_buckets_per_vertex - 1); dir ? b < num_buckets_per_vertex : b >= 0;
                 dir ? ++b : --b) {
                for (int j = 1; j < dim; ++j) {
                    if (bitMap[j][b] == 2) {
                        //need a jump arc
                        if_used = false;
                        for (int b4_i = (dir ? b + 1 : b - 1); dir ? b4_i < num_buckets_per_vertex : b4_i >= 0;
                             dir ? ++b4_i : --b4_i) {
                            if (bitMap[j][b4_i] == 0) {
                                //std::find a going arc
                                std::pair<res_int, int> map1 = (dir
                                                                    ? std::make_pair(b4_i * step_size, j)
                                                                    : std::make_pair((b4_i + 1) * step_size, j));
                                for (int tmp_b = b; dir ? tmp_b < b4_i : tmp_b > b4_i; dir ? ++tmp_b : --tmp_b) {
                                    bitMap[j][tmp_b] = 1;
                                    if constexpr (dir) {
                                        all_forward_buckets[i][tmp_b].jump_arcs.emplace_back(map1);
                                    } else {
                                        all_backward_buckets[i][tmp_b].jump_arcs.emplace_back(map1);
                                    }
                                    ++num_jump_arcs;
                                }
                                if_used = true;
                                break;
                            }
                        }
                        if (!if_used) {
                            for (int tmp_b = b; dir ? tmp_b < num_buckets_per_vertex : tmp_b >= 0;
                                 dir ? ++tmp_b : --tmp_b)
                                bitMap[j][tmp_b] = 1;
                        }
                    }
                }
            }
        }

        (dir ? num_forward_jump_arcs : num_backward_jump_arcs) = num_jump_arcs;
    }
}

#endif // ROUTE_OPT_ARC_ELIMINATION_HPP
