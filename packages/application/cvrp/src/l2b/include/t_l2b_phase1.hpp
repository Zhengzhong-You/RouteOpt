/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_PHASE1_HPP
#define ROUTE_OPT_T_L2B_PHASE1_HPP
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "route_opt_macro.hpp"
#include "l2b_macro.hpp"
#include "l2b_controller.hpp"

namespace RouteOpt::Application::CVRP {
    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::getFeatureDataPhase1(
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Node *node,
        int tree_level,
        const std::vector<double> &lp_solution,
        const std::vector<int> &route_length) {
        if (branching_data_shared.getBranchPair().size() == 1) {
            PRINT_WARNING("no training data phase1: a single candidate");
            return;
        }

        collectEdgeRelatedFeatures(branching_data_shared, tree_level);

        collectStaticFeatures(branching_data_shared, lp_solution, route_length);

        for (const auto &candidate: branching_data_shared.getBranchPair()) {
            const auto &nonzero_idx = getBrConstraintNonzeroIdx(node, candidate);
            collectVariableRelatedFeatures(branching_data_shared, branching_history, candidate, nonzero_idx);
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::collectStaticFeatures(
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        const std::vector<double> &lp_solution,
        const std::vector<int> &route_length) {
        is_in_solution = lp_solution;

        double sum_length = 0;
        int cnt = 0;
        for (int i = 0; i < is_in_solution.size(); ++i) {
            if (std::abs(is_in_solution[i]) > SOL_X_TOLERANCE) {
                sum_length += static_cast<double>(route_length[i]);
                ++cnt;
            }
        }
        sum_length /= cnt;
        average_route_length.second += sum_length;
        ++average_route_length.first;

        auto &branch_pair = branching_data_shared.getBranchPair();
        for (auto &edge: branch_pair) {
            auto &e = edge_tmp_info[edge];
            e.basic_features.emplace_back("cluster_coeff", cluster_coeff);
            e.basic_features.emplace_back("depot_2_center", depot_2_center);
            e.basic_features.emplace_back("average_route_length",
                                          average_route_length.second / average_route_length.first);
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::collectEdgeRelatedFeatures(
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        int tree_level
    ) {
        const auto &branch_pair = branching_data_shared.getBranchPair();
        if (!edge_tmp_info.empty())
            THROW_RUNTIME_ERROR("edge_tmp_info is not empty");

        for (auto &edge: branch_pair) {
            auto &e = edge_tmp_info[edge];
            e.basic_features.emplace_back("tree_level", tree_level);
            e.basic_features.emplace_back("edge_cost",
                                          cost_mat4_vertex_ref.get()[edge.first][edge.second] / max_edge_cost);

            double res = resource_across_arcs_in_forward_sense_ref.get()[edge.first][edge.second].resources[0];
            if (!resource_across_arcs_in_backward_sense_ref.get().empty()) {
                res += resource_across_arcs_in_backward_sense_ref.get()[edge.first][edge.second].resources[0];
                res /= 2;
            }
            e.basic_features.emplace_back("edge_res", res);

            e.basic_features.emplace_back("edge_dis_2_depot",
                                          mid_point_edge_cord_2_depot[edge.first][edge.second]
                                          / max_mid_point_edge_cord_2_depot);

            auto dim_sqr = std::pow(branching_data_shared.getDim(), 2);
            e.basic_features.emplace_back("node_density_in_std_dis_vec_form",
                                          static_cast<double>(std::pow(
                                              node_density_in_std_dis_vec_form[edge.first][edge.second].size(),
                                              2)) / dim_sqr);
            e.basic_features.emplace_back("edge_2_other_convert_dis",
                                          edge_2_other_convert_dis[edge.first][edge.second]);
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::collectVariableRelatedFeatures(
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        const BrCType &candidate,
        const std::vector<int> &nonzero_idx) {
        auto &e = edge_tmp_info[candidate];
        std::vector<double> frac_up;
        frac_up.reserve(nonzero_idx.size());
        for (const auto &col_idx: nonzero_idx) {
            if (std::abs(is_in_solution[col_idx]) > SOL_X_TOLERANCE) {
                frac_up.emplace_back(1 - is_in_solution[col_idx]);
            }
        }
        if (frac_up.empty())
            THROW_RUNTIME_ERROR("frac_up is empty, no nonzero idx");
        e.basic_features.emplace_back("mean_frac_up",
                                      std::accumulate(frac_up.begin(), frac_up.end(), 0.0) / static_cast<int>(frac_up.
                                          size()));
        e.basic_features.emplace_back("min_frac_up", *std::min_element(frac_up.begin(), frac_up.end()));
        e.basic_features.emplace_back("max_frac_up", *std::max_element(frac_up.begin(), frac_up.end()));
        e.basic_features.emplace_back("frac_edge", branching_data_shared.getCandidate().at(candidate));
        auto &exact_improvement_up = branching_history.exact_improvement_up;
        auto &exact_improvement_down = branching_history.exact_improvement_down;
        auto frac_edge_down = branching_data_shared.getCandidate().at(candidate);
        auto frac_edge_up = 1 - frac_edge_down;
        auto if_up = exact_improvement_up.find(candidate) == exact_improvement_up.end() ? false : true;
        auto if_down = exact_improvement_down.find(candidate) == exact_improvement_down.end() ? false : true;
        double improvement_up = if_up
                                    ? (exact_improvement_up.at(candidate).first / exact_improvement_up.at(candidate).
                                       second)
                                    : 0;
        double improvement_down = if_down
                                      ? (exact_improvement_down.at(candidate).first / exact_improvement_down.at(
                                             candidate).
                                         second)
                                      : 0;
        double pseudo_cost_up = std::max(improvement_up * frac_edge_up, 0.0);
        double pseudo_cost_down = std::max(improvement_down * frac_edge_down, 0.0);
        double pseudo_cost_mean = std::sqrt(pseudo_cost_up * pseudo_cost_down);
        e.basic_features.emplace_back("pseudo_cost_geomean_ratio", pseudo_cost_mean);
        e.basic_features.emplace_back("ever_geomean", if_up && if_down);

        auto &lp_testing_improvement_up = branching_history.lp_testing_improvement_up;
        auto &lp_testing_improvement_down = branching_history.lp_testing_improvement_down;
        bool if_find = lp_testing_improvement_down.find(candidate) == lp_testing_improvement_down.end()
                           ? false
                           : true;
        double improvement_lp_up, improvement_lp_down;
        if (if_find) {
            improvement_lp_up = lp_testing_improvement_up.at(candidate).first / lp_testing_improvement_up.at(candidate).
                                second;
            improvement_lp_down = lp_testing_improvement_down.at(candidate).first / lp_testing_improvement_down.
                                  at(candidate).
                                  second;
        } else {
            improvement_lp_up = improvement_lp_down = 0;
        }
        double pseudo_cost_lp_up = std::max(improvement_lp_up * frac_edge_up, 0.0);
        double pseudo_cost_lp_down = std::max(improvement_lp_down * frac_edge_down, 0.0);
        double pseudo_cost_lp_mean = std::sqrt(pseudo_cost_lp_up * pseudo_cost_lp_down);
        e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_up_ratio", pseudo_cost_lp_up);
        e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_down_ratio", pseudo_cost_lp_down);
        e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_geomean_ratio", pseudo_cost_lp_mean);
        e.basic_features.emplace_back(PseudoMark + "improvement_lp_up_ratio", improvement_lp_up);
        e.basic_features.emplace_back(PseudoMark + "improvement_lp_down_ratio", improvement_lp_down);
        e.basic_features.emplace_back(PseudoMark + "ever_lp_find", if_find);
        int times = 0;
        if (branching_history.branch_choice.find(candidate) != branching_history.branch_choice.end()) {
            times = branching_history.branch_choice.at(candidate);
        }
        e.basic_features.emplace_back("branch_times", times);
        auto &lp = edge_long_info[candidate].aver_edge_lp;
        e.basic_features.emplace_back("aver_edge_lp", lp.first / lp.second);
    }
}

#endif // ROUTE_OPT_T_L2B_PHASE1_HPP
