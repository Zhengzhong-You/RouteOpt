/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_PHASE2_HPP
#define ROUTE_OPT_T_L2B_PHASE2_HPP
#include <vector>
#include <iostream>
#include <numeric>
#include "solver.hpp"
#include "route_opt_macro.hpp"
#include "helper_l2b.hpp"
#include "l2b_controller.hpp"
#include "l2b_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::getFeatureDataPhase2(
        Node *node,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing, double current_gap) {
        if (branching_data_shared.getBranchPair().size() == 1) {
            PRINT_REMIND("no data in phase2: a single candidate");
            return;
        }
        collectOneSideEdgeFeatures(branching_data_shared, branching_history, current_gap);
        dual_rc.clear();
        branching_testing.testing(node, branching_history, branching_data_shared,
                                  Branching::CandidateSelector::TestingPhase::LP);
        collectResolvingFeatures(branching_testing.getEdgeInfo());
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::collectOneSideEdgeFeatures(
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        double current_local_gap
    ) {
        auto &lp_his_edge0 = branching_history.lp_testing_improvement_down;
        auto &lp_his_edge1 = branching_history.lp_testing_improvement_up;
        auto &heuristic_his_edge0 = branching_history.heuristic_improvement_down;
        auto &heuristic_his_edge1 = branching_history.heuristic_improvement_up;
        auto &exact_his_edge0 = branching_history.exact_improvement_down;
        auto &exact_his_edge1 = branching_history.exact_improvement_up;

        for (auto &edge: branching_data_shared.getBranchPair()) {
            auto &tmp_edge0_info = edge_tmp_info[edge].extra_features_edge0;
            auto &tmp_edge1_info = edge_tmp_info[edge].extra_features_edge1;
            tmp_edge0_info.emplace_back("current_gap", current_local_gap);
            tmp_edge1_info.emplace_back("current_gap", current_local_gap);

            L2BDetail::addEdgeInfo(tmp_edge0_info, lp_his_edge0, edge, "indicator_lp_change",
                                   "historical_lp_obj_change");
            L2BDetail::addEdgeInfo(tmp_edge1_info, lp_his_edge1, edge, "indicator_lp_change",
                                   "historical_lp_obj_change");

            L2BDetail::addEdgeInfo(tmp_edge0_info,
                                   heuristic_his_edge0,
                                   edge,
                                   "indicator_heuristic_change",
                                   "historical_heuristic_obj_change");
            L2BDetail::addEdgeInfo(tmp_edge1_info,
                                   heuristic_his_edge1,
                                   edge,
                                   "indicator_heuristic_change",
                                   "historical_heuristic_obj_change");

            L2BDetail::addEdgeInfo(tmp_edge0_info, exact_his_edge0, edge, "indicator_exact_change",
                                   "historical_exact_obj_change");
            L2BDetail::addEdgeInfo(tmp_edge1_info, exact_his_edge1, edge, "indicator_exact_change",
                                   "historical_exact_obj_change");
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::collectResolvingDualRC(
        const Solver &solver,
        const BrCType &candidate,
        int BeforeNumRow,
        bool if_left) {
        std::vector<double> duals(3, 0);
        L2BDetail::checkOptimalStatus(solver);
        if (candidate.first) SAFE_SOLVER(solver.getDual(candidate.first - 1, 1, duals.data()))
        SAFE_SOLVER(solver.getDual(candidate.second - 1, 1, duals.data() + 1))
        SAFE_SOLVER(solver.getDual(BeforeNumRow, 1, duals.data() + 2))
        double rc_edge;

        auto dis = cost_mat4_vertex_ref.get()[candidate.first][candidate.second];
        if (equalFloat(dis, 0.)) {
            PRINT_REMIND("edge cost is 0, edge: " + std::to_string(candidate.first) + " " +
                std::to_string(candidate.second));
            rc_edge = -std::numeric_limits<float>::max();
        } else {
            rc_edge = 1 - (0.5 * (duals[0] + duals[1]) + duals[2]) / cost_mat4_vertex_ref.get()[candidate.first][
                          candidate.second];
        }

        if (if_left) {
            dual_rc.emplace_back();
            auto &dual = dual_rc.back();
            dual.candidate = candidate;
            dual.dual1 = duals[2];
            dual.rc1 = rc_edge;
        } else {
            auto &dual = dual_rc.back();
            dual.dual2 = duals[2];
            dual.rc2 = rc_edge;
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::collectResolvingFeatures(
        const std::vector<Branching::CandidateScoreInfo<BrCType> > &edge_info) {
        std::unordered_map<BrCType, size_t, Hasher> edge_to_index;
        for (size_t i = 0; i < edge_info.size(); ++i) {
            edge_to_index[edge_info[i].brc] = i;
        }


        std::sort(dual_rc.begin(), dual_rc.end(),
                  [&edge_to_index](const auto &a, const auto &b) {
                      return edge_to_index.at(a.candidate) < edge_to_index.at(b.candidate);
                  });

        for (int i = 0; i < edge_info.size(); ++i) {
            auto &e_info = edge_info[i];
            auto &d_info = dual_rc[i];
            auto &e = edge_tmp_info[edge_info[i].brc];
            auto &edge = e_info.brc;
            if (e_info.brc != d_info.candidate)
                THROW_RUNTIME_ERROR("edge is not equal, edge_info: " + std::to_string(e_info.brc.first) + " " +
                std::to_string(e_info.brc.second) + ", dual_rc: " + std::to_string(d_info.candidate.first)
                + " " + std::to_string(d_info.candidate.second));
            e.extra_features_edge0.emplace_back("rc_edge", d_info.rc1);
            e.extra_features_edge0.emplace_back("dual", d_info.dual1);
            e.extra_features_edge0.emplace_back("dif", e_info.dif1);
            updateOneSideLPChange(edge, e_info.dif1, true);
            findDiscrepancyResolvingFeatures(edge, true);
            e.extra_features_edge1.emplace_back("rc_edge", d_info.rc2);
            e.extra_features_edge1.emplace_back("dual", d_info.dual2);
            e.extra_features_edge1.emplace_back("dif", e_info.dif2);
            updateOneSideLPChange(edge, e_info.dif2, false);
            findDiscrepancyResolvingFeatures(edge, false);
            double pro = edge_lp_change.at(edge).first * edge_lp_change.at(edge).second;
            e.resolving_lp_features.emplace_back("product", pro);
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::updateOneSideLPChange(const BrCType &candidate, double obj_change,
                                                                       bool if_left) {
        if_left ? (edge_lp_change[candidate].first = obj_change) : (edge_lp_change[candidate].second = obj_change);
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::findDiscrepancyResolvingFeatures(
        const BrCType &candidate, bool if_left) {
        auto &long_ =
                if_left
                    ? edge_long_info[candidate].aver_exact_lp_discrepancy_down
                    : edge_long_info[candidate].aver_exact_lp_discrepancy_up;
        auto &e = if_left
                      ? edge_tmp_info[candidate].extra_features_edge0
                      : edge_tmp_info[candidate].extra_features_edge1;
        auto string1 = if_left ? "aver_exact_lp_discrepancy_down" : "aver_exact_lp_discrepancy_up";
        if (long_.second == 0) {
            e.emplace_back("dis_times", 0);
            e.emplace_back(string1, 0);
        } else {
            e.emplace_back("dis_times", long_.second);
            e.emplace_back(string1, long_.first / long_.second);
        }
    }
}

#endif // ROUTE_OPT_T_L2B_PHASE2_HPP
