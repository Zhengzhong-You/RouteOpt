/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_CANDIDATE_TESTING_HPP
#define ROUTE_OPT_CANDIDATE_TESTING_HPP
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <boost/math/tools/roots.hpp>

#include "route_opt_macro.hpp"
#include "branching_macro.hpp"
#include "candidate_selector_macro.hpp"
#include "candidate_selector_controller.hpp"

namespace RouteOpt::Branching::CandidateSelector {
    namespace CandidateSelectorDetail {
        template<typename BrCType>
        void printScore(const std::vector<CandidateScoreInfo<BrCType> > &edge_info) {
            std::cout << SMALL_PHASE_SEPARATION;
            for (int i = 0; i < static_cast<int>(edge_info.size()); ++i) {
                auto &edge = edge_info[i];
                std::ostringstream oss;
                oss << std::setw(3) << std::right << edge.brc.first
                        << "-"
                        << std::setw(3) << std::left << edge.brc.second
                        << " , "
                        << std::setw(6) << std::right << std::fixed << std::setprecision(2)
                        << edge.score;
                std::cout << std::setw(PRINT_COL_WIDTH) << std::left << oss.str() << "| ";
                if ((i + 1) % PRINT_NUM_COLS == 0) {
                    std::cout << std::endl;
                }
            }
            if (edge_info.size() % PRINT_NUM_COLS != 0) {
                std::cout << std::endl;
            }
            std::cout << SMALL_PHASE_SEPARATION;
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void BranchingTesting<Node, BrCType, Hasher>::testing(
        Node *node, BranchingHistory<BrCType, Hasher> &branching_history,
        BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        TestingPhase phase) {
        auto [marker, num_testing, processTestingFunction] = [&]()
            -> std::tuple<std::string, int &, std::function<void(Node *, const BrCType &, double &, double &)> > {
                    switch (phase) {
                        case TestingPhase::LP:
                            return {"lp", num_phase1, processLPTestingFunction};
                        case TestingPhase::Heuristic:
                            return {"heuristic", num_phase2, processHeurTestingFunction};
                        case TestingPhase::Exact:
                            return {"exact", num_phase3, processExactTestingFunction};
                        default:
                            THROW_RUNTIME_ERROR("invalid phase");
                            std::terminate(); //in case the compiler complains
                    }
                }();


        auto &branch_pair = branching_data_shared.refBranchPair();
        if (branch_pair.size() <= 1) {
            CANDIDATE_SELECTOR_VERBOSE_EXEC(
                PRINT_REMIND("no " + marker + " testing is needed: too few candidates: " + std::to_string(branch_pair.
                    size())))
            return;
        }


        edge_info.assign(branch_pair.size(), CandidateScoreInfo<BrCType>());

        auto eps = TimeSetter::measure([&]() {
            int cnt = 0;
            for (const auto &edge: branch_pair) {
                double dif1, dif2;
                processTestingFunction(node, edge, dif1, dif2);
                edge_info[cnt].dif1 = dif1;
                edge_info[cnt].dif2 = dif2;
                edge_info[cnt++].brc = edge;
            }
        });

        CANDIDATE_SELECTOR_VERBOSE_EXEC(printTimeMessage("testing " + marker, eps);)

        reviseExtremeUnbalancedScore(branching_history, phase);
        std::sort(edge_info.begin(), edge_info.end(), [](const auto &a, const auto &b) {
            return a.score > b.score;
        });

        int size = std::min(num_testing, static_cast<int>(edge_info.size()));
        branch_pair.resize(size);

        std::transform(edge_info.begin(), edge_info.begin() + size, branch_pair.begin(),
                       [](const auto &a) {
                           return a.brc;
                       });

        CANDIDATE_SELECTOR_VERBOSE_EXEC(
            CandidateSelectorDetail::printScore(edge_info);
        )
    }


    template<typename Node, typename BrCType, typename Hasher>
    void BranchingTesting<Node, BrCType, Hasher>::reviseExtremeUnbalancedScore(
        BranchingHistory<BrCType, Hasher> &branching_history,
        TestingPhase phase) {
        double geo_ratio = 0;

        for (auto &edge: edge_info) {
            edge.if_right_max = (edge.dif2 >= edge.dif1);
            if (edge.if_right_max) {
                edge.ratio = edge.dif2 / edge.dif1;
            } else {
                edge.ratio = edge.dif1 / edge.dif2;
            }
            geo_ratio += std::log(edge.ratio);
        }


        geo_ratio = std::exp(geo_ratio / static_cast<double>(edge_info.size())) * ADJUSTMENT_SCORE_RATIO;
        double max_dif = 0;

        for (auto &edge: edge_info) {
            if (edge.ratio > geo_ratio) continue;
            if (edge.if_right_max) {
                max_dif = std::max(max_dif, edge.dif2);
            } else {
                max_dif = std::max(max_dif, edge.dif1);
            }
        }


        for (auto &edge: edge_info) {
            if (edge.ratio <= geo_ratio) continue;
            if (edge.if_right_max) {
                edge.dif2 = std::min(max_dif, edge.dif2);
            } else {
                edge.dif1 = std::min(max_dif, edge.dif1);
            }
        }

        for (auto &edge: edge_info) {
            edge.score = edge.dif1 * edge.dif2;
        }


        auto &testing_improvement_down = [&]() -> auto &{
            if (phase == TestingPhase::LP) {
                return branching_history.lp_testing_improvement_down;
            }
            if (phase == TestingPhase::Heuristic) {
                return branching_history.heuristic_improvement_down;
            }
            if (phase == TestingPhase::Exact) {
                return branching_history.exact_improvement_down;
            }
            THROW_RUNTIME_ERROR("invalid phase");
            std::terminate(); //in case the compiler complains
        }();

        auto &testing_improvement_up = [&]() -> auto &{
            if (phase == TestingPhase::LP) {
                return branching_history.lp_testing_improvement_up;
            }
            if (phase == TestingPhase::Heuristic) {
                return branching_history.heuristic_improvement_up;
            }
            if (phase == TestingPhase::Exact) {
                return branching_history.exact_improvement_up;
            }
            THROW_RUNTIME_ERROR("invalid phase");
            std::terminate(); //in case the compiler complains
        }();

        for (auto &edge: edge_info) {
            testing_improvement_down[edge.brc].first += edge.dif1;
            ++testing_improvement_down[edge.brc].second;
            testing_improvement_up[edge.brc].first += edge.dif2;
            ++testing_improvement_up[edge.brc].second;
        }
    }
}

#endif // ROUTE_OPT_CANDIDATE_TESTING_HPP
