/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_INITIAL_SCREEN_HPP
#define ROUTE_OPT_INITIAL_SCREEN_HPP
#include <vector>
#include <unordered_map>
#include <iostream>
#include "route_opt_macro.hpp"
#include "branching_macro.hpp"
#include "candidate_selector_macro.hpp"

namespace RouteOpt::Branching {
    namespace HistoryDetail {
        inline void getGeomean(double &old_v, int &n, double new_v) {
            if (new_v < TOLERANCE) return;
            if (n == 0) {
                old_v = new_v;
                ++n;
                return;
            }
            old_v = std::exp((n * std::log(old_v) + std::log(new_v)) / (n + 1));
            ++n;
        }
    }

    template<typename BrCType, typename Hasher>
    double checkImprovements(
        const BranchingHistory<BrCType, Hasher> &branching_history,
        const std::pair<const BrCType, double> &candidate,
        bool if_only_check_exact = false) {
        // Arrays of pointers to the corresponding 'up' and 'down' maps
        std::array<const std::unordered_map<BrCType, std::pair<double, int>, Hasher> *, 3> improvement_ups = {
            &branching_history.exact_improvement_up,
            &branching_history.heuristic_improvement_up,
            &branching_history.lp_testing_improvement_up
        };

        std::array<const std::unordered_map<BrCType, std::pair<double, int>, Hasher> *, 3> improvement_downs = {
            &branching_history.exact_improvement_down,
            &branching_history.heuristic_improvement_down,
            &branching_history.lp_testing_improvement_down
        };

        // Iterating over the arrays
        for (size_t i = 0; i < improvement_ups.size(); ++i) {
            auto &up_map = improvement_ups[i];
            auto &down_map = improvement_downs[i]; // Select the matching down std::map

            auto up_iter = up_map->find(candidate.first);
            if (up_iter != up_map->end()) {
                auto down_iter = down_map->find(candidate.first);
                if (down_iter != down_map->end()) {
                    double frac_down = candidate.second;
                    double frac_up = 1 - frac_down;
                    double up = up_iter->second.first / up_iter->second.second;
                    double down = down_iter->second.first / down_iter->second.second;
                    return up * frac_up * down * frac_down;
                }
            }
            if (if_only_check_exact) break;
        }

        return INVALID_BR_SCORE;
    }

    template<typename BrCType, typename Hasher>
    void BranchingHistory<BrCType, Hasher>::initialScreen(
        BranchingDataShared<BrCType, Hasher> &branching_data_shared, int num) {
        std::vector<std::pair<BrCType, double> > fracEdges;
        std::vector<std::pair<BrCType, double> > OldBranch;
        int num_all_frac_edge = 0;
        for (const auto &edge: branching_data_shared.getCandidate()) {
            if (checkFrac(edge.second, CANDIDATE_TOLERANCE)) {
                ++num_all_frac_edge;
                if (auto res = checkImprovements(*this, edge); !equalFloat(res, INVALID_BR_SCORE)) {
                    OldBranch.emplace_back(edge.first, res);
                } else {
                    fracEdges.emplace_back(edge.first, std::abs(edge.second - 0.5));
                }
            }
        }

        std::sort(OldBranch.begin(), OldBranch.end(), [](const auto &a, const auto &b) {
            return a.second > b.second;
        });
        std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
            return a.second < b.second;
        });

        int all_branch_phase0 = std::min(num, num_all_frac_edge);
        int sudo_cap = std::min(static_cast<int>(OldBranch.size()), static_cast<int>(all_branch_phase0 * PSEUDO_FRAC));
        int frac_cap = std::min(static_cast<int>(fracEdges.size()), all_branch_phase0 - sudo_cap);
        if (frac_cap < all_branch_phase0 - sudo_cap) {
            sudo_cap = std::min(static_cast<int>(OldBranch.size()), all_branch_phase0 - frac_cap);
        }
        all_branch_phase0 = sudo_cap + frac_cap;

        auto &branch_pair = branching_data_shared.refBranchPair();
        branch_pair.resize(all_branch_phase0);
        std::transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, branch_pair.begin(), [](const auto &a) {
            return a.first;
        });

        std::transform(fracEdges.begin(), fracEdges.begin() + frac_cap, branch_pair.begin() + sudo_cap,
                       [](const auto &a) {
                           return a.first;
                       });


        BRANCHING_VERBOSE_EXEC(
            std::cout<<"pseudo= "<<sudo_cap<<" frac= "<<frac_cap<<std::endl;)
    };

    template<typename BrCType, typename Hasher>
    bool BranchingHistory<BrCType, Hasher>::isRecordedCandidate(const BrCType &candidate) const {
        return checkImprovements(*this, std::pair<const BrCType, double>(candidate, 0.)) != INVALID_BR_SCORE;
    }


    template<typename BrCType, typename Hasher>
    bool BranchingHistory<BrCType, Hasher>::isOnceBranched(const BrCType &candidate) const {
        return checkImprovements(*this, std::pair<const BrCType, double>(candidate, 0.), true) != INVALID_BR_SCORE;
    }


    template<typename BrCType, typename Hasher>
    void BranchingHistory<BrCType, Hasher>::recordExactPerScore(const BrCType &edge, double old_val, double now_val,
                                                                bool dir,
                                                                int tree_level) {
        auto dif = now_val - old_val;
        if (dif < RC_TOLERANCE * now_val) {
            THROW_RUNTIME_ERROR("branching leads to a decreasing in lb: from " + std::to_string(old_val) + " to " +
                std::to_string(now_val));
        }
        dif = std::max(dif, TOLERANCE);
        if (dir) {
            exact_improvement_up[edge].first += dif;
            ++exact_improvement_up[edge].second;
        } else {
            exact_improvement_down[edge].first += dif;
            ++exact_improvement_down[edge].second;
        }

        if (increase_depth.size() <= tree_level) increase_depth.resize(tree_level + 1);
        auto &r_star = increase_depth[tree_level];

        auto &recordings = dir ? r_star.second.second : r_star.first.second;
        auto &increase = dir ? r_star.second.first : r_star.first.first;

        HistoryDetail::getGeomean(increase, recordings, dif);
    }

    template<typename BrCType, typename Hasher>
    double BranchingHistory<BrCType, Hasher>::tellBranchingIncreaseVal(int tree_level) {
        if (increase_depth.size() <= tree_level) {
            THROW_RUNTIME_ERROR("the tree level is too deep!");
        }
        auto &left_increase = increase_depth[tree_level].first;
        auto &right_increase = increase_depth[tree_level].second;
        double r_star;
        if (left_increase.second == 0 || right_increase.second == 0) {
            r_star = std::numeric_limits<float>::max();
            for (auto &[fst, snd]: increase_depth) {
                if (fst.second == 0 || snd.second == 0) continue;
                r_star = std::min(r_star, std::sqrt(fst.first * snd.first));
            }
            if (r_star == std::numeric_limits<float>::max()) {
                const double val = (left_increase.second == 0) ? right_increase.first : left_increase.first;
                r_star = val * INITIAL_CUTTING_BRANCHING_RATIO;
            }
        } else {
            r_star = std::sqrt(left_increase.first * right_increase.first);
        }
        return r_star * R_DISCOUNT;
    }
}

#endif // ROUTE_OPT_INITIAL_SCREEN_HPP
