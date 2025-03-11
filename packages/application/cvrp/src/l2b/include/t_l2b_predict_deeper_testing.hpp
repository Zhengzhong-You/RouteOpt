/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_PREDICT_DEEPER_TESTING_HPP
#define ROUTE_OPT_T_L2B_PREDICT_DEEPER_TESTING_HPP
#include <iostream>
#include <fstream>
#include "l2b_predict.hpp"
#include "l2b_macro.hpp"
#include "route_opt_macro.hpp"


namespace RouteOpt::Application::CVRP {
    namespace L2BSecondStageDetail {
        inline void findRho(const std::vector<std::pair<double, int> > &his_recording, int d1, int d2, double &rho1,
                            double &rho2) {
            rho1 = std::numeric_limits<double>::max();
            rho2 = -std::numeric_limits<double>::max();

            // Find rho1: the smallest rho_k such that h(rho_k) > D1
            if (d1 == BEST_PRE) {
                rho1 = 1.;
            } else {
                for (const auto &[fst, snd]: his_recording) {
                    if (snd > d1 && fst < rho1) {
                        rho1 = fst;
                    }
                }
            }


            // Find rho2: the largest rho_k such that h(rho_k) < D2
            if (d2 == BEST_PRE) {
                rho2 = 1.;
            } else {
                for (const auto &[fst, snd]: his_recording) {
                    if (snd < d2 && fst > rho2) {
                        rho2 = fst;
                    }
                }
            }
            // If no valid rho1 or rho2 is found, report an error
            if (rho1 == std::numeric_limits<double>::max()) {
                THROW_RUNTIME_ERROR("rho1 is not found, please check the data");
            }
            if (rho2 == -std::numeric_limits<double>::max()) {
                THROW_RUNTIME_ERROR("rho2 is not found, please check the data");
            }
        }

        inline void getADD(int d1, int d2, double delta, double d_threshold, double rho1, double rho2,
                           double &a_result) {
            if (rho2 == 0.) {
                THROW_RUNTIME_ERROR("rho2 cannot be zero to avoid division error");
            }

            double ratio = rho1 / rho2;
            if (d1 == d2 && d1 <= d_threshold) {
                a_result = std::min(delta, ratio);
            } else if (d1 < d2) {
                a_result = std::min(1.0, ratio);
            } else {
                a_result = ratio;
            }
        }

        template<typename BrCType, typename Hasher>
        void getDeltaA(const BrCType &candidate1, const BrCType &candidate2,
                       const std::unordered_map<BrCType,
                           std::pair<LPNPreInfo, LPNPreInfo>,
                           Hasher> &edge_lp_pre,
                       const std::vector<std::pair<double, int> > &his_recording,
                       double delta, double d_threshold,
                       double &delta_a) {
            auto &e1_left = edge_lp_pre.at(candidate1).first.pre;
            auto &e1_right = edge_lp_pre.at(candidate1).second.pre;
            auto &e2_left = edge_lp_pre.at(candidate2).first.pre;
            auto &e2_right = edge_lp_pre.at(candidate2).second.pre;

            double rho1, rho2, a1, a2, a3, a4;

            findRho(his_recording, e1_left, e2_left, rho1, rho2);
            getADD(e1_left, e2_left, delta, d_threshold, rho1, rho2, a1);

            findRho(his_recording, e1_right, e2_right, rho1, rho2);
            getADD(e1_right, e2_right, delta, d_threshold, rho1, rho2, a2);

            findRho(his_recording, e1_left, e2_right, rho1, rho2);
            getADD(e1_left, e2_right, delta, d_threshold, rho1, rho2, a3);

            findRho(his_recording, e1_right, e2_left, rho1, rho2);
            getADD(e1_right, e2_left, delta, d_threshold, rho1, rho2, a4);
            delta_a = std::min({a1 * a2, a3 * a4});
        }

        template<typename BrCType, typename Hasher>
        bool isWin(const BrCType &candidate1, const BrCType &candidate2,
                   std::unordered_map<BrCType,
                       std::pair<LPNPreInfo, LPNPreInfo>,
                       Hasher> &edge_lp_pre,
                   const std::vector<std::pair<double, int> > &his_recording,
                   double delta, double d_threshold) {
            double delta_a;
            getDeltaA(candidate1, candidate2, edge_lp_pre, his_recording, delta, d_threshold, delta_a);

            auto &S_i_left = edge_lp_pre.at(candidate1).first.lp;
            auto &S_i_right = edge_lp_pre.at(candidate1).second.lp;
            auto &S_j_left = edge_lp_pre.at(candidate2).first.lp;
            auto &S_j_right = edge_lp_pre.at(candidate2).second.lp;

            double S_i_lr = S_i_left * S_i_right;
            double S_j_lr = S_j_left * S_j_right;


            if (S_i_lr >= S_j_lr * delta_a) {
                edge_lp_pre.erase(candidate2);
                return true;
            }

            return false;
        }

        template<typename BrCType, typename Hasher>
        void findBest_N_Competitor(const std::unordered_map<BrCType,
                                       std::pair<LPNPreInfo, LPNPreInfo>,
                                       Hasher> &edge_lp_pre, BrCType &best, BrCType &competitor) {
            double best_value = -std::numeric_limits<double>::max();
            double second_best_value = -std::numeric_limits<double>::max();
            int best_pre = std::numeric_limits<int>::max(); // For tie-breaking
            int second_best_pre = std::numeric_limits<int>::max();

            for (const auto &[key, value]: edge_lp_pre) {
                double S_i_lr = value.first.lp * value.second.lp;
                int max_pre = std::max(value.first.pre, value.second.pre);

                if (S_i_lr > best_value || (S_i_lr == best_value && max_pre < best_pre)) {
                    competitor = best;
                    second_best_value = best_value;
                    second_best_pre = best_pre;

                    best = key;
                    best_value = S_i_lr;
                    best_pre = max_pre;
                } else if (S_i_lr > second_best_value || (S_i_lr == second_best_value && max_pre < second_best_pre)) {
                    competitor = key;
                    second_best_value = S_i_lr;
                    second_best_pre = max_pre;
                }
            }
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::deeperTesting(
        Node *node,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history) {
        auto &branch_pair = branching_data_shared.refBranchPair();

        num_deep_testing = 0;
        while (true) {
            if (edge_lp_pre.size() == 1) break;
            BrCType best, competitor;
            L2BSecondStageDetail::findBest_N_Competitor(edge_lp_pre, best, competitor);
            pickOneNodeCGTesting(node, best, competitor, branching_history);
        }


        branch_pair = {edge_lp_pre.begin()->first};

        CANDIDATE_SELECTOR_VERBOSE_EXEC(std::cout<<"brc= "<<branch_pair.front().first<<"-"<<branch_pair.front().second
            <<std::endl);
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::pickOneNodeCGTesting(Node *node, const BrCType &best,
                                                              const BrCType &competitor,
                                                              Branching::BranchingHistory<BrCType, Hasher> &
                                                              branching_history) {
        if (L2BSecondStageDetail::isWin(best, competitor, edge_lp_pre, his_recording,
                                        THRESHOLD_DELTA, TRUST_PRE_BAR)) {
            return;
        }

        auto &best_info = edge_lp_pre.at(best);
        auto &competitor_info = edge_lp_pre.at(competitor);

        int D_max_best = std::max(best_info.first.pre, best_info.second.pre);
        int D_min_best = std::min(best_info.first.pre, best_info.second.pre);
        bool if_right_max_best = (best_info.second.pre > best_info.first.pre);

        int D_max_comp = std::max(competitor_info.first.pre, competitor_info.second.pre);
        int D_min_comp = std::min(competitor_info.first.pre, competitor_info.second.pre);


        bool condition3 = (D_max_best < D_max_comp) || (D_max_best == D_max_comp && D_max_best <= TRUST_PRE_BAR);

        if (!condition3) {
            testCGOneSide(node, best, if_right_max_best, branching_history);
            return;
        }

        bool condition4 = (D_min_best < D_min_comp) || (D_min_best == D_min_comp && D_min_best <= TRUST_PRE_BAR);

        if (!condition4) {
            testCGOneSide(node, best, !if_right_max_best, branching_history);
            return;
        }

        THROW_RUNTIME_ERROR("conditions 3 and 4 are both satisfied, please check the data");
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::testCGOneSide(Node *node,
                                                       const BrCType &candidate, bool dir,
                                                       Branching::BranchingHistory<BrCType, Hasher> &branching_history
    ) {
        double dif1;
        std::cout << "candidate= " << candidate.first << "-" << candidate.second << " : " << (dir ? "+" : "-") <<
                std::endl;
        std::cout << SMALL_PHASE_SEPARATION;
        ++num_deep_testing;
        getCGOneSideScore(node, candidate, dif1, dir);
        l2b_controller_ref.get().recordDiscrepancyLongInfo(candidate, dif1, !dir);

        auto &e_info = dir ? edge_lp_pre[candidate].second : edge_lp_pre[candidate].first;
        his_recording.emplace_back(e_info.lp / dif1, e_info.pre);
        e_info.lp = dif1;
        e_info.pre = BEST_PRE;

        auto &testing_improvement = dir
                                        ? branching_history.heuristic_improvement_up
                                        : branching_history.heuristic_improvement_down;

        testing_improvement[candidate].first += dif1;
        ++testing_improvement[candidate].second;
    }
}

#endif // ROUTE_OPT_T_L2B_PREDICT_DEEPER_TESTING_HPP
