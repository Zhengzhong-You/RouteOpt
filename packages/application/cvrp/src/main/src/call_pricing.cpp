/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "cvrp_macro.hpp"

namespace RouteOpt::Application::CVRP {
    namespace {
        bool tellIfRouteValueIntegral(const std::vector<double> &x, const std::vector<SequenceInfo> &cols) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                const auto val = x[i];
                if (val < TOLERANCE) continue;
                if (cols[i].forward_concatenate_pos < -1 || std::abs(val - 1) > TOLERANCE) {
                    return false;
                }
            }
            return true;
        }

        bool tellIfEdgeValueIntegral(
            const std::unordered_map<std::pair<int, int>, double, PairHasher> &edge_map) {
            for (const auto &[edge, val]: edge_map) {
                if (checkFrac(val, Branching::CANDIDATE_TOLERANCE)) {
                    return false;
                }
            }
            return true;
        }

        bool tryBuildSpecialEdgeIntegralRoute(const SequenceInfo &lhs,
                                             const SequenceInfo &rhs,
                                             SequenceInfo &merged) {
            std::unordered_map<int, int> lhs_count;
            std::unordered_map<int, int> rhs_count;
            for (auto node: lhs.col_seq) ++lhs_count[node];
            for (auto node: rhs.col_seq) ++rhs_count[node];

            int shared_vertex = -1;
            for (const auto &[node, count]: lhs_count) {
                auto it = rhs_count.find(node);
                if (count == 1 && it != rhs_count.end() && it->second == 1) {
                    if (shared_vertex != -1) return false;
                    shared_vertex = node;
                }
            }
            if (shared_vertex == -1) return false;

            int lhs_pos = -1;
            int rhs_pos = -1;
            for (int i = 0; i < static_cast<int>(lhs.col_seq.size()); ++i) {
                if (lhs.col_seq[i] == shared_vertex) {
                    lhs_pos = i;
                    break;
                }
            }
            for (int i = 0; i < static_cast<int>(rhs.col_seq.size()); ++i) {
                if (rhs.col_seq[i] == shared_vertex) {
                    rhs_pos = i;
                    break;
                }
            }
            if (lhs_pos == -1 || rhs_pos == -1) return false;
            if (rhs_pos + 1 >= static_cast<int>(rhs.col_seq.size())) return false;

            SequenceInfo candidate;
            candidate.col_seq.reserve(lhs_pos + 1 + static_cast<int>(rhs.col_seq.size()) - rhs_pos - 1);
            for (int i = 0; i <= lhs_pos; ++i) candidate.col_seq.emplace_back(lhs.col_seq[i]);
            for (int i = rhs_pos + 1; i < static_cast<int>(rhs.col_seq.size()); ++i) {
                candidate.col_seq.emplace_back(rhs.col_seq[i]);
            }
            if (candidate.col_seq.empty()) return false;

            std::vector<int> visit_count(MAX_NUM_CUSTOMERS, 0);
            for (auto node: candidate.col_seq) {
                if (++visit_count[node] != 1) return false;
            }
            candidate.forward_concatenate_pos = static_cast<int>(candidate.col_seq.size()) - 1;
            merged = std::move(candidate);
            return true;
        }

        // Special case repair for edge-integral but route-fractional solutions such as
        // 0-18-35-1-35-18-0 and 0-16-14-1-14-16-0. The two 0.5-routes share exactly one
        // singleton customer, so we concatenate the prefix of one route with the suffix of
        // the other to obtain an integer route like 0-18-35-1-14-16-0.
        bool repairSpecialEdgeIntegralRouteFractionality(const std::vector<double> &x,
                                                         const std::vector<SequenceInfo> &cols,
                                                         int dim,
                                                         const std::vector<double> &demand,
                                                         double cap,
                                                         std::vector<SequenceInfo> &repaired_routes) {
            struct FractionalRoute {
                double value{};
                SequenceInfo route{};
            };

            repaired_routes.clear();
            repaired_routes.reserve(cols.size());
            std::vector<FractionalRoute> fractional_routes;
            fractional_routes.reserve(cols.size());

            for (int i = 0; i < static_cast<int>(x.size()); ++i) {
                if (x[i] <= SOL_X_TOLERANCE) continue;
                if (std::abs(x[i] - 1.0) <= TOLERANCE && cols[i].forward_concatenate_pos >= -1) {
                    repaired_routes.emplace_back(cols[i]);
                } else {
                    fractional_routes.emplace_back(FractionalRoute{x[i], cols[i]});
                }
            }
            if (fractional_routes.empty()) return false;

            std::vector<char> used(fractional_routes.size(), false);
            for (int i = 0; i < static_cast<int>(fractional_routes.size()); ++i) {
                if (used[i]) continue;
                if (std::abs(fractional_routes[i].value - 0.5) > TOLERANCE) return false;

                bool matched = false;
                for (int j = i + 1; j < static_cast<int>(fractional_routes.size()); ++j) {
                    if (used[j]) continue;
                    if (std::abs(fractional_routes[i].value - fractional_routes[j].value) > TOLERANCE) continue;

                    SequenceInfo merged;
                    if (tryBuildSpecialEdgeIntegralRoute(fractional_routes[i].route, fractional_routes[j].route,
                                                         merged) ||
                        tryBuildSpecialEdgeIntegralRoute(fractional_routes[j].route, fractional_routes[i].route,
                                                         merged)) {
                        repaired_routes.emplace_back(std::move(merged));
                        used[i] = true;
                        used[j] = true;
                        matched = true;
                        break;
                    }
                }
                if (!matched) return false;
            }

            std::vector<int> coverage(dim, 0);
            for (const auto &route: repaired_routes) {
                double local_cap = 0;
                for (auto node: route.col_seq) {
                    if (node <= 0 || node >= dim) return false;
                    ++coverage[node];
                    local_cap += demand[node];
                }
                if (local_cap > cap + TOLERANCE) return false;
            }
            for (int i = 1; i < dim; ++i) {
                if (coverage[i] != 1) return false;
            }
            return true;
        }
    }

    void CVRPSolver::callPricingAtBeg(BbNode *node) {
        setEnv(node);
        if constexpr (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2) {
            GetTrainingData<BbNode, std::pair<int, int>,
                PairHasher>::checkIfStopGeneratingData(node->getTreeSize(), node->refIfTerminate());
            if (node->getIfTerminate()) return;
        }
        double eps;
        callPricing(node, std::numeric_limits<float>::max(), eps);
        if (node->getIfRootNode() && app_type == APPLICATION_TYPE::VRPTW)
            augmentNGRound(
                node, pricing_controller.refNG());
    }

    void CVRPSolver::callPricing(BbNode *node, double labeling_time_limit, double &time_4_pure_pricing) {
        if (node->getIfInEnumState()) {
        ENU:
            callInspection(node, time_4_pure_pricing);
            BbNode::regenerateEnumMat(node, nullptr, false, optimal_dual_vector);
        } else {
            callLabeling(node, labeling_time_limit, time_4_pure_pricing);
            if (node->getIfInEnumState()) {
                goto ENU;
            }
            if constexpr (ml_type != ML_TYPE::ML_GET_DATA_1 || ml_type != ML_TYPE::ML_GET_DATA_2) {
                if (pricing_controller.getIfCompleteCG()) {
                    std::vector<int> cstr_index;
                    node->findNonActiveCuts(optimal_dual_vector, cstr_index);
                    SAFE_SOLVER(node->refSolver().reoptimize(SOLVER_DUAL_SIMPLEX));
                }
            }
        }

        if (glob_timer.getTime() > TIME_LIMIT) {
            PRINT_REMIND("time limit reached!");
            node->refIfTerminate() = true;
        }

        if constexpr (ml_type != ML_TYPE::ML_NO_USE) {
            if (!node->getIfTerminate()) l2b_controller.recordEdgeLongInfo(BbNode::obtainSolEdgeMap(node));
        }

        if (!node->getIfTerminate() && pricing_controller.getIfCompleteCG()) {
            std::vector<double> x(node->getCols().size());
            if (!x.empty()) {
                SAFE_SOLVER(node->refSolver().getX(0, static_cast<int>(x.size()), x.data()))
            }
            if (!tellIfRouteValueIntegral(x, node->getCols())) {
                const auto edge_map = BbNode::obtainSolEdgeMap(node);
                if (tellIfEdgeValueIntegral(edge_map)) {
                    double current_lp;
                    SAFE_SOLVER(node->refSolver().getObjVal(&current_lp))
                    std::vector<SequenceInfo> repaired_routes;
                    if (!repairSpecialEdgeIntegralRouteFractionality(x, node->getCols(), dim, demand,
                                                                     cap, repaired_routes)) {
                        THROW_RUNTIME_ERROR(
                            "BUG: edge-integral but route-fractional LP solution could not be repaired by the special-case concatenation logic");
                    }

                    std::vector<double> repaired_x(repaired_routes.size(), 1.0);
                    bool if_feasible;
                    checkSolutionFeasibility(repaired_x, repaired_routes, if_feasible);
                    if (!if_feasible) {
                        THROW_RUNTIME_ERROR(
                            "BUG: repaired edge-integral route-fractional LP solution failed feasibility verification");
                    }

                    const auto repaired_obj = current_lp;

                    if (repaired_obj + TOLERANCE < ub) {
                        ub = repaired_obj;
                        ip_opt_sol.clear();
                        ip_opt_sol.reserve(repaired_routes.size());
                        for (const auto &route: repaired_routes) {
                            ip_opt_sol.emplace_back(route.col_seq);
                        }
                        std::cout << "\x1b[36mUpdated UB: " << ub << "\x1b[0m" << std::endl;
                    }
                    node->refValue() = repaired_obj;
                    node->refIfTerminate() = true;
                    PRINT_REMIND(
                        "edge-integral but route-fractional LP solution repaired into a feasible integer solution");
                }
            }
        }
    }

    void CVRPSolver::callInspection(BbNode *node, double &time_4_pure_pricing) {
        constexpr bool if_update_column_pool = true;
        constexpr bool if_allow_delete_col = true;
        time_4_pure_pricing = TimeSetter::measure([&]() {
            solveLPByInspection(node, if_update_column_pool, if_allow_delete_col);
        });
    };

    void CVRPSolver::callLabeling(BbNode *node, double labeling_time_limit, double &time_4_pure_pricing) {
        constexpr bool if_open_heur = true;
        constexpr bool if_open_exact = true;
        constexpr bool if_update_node_val = true;
        constexpr bool if_possible_terminate_early = false;
        constexpr bool if_fix_row = false;
        constexpr bool if_allow_delete_col = true;

        bool if_fix_meet_point = !node->getIfRootNode();
        bool if_consider_regenerate_bucket_graph = node->getIfRootNode();

        time_4_pure_pricing = TimeSetter::measure([&]() {
            solveLPInLabeling(node, if_open_heur, if_open_exact, if_update_node_val,
                              if_consider_regenerate_bucket_graph, if_possible_terminate_early,
                              if_fix_row, if_fix_meet_point, if_allow_delete_col, labeling_time_limit);
        });


        if (node->getIfTerminate() || !pricing_controller.getIfCompleteCG()) return;
        if (!node->getIfRootNode() && (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2)) return;
        pricing_controller.eliminateArcs<!IF_SYMMETRY_PROHIBIT>(node->getRCCs(), node->getR1Cs(),
                                                                node->getBrCs(), optimal_dual_vector, ub,
                                                                node->calculateOptimalGap(ub),
                                                                node->refLastGap(),
                                                                node->refNumForwardBucketArcs(),
                                                                node->refNumBackwardBucketArcs(),
                                                                node->refNumForwardJumpArcs(),
                                                                node->refNumBackwardJumpArcs());

        if constexpr (ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type == ML_TYPE::ML_GET_DATA_2) return;
        if (!node->getRCCs().empty() || !node->getR1Cs().empty()) callEnumeration(node);
    };
}
