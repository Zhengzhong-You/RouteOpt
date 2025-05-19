/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_PREDICT_HPP
#define ROUTE_OPT_T_L2B_PREDICT_HPP
#include <iostream>
#include <fstream>
#include "l2b_predict.hpp"
#include "l2b_macro.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Application::CVRP {
    namespace L2BPredictDetail {
        template<typename BrCType, typename Hasher>
        bool checkIfCGRecorded(
            const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            const BrCType &candidate) {
            // Arrays of pointers to the corresponding 'up' and 'down' maps
            std::array<const std::unordered_map<BrCType, std::pair<double, int>, Hasher> *, 2> improvement_ups = {
                &branching_history.exact_improvement_up,
                &branching_history.heuristic_improvement_up,
            };

            std::array<const std::unordered_map<BrCType, std::pair<double, int>, Hasher> *, 2> improvement_downs = {
                &branching_history.exact_improvement_down,
                &branching_history.heuristic_improvement_down,
            };

            // Iterating over the arrays
            for (size_t i = 0; i < improvement_ups.size(); ++i) {
                auto &up_map = improvement_ups[i];
                auto &down_map = improvement_downs[i]; // Select the matching down std::map

                auto up_iter = up_map->find(candidate);
                if (up_iter != up_map->end()) {
                    auto down_iter = down_map->find(candidate);
                    if (down_iter != down_map->end()) {
                        return true;
                    }
                }
            }
            return false;
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::predict_model_2(
        const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared, bool if_trust) {
        const auto &branch_pair = branching_data_shared.getBranchPair();
        edge_lp_pre.clear();
        his_recording.clear();
        his_recording.emplace_back(std::numeric_limits<float>::max(),
                                   std::numeric_limits<int>::max());
        his_recording.emplace_back(1., -1);
        if (branch_pair.size() == 1) {
            edge_lp_pre[branch_pair[0]].first.lp = 0.; //place holder
            return;
        }
        for (auto &edge: branch_pair) {
            edge_lp_pre[edge].first.lp = std::min(l2b_controller_ref.get().getEdgeLPChange().at(edge).first,
                                                  TRUST_SCORE);
            edge_lp_pre[edge].second.lp = std::min(l2b_controller_ref.get().getEdgeLPChange().at(edge).second,
                                                   TRUST_SCORE);
        }

        DMatrixHandle test;
        bst_ulong output_length;
        const float *output_result;
        output_result = nullptr;

        auto &edge = l2b_controller_ref.get().edge_tmp_info[branch_pair[0]];
        auto numFeatures = static_cast<int>(edge.basic_features.size() + edge.extra_features_edge0.size() + edge.
                                            resolving_lp_features.
                                            size());
        BoosterHandle booster = booster_2;
        int num_row = 2 * static_cast<int>(branch_pair.size());
        auto data = new float[num_row * numFeatures];

        for (int i = 0; i < num_row; i++) {
            auto &tmp_edge = l2b_controller_ref.get().edge_tmp_info[branch_pair[i / 2]];
            int j = 0;
            DEBUG_INPUT_DATA_CALL(
                std::cout<<"edge= " << branch_pair[i / 2].first << " " << branch_pair[i / 2].second << std::endl;)
            for (auto &fs: tmp_edge.basic_features) {
                DEBUG_INPUT_DATA_CALL(L2BDetail::debugInputData(fs))
                DEBUG_INPUT_DATA_CALL(L2BDetail::printAllInputData(fs))
                data[i * numFeatures + j] = static_cast<float>(fs.second);
                ++j;
            }
            auto &extra = i % 2 == 0 ? tmp_edge.extra_features_edge0 : tmp_edge.extra_features_edge1;
            for (auto &fs: extra) {
                DEBUG_INPUT_DATA_CALL(L2BDetail::debugInputData(fs))
                DEBUG_INPUT_DATA_CALL(L2BDetail::printAllInputData(fs))
                data[i * numFeatures + j] = static_cast<float>(fs.second);
                ++j;
            }
            for (auto &fs: tmp_edge.resolving_lp_features) {
                DEBUG_INPUT_DATA_CALL(L2BDetail::debugInputData(fs))
                DEBUG_INPUT_DATA_CALL(L2BDetail::printAllInputData(fs))
                data[i * numFeatures + j] = static_cast<float>(fs.second);
                ++j;
            }
            DEBUG_INPUT_DATA_CALL(std::cout<<SMALL_PHASE_SEPARATION)
        }

        DEBUG_INPUT_DATA_CALL(L2BDetail::checkMLInputData(num_row, data, numFeatures))

        SAFE_XGBOOST_CALL(XGDMatrixCreateFromMat(data,
            num_row,
            numFeatures,
            std::numeric_limits<float>::quiet_NaN(),
            &test))

        SAFE_XGBOOST_CALL(XGBoosterPredict(booster, test, 1, 0, 0, &output_length, &output_result))

        float min_val = *std::min_element(output_result, output_result + output_length);
        float max_val = *std::max_element(output_result, output_result + output_length);
        for (unsigned int i = 0; i < output_length; i++) {
            auto &e = edge_lp_pre[branch_pair[i / 2]];

            auto val = std::max(
                static_cast<int>((output_result[i] - min_val) / (max_val - min_val) * (MAX_R_SECOND_STAGE + 1) -
                                 TOLERANCE),
                0);
            auto &lp = i % 2 == 0 ? e.first.lp : e.second.lp;
            if (lp >= TRUST_SCORE && val == 0) {
                val = BEST_PRE;
                lp = TRUST_SCORE + 1; // tiebreaker
            }
            if (!L2BPredictDetail::checkIfCGRecorded<BrCType, Hasher>(branching_history, branch_pair[i / 2]) &&
                !if_trust) {
                val = MAX_R_SECOND_STAGE;
            }
            i % 2 == 0 ? e.first.pre = val : e.second.pre = val;
        }

        //print edge_lp_pre
        for (auto &pr_: edge_lp_pre) {
            auto &pr = pr_.second.first;
            auto &pr2 = pr_.second.second;
            CANDIDATE_SELECTOR_VERBOSE_EXEC(
                std::cout << "edge= " << pr_.first.first << " " << pr_.first.second << " lp= " << pr.lp << " pre= "
                << pr.pre << " lp2= " << pr2.lp << " pre2= " << pr2.pre << " score= " << pr.lp * pr2.lp <<
                std::endl;
                std::cout<<SMALL_PHASE_SEPARATION;
            )
        }

        SAFE_XGBOOST_CALL(XGDMatrixFree(test))
        delete[] data;
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::useModelPhase1(
        Node *node,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
        int tree_level,
        const std::vector<double> &lp_solution,
        const std::vector<int> &route_length) {
        if (branching_data_shared.refBranchPair().size() <= 1) {
            CANDIDATE_SELECTOR_VERBOSE_EXEC(
                PRINT_WARNING("too few candidates in phase1: " + std::to_string(branching_data_shared.refBranchPair().
                        size()) +
                    " no further selection is necessary"));
            return;
        }
        auto num = std::min(branching_testing.getNumPhase0(),
                            static_cast<int>(branching_data_shared.refBranchPair().size()));

        l2b_controller_ref.get().getFeatureDataPhase1(branching_data_shared, branching_history, node,
                                                      tree_level,
                                                      lp_solution, route_length);
        l2b_controller_ref.get().distinguishBranchPairSource(branching_data_shared, branching_history);

        auto local_ratio_pseudo =
                static_cast<double>(l2b_controller_ref.get().branch_pair_from_pseudo.size())
                / static_cast<double>(l2b_controller_ref.get().branch_pair_from_pseudo.size()
                                      + l2b_controller_ref.get().branch_pair_from_fractional.size());

        int pseudo_size = std::min(static_cast<int>(num * local_ratio_pseudo),
                                   static_cast<int>(l2b_controller_ref.get().branch_pair_from_pseudo.size()));

        int frac_size = std::min(num - pseudo_size,
                                 static_cast<int>(l2b_controller_ref.get().branch_pair_from_fractional.size()));
        if (pseudo_size + frac_size < num) {
            pseudo_size = std::min(static_cast<int>(l2b_controller_ref.get().branch_pair_from_pseudo.size()),
                                   num - frac_size);
        }

        num = pseudo_size + frac_size;

        auto &branch_pair = branching_data_shared.refBranchPair();

        auto from_pseudo = l2b_controller_ref.get().branch_pair_from_pseudo;
        predict_model_1(from_pseudo);
        std::transform(from_pseudo.begin(),
                       from_pseudo.begin() + pseudo_size,
                       branch_pair.begin(),
                       [](const auto &a) {
                           return a;
                       });
        auto from_frac = l2b_controller_ref.get().branch_pair_from_fractional;
        predict_model_1(from_frac);
        std::transform(from_frac.begin(),
                       from_frac.begin() + frac_size,
                       branch_pair.begin() + pseudo_size,
                       [](const auto &a) {
                           return a;
                       });
        branch_pair.resize(num);
        CANDIDATE_SELECTOR_VERBOSE_EXEC(
            std::cout << "#candidates in phase1= " << num<<
            std::endl;
        )
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::predict_model_1(std::vector<BrCType> &candidate_vector) {
        if (candidate_vector.empty()) return;

        DMatrixHandle test;
        bst_ulong output_length;
        const float *output_result;
        output_result = nullptr;
        auto handle = booster_1;

        auto &edge = l2b_controller_ref.get().edge_tmp_info[candidate_vector.front()];
        int numFeatures = static_cast<int>(edge.basic_features.size());
        auto data = new float[candidate_vector.size() * numFeatures];

        for (int i = 0; i < candidate_vector.size(); i++) {
            auto &tmp_edge = l2b_controller_ref.get().edge_tmp_info[candidate_vector[i]];
            int j = 0;
            DEBUG_INPUT_DATA_CALL(
                std::cout<<"edge= " << candidate_vector[i].first << " " << candidate_vector[i].second << std::endl;)
            for (auto &fs: tmp_edge.basic_features) {
                DEBUG_INPUT_DATA_CALL(L2BDetail::debugInputData(fs))
                DEBUG_INPUT_DATA_CALL(L2BDetail::printAllInputData(fs))
                data[i * numFeatures + j] = static_cast<float>(fs.second);
                ++j;
            }
            DEBUG_INPUT_DATA_CALL(std::cout<<SMALL_PHASE_SEPARATION)
        }

        DEBUG_INPUT_DATA_CALL(L2BDetail::checkMLInputData(static_cast<int>(candidate_vector.size()), data, numFeatures))


        SAFE_XGBOOST_CALL(XGDMatrixCreateFromMat(data,
            static_cast<int>(candidate_vector.size()),
            numFeatures,
            std::numeric_limits<float>::quiet_NaN(),
            &test))

        SAFE_XGBOOST_CALL(XGBoosterPredict (handle, test, 1, 0, 0, &output_length, &output_result))


        std::vector<std::pair<BrCType, double> > combined(candidate_vector.size());
        std::transform(candidate_vector.begin(), candidate_vector.end(), output_result, combined.begin(),
                       [](const auto &val, const auto &score) { return std::make_pair(val, score); });

        std::sort(combined.begin(), combined.end(), [](const auto &a, const auto &b) {
            return a.second > b.second;
        });

        std::transform(combined.begin(), combined.end(), candidate_vector.begin(), [](const auto &a) {
            return a.first;
        });

        SAFE_XGBOOST_CALL(XGDMatrixFree(test))

        delete[] data;
    }


    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::loadModel(int phase, const std::string &model_path) {
        if (phase != 1 && phase != 2) {
            THROW_RUNTIME_ERROR("phase must be 1 or 2, but got " + std::to_string(phase));
        }
        if (model_path.empty()) {
            PRINT_WARNING("model " + std::to_string(phase) + " is not loaded, because the model path is empty");
            return;
        }
        auto &booster = phase == 1 ? booster_1 : booster_2;
        SAFE_XGBOOST_CALL(XGBoosterCreate(nullptr, 0, &booster))
        SAFE_XGBOOST_CALL(XGBoosterLoadModel(booster, model_path.c_str()))
    }


    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::useModelPhase2(
        Node *node,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
        double local_gap, int tree_level) {
        branching_testing.refLPTimeCnt().first = TimeSetter::measure([&]() {
            l2b_controller_ref.get().getFeatureDataPhase2(node, branching_data_shared, branching_history,
                                                          branching_testing,
                                                          local_gap);
        });
        branching_testing.refLPTimeCnt().second = branching_testing.getNumPhase0() == 1
                                                      ? 0
                                                      : 2 * branching_testing.getNumPhase0();

        predict_model_2(branching_history, branching_data_shared,
                        branching_testing.getNumPhase0() < CONSERVATIVE_LP_TESTING_THRESHOLD);

        branching_testing.refHeuristicTimeCnt().first = TimeSetter::measure([&]() {
            deeperTesting(node, branching_data_shared, branching_history);
        });
        branching_testing.refHeuristicTimeCnt().second = num_deep_testing;
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Predict<Node, BrCType, Hasher>::useMLInGeneralFramework(
        Node *node,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
        int tree_level,
        double local_gap,
        const std::vector<double> &lp_solution,
        const std::vector<int> &route_length) {
        branching_history.initialScreen(branching_data_shared, L2B_PHASE0);
        useModelPhase1(node, branching_history, branching_data_shared, branching_testing, tree_level,
                       lp_solution, route_length);
        useModelPhase2(node, branching_history, branching_data_shared, branching_testing, local_gap, tree_level);
    }
}

#endif // ROUTE_OPT_T_L2B_PREDICT_HPP
