/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_TRAIN_HPP
#define ROUTE_OPT_T_L2B_TRAIN_HPP
#include <iostream>
#include <fstream>
#include "l2b_train.hpp"
#include "l2b_predict.hpp"
#include "l2b_macro.hpp"
#include "route_opt_macro.hpp"


namespace RouteOpt::Application::CVRP {
    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::initOutputPath() {
        mkDir(TRAIN_FOLDER_LP.data());
        mkDir(TRAIN_FOLDER_EXACT.data());
        lp_output_path = std::string(TRAIN_FOLDER_LP) + "/" + file_name_ref.get() + ".txt";
        exact_output_path = std::string(TRAIN_FOLDER_EXACT) + "/" + file_name_ref.get() + ".txt";
    }

    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::generateModelPhase1(
        Node *node,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
        int tree_level,
        const std::vector<double> &lp_solution,
        const std::vector<int> &route_length) {
        branching_history.initialScreen(branching_data_shared, L2B_PHASE0);
        l2b_controller_ref.get().getFeatureDataPhase1(branching_data_shared, branching_history, node, tree_level,
                                                      lp_solution, route_length);
        l2b_controller_ref.get().distinguishBranchPairSource(branching_data_shared, branching_history);
        simulateWriteLPPseudoCost(node, branching_history, branching_testing, branching_data_shared);

        branching_testing.setNumPhase3(1);
        branching_testing.testing(node, branching_history, branching_data_shared,
                                  Branching::CandidateSelector::TestingPhase::Exact);
        writeTrainingLPFile(branching_testing);
    }

    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::simulateWriteLPPseudoCost(
        Node *node,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared) {
        auto tmp_pair = branching_data_shared.getBranchPair();

        auto ratio_pseudo =
                static_cast<double>(l2b_controller_ref.get().branch_pair_from_pseudo.size())
                / static_cast<double>(l2b_controller_ref.get().branch_pair_from_pseudo.size() + l2b_controller_ref.get()
                                      .branch_pair_from_fractional.size());


        auto num = std::min(L2B_SIMULATE_NUM, static_cast<int>(tmp_pair.size()));
        int pseudo_size = static_cast<int>(num * ratio_pseudo);
        int frac_size = num - pseudo_size;

        auto &branch_pair = branching_data_shared.refBranchPair();
        branch_pair.clear();

        std::transform(l2b_controller_ref.get().branch_pair_from_pseudo.begin(),
                       l2b_controller_ref.get().branch_pair_from_pseudo.begin() + pseudo_size,
                       std::back_inserter(branch_pair),
                       [](const BrCType &edge) {
                           return edge;
                       });


        std::transform(l2b_controller_ref.get().branch_pair_from_fractional.begin(),
                       l2b_controller_ref.get().branch_pair_from_fractional.begin() + frac_size,
                       std::back_inserter(branch_pair),
                       [](const BrCType &edge) {
                           return edge;
                       });

        branching_testing.setNumPhase1(0);
        branching_testing.testing(node, branching_history, branching_data_shared,
                                  Branching::CandidateSelector::TestingPhase::LP);

        //recover
        branch_pair = tmp_pair;
    }

    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::writeTrainingLPFile(
        const Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing) {
        const auto &edge_info = branching_testing.getEdgeInfo();
        for (auto &edge: edge_info) {
            l2b_controller_ref.get().edge_tmp_info[edge.brc].sb_scores = edge.score;
        }
        const auto &path = lp_output_path;
        std::ofstream trainingData;
        trainingData.open(path, std::ios::app);
        std::vector<std::pair<int, int> > record;
        record.reserve(l2b_controller_ref.get().edge_tmp_info.size());
        for (auto &tmp_info: l2b_controller_ref.get().edge_tmp_info) {
            auto find_iter = std::find_if(tmp_info.second.basic_features.begin(), tmp_info.second.basic_features.end(),
                                          [](const std::pair<std::string, double> &p) {
                                              return p.first == PseudoMark + "ever_lp_find";
                                          });
            if (std::abs(find_iter->second - 1) < TOLERANCE) {
                trainingData << tmp_info.second.sb_scores;
                trainingData << " qid:" << qid;
                int cnt = 0;
                for (auto &feature: tmp_info.second.basic_features) {
                    trainingData << " " << cnt << ":" << static_cast<float>(feature.second);
                    DEBUG_INPUT_DATA_CALL(trainingData << " (" << feature.first << ") ";
                    )
                    ++cnt;
                }
                trainingData << std::endl;
            } else {
                record.emplace_back(tmp_info.first);
            }
        }

        ++qid;

        for (auto &pr: record) {
            auto &tmp_info = l2b_controller_ref.get().edge_tmp_info[pr];
            trainingData << tmp_info.sb_scores;
            trainingData << " qid:" << qid;
            int cnt = 0;
            for (auto &feature: tmp_info.basic_features) {
                trainingData << " " << cnt << ":" << static_cast<float>(feature.second);
                DEBUG_INPUT_DATA_CALL(trainingData << " (" << feature.first << ") ";
                )
                ++cnt;
            }
            trainingData << std::endl;
        }

        ++qid;

        trainingData.close();
    }

    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::generateModelPhase2(
        Node *node,
        Predict<Node, BrCType, Hasher> &predict,
        Branching::BranchingHistory<BrCType, Hasher> &branching_history,
        Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
        int tree_level,
        double local_gap,
        const std::vector<double> &lp_solution,
        const std::vector<int> &route_length) {
        branching_history.initialScreen(branching_data_shared, L2B_PHASE0);
        predict.useModelPhase1(node, branching_history, branching_data_shared, branching_testing, tree_level,
                               lp_solution, route_length);

        l2b_controller_ref.get().getFeatureDataPhase2(node, branching_data_shared, branching_history, branching_testing,
                                                      local_gap);

        branching_testing.setNumPhase3(1);
        branching_testing.testing(node, branching_history, branching_data_shared,
                                  Branching::CandidateSelector::TestingPhase::Exact);
        for (auto &edge: branching_testing.getEdgeInfo()) {
            l2b_controller_ref.get().recordDiscrepancyLongInfo(edge.brc, edge.dif1, true);
            l2b_controller_ref.get().recordDiscrepancyLongInfo(edge.brc, edge.dif2, false);
            findCGScore4OneEdge(edge.brc, edge.dif1, true);
            findCGScore4OneEdge(edge.brc, edge.dif2, false);
        }
        writeTrainingExactFile(branching_testing);
    }


    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::writeOneEdgeInfo(const BrCType &edge,
                                                                  std::ofstream &trainingData, bool if_left) {
        auto &tmp_info = l2b_controller_ref.get().edge_tmp_info[edge];
        if (tmp_info.resolving_lp_features.empty())return;
        if (tmp_info.resolving_lp_features.back().first != "product") {
            throw std::runtime_error("Error: the last feature is not product but "
                                     + tmp_info.resolving_lp_features.back().first);
        }
        trainingData << (if_left
                             ? (l2b_controller_ref.get().getEdgeLPChange().at(edge).first / edge_cg_change.at(edge).
                                first)
                             : (l2b_controller_ref.get().getEdgeLPChange().at(edge).second
                                / edge_cg_change.at(edge).second));
        //this is gap! (this time, we try if ratio would work!)
        trainingData << " qid:" << qid;
        int cnt = 0;
        for (auto &feature: tmp_info.basic_features) {
            trainingData << " " << cnt << ":" << static_cast<float>(feature.second);
            DEBUG_INPUT_DATA_CALL(std::cout << " " << cnt << " " << feature.first << " " << feature.second << " ";
            )
            ++cnt;
        }
        auto &extra = if_left ? tmp_info.extra_features_edge0 : tmp_info.extra_features_edge1;
        for (auto &feature: extra) {
            trainingData << " " << cnt << ":" << static_cast<float>(feature.second);
            DEBUG_INPUT_DATA_CALL(std::cout << " " << cnt << " " << feature.first << " " << feature.second << " ";
            )
            ++cnt;
        }
        for (auto &feature: tmp_info.resolving_lp_features) {
            trainingData << " " << cnt << ":" << static_cast<float>(feature.second);
            DEBUG_INPUT_DATA_CALL(std::cout << " " << cnt << " " << feature.first << " " << feature.second << " ";
            )
            ++cnt;
        }
        trainingData << std::endl;
        DEBUG_INPUT_DATA_CALL(std::cout << std::endl;
        )
    }


    template<typename Node, typename BrCType, typename Hasher>
    void GetTrainingData<Node, BrCType, Hasher>::writeTrainingExactFile(
        const Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing) {
        const auto &edge_info = branching_testing.getEdgeInfo();
        for (auto &edge: edge_info) {
            l2b_controller_ref.get().edge_tmp_info[edge.brc].sb_scores = edge.score;
        }
        auto &path = exact_output_path;
        std::ofstream trainingData;
        trainingData.open(path, std::ios::app);
        for (auto &edge: edge_info) {
            writeOneEdgeInfo(edge.brc, trainingData, true);
            writeOneEdgeInfo(edge.brc, trainingData, false);
        }

        ++qid;

        trainingData.close();
        edge_cg_change.clear();
        DEBUG_INPUT_DATA_CALL(L2BDetail::printFeatures(l2b_controller_ref.get().edge_tmp_info));
    }
}


#include "t_l2b_train.hpp"

#endif // ROUTE_OPT_T_L2B_TRAIN_HPP
