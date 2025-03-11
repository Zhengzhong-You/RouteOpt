/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_L2B_TRAIN_HPP
#define ROUTE_OPT_L2B_TRAIN_HPP
#include "l2b_controller.hpp"
#include "l2b_predict.hpp"
#include "candidate_selector_controller.hpp"
#include "l2b_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<typename BrCType>
    struct CandidateSourceInfo {
        std::vector<BrCType> branch_pair_from_pseudo{};
        std::vector<BrCType> branch_pair_from_fractional{};
    };

    template<typename Node, typename BrCType, typename Hasher>
    class GetTrainingData {
    public:
        explicit GetTrainingData(
            Learning2Branch<Node, BrCType, Hasher> &l2b_controller,
            const std::string &file_name,
            bool if_init
        )
            : l2b_controller_ref(l2b_controller),
              file_name_ref(file_name) {
            if (if_init) initOutputPath();
        }

        static void checkIfStopGeneratingData(int tree_level, bool &if_terminate) {
            if (if_terminate) return;
            if (tree_level > MAX_TREE_LEVEL) {
                if_terminate = true;
                PRINT_REMIND("maximum tree level= "+std::to_string(MAX_TREE_LEVEL)+" reached, stop generating data");
            }
        }

        void generateModelPhase1(
            Node *node,
            Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
            int tree_level,
            const std::vector<double> &lp_solution,
            const std::vector<int> &route_length);

        void generateModelPhase2(
            Node *node,
            Predict<Node, BrCType, Hasher> &predict,
            Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
            int tree_level,
            double local_gap,
            const std::vector<double> &lp_solution,
            const std::vector<int> &route_length);

        void findCGScore4OneEdge(const std::pair<int, int> &edge, double dif, bool if_left) {
            if_left ? edge_cg_change[edge].first = dif : edge_cg_change[edge].second = dif;
        }

        GetTrainingData() = delete;

        ~GetTrainingData() = default;

    private:
        //refer
        std::reference_wrapper<Learning2Branch<Node, BrCType, Hasher> > l2b_controller_ref;
        std::reference_wrapper<const std::string> file_name_ref;
        //private
        int qid{};
        std::string lp_output_path{};
        std::string exact_output_path{};
        std::unordered_map<BrCType, std::pair<double, double>, Hasher> edge_cg_change{};

        void writeOneEdgeInfo(const BrCType &edge,
                              std::ofstream &trainingData, bool if_left);

        void initOutputPath();

        void simulateWriteLPPseudoCost(
            Node *node,
            Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
            Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared);

        void writeTrainingLPFile(
            const Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing);

        void writeTrainingExactFile(
            const Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing);
    };
}

#include "t_l2b_train.hpp"
#endif // ROUTE_OPT_L2B_TRAIN_HPP
