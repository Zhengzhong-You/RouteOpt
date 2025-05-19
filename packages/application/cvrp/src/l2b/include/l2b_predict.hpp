/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_L2B_PREDICT_HPP
#define ROUTE_OPT_L2B_PREDICT_HPP
#include <xgboost/c_api.h>
#include "l2b_controller.hpp"
#include "candidate_selector_controller.hpp"
#include "l2b_macro.hpp"

namespace RouteOpt::Application::CVRP {
    struct LPNPreInfo {
        double lp;
        int pre;
    };

    template<typename Node, typename BrCType, typename Hasher>
    class Predict {
    public:
        explicit Predict(
            Learning2Branch<Node, BrCType, Hasher> &l2b_controller,
            const std::function<void(Node *, const BrCType &, double &, bool)> &getCGOneSideScore,
            bool if_init,
            bool if_load_model_2 = true
        )
            : l2b_controller_ref(l2b_controller),
              getCGOneSideScore(getCGOneSideScore) {
            if (if_init) {
                loadModel(1, model_1.data());
                if (if_load_model_2)loadModel(2, model_2.data());
            }
        }

        void useModelPhase1(
            Node *node, Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
            int tree_level,
            const std::vector<double> &lp_solution,
            const std::vector<int> &route_length);

        void useMLInGeneralFramework(Node *node,
                                     Branching::BranchingHistory<BrCType, Hasher> &branching_history,
                                     Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
                                     Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &
                                     branching_testing,
                                     int tree_level,
                                     double local_gap,
                                     const std::vector<double> &lp_solution,
                                     const std::vector<int> &route_length);

        auto getNumDeepTesting() const {
            return num_deep_testing;
        }

        Predict() = delete;

        ~Predict() {
            if (booster_1) {
                XGBoosterFree(booster_1);
                booster_1 = nullptr;
            }
            if (booster_2) {
                XGBoosterFree(booster_2);
                booster_2 = nullptr;
            }
        }

    private:
        //refer
        std::reference_wrapper<Learning2Branch<Node, BrCType, Hasher> > l2b_controller_ref;
        std::function<void(Node *, const BrCType &, double &, bool)> getCGOneSideScore;

        //private
        int num_deep_testing{};
        BoosterHandle booster_1{};
        BoosterHandle booster_2{};
        std::unordered_map<BrCType,
            std::pair<LPNPreInfo, LPNPreInfo>,
            Hasher> edge_lp_pre{}; //std::pair<double, double>, first is lp, second is pre
        std::vector<std::pair<double, int> > his_recording{}; //first is the ratio, second is the pre


        void loadModel(int phase, const std::string &model_path);

        void predict_model_1(std::vector<BrCType> &candidate_vector);

        void predict_model_2(const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
                             const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
                             bool if_trust);

        void testCGOneSide(Node *node, const BrCType &candidate, bool dir,
                           Branching::BranchingHistory<BrCType, Hasher> &branching_history);

        void deeperTesting(Node *node, Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
                           Branching::BranchingHistory<BrCType, Hasher> &branching_history);

        void pickOneNodeCGTesting(Node *node, const BrCType &best, const BrCType &competitor,
                                  Branching::BranchingHistory<BrCType, Hasher> &branching_history);

        void useModelPhase2(Node *node,
                            Branching::BranchingHistory<BrCType, Hasher> &branching_history,
                            Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
                            Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
                            double local_gap, int tree_level);
    };
}

#include "t_l2b_predict.hpp"
#include "t_l2b_predict_deeper_testing.hpp"

#endif // ROUTE_OPT_L2B_PREDICT_HPP
