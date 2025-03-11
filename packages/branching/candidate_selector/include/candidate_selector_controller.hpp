/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file candidate_selector_controller.hpp
 * @brief CandidateSelectorController class for managing candidate selection in RouteOpt.
 *
 * This header defines the BranchingTesting class, which encapsulates the candidate selection process
 * during branching in the RouteOpt framework. It manages different testing phases (LP, Heuristic, Exact)
 * by applying user-defined processing functions on nodes, and measuring the corresponding execution times.
 */

#ifndef ROUTE_OPT_CANDIDATE_SELECTOR_CONTROLLER_HPP
#define ROUTE_OPT_CANDIDATE_SELECTOR_CONTROLLER_HPP

#include <unordered_map>
#include <vector>
#include <functional>
#include "route_opt_macro.hpp"
#include "candidate_selector_macro.hpp"
#include "branching_macro.hpp"

namespace RouteOpt::Branching::CandidateSelector {
    /**
     * @brief The BranchingTesting class encapsulates the candidate selection process
     *        during branching in the RouteOpt framework.
     *
     * This templated class manages different testing phases (LP, Heuristic, Exact)
     * by applying user-defined processing functions on nodes, and measuring the
     * corresponding execution times. It also maintains phase-specific counters
     * and provides an interface to retrieve the best candidate branching decision.
     *
     * @tparam Node    Type representing a node in the branching tree.
     * @tparam BrCType Type representing a branching candidate.
     * @tparam Hasher  Hash function for candidate keys.
     */
    template<typename Node, typename BrCType, typename Hasher>
    class BranchingTesting {
    public:
        /**
         * @brief Constructs a BranchingTesting object with phase parameters and testing functions.
         *
         * @param num_phase0             Number of tests to perform in phase 0 (LP phase).
         * @param num_phase1             Number of tests to perform in phase 1 (Heuristic phase).
         * @param num_phase2             Number of tests to perform in phase 2 (Exact phase).
         * @param num_phase3             Number of tests to perform in phase 3.
         * @param processLPTestingFunction    Function to process LP testing.
         * @param processHeurTestingFunction  Function to process heuristic testing.
         * @param processExactTestingFunction Function to process exact testing.
         */
        BranchingTesting(int num_phase0, int num_phase1, int num_phase2, int num_phase3,
                         const std::function<void(Node *, const BrCType &, double &, double &)> &
                         processLPTestingFunction,
                         const std::function<void(Node *, const BrCType &, double &, double &)> &
                         processHeurTestingFunction,
                         const std::function<void(Node *, const BrCType &, double &, double &)> &
                         processExactTestingFunction) {
            setNumPhase0(num_phase0);
            setNumPhase1(num_phase1);
            setNumPhase2(num_phase2);
            setNumPhase3(num_phase3);
            setProcessLPTestingFunction(processLPTestingFunction);
            setProcessHeurTestingFunction(processHeurTestingFunction);
            setProcessExactTestingFunction(processExactTestingFunction);
        }

        /**
         * @brief Sets the number of tests for phase 0 (LP phase).
         * @param num Number of tests.
         */
        void setNumPhase0(int num) {
            num_phase0 = num;
        }

        /**
         * @brief Sets the number of tests for phase 1 (Heuristic phase).
         * @param num Number of tests.
         */
        void setNumPhase1(int num) {
            num_phase1 = num;
        }

        /**
         * @brief Sets the number of tests for phase 2 (Exact phase).
         * @param num Number of tests.
         */
        void setNumPhase2(int num) {
            num_phase2 = num;
        }

        /**
         * @brief Sets the number of tests for phase 3.
         * @param num Number of tests.
         */
        void setNumPhase3(int num) {
            num_phase3 = num;
        }

        /**
         * @brief Sets the function used for processing LP testing.
         * @param func Function to process LP testing.
         */
        void setProcessLPTestingFunction(const std::function<void(Node *, const BrCType &, double &, double &)> &func) {
            processLPTestingFunction = func;
        }

        /**
         * @brief Sets the function used for processing heuristic testing.
         * @param func Function to process heuristic testing.
         */
        void setProcessHeurTestingFunction(
            const std::function<void(Node *, const BrCType &, double &, double &)> &func) {
            processHeurTestingFunction = func;
        }

        /**
         * @brief Sets the function used for processing exact testing.
         * @param func Function to process exact testing.
         */
        void setProcessExactTestingFunction(
            const std::function<void(Node *, const BrCType &, double &, double &)> &func) {
            processExactTestingFunction = func;
        }

        /**
         * @brief Retrieves the edge score information generated during testing.
         * @return Constant reference to the vector of edge score information.
         */
        const std::vector<CandidateScoreInfo<BrCType> > &getEdgeInfo() const {
            return edge_info;
        }

        /**
         * @brief Gets the number of tests configured for phase 0.
         * @return Number of tests in phase 0.
         */
        [[nodiscard]] int getNumPhase0() const {
            return num_phase0;
        }

        /**
         * @brief Gets the number of tests configured for phase 1.
         * @return Number of tests in phase 1.
         */
        [[nodiscard]] int getNumPhase1() const {
            return num_phase1;
        }

        /**
         * @brief Gets the number of tests configured for phase 2.
         * @return Number of tests in phase 2.
         */
        [[nodiscard]] int getNumPhase2() const {
            return num_phase2;
        }

        /**
         * @brief Gets the number of tests configured for phase 3.
         * @return Number of tests in phase 3.
         */
        [[nodiscard]] int getNumPhase3() const {
            return num_phase3;
        }

        /**
         * @brief Provides a mutable reference to the LP testing time counter.
         * @return Reference to the LP time counter (pair: time and count).
         */
        auto &refLPTimeCnt() {
            return lp_time_cnt;
        }

        /**
         * @brief Provides a mutable reference to the heuristic testing time counter.
         * @return Reference to the heuristic time counter (pair: time and count).
         */
        auto &refHeuristicTimeCnt() {
            return heuristic_time_cnt;
        }

        /**
         * @brief Provides a mutable reference to the exact testing time counter.
         * @return Reference to the exact time counter (pair: time and count).
         */
        auto &refExactTimeCnt() {
            return exact_time_cnt;
        }

        /**
         * @brief Executes the testing phases and returns the best branching candidate.
         *
         * This function applies the LP, Heuristic, and Exact testing functions on the given node.
         * It measures the execution time for each phase and uses the branching history and shared
         * data to screen and select the best candidate.
         *
         * @param node                 Pointer to the current node.
         * @param branching_history    Reference to the branching history.
         * @param branching_data_shared Reference to the shared branching data.
         * @param candidate_map        Map of candidate branching decisions and their scores.
         * @return The best candidate branching decision.
         */
        const BrCType &getBestCandidate(Node *node,
                                        BranchingHistory<BrCType, Hasher> &branching_history,
                                        BranchingDataShared<BrCType, Hasher> &branching_data_shared,
                                        const std::unordered_map<BrCType, double, Hasher> &candidate_map) {
            // Update the candidate map in the shared data.
            branching_data_shared.refCandidateMap() = candidate_map;
            // Perform initial screening based on the LP phase.
            branching_history.initialScreen(branching_data_shared, num_phase0);
            // Measure LP testing time.
            lp_time_cnt.first = TimeSetter::measure([&]() {
                testing(node, branching_history, branching_data_shared, TestingPhase::LP);
            });
            lp_time_cnt.second = num_phase0 == 1 ? 0 : 2 * num_phase0;
            // Measure heuristic testing time.
            heuristic_time_cnt.first = TimeSetter::measure([&]() {
                testing(node, branching_history, branching_data_shared, TestingPhase::Heuristic);
            });
            heuristic_time_cnt.second = num_phase1 == 1 ? 0 : 2 * num_phase1;
            // Measure exact testing time.
            exact_time_cnt.first = TimeSetter::measure([&]() {
                testing(node, branching_history, branching_data_shared, TestingPhase::Exact);
            });
            exact_time_cnt.second = num_phase2 == 1 ? 0 : 2 * num_phase2;
            // Return the best candidate as determined by the initial screening.
            return branching_data_shared.refBranchPair().front();
        }

        /**
         * @brief Performs testing on the given node for a specified testing phase.
         *
         * This function is responsible for executing the appropriate testing function
         * (LP, Heuristic, or Exact) based on the provided phase parameter.
         *
         * @param node                 Pointer to the node being tested.
         * @param branching_history    Reference to the branching history.
         * @param branching_data_shared Reference to the shared branching data.
         * @param phase                The testing phase to execute.
         */
        void testing(
            Node *node,
            BranchingHistory<BrCType, Hasher> &branching_history,
            BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            TestingPhase phase);

        /**
         * @brief Updates the BKF controllers with the measured testing times.
         *
         * This function iterates over the provided BKF controllers and sets their testing
         * and node processing times based on the measured time counters from different testing phases.
         *
         * @param eps             Additional elapsed time.
         * @param bkf_controllers Vector of BKFController objects to update.
         */
        void updateBKFtime(double eps, std::vector<BKF::BKFController> &bkf_controllers) const {
            for (int i = 0; i < bkf_controllers.size(); ++i) {
                auto &bkf = bkf_controllers[i];
                if (i == 0) {
                    bkf.setTestingTime(lp_time_cnt.first, lp_time_cnt.second);
                    bkf.setNodeTime(eps + (heuristic_time_cnt.first + exact_time_cnt.first) / 2);
                } else if (i == 1) {
                    bkf.setTestingTime(heuristic_time_cnt.first, heuristic_time_cnt.second);
                    bkf.setNodeTime(eps + (exact_time_cnt.first + lp_time_cnt.first) / 2);
                } else if (i == 2) {
                    bkf.setTestingTime(exact_time_cnt.first, exact_time_cnt.second);
                    bkf.setNodeTime(eps + (lp_time_cnt.first + heuristic_time_cnt.first) / 2);
                } else
                    THROW_RUNTIME_ERROR("BKFController only supports 3 phases, but got " + std::to_string(i));
            }
        }

        // Delete the default constructor to enforce proper initialization.
        BranchingTesting() = delete;

        ~BranchingTesting() = default;

    private:
        // Number of testing phases for LP, Heuristic, Exact, and an additional phase.
        int num_phase0{};
        int num_phase1{};
        int num_phase2{};
        int num_phase3{};
        // Time counters for each testing phase (pair: total time, count of tests).
        std::pair<double, int> lp_time_cnt{}; // LP testing time and count.
        std::pair<double, int> heuristic_time_cnt{}; // Heuristic testing time and count.
        std::pair<double, int> exact_time_cnt{}; // Exact testing time and count.
        // Edge score information collected during testing.
        std::vector<CandidateScoreInfo<BrCType> > edge_info{};
        // Function objects for processing each type of testing.
        std::function<void(Node *, const BrCType &, double &, double &)> processLPTestingFunction{};
        std::function<void(Node *, const BrCType &, double &, double &)> processHeurTestingFunction{};
        std::function<void(Node *, const BrCType &, double &, double &)> processExactTestingFunction{};

        /**
         * @brief Revises the score for extremely unbalanced candidates.
         *
         * This function adjusts the candidate scores for a
         * given testing phase if the scores are extremely unbalanced.
         *
         * @param branching_history Reference to the branching history.
         * @param phase             The testing phase for which to revise scores.
         */
        void reviseExtremeUnbalancedScore(BranchingHistory<BrCType, Hasher> &branching_history,
                                          TestingPhase phase);
    };
}

#include "candidate_testing.hpp"
#include "initial_screen.hpp"
#endif // ROUTE_OPT_CANDIDATE_SELECTOR_CONTROLLER_HPP
