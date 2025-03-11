/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file bkf_controller.hpp
 * @brief BKFController class definition for managing BKF branching operations.
 *
 * This header file defines the BKFController class, which is responsible for managing
 * the BKF (Best K Formula) branching operations, including parameter estimation,
 * time measurement, and candidate selection.
 */

#ifndef ROUTE_OPT_BKF_CONTROLLER_HPP
#define ROUTE_OPT_BKF_CONTROLLER_HPP

#include <vector>
#include <unordered_map>
#include <functional>
#include <iostream>
#include "route_opt_macro.hpp"
#include "helper_bkf.hpp"

namespace RouteOpt::Branching::BKF {
    class BKFController; // Forward declaration for BKFController

    /**
     * @brief Shared data class for BKF branching operations.
     *
     * BKFDataShared holds common data used by BKF controllers, such as node counters,
     * the current best "f" value, and information regarding r_star at different depths.
     */
    class BKFDataShared {
    public:
        /**
         * @brief Get the current best "f" value.
         * @return The current f value.
         */
        [[nodiscard]] double getF() const { return f; }

        /**
         * @brief Update the f value if the new value is greater.
         * @param f The new f value.
         */
        void updateF(double f) {
            if (this->f < f) {
                std::cout << "update f from " << this->f << " to " << f << std::endl;
                this->f = f;
            }
        }

        /**
         * @brief Get the number of nodes processed before enumeration state (B4).
         * @return The number of nodes in the before enumeration state.
         */
        [[nodiscard]] int getNodeB4() const {
            return num_node_b4;
        }

        /**
         * @brief Get the number of nodes processed after enumeration state (Af).
         * @return The number of nodes in the after enumeration state.
         */
        [[nodiscard]] int getNodeAf() const {
            return num_node_af;
        }

        /**
         * @brief Increment the counter for nodes in the "B4" state.
         */
        void increaseNodeB4() {
            ++num_node_b4;
        }

        /**
         * @brief Increment the counter for nodes in the "Af" state.
         */
        void increaseNodeAf() {
            ++num_node_af;
        }

        /**
         * @brief Get the current best r value.
         * @return The current r_best value.
         */
        [[nodiscard]] double getCurrentRBest() const { return current_r_best; }

        /**
         * @brief Calculate the r_star value based on lifting and tree level.
         *
         * This function uses the provided parameters along with the BKFController
         * to update the r_star value at the current depth.
         *
         * @param lift      The lift value from branching.
         * @param tree_level The current tree level.
         * @param dir       The branching direction.
         * @param node_idx  The index of the node.
         * @param controller Reference to a BKFController instance.
         */
        void calculateRStar(double lift, int tree_level, bool dir, int node_idx,
                            BKFController &controller);

        /**
         * @brief Get the vector storing r_star depth information.
         * @return Const reference to the r_star_depth vector.
         */
        [[nodiscard]] const auto &getRStarDepth() const { return r_star_depth; }

        /**
         * @brief Get a mutable reference to the r_star_depth vector.
         * @return Reference to the r_star_depth vector.
         */
        [[nodiscard]] auto &refRStarDepth() { return r_star_depth; }

        BKFDataShared() = default;

        ~BKFDataShared() = default;

    private:
        // Counters for nodes before (B4) and after (Af) enumeration state.
        int num_node_b4{};
        int num_node_af{};
        // Current best f value.
        double f{};
        // Current best r value.
        double current_r_best{};
        // Stores r_star information at various depths.
        std::vector<std::pair<std::pair<double, int>, std::pair<double, int> > > r_star_depth{};
    };

    /**
     * @brief Controller class for BKF branching operations.
     *
     * BKFController manages various aspects of the branching process, including
     * parameter estimation, time measurement, and best k determination.
     */
    class BKFController {
    public:
        /**
         * @brief Set the estimated m and n values.
         *
         * @param m Estimated parameter m.
         * @param n Total number of elements (n).
         *
         * This function also calls evaluateM1() to update internal parameters.
         */
        void setMN(int m, int n) {
            all_n = n;
            est_m = m;
            evaluateM1();
        }

        /**
         * @brief Set the temporary flag indicating whether the current state is "B4".
         *
         * @param if_b4 Boolean flag: true if in "B4" state, false otherwise.
         */
        void setIfB4(bool if_b4) {
            tmp_if_b4 = if_b4;
        }

        /**
         * @brief Set the testing time for a given number of tests.
         *
         * @param t Total testing time.
         * @param n Number of tests.
         *
         * If n is zero, the function assigns the entire time to the corresponding testing time.
         * Otherwise, it averages the time over n tests.
         */
        void setTestingTime(double t, int n) {
            if (n == 0)
                tmp_single_testing_time = tmp_if_b4 ? t_for_one_testing_b4 : t_for_one_testing_af;
            else
                tmp_single_testing_time = t / n;
        }

        /**
         * @brief Set the time taken for processing a single node.
         *
         * @param t Time for a single node.
         */
        void setNodeTime(double t) {
            tmp_single_node_time = t;
        }

        /**
         * @brief Update time measurement statistics.
         *
         * @param sharedData Reference to shared BKF data.
         *
         * This function updates the average testing and node processing times
         * based on the current state and number of nodes processed.
         */
        void updateTimeMeasure(BKFDataShared &sharedData) {
            double &t = tmp_if_b4 ? t_for_one_testing_b4 : t_for_one_testing_af;
            double &c = tmp_if_b4 ? c_for_node_b4 : c_for_node_af;
            int n = tmp_if_b4 ? sharedData.getNodeB4() : sharedData.getNodeAf();
            updateStateAverage(tmp_single_testing_time, t, n);
            updateStateAverage(tmp_single_node_time, c, n);
            std::cout << "t= " << t << " c= " << c << " n= " << n << " if_b4= " << tmp_if_b4 << std::endl;
        }

        /**
         * @brief Get the best candidate parameter k based on shared data.
         *
         * @param sharedData Shared BKF data.
         * @param ub Upper bound.
         * @param lb Lower bound.
         * @return The best candidate k value.
         */
        int getBestK(const BKFDataShared &sharedData, double ub, double lb);

        /**
         * @brief Update the best k value for a given node.
         *
         * @param k The candidate k value.
         * @param node_idx Index of the node.
         */
        void updateOptK(int k, int node_idx) {
            opt_k[node_idx] = k;
        }

        /**
         * @brief Retrieve the best k value for a given node.
         *
         * @param node_idx Index of the node.
         * @return The best k value if found, otherwise returns est_m.
         *
         * Once the best k is accessed, it is removed from the internal map.
         * If no value is found, a reminder is printed and est_m is returned.
         */
        double getOptK(int node_idx) {
            auto it = opt_k.find(node_idx);
            if (it == opt_k.end()) {
                PRINT_REMIND("opt_k not found for node " + std::to_string(node_idx) + ", use m directly");
                return est_m;
            }
            auto res = it->second;
            opt_k.erase(it);
            return res;
        }

        /**
         * @brief Get the current alpha value.
         * @return The alpha value.
         */
        [[nodiscard]] double getAlpha() const { return alpha; }

        BKFController() = default;

        ~BKFController() = default;

    private:
        // Initialization flag.
        bool if_init{};
        // Temporary testing time for a single test.
        double tmp_single_testing_time{};
        // Temporary node processing time.
        double tmp_single_node_time{};
        // Temporary flag indicating state B4.
        bool tmp_if_b4{};
        // Alpha parameter used in BKF calculations.
        double alpha{};
        // Total number of elements (n) and estimated m parameter.
        double all_n{};
        double est_m{};
        // Map storing the best k value for nodes; once retrieved, the entry is removed.
        std::unordered_map<int, double> opt_k{};
        // Timing variables for testing in the "B4" state.
        double t_for_one_testing_b4{};
        double c_for_node_b4{};
        // Timing variables for testing in the "Af" state.
        double t_for_one_testing_af{};
        double c_for_node_af{};

        /**
         * @brief Evaluate and update internal parameter M1 based on current settings.
         */
        void evaluateM1();
    };
}

#endif // ROUTE_OPT_BKF_CONTROLLER_HPP
