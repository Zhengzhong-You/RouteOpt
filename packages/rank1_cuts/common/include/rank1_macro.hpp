/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rank1_macro.hpp
 * @brief Definitions for Rank-1 Cuts and related data structures in RouteOpt.
 *
 * This header defines constants, structures, and a shared data class used for managing
 * Rank-1 cuts (R1Cs) in the pricing process. These definitions include parameters
 * for tolerances, maximum limits, and an enumeration for pricing difficulty levels.
 */

#ifndef ROUTE_OPT_RANK1_MACRO_HPP
#define ROUTE_OPT_RANK1_MACRO_HPP

#include <vector>
#include <unordered_map>
#include <bitset>

#include "rank1_separation_macro.hpp"

namespace RouteOpt::Rank1Cuts {
    // Maximum number of Rank-1 cuts (R1Cs) considered during the pricing process.
    constexpr int MAX_NUM_R1CS_IN_PRICING = 2048;
    // Maximum possible number of R1Cs that can be associated with a single vertex i, i.e., i in C or i in M.
    constexpr int MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX = 128;
    // Tolerance used for numerical comparisons of this package.
    constexpr double RANK1_TOLERANCE = 1e-6;
    // Indicator for an invalid Rank-1 cut.
    constexpr int RANK1_INVALID = -1;
    // Initial index value for a Rank-1 cut.
    constexpr int INITIAL_IDX_R1C = -1;
    // Maximum |C|
    constexpr int MAX_RANK_ROW = 9;

    // Bitset type used to index or track R1Cs during the pricing process.
    using r1cIndex = std::bitset<MAX_NUM_R1CS_IN_PRICING>;

    // Enumeration for pricing difficulty levels.
    // This helps in categorizing or filtering cuts based on pricing complexity.
    // For now, EXTREMELY_HARD level is treated as the same as HARD level.
    enum class PRICING_HARD_LEVEL {
        EASY = 0, // Easy level of difficulty
        HARD = 1, // Hard level of difficulty
        EXTREMELY_HARD = 2 // Extremely hard level of difficulty
    };

    /**
     * @brief Structure representing a Rank-1 cut (R1c).
     *
     * The structure contains information about a cut and its associated plan, along with
     * an index, a right-hand side value, and additional memory for storing arc-related data.
     */
    struct R1c {
        // Pair containing the cut information (vector of integers) and a plan index.
        std::pair<std::vector<int>, int> info_r1c{};
        // Index of the current Rank-1 cut, initialized to INITIAL_IDX_R1C.
        int idx_r1c{INITIAL_IDX_R1C};
        // The right-hand side (rhs) value of the cut;
        int rhs{};
        // General arc memory representation;
        std::vector<std::pair<int, int> > arc_mem{};

        bool tellIfNodeMemory() const {
            routeOptLong c = 0;
            for (auto &i: info_r1c.first) {
                c.set(i);
            }
            for (auto &i: arc_mem) {
                c.set(i.first);
                c.set(i.second);
            }
            auto size = c.count();
            if (size * (size - 1) / 2 == arc_mem.size()) {
                std::cout << "this is node memory cut, size: " << size << std::endl;
                return true;
            }
            std::cout << "this is arc memory cut, size: " << size << std::endl;
            return false;
        }
    };

    /**
     * @brief Shared data class for Rank-1 cuts multipliers.
     *
     * This class pre-computes and stores optimal multipliers for various cut dimensions,
     * which are used during separation and the pricing process to quickly retrieve the multipliers, denominator,
     * and right-hand side for a given cut.
     */
    class Rank1CutsDataShared {
    public:
        /**
         * @brief Constructor.
         *
         * Initializes the shared data with a specified customer dimension and generates
         * the optimal multiplier mapping.
         *
         * @param dim Instance Dimension parameter (Number of customers + Number of depots)
         */
        explicit Rank1CutsDataShared(int dim): dim(dim) {
            generateOptimalMultiplier();
        }

        /**
         * @brief Retrieves the multiplier plan information for a given cut dimension and plan index.
         *
         * The plan information includes the states vector, the denominator, and the right-hand side (rhs).
         *
         * @param multiplier [out] Vector of integers representing the enumerator of the multiplier.
         * @param denominator [out] The denominator used in the multiplier.
         * @param rhs [out] The right-hand side value of the cut.
         * @param cut_dim The dimension of the cut used as the key in the multiplier mapping.
         * @param plan_idx The index of the specific plan within the given cut dimension.
         */
        void getPlanInfo(std::vector<int> &multiplier, int &denominator, int &rhs,
                         int cut_dim,
                         int plan_idx) const {
            const auto &plan = map_rank1_multiplier.at(cut_dim).at(plan_idx);
            multiplier = std::get<0>(plan);
            denominator = std::get<1>(plan);
            rhs = std::get<2>(plan);
        }

        /**
         * @brief Retrieves the multiplier vector (only enumerator) for a given cut dimension and plan index.
         *
         * @param cut_dim The dimension of the cut.
         * @param plan_idx The plan index for which the states vector is needed.
         * @return const std::vector<int>& The states vector.
         */
        [[nodiscard]] const std::vector<int> &getMultiplier(int cut_dim, int plan_idx) const {
            return std::get<0>(map_rank1_multiplier.at(cut_dim).at(plan_idx));
        }

        /**
         * @brief Retrieves the denominator for a given cut dimension and plan index.
         *
         * @param cut_dim The dimension of the cut.
         * @param plan_idx The plan index for which the denominator is needed.
         * @return int The denominator value.
         */
        [[nodiscard]] int getDenominator(int cut_dim, int plan_idx) const {
            return std::get<1>(map_rank1_multiplier.at(cut_dim).at(plan_idx));
        }

        /**
         * @brief Retrieves the right-hand side (rhs) value for a given cut dimension and plan index.
         *
         * @param cut_dim The dimension of the cut.
         * @param plan_idx The plan index for which the rhs is needed.
         * @return int The rhs value.
         */
        [[nodiscard]] int getRhs(int cut_dim, int plan_idx) const {
            return std::get<2>(map_rank1_multiplier.at(cut_dim).at(plan_idx));
        }

        /**
         * @brief Returns the dimension used in the shared data.
         *
         * @return int The current instance dimension.
         */
        [[nodiscard]] int getDim() const {
            return dim;
        }

    private:
        // Instance dimension parameter used in generating the optimal multiplier.
        int dim{};
        // Mapping that stores the optimal multiplier plans.
        // Key: cut dimension (int).
        // Value: vector of tuples, where each tuple contains:
        //        - std::vector<int>: enumerator vector
        //        - int: denominator value
        //        - int: right-hand side (rhs) value
        std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, int> > > map_rank1_multiplier{};

        /**
         * @brief Generates the optimal multiplier mapping.
         *
         * This private method pre-computes and fills the mapping with the best multiplier plans
         * for each cut dimension. The mapping is then used to quickly retrieve multiplier information
         * during the pricing process.
         */
        void generateOptimalMultiplier();
    };
}

#endif // ROUTE_OPT_RANK1_MACRO_HPP
