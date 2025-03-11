/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file branching_macro.hpp
 * @brief Macros, constants, and utility functions for the RouteOpt branching module.
 *
 * This header provides definitions for constants (such as invalid branch scores and tolerances),
 * utility functions (e.g., for verbose logging), and custom data structures (e.g., BranchingHistory,
 * EdgeScoreInfo) used in the RouteOpt branching module.
 */

#ifndef ROUTE_OPT_BRANCHING_MACRO_HPP
#define ROUTE_OPT_BRANCHING_MACRO_HPP

#include <unordered_map>
#include <functional>

namespace RouteOpt::Branching {
    // Constant representing an invalid branch score.
    constexpr double INVALID_BR_SCORE = -1.;

    // Tolerance for edge comparisons during branching decisions.
    constexpr double CANDIDATE_TOLERANCE = 1e-4;

    // Pseudo-fraction value used for fractional branching decisions.
    constexpr double PSEUDO_FRAC = 0.5;

    // Initial ratio used when performing cutting branching.
    constexpr double INITIAL_CUTTING_BRANCHING_RATIO = 0.2;

    // Discount factor used in the branching process. For the TW benchmark, a value of 0.72 is used.
    constexpr double R_DISCOUNT = 0.72; // for tw benchmark 0.5 will be much better!

    // Verbose macro to control debugging output for branching.
#define BRANCHING_VERBOSE

#ifdef BRANCHING_VERBOSE
    // If verbose mode is enabled, execute the code inside the macro.
#define BRANCHING_VERBOSE_EXEC(...) __VA_ARGS__;
#else
        // Otherwise, do nothing.
        BRANCHING_VERBOSE_EXEC();
#endif

    // Forward declaration for the templated shared branching data class.
    template<typename BrCType, typename Hasher>
    class BranchingDataShared;

    /**
     * @brief Templated structure to record branching history.
     *
     * Stores various maps to record improvements in branching decisions for different
     * methods (exact, heuristic, LP testing) as well as the count of branch choices.
     * Also tracks depth increases during the branching process.
     *
     * @tparam BrCType Type representing a branching candidate (e.g., an edge).
     * @tparam Hasher  Hash function used for the candidate type.
     */
    template<typename BrCType, typename Hasher>
    struct BranchingHistory {
        // Maps to record improvement in branch score for exact testing (up and down directions).
        std::unordered_map<BrCType, std::pair<double, int>, Hasher> exact_improvement_up{};
        std::unordered_map<BrCType, std::pair<double, int>, Hasher> exact_improvement_down{};

        // Maps to record improvement in branch score for heuristic testing (up and down directions).
        std::unordered_map<BrCType, std::pair<double, int>, Hasher> heuristic_improvement_up{};
        std::unordered_map<BrCType, std::pair<double, int>, Hasher> heuristic_improvement_down{};

        // Maps to record improvement in branch score for LP testing (up and down directions).
        std::unordered_map<BrCType, std::pair<double, int>, Hasher> lp_testing_improvement_up{};
        std::unordered_map<BrCType, std::pair<double, int>, Hasher> lp_testing_improvement_down{};

        // Map to count the number of times a branch candidate is chosen.
        std::unordered_map<BrCType, int, Hasher> branch_choice{};

        // Vector to handle increases in branching depth; each element stores paired values.
        std::vector<std::pair<std::pair<double, int>, std::pair<double, int> > > increase_depth{};

        // Record an exact improvement score for a candidate edge.
        void recordExactPerScore(const BrCType &edge, double old_val, double now_val, bool dir,
                                 int tree_level);

        // Return a value representing the increase in branching score at a given tree level.
        double tellBranchingIncreaseVal(int tree_level);

        // Perform an initial screening of candidates using shared branching data.
        void initialScreen(BranchingDataShared<BrCType, Hasher> &branching_data_shared, int num);

        // Check if a candidate edge has been recorded previously.
        bool isRecordedCandidate(const BrCType &edge) const;

        // Check if a candidate edge has been branched at least once.
        bool isOnceBranched(const BrCType &edge) const;
    };

    /**
     * @brief Structure to store edge score information.
     *
     * Contains an edge, differences for left and right score components, the overall score,
     * a ratio (e.g., max/min score), and a flag indicating if the right side has the maximum score.
     *
     * @tparam BrCType Type representing a branching candidate (e.g., an edge).
     */
    template<typename BrCType>
    struct CandidateScoreInfo {
        BrCType brc{}; // The candidate.
        double dif1{}, dif2{}; // Score differences for left and right parts.
        double score{}; // Overall computed score.
        double ratio{1.}; // Ratio of maximum to minimum score.
        bool if_right_max{}; // Flag indicating if the right side is the maximum.
    };

    /**
     * @brief Structure to represent branching bounds.
     *
     * Stores the upper and lower bounds for a branch.
     */
    struct BranchingBound {
        double ub{}; // Upper bound.
        double lb{}; // Lower bound.
    };

    /**
     * @brief Templated class for storing shared branching data.
     *
     * Maintains the problem dimension, a vector of branch pairs, and a candidate map
     * mapping an edge (of type BrCType) to its associated value.
     *
     * @tparam BrCType Type representing a branching candidate.
     * @tparam Hasher  Hash function used for the candidate.
     */
    template<typename BrCType, typename Hasher>
    class BranchingDataShared {
    public:
        /**
         * @brief Constructor initializing the shared data with a given dimension.
         * @param dim The dimension (e.g., number of customers + depot).
         */
        explicit BranchingDataShared(int dim) : dim(dim) {
        }

        /**
         * @brief Get the dimension of the instance.
         * @return The dimension.
         */
        [[nodiscard]] int getDim() const {
            return dim;
        }

        /**
         * @brief Reference to the candidate map.
         * @return A reference to the map linking candidates to their values.
         */
        auto &refCandidateMap() {
            return candidate_map;
        }

        /**
         * @brief Get the candidate map.
         * @return A const reference to the candidate map.
         */
        const auto &getCandidate() const {
            return candidate_map;
        }

        /**
         * @brief Get the branch pair vector.
         * @return A const reference to the branch pair vector.
         */
        const auto &getBranchPair() const {
            return branch_pair;
        }

        /**
         * @brief Reference to the branch pair vector.
         * @return A reference to the branch pair vector.
         */
        auto &refBranchPair() {
            return branch_pair;
        }

        // Delete the default constructor to enforce initialization with a dimension.
        BranchingDataShared() = delete;

        ~BranchingDataShared() = default;

    private:
        int dim{}; // Problem dimension (e.g., number of customers + depot).
        // Vector storing branch pairs, which are used for candidate selection.
        std::vector<BrCType> branch_pair{};
        // Map storing candidate and their associated nonzero values.
        std::unordered_map<BrCType, double, Hasher> candidate_map{};
    };
} // namespace RouteOpt::Branching

#endif // ROUTE_OPT_BRANCHING_MACRO_HPP
