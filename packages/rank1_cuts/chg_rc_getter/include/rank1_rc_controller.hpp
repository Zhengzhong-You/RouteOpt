/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rank1_rc_controller.hpp
 * @brief Controller and utilities for managing Rank-1 cut pricing statistics and states.
 *
 * This header defines structures and a controller class for handling the retrieval,
 * update, and manipulation of Rank-1 cut related data in the context of pricing.
 * It includes operations for copying pricing statistics, updating states,
 * performing dominance checks, concatenating state information, and assigning label memory.
 */

#ifndef ROUTE_OPT_RANK1_RC_CONTROLLER_HPP
#define ROUTE_OPT_RANK1_RC_CONTROLLER_HPP

#include <algorithm>
#include "rank1_macro.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Rank1Cuts::RCGetter {
    /**
     * @brief Structure holding Rank-1 cut pricing statistics.
     *
     * This structure stores the number of valid cuts and pointers to arrays containing
     * indices of valid cuts and a complete cut map. It also provides a method to copy data from another instance.
     */
    struct R1CPricingStat {
        int num{}; ///< Number of valid cuts.
        int *valid_cut_idx{}; ///< Pointer to an array of valid cut indices.
        int *cut_map{}; ///< Pointer to the cut state.

        /**
         * @brief Copies pricing statistics from another R1CPricingStat instance.
         *
         * @param other The source R1CPricingStat instance.
         */
        void copyFrom(const R1CPricingStat &other) {
            num = other.num;
            std::copy_n(other.valid_cut_idx, num, valid_cut_idx);
            std::copy_n(other.cut_map, MAX_NUM_R1CS_IN_PRICING, cut_map);
        }
    };

    /**
     * @brief Structure representing the usage states for Rank-1 cuts.
     *
     * This structure contains various state representations
     */
    struct R1CUseStates {
        std::vector<int> v_union_mem; ///< Vector representing union (C+M) memory
        std::pair<std::vector<int>, std::vector<int> > add_vec{};
        ///fst cut dimension, snd sparse representation of cuts.

        R1CUseStates(int cut_dim) {
            add_vec.first.assign(cut_dim, 0);
        }

        R1CUseStates() = delete;
    };

    /**
     * @brief Controller for Rank-1 cut pricing and state management.
     *
     * The Rank1RCController class provides functionality to:
     * - Initialize label memory for pricing statistics.
     * - Retrieve dual values from the column generation (CG) process.
     * - Update Rank-1 cut states.
     * - Apply dominance rules and state concatenation.
     * - Manage dynamic memory for label assignments.
     */
    class Rank1RCController {
    public:
        /**
         * @brief Constructs a Rank1RCController with a reference to shared Rank-1 cut data.
         *
         * @param data_shared Reference to the shared Rank-1 cuts data (used for retrieving multiplier plans and related data).
         */
        explicit Rank1RCController(Rank1CutsDataShared &data_shared)
            : data_shared_ref(std::ref(data_shared)) {
        }

        /**
         * @brief Initializes the label fields of a given pricing statistic.
         *
         * This static method fills the entire cut_map array with zeros.
         *
         * @param r1c The R1CPricingStat instance to initialize.
         */
        static void initializeLabel(const R1CPricingStat &r1c) {
            std::fill_n(r1c.cut_map, MAX_NUM_R1CS_IN_PRICING, 0);
        }

        /**
         * @brief Retrieves Rank-1 dual values.
         *
         * This method updates the dual map for the given set of cuts using the provided π-vector.
         *
         * @param cuts Vector of Rank-1 cuts.
         * @param pi_vector Vector of dual prices.
         */
        void getRank1DualsInCG(
            const std::vector<R1c> &cuts,
            const std::vector<double> &pi_vector);

        /**
         * @brief Updates the Rank-1 cut states.
         *
         * This method updates the state (rc) and outputs new pricing statistics based on input states,
         * transferring information from the customer "from" to "to".
         *
         * @param rc [in,out] Reference to the state value.
         * @param out_states [out] Output pricing statistics after update.
         * @param in_states Input pricing statistics.
         * @param from Starting customer for the update.
         * @param to Ending customer for the update.
         */
        void updateR1CStates(double &rc, R1CPricingStat &out_states,
                             const R1CPricingStat &in_states,
                             int from,
                             int to);

        /**
         * @brief Checks if dominance holds between two sets of pricing statistics.
         *
         * This method determines whether the output states dominate the input states by comparing
         * a computed gap.
         *
         * @param gap [in,out] The gap value used for the dominance test.
         * @param out_states The output pricing statistics.
         * @param in_states The input pricing statistics.
         * @return true if dominance condition is satisfied (`out` is better than `in`); false otherwise.
         */
        bool doR1CDominance(double &gap,
                            const R1CPricingStat &out_states,
                            const R1CPricingStat &in_states);

        /**
         * @brief Concatenates Rank-1 cut states.
         *
         * This method concatenates the input states into the output states based on the provided requirement.
         *
         * @param rc [in,out] Reference to the current rc.
         * @param req The required value for concatenation.
         * @param out_states [out] Output pricing statistics after concatenation.
         * @param in_states Input pricing statistics.
         * @param out `out` customer index for the output state.
         * @param in `in` customer index for the input state.
         * @return true if concatenation is successful; false otherwise.
         */
        bool concatenateR1CStates(double &rc, double req, R1CPricingStat &out_states,
                                  const R1CPricingStat &in_states, int out, int in);

        /**
         * @brief Assigns label memory to an array of label objects.
         *
         * This templated method assigns contiguous memory space for the cut_map and valid_cut_idx
         * fields for each label object in the provided array.
         *
         * @tparam T Type of label objects.
         * @tparam R1cPtr Pointer-to-member type pointing to the R1c field within label objects.
         * @param all_label Pointer to the array of label objects.
         * @param label_assign Number of labels to assign memory for.
         * @param r1c_ptr Pointer-to-member indicating the R1c field in the label objects.
         *
         * @note Throws a runtime error if all_label is a null pointer.
         */
        template<typename T, typename R1cPtr>
        void assignLabelMem(T *all_label, size_t label_assign, R1cPtr r1c_ptr) {
            if (all_label == nullptr)
                THROW_RUNTIME_ERROR("all_label is nullptr");
            delete[] label_int_space;
            size_t constexpr label_int_space_len = MAX_NUM_R1CS_IN_PRICING + MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX;
            label_int_space = new int[label_int_space_len * label_assign];
            size_t cut_map_beg = 0;
            for (size_t i = 0; i < label_assign; ++i) {
                auto &r1c = all_label[i].*r1c_ptr;
                r1c.cut_map = label_int_space + cut_map_beg;
                cut_map_beg += MAX_NUM_R1CS_IN_PRICING;
                r1c.valid_cut_idx = label_int_space + cut_map_beg;
                cut_map_beg += MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX;
            }
        }

        // Delete default constructor to enforce proper initialization.
        Rank1RCController() = delete;

        /**
         * @brief Destructor.
         *
         * Releases dynamically allocated memory for label_int_space.
         */
        ~Rank1RCController() {
            delete[] label_int_space;
        }

    private:
        int *label_int_space{}; ///< Dynamically allocated memory space for label assignments.
        std::vector<R1CUseStates> cg_v_cut_map{}; ///< Vector of cut's state usage for column generation.
        std::vector<std::vector<r1cIndex> > cg_v_v_use_states{};
        std::vector<int> cg_r1c_denominator{}; ///< Vector storing denominators for Rank-1 cuts.
        std::vector<double> rank1_dual{}; ///< Vector storing nonzero dual values.
        std::vector<double> revised_rank1_dual{}; ///< Vector storing revised dual values for dominance rules.

        // Reference to shared Rank-1 cuts data.
        const std::reference_wrapper<Rank1CutsDataShared> data_shared_ref;
    };
} // namespace RouteOpt::Rank1Cuts::RCGetter

#endif // ROUTE_OPT_RANK1_RC_CONTROLLER_HPP
