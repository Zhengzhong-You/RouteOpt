/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rank1_coefficient_controller.hpp
 * @brief Rank-1 Coefficient Getter for calculating coefficients of Rank-1 cuts.
 *
 * This header defines the Rank1CoefficientGetter class, which provides methods to
 * compute the coefficients of Rank-1 cuts in the context of sequence information
 * and solver data. It supports both full memory and limited memory approaches
 * for coefficient calculation.
 */

#ifndef ROUTE_OPT_RANK1_COEFFICIENT_CONTROLLER_HPP
#define ROUTE_OPT_RANK1_COEFFICIENT_CONTROLLER_HPP

#include <vector>
#include <unordered_map>
#include <Eigen/Sparse>
#include "rank1_macro.hpp"
#include "route_opt_macro.hpp"
#include "solver.hpp"

namespace RouteOpt::Rank1Cuts::CoefficientGetter {
    /**
     * @brief Class responsible for calculating coefficients of Rank-1 cuts.
     *
     * The Rank1CoefficientGetter class provides methods to compute the coefficients
     * of Rank-1 cuts based on sequence information and solver data. It supports both
     * full memory and limited memory approaches for coefficient calculation.
     */
    class Rank1CoefficientGetter {
    public:
        /**
         * @brief Constructor initializing with shared Rank-1 cuts data.
         *
         * @param data_shared Reference to the shared Rank-1 cuts data.
         */
        explicit Rank1CoefficientGetter(Rank1CutsDataShared &data_shared)
            : data_shared_ref(std::ref(data_shared)) {
        }

        /**
         * @brief Computes the coefficients of Rank-1 cuts.
         *
         * This templated method calculates the coefficients of Rank-1 cuts based on
         * the provided sequence information, cuts, solver data, and memory limitation flag.
         * The results are stored in the provided matrix.
         *
         * @tparam MatrixType Type of the matrix to store the coefficients.
         * @param seq_info Vector of SequenceInfo objects representing sequence data.
         * @param cuts Vector of Rank-1 cuts.
         * @param solver Pointer to the Solver instance.
         * @param if_limit_mem Flag indicating whether to use limited memory approach.
         * @param mat Matrix to store the computed coefficients.
         */
        template<typename MatrixType>
        void getR1CCoeffs(
            const std::vector<SequenceInfo> &seq_info,
            const std::vector<R1c> &cuts,
            const Solver *solver,
            bool if_limit_mem,
            MatrixType &mat);

        /**
         * @brief Builds the enumeration matrix for Rank-1 cuts.
         *
         * This method constructs the enumeration matrix for Rank-1 cuts based on the
         * provided sparse matrix, cuts, starting index, and triplets.
         *
         * @param mat Sparse matrix representing the coefficients.
         * @param r1cs Vector of Rank-1 cuts.
         * @param start Starting row index in LP model for processing.
         * @param triplets Vector of Eigen::Triplet to store the non-zero entries.
         */
        void buildR1CEnuMatrix(
            const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat,
            const std::vector<R1c> &r1cs,
            int start,
            std::vector<Eigen::Triplet<double> > &triplets
        ) const;

        /**
         * @brief Recovers Rank-1 cuts to full memory (used in enumeration state).
         *
         * This method recovers the Rank-1 cuts to full memory. It updates the input vector with the
         * recovered cuts and may adjust internal state using the solver data to ensure
         * consistency between the sequence information and the LP model.
         *
         * @param r1cs Vector of Rank-1 cuts to be updated with the recovered data.
         * @param cols Vector of SequenceInfo objects that represent the columns to be processed.
         * @param solver Reference to the Solver instance used for recovering constraints coefficients of Rank-1 cuts.
         */

        void recoverR1CsInEnum(std::vector<R1c> &r1cs,
                               const std::vector<SequenceInfo> &cols,
                               Solver &solver);

        // Delete default constructor to enforce proper initialization.
        Rank1CoefficientGetter() = delete;

        // Default destructor.
        ~Rank1CoefficientGetter() = default;

    private:
        std::vector<std::pair<r1cIndex, std::vector<int> > > lp_v_cut_map{}; ///< Vector storing cut maps for LP.
        std::vector<std::vector<std::vector<int> > > lp_v_v_use_states{}; ///< 3D vector representing state usage.
        std::vector<int> lp_r1c_denominator{}; ///< Vector storing denominators for Rank-1 cuts.

        // Reference to shared Rank-1 cuts data.
        const std::reference_wrapper<Rank1CutsDataShared> data_shared_ref;

        /**
         * @brief Computes full memory coefficients of Rank-1 cuts.
         *
         * This templated method calculates the coefficients of Rank-1 cuts using the full
         * memory approach based on the provided sequence information, cuts, and matrix.
         *
         * @tparam MatrixType Type of the matrix to store the coefficients.
         * @param seq_info Vector of SequenceInfo objects representing sequence data.
         * @param cuts Vector of Rank-1 cuts.
         * @param mat Matrix to store the computed coefficients.
         */
        template<typename MatrixType>
        void getFullMemR1CCoeffs(
            const std::vector<SequenceInfo> &seq_info,
            const std::vector<R1c> &cuts,
            MatrixType &mat);

        /**
         * @brief Computes full memory coefficients of Rank-1 cuts using solver data.
         *
         * This templated method calculates the coefficients of Rank-1 cuts using the full
         * memory approach based on the provided LP solver data, cuts, and matrix.
         *
         * @tparam MatrixType Type of the matrix to store the coefficients.
         * @param solver Pointer to the Solver instance.
         * @param cuts Vector of Rank-1 cuts.
         * @param mat Matrix to store the computed coefficients.
         */
        template<typename MatrixType>
        void getFullMemR1CCoeffs(
            const Solver *solver,
            const std::vector<R1c> &cuts,
            MatrixType &mat);

        /**
         * @brief Extends the coefficient states for Rank-1 cuts.
         *
         * This method updates the states, sparse representation, count map, and valid sparse
         * number based on the transition from one customer to another.
         *
         * @param states Vector representing the current states.
         * @param sparse_rep Vector representing the sparse representation.
         * @param cnt Unordered map counting occurrences of each state.
         * @param valid_sparse_num Reference to the number of valid sparse entries.
         * @param from Starting customer index.
         * @param to Ending customer index.
         */
        void getCoefficientExtendR1C(std::vector<int> &states,
                                     std::vector<int> &sparse_rep,
                                     std::unordered_map<int, int> &cnt,
                                     int &valid_sparse_num,
                                     int from,
                                     int to);

        /**
         * @brief Computes limited memory coefficients of Rank-1 cuts.
         *
         * This templated method calculates the coefficients of Rank-1 cuts using the limited
         * memory approach based on the provided sequence information and matrix.
         *
         * @tparam MatrixType Type of the matrix to store the coefficients.
         * @param seq_info Vector of SequenceInfo objects representing sequence data.
         * @param mat Matrix to store the computed coefficients.
         */
        template<typename MatrixType>
        void getLimitedR1CCoeffs(const std::vector<SequenceInfo> &seq_info,
                                 MatrixType &mat);

        /**
        * @brief Prepares internal data structures for limited memory Rank-1 cut coefficient computation.
        *
        * This method initializes and sets up the necessary internal data structures to facilitate
        * the computation of Rank-1 cut coefficients in a limited memory environment. It processes
        * the provided cuts to extract and organize information required for efficient coefficient
        * calculation.
        *
        * @param cuts Vector of Rank-1 cuts to be processed.
        */
        void getLimitedR1CPre(const std::vector<R1c> &cuts);
    };
}

#include "r1c_coeffs_getter.hpp"
#include "full_r1c_coeffs_getter.hpp"
#endif // ROUTE_OPT_RANK1_COEFFICIENT_CONTROLLER_HPP
