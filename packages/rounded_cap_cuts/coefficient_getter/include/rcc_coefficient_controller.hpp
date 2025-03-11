/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rcc_coefficient_controller.hpp
 * @brief Controller for RCC coefficient generation in RouteOpt.
 *
 * This header defines the RCCCoefficientController class, which provides methods
 * to compute coefficient matrices and build enumeration matrices for Rounded Cap Cuts (RCCs).
 * The coefficients are computed based on sequence information and RCC, and the
 * enumeration matrix is constructed using a sparse matrix representation.
 */

#ifndef ROUTE_OPT_RCC_COEFFICIENT_CONTROLLER_HPP
#define ROUTE_OPT_RCC_COEFFICIENT_CONTROLLER_HPP


#include <Eigen/Sparse>
#include "route_opt_macro.hpp"
#include "rcc_macro.hpp"
#include "solver.hpp"


namespace RouteOpt::RCCs::CoefficientGetter {
    /**
     * @brief Controller for generating RCC coefficient matrices.
     *
     * The RCCCoefficientController class provides static methods to compute coefficient
     * matrices for RCCs based on input sequence information and cut data. Additionally,
     * it offers functionality to build an enumeration matrix using Eigen's sparse matrix
     * representation.
     */
    class RCCCoefficientController {
    public:
        /**
         * @brief Computes the coefficient matrix for RCCs.
         *
         * This templated static method calculates the coefficient matrix based on the provided
         * sequence information and RCC cuts. The method supports an option to compute an
         * elementary coefficient matrix if required.
         *
         * @tparam MatrixType The type of the matrix to be computed (e.g., Eigen::MatrixXd).
         * @param seq_info Vector of sequence information used in the computation.
         * @param cuts Vector of RCC cuts.
         * @param if_elementary Flag indicating whether to achieve coefficient matrix using form 3 of RCC.
         * @param mat [in,out] The matrix to be computed and updated with coefficients.
         */
        template<typename MatrixType>
        static void getCoefficientRCC(const std::vector<SequenceInfo> &seq_info,
                                      const std::vector<Rcc> &cuts,
                                      bool if_elementary,
                                      MatrixType &mat);

        /**
         * @brief Builds the RCC enumeration matrix.
         *
         * This static method constructs an enumeration matrix (as a list of triplets) from a given
         * sparse matrix (which contains dim-1 rows) and the RCC cuts. The starting index is used
         * to offset the enumeration.
         *
         * @param mat Sparse matrix (with dim-1 rows) used as the base for building the enumeration matrix. (set-partitioning constraints)
         * @param rccs Vector of RCC cuts used in the enumeration process.
         * @param start Starting row index.
         * @param triplets [in,out] Vector of Eigen::Triplet<double> used to build the sparse enumeration matrix.
         */
        static void buildRCCEnuMatrix(
            const Eigen::SparseMatrix<double> &mat, // Matrix with dim-1 rows.
            const std::vector<Rcc> &rccs,
            int start,
            std::vector<Eigen::Triplet<double> > &triplets
        );


        /**
         * @brief Recovers RCCs after entering enumeration state.
         *
         * This static method converts RCCs to form 3, and update the constraint matrix with newly calculated coefficients.
         *
         * @param rccs [in,out] Vector of RCC cuts that will be updated with recovered strengthened form, i.e., form 3.
         * @param cols Vector of SequenceInfo objects representing the columns of the LP.
         * @param solver Reference to the solver instance (LP model).
         */
        static void recoverRCCsInEnum(std::vector<Rcc> &rccs,
                                      const std::vector<SequenceInfo> &cols,
                                      Solver &solver);


        // Default constructor.
        RCCCoefficientController() = default;

        // Default destructor.
        ~RCCCoefficientController() = default;
    };
} // namespace RouteOpt::RCCs::CoefficientGetter

#include "get_rcc_coefficient.hpp"  // Include template implementation.

#endif // ROUTE_OPT_RCC_COEFFICIENT_CONTROLLER_HPP
