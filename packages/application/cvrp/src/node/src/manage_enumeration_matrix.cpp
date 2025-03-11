/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "node.hpp"
#include "node_macro.hpp"

namespace RouteOpt::Application::CVRP {
    void BbNode::changeEnumMatByCuts(
        const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter) {
        if (index_columns_in_enumeration_column_pool.size() == 0) return;

        int oldNum = 0;
        for (auto &it: matrix_in_enumeration) oldNum += static_cast<int>(it.rows());

        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        auto size_enumeration_col_pool = static_cast<int>(index_columns_in_enumeration_column_pool.size());

        sparseRowMatrixXd mat(num_row - oldNum, size_enumeration_col_pool);
        std::vector<Eigen::Triplet<double> >
                triplets;
        triplets.reserve(
            static_cast<size_t>(static_cast<double>(num_row) * static_cast<double>(size_enumeration_col_pool) *
                                NodeLPDensityEstimation));

        buildRCCInEnuMatrix(triplets, oldNum);
        buildAllR1CInEnuMatrix(triplets, rank1_coefficient_getter, oldNum);

        mat.setFromTriplets(triplets.begin(), triplets.end());
        matrix_in_enumeration.push_back(std::move(mat));
    }
}
