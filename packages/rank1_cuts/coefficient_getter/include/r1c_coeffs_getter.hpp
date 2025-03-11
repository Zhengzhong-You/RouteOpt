/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_R1C_COEFFS_GETTER_HPP
#define ROUTE_OPT_R1C_COEFFS_GETTER_HPP
#include "rank1_coefficient_controller.hpp"
#include "limited_r1c_coeffs_getter.hpp"


namespace RouteOpt::Rank1Cuts::CoefficientGetter {
    template<typename MatrixType>
    void Rank1CoefficientGetter::getR1CCoeffs(
        const std::vector<SequenceInfo> &seq_info,
        const std::vector<R1c> &cuts,
        const Solver *solver,
        bool if_limit_mem,
        MatrixType &mat) {
        if (if_limit_mem) {
            getLimitedR1CPre(cuts);
            getLimitedR1CCoeffs(seq_info, mat);
        } else {
            if (solver != nullptr) {
                getFullMemR1CCoeffs(solver, cuts, mat);
            } else {
                getFullMemR1CCoeffs(seq_info, cuts, mat);
            }
        }
    }
}


#endif // ROUTE_OPT_R1C_COEFFS_GETTER_HPP
