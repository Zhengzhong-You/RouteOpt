/*
* Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include "rcc_coefficient_controller.hpp"
#include "solver.hpp"

namespace RouteOpt::RCCs::CoefficientGetter {
    void RCCCoefficientController::recoverRCCsInEnum(std::vector<Rcc> &rccs,
                                                     const std::vector<SequenceInfo> &cols,
                                                     Solver &solver) {
        if (rccs.empty()) return;

        for (auto &rcc: rccs) {
            auto idx = rcc.idx_rcc;
            double k;
            if (rcc.form_rcc == static_cast<int>(RCCs::RCCForm::RCC_FORM_1)) {
                k = static_cast<int>(rcc.info_rcc_customer.size()) - rcc.rhs;
            } else if (rcc.form_rcc == static_cast<int>(RCCs::RCCForm::RCC_FORM_2)) {
                k = static_cast<int>(rcc.info_rcc_outside_customer.size()) - rcc.rhs;
            } else
                THROW_RUNTIME_ERROR("rcc form does not match");
            auto sense = SOLVER_GREATER_EQUAL;
            SAFE_SOLVER(solver.setRhs(idx, 1, &sense, &k))
            rcc.form_rcc = static_cast<int>(RCCs::RCCForm::RCC_FORM_3);
            rcc.rhs = k;
        }

        Eigen::SparseMatrix<double, Eigen::RowMajor> mat;
        getCoefficientRCC(
            cols, rccs, true, mat);

        auto size = cols.size() * rccs.size();
        std::vector<double> val(size, 0);
        for (auto i = 0; i < mat.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, i); it; ++it) {
                val[i * cols.size() + it.col()] = it.value();
            }
        }

        //the 1st column;
        for (auto i = 0; i < rccs.size(); ++i) {
            val[i * cols.size()] = rccs[i].rhs;
        }

        std::vector<int> cind(size), vind(size);
        int cols_size = static_cast<int>(cols.size());
        for (auto i = 0; i < rccs.size(); ++i) {
            std::fill(cind.begin() + i * cols_size, cind.begin() + (i + 1) * cols_size, rccs[i].idx_rcc);
            std::iota(vind.begin() + i * cols_size, vind.begin() + (i + 1) * cols_size, 0);
        }

        SAFE_SOLVER(solver.XchangeCoeffs(cind.size(), cind.data(), vind.data(), val.data()))
        SAFE_SOLVER(solver.reoptimize())
        double lp_val;
        SAFE_SOLVER(solver.getObjVal(&lp_val))
        std::cout << "recover rcc lp val= " << lp_val << std::endl;
    }
}
