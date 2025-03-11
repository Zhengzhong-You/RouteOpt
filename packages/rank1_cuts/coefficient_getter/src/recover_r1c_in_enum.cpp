/*
* Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <vector>
#include <unordered_map>
#include "rank1_coefficient_controller.hpp"


namespace RouteOpt::Rank1Cuts::CoefficientGetter {
    void Rank1CoefficientGetter::recoverR1CsInEnum(std::vector<R1c> &r1cs,
                                                   const std::vector<SequenceInfo> &cols,
                                                   Solver &solver) {
        if (r1cs.empty()) return;
        for (auto &r1c: r1cs) r1c.arc_mem.clear();

        Eigen::SparseMatrix<double, Eigen::ColMajor> mat;
        getR1CCoeffs(cols, r1cs, &solver, false, mat);

        auto nz = static_cast<int>(mat.nonZeros());
        std::vector<int> cind(nz), vind(nz);
        std::vector<double> val(nz);

        nz = 0;
        //since the first column is the rhs, and does not need to be changed
        for (auto i = 1; i < mat.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
                cind[nz] = r1cs[it.row()].idx_r1c;
                vind[nz] = i;
                val[nz] = it.value();
                ++nz;
            }
        }

        cind.resize(nz);
        vind.resize(nz);
        val.resize(nz);

        SAFE_SOLVER(
            solver.XchangeCoeffs(
                cind.size(),
                cind.data(),
                vind.data(),
                val.data()
            ))

        SAFE_SOLVER(solver.reoptimize())
        double lp_val;
        SAFE_SOLVER(solver.getObjVal(&lp_val))
        std::cout << "recover rank1c lp val= " << lp_val << std::endl;
    }
}
