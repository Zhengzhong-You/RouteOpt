/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_FULL_R1C_COEFFS_GETTER_HPP
#define ROUTE_OPT_FULL_R1C_COEFFS_GETTER_HPP
#include <Eigen/Sparse>
#include "rank1_coefficient_controller.hpp"


namespace RouteOpt::Rank1Cuts::CoefficientGetter {
    namespace FullMemeoryRank1GetterDetail {
        inline void helpMapVisCustomerGetter(
            const Rank1CutsDataShared &rank1CutsDataShared,
            const std::vector<SequenceInfo> &seq_info,
            std::vector<Eigen::SparseVector<int> > &vis_customer
        ) {
            auto dim = rank1CutsDataShared.getDim();
            vis_customer.assign(dim, Eigen::SparseVector<int>(seq_info.size()));

            for (int i = 0; i < seq_info.size(); ++i) {
                const auto &col = seq_info[i].col_seq;
                std::unordered_map<int, int> mp;
                mp.reserve(dim);
                for (const auto &j: col) {
                    ++mp[j];
                }
                for (const auto &pr: mp) {
                    vis_customer[pr.first].insert(i) = pr.second;
                }
            }
        }
    }

    template<typename MatrixType>
    void Rank1CoefficientGetter::getFullMemR1CCoeffs(
        const Solver *solver,
        const std::vector<R1c> &cuts,
        MatrixType &mat) {
        std::vector<Eigen::Triplet<int> > triplets;
        int num_col;
        SAFE_SOLVER(solver->getNumCol(&num_col))
        triplets.reserve(cuts.size() * num_col);
        const auto &rank1CutsDataShared = data_shared_ref.get();

        size_t numnzP;
        std::vector<int> ai_col(num_col);
        std::vector<size_t> solver_beg(2);
        std::vector<int> solver_ind(num_col);
        std::vector<double> solver_val(num_col);

        for (int row = 0; row < cuts.size(); ++row) {
            auto &cut = cuts[row].info_r1c;
            std::fill(ai_col.begin(), ai_col.end(), 0);
            int cut_size = static_cast<int>(cut.first.size());
            const auto &multi = rank1CutsDataShared.getMultiplier(cut_size, cut.second);
            auto denominator = rank1CutsDataShared.getDenominator(cut_size, cut.second);
            auto rhs = rank1CutsDataShared.getRhs(cut_size, cut.second);
            int count = 0;
            for (auto i: cut.first) {
                numnzP = num_col;
                SAFE_SOLVER(solver->XgetConstraints(&numnzP,
                    solver_beg.data(),
                    solver_ind.data(),
                    solver_val.data(),
                    i - 1,
                    1))
                for (size_t j = 0; j < numnzP; ++j) {
                    ai_col[solver_ind[j]] += multi[count];
                }
                ++count;
            }
            std::vector<Eigen::Triplet<int> > tmp_triplets(num_col);
            int cnt = 0;
            if (rhs > 0.1) tmp_triplets[cnt++] = {row, 0, rhs};
            for (int i = 0; i < num_col; ++i) {
                auto tmp = static_cast<int>(ai_col[i] / denominator);
                if (tmp) tmp_triplets[cnt++] = {row, i, tmp};
            }
            triplets.insert(triplets.end(), tmp_triplets.begin(), tmp_triplets.begin() + cnt);
        }
        mat.resize(cuts.size(), num_col);
        mat.setZero();
        mat.setFromTriplets(triplets.begin(), triplets.end());
    }

    template<typename MatrixType>
    void Rank1CoefficientGetter::getFullMemR1CCoeffs(
        const std::vector<SequenceInfo> &seq_info,
        const std::vector<R1c> &cuts,
        MatrixType &mat) {
        std::vector<Eigen::Triplet<int> > triplets;
        auto num_col = static_cast<int>(seq_info.size());
        auto num_cut = static_cast<int>(cuts.size());
        triplets.reserve(num_cut * num_col);

        std::vector<Eigen::SparseVector<int> > vis_customer;
        //vis sequence
        const auto &rank1CutsDataShared = data_shared_ref.get();
        auto dim = rank1CutsDataShared.getDim();
        FullMemeoryRank1GetterDetail::helpMapVisCustomerGetter(rank1CutsDataShared, seq_info, vis_customer);

        //vis cuts
        for (int i = 0; i < cuts.size(); ++i) {
            const auto &c = cuts[i].info_r1c;
            auto c_size = static_cast<int>(c.first.size());
            Eigen::SparseVector<int> tmp_sparse_vector(num_col);
            tmp_sparse_vector.reserve(num_col);

            std::vector<int> multi;
            int denominator, rhs;
            data_shared_ref.get().getPlanInfo(multi, denominator, rhs, c_size, c.second);

            for (int j = 0; j < c_size; ++j) {
                auto customer = c.first[j];
                tmp_sparse_vector += vis_customer[customer] * multi[j];
            }
            tmp_sparse_vector /= denominator;
            for (int k = 0; k < tmp_sparse_vector.outerSize(); ++k) {
                for (Eigen::SparseVector<int>::InnerIterator it(tmp_sparse_vector, k); it; ++it) {
                    if (it.value() == 0) continue;
                    triplets.emplace_back(i, it.index(), it.value());
                }
            }
        }
        mat.resize(num_cut, num_col);
        mat.setZero();
        mat.setFromTriplets(triplets.begin(), triplets.end());
    }


    inline void Rank1CoefficientGetter::buildR1CEnuMatrix(
        const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat,
        const std::vector<R1c> &r1cs,
        int start,
        std::vector<Eigen::Triplet<double> > &triplets) const {
        if (r1cs.empty()) return;
        auto num_col = static_cast<int>(mat.cols());
        Eigen::SparseMatrix<double, Eigen::RowMajor> tmp(1, num_col);
        for (auto &r1c: r1cs) {
            if (r1c.idx_r1c < start)continue;
            tmp.setZero();
            auto &info = r1c.info_r1c;
            std::vector<int> multi;
            int denominator, rhs;
            data_shared_ref.get().getPlanInfo(multi, denominator, rhs, static_cast<int>(info.first.size()),
                                              info.second);
            int count = 0;
            for (auto &i: info.first) {
                tmp += mat.row(i - 1) * multi[count++];
            }
            tmp /= denominator;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
                auto val_ = static_cast<int>(it.value() + TOLERANCE);
                if (val_)triplets.emplace_back(r1c.idx_r1c - start, it.col(), val_);
            }
        }
    }
}


#endif // ROUTE_OPT_FULL_R1C_COEFFS_GETTER_HPP
