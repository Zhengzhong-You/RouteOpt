/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_GET_RCC_COEFFICIENT_HPP
#define ROUTE_OPT_GET_RCC_COEFFICIENT_HPP
#include <Eigen/Sparse>
#include <vector>
#include "rcc_coefficient_controller.hpp"


namespace RouteOpt::RCCs::CoefficientGetter {
    namespace CoefficientGetterDetail {
        inline void getMap(const std::vector<SequenceInfo> &cols,
                           std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int> >, PairHasher> &
                           map_idx) {
            map_idx.clear();
            for (int i = 0; i < static_cast<int>(cols.size()); ++i) {
                const auto &c = cols[i].col_seq;
                int b4 = 0;
                for (int j: c) {
                    std::pair<int, int> key = (b4 < j) ? std::make_pair(b4, j) : std::make_pair(j, b4);
                    auto &vec = map_idx[key];
                    if (!vec.empty() && vec.back().first == i) {
                        ++vec.back().second;
                    } else {
                        vec.emplace_back(i, 1);
                    }
                    b4 = j;
                }
                auto &vec = map_idx[{0, b4}];
                if (!vec.empty() && vec.back().first == i) {
                    ++vec.back().second;
                } else {
                    vec.emplace_back(i, 1);
                }
            }
        }
    }

    template<typename MatrixType>
    void RCCCoefficientController::getCoefficientRCC(const std::vector<SequenceInfo> &seq_info,
                                                     const std::vector<Rcc> &cuts,
                                                     bool if_elementary,
                                                     MatrixType &mat) {
        mat.resize(cuts.size(), seq_info.size());
        mat.setZero();
        std::vector<Eigen::Triplet<double> > triplets;
        triplets.reserve(cuts.size() * seq_info.size());
        if (if_elementary) {
            for (const auto &rcc: cuts) {
                if (rcc.form_rcc != static_cast<int>(RCCForm::RCC_FORM_3))
                    THROW_RUNTIME_ERROR("RCC form does not match for the elementary cols");
            }
            std::unordered_map<int, Eigen::RowVectorXi> v_col_map;
            Eigen::RowVectorXi tmp;
            for (int col = 0; col < seq_info.size(); ++col) {
                for (auto j: seq_info[col].col_seq) {
                    if (v_col_map.find(j) == v_col_map.end()) {
                        v_col_map[j] = Eigen::RowVectorXi::Zero(static_cast<int>(seq_info.size()));
                    }
                    v_col_map[j][col] = 1;
                }
            }
            tmp.resize(static_cast<int>(seq_info.size()));
            for (int idx = 0; idx < cuts.size(); ++idx) {
                const auto &rcc = cuts[idx];
                tmp.setZero();
                auto &customerInfo = rcc.info_rcc_customer;
                for (int j: customerInfo) {
                    tmp += v_col_map[j];
                }
                for (int j = 0; j < seq_info.size(); ++j) {
                    if (tmp[j] != 0) {
                        triplets.emplace_back(idx, j, 1);
                    }
                }
            }
            mat.setFromTriplets(triplets.begin(), triplets.end());
        } else {
            std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int> >, PairHasher> map_idx;
            CoefficientGetterDetail::getMap(seq_info, map_idx);
            std::vector<double> sup(seq_info.size());
            for (int idx = 0; idx < cuts.size(); ++idx) {
                const auto &rcc = cuts[idx];
                std::fill_n(sup.begin(), seq_info.size(), 0);
                if (rcc.form_rcc == static_cast<int>(RCCForm::RCC_FORM_1)) {
                    auto &customer_info = rcc.info_rcc_customer;
                    for (auto i = customer_info.begin(); i != customer_info.end(); ++i) {
                        auto j = i;
                        ++j;
                        for (; j != customer_info.end(); ++j) {
                            auto pr = *i < *j ? std::make_pair(*i, *j) : std::make_pair(*j, *i);
                            for (auto &col: map_idx[pr]) {
                                sup[col.first] += col.second;
                            }
                        }
                    }
                } else if (rcc.form_rcc == static_cast<int>(RCCForm::RCC_FORM_2)) {
                    auto &customer_info = rcc.info_rcc_customer;
                    auto &outside_customer_info = rcc.info_rcc_outside_customer;
                    for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
                        auto j = i;
                        ++j;
                        for (; j != outside_customer_info.end(); ++j) {
                            auto pr = *i < *j ? std::make_pair(*i, *j) : std::make_pair(*j, *i);
                            for (auto &col: map_idx[pr]) {
                                sup[col.first] += col.second;
                            }
                        }
                    }

                    for (auto customer_it: outside_customer_info) {
                        for (auto &col: map_idx[{0, customer_it}]) {
                            sup[col.first] += 0.5 * col.second;
                        }
                    }
                    for (auto customer_it: customer_info) {
                        for (auto &col: map_idx[{0, customer_it}]) {
                            sup[col.first] -= 0.5 * col.second;
                        }
                    }
                } else
                    THROW_RUNTIME_ERROR("RCC form does not match");

                for (int i = 0; i < seq_info.size(); ++i) {
                    if (equalFloat(sup[i], 0)) continue;
                    triplets.emplace_back(idx, i, sup[i]);
                }
            }
            mat.setFromTriplets(triplets.begin(), triplets.end());
        }
    }

    inline void RCCCoefficientController::buildRCCEnuMatrix(
        const Eigen::SparseMatrix<double> &mat, //mat only has dim-1 rows
        const std::vector<Rcc> &rccs,
        int start,
        std::vector<Eigen::Triplet<double> > &triplets
    ) {
        if (rccs.empty()) return;
        std::vector<Eigen::Triplet<double> > another_triplets;
        auto size_enumeration_col_pool = static_cast<int>(mat.cols());
        another_triplets.reserve(size_enumeration_col_pool * rccs.size());
        std::unordered_map<int, int> new_idx_map;
        new_idx_map.reserve(rccs.size());
        int cnt = 0;
        for (int c = 0; c < rccs.size(); ++c) {
            auto &rcc = rccs[c];
            if (rcc.idx_rcc < start)continue;
            new_idx_map[cnt] = c;
            for (auto &i: rcc.info_rcc_customer) {
                another_triplets.emplace_back(cnt, i - 1, 1);
            }
            ++cnt;
        }
        Eigen::SparseMatrix<double, Eigen::RowMajor> tmp(cnt, mat.rows());
        tmp.setFromTriplets(another_triplets.begin(), another_triplets.end());
        tmp = tmp * mat;

        for (int c = 0; c < cnt; ++c) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, c); it; ++it) {
                triplets.emplace_back(rccs[new_idx_map[c]].idx_rcc - start, it.col(), 1);
            }
        }
    }
}

#include "get_rcc_coefficient.hpp"

#endif // ROUTE_OPT_GET_RCC_COEFFICIENT_HPP
