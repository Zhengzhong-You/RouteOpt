/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "add_columns_macro.hpp"
#include "add_column_controller.hpp"
#include "rcc_coefficient_controller.hpp"
#include "node_macro.hpp"

namespace RouteOpt::Application::CVRP {
    namespace AddColumnDetail {
        void checkRC(
            const std::vector<R1c> &r1cs,
            const std::vector<std::tuple<Label *, Label *, double> > &negative_rc_label_tuple,
            const std::vector<SequenceInfo> &new_cols,
            const RowVectorXd &rc,
            const RowVectorXd &pi4_labeling,
            const sparseColMatrixXd &mat,
            int num_row) {
            bool if_wrong_pricing_column = false;
            for (int i = 0; i < negative_rc_label_tuple.size(); ++i) {
                auto &label_tuple = negative_rc_label_tuple[i];
                auto rc_n = std::get<2>(label_tuple);
                if (std::abs(rc_n - rc(i)) > TOLERANCE) {
                    std::cout << "rc_n= " << rc_n << " rc= " << rc(i) << std::endl;
                    auto rc_dif = rc_n - rc(i);
                    for (int j = 0; j < num_row; ++j) {
                        if (std::abs(pi4_labeling[j] - rc_dif) < TOLERANCE) {
                            std::cout << "pi4_labeling[j]= " << pi4_labeling[j] << " rc_dif= " << rc_dif << std::endl;
                            std::cout << "j= " << j << std::endl;
                            for (auto &r1c: r1cs) {
                                if (r1c.idx_r1c == j) {
                                    std::cout << " arc_mem= ";
                                    for (auto &it: r1c.arc_mem) {
                                        // for (auto &it2: it.first) std::cout << it2 << " ";
                                        std::cout << " to " << it.second << std::endl;
                                    }
                                    std::cout << " info= ";
                                    for (auto &it: r1c.info_r1c.first)std::cout << it << " ";
                                    std::cout << " | snd= " << r1c.info_r1c.second << std::endl;

                                    std::cout << std::endl;
                                    std::cout << " col seq= ";
                                    for (auto &it: new_cols[i].col_seq)std::cout << it << " ";
                                    std::cout << std::endl;
                                    std::cout << "col forward_concatenate_pos= " << new_cols[i].forward_concatenate_pos
                                            <<
                                            std::endl;
                                    std::cout << " rc col coeff= " << mat.coeff(j, i) << std::endl;
                                    if_wrong_pricing_column = true;
                                }
                            }
                        }
                    }
                }
            }

            if (if_wrong_pricing_column) {
                std::cout << "col size= " << rc.size() << std::endl;
                THROW_RUNTIME_ERROR("wrong pricing column");
            }
        }
    }

    void CVRP_AddColumnController::calculateColumnCoefficientsB4Enumeration(bool if_enu,
                                                                            sparseColMatrixXd &mat,
                                                                            Eigen::RowVectorXd &cost,
                                                                            std::unordered_map<std::pair<int, int>,
                                                                                std::vector<int>
                                                                                , PairHasher> &
                                                                            edge_map) {
        int ccnt_cnt;
        int past_node;
        double cost_sum;

        int ccnt = static_cast<int>(new_cols_ref.get().size());
        edge_map.clear();
        edge_map.reserve(dim * dim);
        mat.resize(num_row, ccnt);
        mat.setZero();
        cost.resize(ccnt);

        std::unordered_map<std::pair<int, int>, double, PairHasher> coeff_map;
        coeff_map.reserve(dim * dim);

        std::unordered_map<int, int> single_size_route; //customer, col idx

        ccnt_cnt = 0;
        for (auto &c: new_cols_ref.get()) {
            auto &col = c.col_seq;
            if (col.size() == 1) single_size_route[col[0]] = ccnt_cnt;
            past_node = 0;
            cost_sum = 0;
            for (int curr_node: col) {
                cost_sum += cost_mat4_vertex_ref.get()[past_node][curr_node];
                ++coeff_map[{curr_node - 1, ccnt_cnt}];
                auto pr = past_node < curr_node
                              ? std::make_pair(past_node, curr_node)
                              : std::make_pair(curr_node, past_node);
                edge_map[pr].emplace_back(ccnt_cnt);
                past_node = curr_node;
            }
            edge_map[std::make_pair(0, past_node)].emplace_back(ccnt_cnt);
            cost_sum += cost_mat4_vertex_ref.get()[past_node][0];
            cost(ccnt_cnt) = cost_sum;
            ++ccnt_cnt;
        }


        sparseColMatrixXd mat_rcc;
        RCCs::CoefficientGetter::RCCCoefficientController::getCoefficientRCC(
            new_cols_ref.get(), *rccs_ptr, if_enu, mat_rcc);

        for (auto &br: *brcs_ptr) {
            int idx = br.idx_brc;
            if (idx == INVALID_BRC_INDEX) continue;
            int ai = br.edge.first;
            int aj = br.edge.second;
            auto pr = ai < aj ? std::make_pair(ai, aj) : std::make_pair(aj, ai);
            for (auto it_map: edge_map[pr]) ++coeff_map[{idx, it_map}];
            if (ai == 0 && single_size_route.find(aj) != single_size_route.end()) {
                coeff_map[{idx, single_size_route[aj]}] = 1; //force to be 1
            }
        }

        sparseColMatrixXd mat_r1c;
        rank1_coefficient_getter_ref.get().getR1CCoeffs(new_cols_ref.get(), *r1cs_ptr, nullptr, !if_enu, mat_r1c);

        size_t n = mat_rcc.nonZeros() + coeff_map.size() + mat_r1c.nonZeros() + ccnt;
        std::vector<Eigen::Triplet<double> > triplet(n);
        n = 0;

        for (auto &it: coeff_map) {
            triplet[n++] = {it.first.first, it.first.second, it.second};
        }

        std::vector<int> lp_r1c_map(r1cs_ptr->size()), lp_rcc_map(rccs_ptr->size());
        std::transform(r1cs_ptr->begin(), r1cs_ptr->end(), lp_r1c_map.begin(),
                       [](const R1c &r1c) { return r1c.idx_r1c; });
        std::transform(rccs_ptr->begin(), rccs_ptr->end(), lp_rcc_map.begin(),
                       [](const Rcc &rcc) { return rcc.idx_rcc; });

        for (int i = 0; i < mat_r1c.outerSize(); ++i) {
            for (sparseColMatrixXd::InnerIterator it(mat_r1c, i); it; ++it) {
                triplet[n++] = {lp_r1c_map[it.row()], static_cast<int>(it.col()), it.value()};
            }
            for (sparseColMatrixXd::InnerIterator it(mat_rcc, i); it; ++it) {
                triplet[n++] = {lp_rcc_map[it.row()], static_cast<int>(it.col()), it.value()};
            }
        }

        //  remove_vehicle_constraint_call()
        auto real_dim = dim - 1;
        for (int i = 0; i < ccnt; ++i) {
            triplet[n++] = {real_dim, i, 1};
        }
        triplet.resize(n);

        SAFE_EIGEN(mat.setFromTriplets(triplet.begin(), triplet.end());)
    }

    void CVRP_AddColumnController::addColumns(int &ccnt, const std::vector<double> &pi4_labeling, bool if_check_rc,
                                              bool if_enu) {
        if (ccnt == 0) return;

        SAFE_SOLVER(solver_ptr->updateModel())
        SAFE_SOLVER(solver_ptr->getNumRow(&num_row))

        sparseColMatrixXd mat;
        Eigen::RowVectorXd cost;
        std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> edge_map;

        calculateColumnCoefficientsB4Enumeration(if_enu, mat, cost, edge_map);

        RowVectorXd rc;
        if (if_check_rc) {
            RowVectorXd local_pi(num_row);
            std::transform(pi4_labeling.begin(), pi4_labeling.end(), local_pi.data(),
                           [](double pi) { return pi; });
            SAFE_EIGEN(rc = cost - local_pi * mat)
            if constexpr (CHECK_RC_EVERY_COLUMN) {
                PRINT_DEBUG("check RC for every column");
                AddColumnDetail::checkRC(*r1cs_ptr, negative_rc_label_tuple_ref.get(), new_cols_ref.get(), rc, local_pi,
                                         mat, num_row);
            }
        } else {
            rc.resize(ccnt);
            for (int i = 0; i < ccnt; ++i) rc(i) = -1;
        }
        int ccnt_cnt = 0;
        size_t nzcnt = 0;
        std::vector<size_t> solver_beg(ccnt + 1);
        std::vector<int> solver_ind(num_row * ccnt);
        std::vector<double> solver_val(num_row * ccnt);
        std::vector<double> solver_obj(ccnt);
        for (int i = 0; i < ccnt; ++i) {
            if (rc(i) < RC_TOLERANCE) {
                solver_beg[ccnt_cnt] = nzcnt;
                solver_obj[ccnt_cnt] = cost(i);
                ++ccnt_cnt;
                for (sparseColMatrixXd::InnerIterator it(mat, i); it; ++it) {
                    solver_ind[nzcnt] = static_cast<int>(it.row());
                    solver_val[nzcnt++] = it.value();
                }
                existing_cols_ptr->emplace_back(new_cols_ref.get()[i]);
            } else if (rc(i) > -Rank1Cuts::MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX * dim * RC_TOLERANCE) {
                auto &col = new_cols_ref.get()[i];
                for (auto j: col.col_seq) std::cout << j << " ";
                std::cout << " | " << col.forward_concatenate_pos << std::endl;
                std::cout << "col " << i << " is not allowed! The rc= " << rc(i) << std::endl;
                THROW_RUNTIME_ERROR("wrong pricing column");
            }
            if constexpr (CHECK_PRICING_LABELS) {
                if (!equalFloat(seq_rc_map_ref.get().at(new_cols_ref.get()[i].col_seq), rc(i))) {
                    std::cout << "seq_rc= " << seq_rc_map_ref.get().at(new_cols_ref.get()[i].col_seq) << " rc= " <<
                            rc(i)
                            <<
                            std::endl;
                    THROW_RUNTIME_ERROR("seq_rc is not equal to rc");
                }
            }
        }

        if constexpr (CHECK_PRICING_LABELS) {
            seq_rc_map_ref.get().clear();
        }

        solver_beg[ccnt_cnt] = nzcnt;

        ccnt = ccnt_cnt;
        if (!ccnt) return;

        SAFE_SOLVER(solver_ptr->XaddVars(ccnt_cnt,
            nzcnt,
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            solver_obj.data(),
            nullptr,
            nullptr,
            nullptr,
            nullptr))
        SAFE_SOLVER(solver_ptr->updateModel())
    }


    void CVRP_AddColumnController::addColumnsByInspection(const std::vector<int> &Col_added) {
        if (Col_added.empty()) return;
        auto &ptr = *idx_col_pool_ptr;
        auto &ptr_cost = *cost_col_pool_ptr;
        auto &mat = *matrix_col_pool_ptr;
        auto &cols = *existing_cols_ptr;
        auto size_enumeration_col_pool = static_cast<int>(ptr.size());

        // int idx = num_col;
        auto idx = cols.size();
        cols.resize(cols.size() + Col_added.size());
        for (auto &col: Col_added) {
            auto &seq = cols[idx++].col_seq; //no other information is needed!
            for (auto j = ptr[col] + 1;; ++j) {
                int curr_node = col_pool4_pricing[j];
                if (!curr_node)break;
                seq.emplace_back(curr_node);
            }
        }

        std::vector<int> if_col_added(size_enumeration_col_pool, -1);
        int cnt = 0;
        for (auto col: Col_added) {
            if_col_added[col] = cnt++;
        }

        SAFE_SOLVER(solver_ptr->getNumRow(&num_row))
        Eigen::SparseMatrix<double, Eigen::ColMajor> tmp_mat(num_row, (int) Col_added.size());
        std::vector<Eigen::Triplet<double> > triplets;
        triplets.reserve(
            static_cast<size_t>(static_cast<double>(Col_added.size()) * num_row * NodeLPDensityEstimation));
        size_t num = 0;
        for (auto &it: mat) {
            for (int i = 0; i < it.rows(); ++i) {
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_it(it, i); inner_it; ++
                     inner_it) {
                    if (if_col_added[inner_it.col()] != -1) {
                        triplets.emplace_back(num, if_col_added[inner_it.col()], inner_it.value());
                    }
                }
                ++num;
            }
        }

        SAFE_EIGEN(tmp_mat.setFromTriplets(triplets.begin(), triplets.end());)
        auto nonZeros = tmp_mat.nonZeros();
        std::vector<double> solver_obj(Col_added.size()), solver_val(nonZeros);
        std::vector<size_t> solver_beg(Col_added.size() + 1);
        std::vector<int> solver_ind(nonZeros);
        int ccnt = 0;
        size_t nz = 0;
        for (auto col: Col_added) {
            solver_obj[ccnt] = ptr_cost[col];
            solver_beg[ccnt] = nz;
            for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(tmp_mat, ccnt); it; ++it) {
                solver_ind[nz] = static_cast<int>(it.row());
                solver_val[nz++] = it.value();
            }
            ++ccnt;
        }
        solver_beg[ccnt] = nz;
        SAFE_SOLVER(solver_ptr->XaddVars(ccnt,
            solver_ind.size(),
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            solver_obj.data(),
            nullptr,
            nullptr,
            nullptr,
            nullptr))
        SAFE_SOLVER(solver_ptr->updateModel())
    }
}
