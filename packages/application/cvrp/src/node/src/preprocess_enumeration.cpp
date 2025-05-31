/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include "node.hpp"
#include "node_macro.hpp"


namespace RouteOpt::Application::CVRP {
    void BbNode::preprocessEnumeration(
        const int *col_pool4_pricing,
        Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter,
        const std::vector<routeOptLong> &ng_mem4_vertex, const std::vector<double> &optional_demand,
        double optional_cap) {
        printHeadLines("Process Enumeration Matrix");
        deleteColumnByNGMemory(1, ng_mem4_vertex, true);
        for (auto &rcc: rccs)rcc.if_keep = false;

        if (!optional_demand.empty())(cleanColumnsCapInfeasible(optional_demand, optional_cap));

        deleteBranchCutsAndR1C1s();
        rank1_coefficient_getter.recoverR1CsInEnum(r1cs, cols, solver);
        RCCs::CoefficientGetter::RCCCoefficientController::recoverRCCsInEnum(rccs, cols, solver);
        generateVertex2IndexColsAndEdge2IndexCols(col_pool4_pricing, rank1_coefficient_getter);
    }

    void BbNode::deleteBranchCutsAndR1C1s() {
        if (brcs.empty() && r1cs.empty()) return;
        int num_col;
        SAFE_SOLVER(solver.getNumCol(&num_col))

        if (!brcs.empty()) {
            std::vector<int> solver_ind(num_col);
            std::vector<size_t> solver_beg(num_col + 1);
            std::vector<double> solver_val(num_col);
            std::set<int> delete_col;
            std::vector<int> ai_col(num_col);
            std::vector<int> aj_col(num_col);

            size_t numnzP;
            for (auto &brc: brcs) {
                auto ai = brc.edge.first;
                auto aj = brc.edge.second;
                if (brc.br_dir) {
                    for (int i = 0; i < num_col; ++i) {
                        ai_col[i] = 0;
                        aj_col[i] = 0;
                    }
                    if (ai) {
                        numnzP = num_col;
                        SAFE_SOLVER(solver.XgetConstraints(&numnzP,
                            solver_beg.data(),
                            solver_ind.data(),
                            solver_val.data(),
                            ai - 1,
                            1))
                        for (size_t i = 0; i < numnzP; ++i)ai_col[solver_ind[i]] = 1;
                    }
                    numnzP = num_col;
                    SAFE_SOLVER(solver.XgetConstraints(&numnzP,
                        solver_beg.data(),
                        solver_ind.data(),
                        solver_val.data(),
                        aj - 1,
                        1))
                    for (size_t i = 0; i < numnzP; ++i)aj_col[solver_ind[i]] = 1;
                    numnzP = num_col;
                    SAFE_SOLVER(solver.XgetConstraints(&numnzP,
                        solver_beg.data(),
                        solver_ind.data(),
                        solver_val.data(),
                        brc.idx_brc,
                        1))
                    for (int j = 0; j < solver_ind[0]; ++j)if (aj_col[j] || ai_col[j]) delete_col.insert(j);
                    for (size_t i = 1; i < numnzP; ++i)
                        for (int j = solver_ind[i - 1] + 1; j < solver_ind[i]; ++j)
                            if (aj_col[j] || ai_col[j])delete_col.insert(j);
                    for (int j = solver_ind[numnzP - 1] + 1; j < num_col; ++j)
                        if (aj_col[j] || ai_col[j])
                            delete_col.
                                    insert(j);
                } else {
                    //this is empty since such constraint is not imposed by adding constraint
                }
            }
            delete_col.erase(0);
            std::vector<int> col_idx(delete_col.begin(), delete_col.end());
            rmLPCols(col_idx);
            SAFE_SOLVER(solver.reoptimize(SOLVER_PRIMAL_SIMPLEX))
        }

        SAFE_SOLVER(solver.getNumCol(&num_col)) //need retrieve again
        std::vector<int> solver_ind(num_col);
        int len = 0, keep;
        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        std::vector<int> local_cstr_index(num_row);
        std::vector<int> deleted_cstrs;
        deleted_cstrs.reserve(num_row);
        std::iota(local_cstr_index.data(), local_cstr_index.data() + num_row, 0);
        for (auto &brc: brcs) {
            if (brc.idx_brc == INVALID_BRC_INDEX) continue;
            keep = brc.idx_brc;
            solver_ind[len++] = keep;
            local_cstr_index[keep] = INVALID_ROW_INDEX;
            deleted_cstrs.emplace_back(keep);
        }
        for (auto &r1c: r1cs) {
            if (r1c.info_r1c.first.size() == 1) {
                keep = r1c.idx_r1c;
                solver_ind[len++] = keep;
                local_cstr_index[keep] = INVALID_ROW_INDEX;
                deleted_cstrs.emplace_back(keep);
            }
        }
        if (deleted_cstrs.empty()) {
            return;
        }
        std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
        int delta = 0;
        auto stop_sign = deleted_cstrs.end() - 1;
        for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
            ++delta;
            for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
        }
        ++delta;
        for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;

        for (auto &rcc: rccs) rcc.idx_rcc = local_cstr_index[rcc.idx_rcc];

        for (auto i = r1cs.begin(); i < r1cs.end();) {
            if (local_cstr_index[i->idx_r1c] == INVALID_ROW_INDEX) {
                i = r1cs.erase(i);
            } else {
                i->idx_r1c = local_cstr_index[i->idx_r1c];
                ++i;
            }
        }

        for (auto &brc: brcs) brc.idx_brc = INVALID_BRC_INDEX;

        SAFE_SOLVER(solver.delConstraints(len, solver_ind.data()))
        SAFE_SOLVER(solver.reoptimize())
    }


    void BbNode::generateVertex2IndexColsAndEdge2IndexCols(const int *col_pool4_pricing,
                                                           const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter
                                                           &rank1_coefficient_getter,
                                                           bool if_rm_cols) {
        auto size_pool = static_cast<int>(index_columns_in_enumeration_column_pool.size());
        if (size_pool == 0) return;

        auto &ptr = index_columns_in_enumeration_column_pool;
        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        sparseRowMatrixXd mat(num_row, size_pool);

        std::vector<Eigen::Triplet<double> > triplets;
        triplets.reserve(
            static_cast<size_t>(static_cast<double>(num_row) * static_cast<double>(size_pool) *
                                NodeLPDensityEstimation));

        deleted_columns_in_enumeration_pool.resize(size_pool, false);

        for (int i = 0; i < size_pool; ++i) {
            for (auto j = ptr[i] + 1;; ++j) {
                int curr_node = col_pool4_pricing[j];
                if (!curr_node) break;
                triplets.emplace_back(curr_node - 1, i, 1);
            }
        }

        basic_matrix.resize(dim - 1, size_pool);
        basic_matrix.setFromTriplets(triplets.begin(), triplets.end());
        row_basic_matrix.resize(dim - 1, size_pool);
        row_basic_matrix.setFromTriplets(triplets.begin(), triplets.end());

        for (int i = 0; i < size_pool; ++i) triplets.emplace_back(dim - 1, i, 1);


        auto eps = TimeSetter::measure([&] {
            buildRCCInEnuMatrix(triplets);
        });
        printTimeMessage("obtain rcc coefficient", eps);


        eps = TimeSetter::measure([&] {
            buildAllR1CInEnuMatrix(triplets, rank1_coefficient_getter);
        });
        printTimeMessage("obtain r1c coefficient", eps);

        eps = TimeSetter::measure([&] {
            mat.setFromTriplets(triplets.begin(), triplets.end());
            matrix_in_enumeration.push_back(std::move(mat));
        });
        printTimeMessage("obtain matrix from triplets", eps);

        if (if_rm_cols) {
            eps = TimeSetter::measure([&] {
                rmColByBranchInEnuMatrix(deleted_columns_in_enumeration_pool, false, brcs, col_pool4_pricing);
            });
            printTimeMessage("delete columns by branch", eps);

            eps = TimeSetter::measure([&] {
                std::vector<double> duals(num_row, -1); //to prevent any constraints are deleted in this step;
                regenerateEnumMat(this, nullptr, false, duals);
            });
            printTimeMessage("regenerate enumeration matrix", eps);
            auto size_enumeration_col_pool = static_cast<int>(index_columns_in_enumeration_column_pool.size());
            if (size_pool != size_enumeration_col_pool)
                std::cout << "#routes = " << size_enumeration_col_pool << std::endl;
        }
    }

    void BbNode::buildRCCInEnuMatrix(
        std::vector<Eigen::Triplet<double> > &triplets, int old_num) {
        RCCs::CoefficientGetter::RCCCoefficientController::buildRCCEnuMatrix(basic_matrix, rccs, old_num, triplets);
    }

    void BbNode::buildAllR1CInEnuMatrix(
        std::vector<Eigen::Triplet<double> > &triplets,
        const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter,
        int old_num
    ) {
        auto &mat0 = matrix_in_enumeration.empty() ? row_basic_matrix : matrix_in_enumeration.front();
        rank1_coefficient_getter.buildR1CEnuMatrix(mat0, r1cs, old_num, triplets);
        if (matrix_in_enumeration.empty()) row_basic_matrix.resize(0, 0);
    }

    void BbNode::rmColByBranchInEnuMatrix(
        std::vector<bool> &deleted_columns_in_enumeration_pool,
        bool if_not_use_by_impose_br,
        const std::vector<Brc> &brcs,
        const int *col_pool4_pricing) {
        auto &mat = matrix_in_enumeration.front();
        int size_pool = static_cast<int>(index_columns_in_enumeration_column_pool.size());
        sparseRowMatrixXd tmp(1, size_pool);

        std::vector<int> must_use(brcs.size()), cannot_use(brcs.size());
        int cnt1 = 0, cnt2 = 0;
        for (int c = 0; c < brcs.size(); ++c) {
            if (brcs[c].br_dir) {
                must_use[cnt1++] = c;
            } else cannot_use[cnt2++] = c;
        }
        must_use.resize(cnt1);
        cannot_use.resize(cnt2);

        for (auto i: must_use) {
            auto &brc = brcs[i];
            int ai = brc.edge.first, aj = brc.edge.second;
            if (ai) {
                tmp = mat.row(ai - 1) + mat.row(aj - 1);
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
                    if (it.value() > 0.5) {
                        if (deleted_columns_in_enumeration_pool[it.col()]) continue;
                        if (it.value() < 1.5) {
                            deleted_columns_in_enumeration_pool[it.col()] = true;
                        } else {
                            for (auto j = index_columns_in_enumeration_column_pool[it.col()] + 1;; ++j) {
                                int current_node = col_pool4_pricing[j];
                                if (!current_node) break;
                                if (current_node == ai) {
                                    if (col_pool4_pricing[j + 1] != aj && col_pool4_pricing[j - 1] != aj)
                                        deleted_columns_in_enumeration_pool[it.col()] = true;
                                }
                            }
                        }
                    }
                }
            } else {
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, aj - 1); it; ++it) {
                    if (it.value() > 0.5 && !deleted_columns_in_enumeration_pool[it.col()]) {
                        auto j = index_columns_in_enumeration_column_pool[it.col()] + 1;
                        if (col_pool4_pricing[j] == aj) {
                            continue;
                        }
                        for (;; ++j) {
                            int current_node = col_pool4_pricing[j];
                            if (!current_node) break;
                        }
                        if (col_pool4_pricing[j - 1] != aj) {
                            deleted_columns_in_enumeration_pool[it.col()] = true;
                        }
                    }
                }
            }
        }

        if (if_not_use_by_impose_br) {
            for (auto i: cannot_use) {
                auto &brc = brcs[i];
                int ai = brc.edge.first, aj = brc.edge.second;
                if (ai) {
                    tmp = mat.row(ai - 1) + mat.row(aj - 1);
                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
                        if (it.value() > 1.5) {
                            if (deleted_columns_in_enumeration_pool[it.col()]) continue;
                            for (auto j = index_columns_in_enumeration_column_pool[it.col()] + 1;; ++j) {
                                int current_node = col_pool4_pricing[j];
                                if (!current_node) break;
                                if (current_node == ai) {
                                    if (col_pool4_pricing[j + 1] == aj || col_pool4_pricing[j - 1] == aj) {
                                        deleted_columns_in_enumeration_pool[it.col()] = true;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, aj - 1); it; ++it) {
                        if (it.value() > 0.5 && !deleted_columns_in_enumeration_pool[it.col()]) {
                            auto j = index_columns_in_enumeration_column_pool[it.col()] + 1;
                            if (col_pool4_pricing[j] == aj) {
                                deleted_columns_in_enumeration_pool[it.col()] = true;
                                continue;
                            }
                            for (;; ++j) {
                                int current_node = col_pool4_pricing[j];
                                if (!current_node) break;
                            }
                            if (col_pool4_pricing[j - 1] == aj) {
                                deleted_columns_in_enumeration_pool[it.col()] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    void BbNode::regenerateEnumMat(BbNode *node, BbNode *node2, bool if_force, std::vector<double> &duals) {
        if (node == nullptr)std::terminate();
        if (node->if_terminate) return;
        auto &deleted_columns_in_enumeration_pool =
                node2 ? node2->deleted_columns_in_enumeration_pool : node->deleted_columns_in_enumeration_pool;
        auto size_pool = static_cast<int>(node->index_columns_in_enumeration_column_pool.size());
        if (!size_pool) {
            node->deleted_columns_in_enumeration_pool.clear();
            return;
        }
        BbNode *out_node;
        int del_size = 0;
        if (!node2) {
            out_node = node;
            for (int i = 0; i < size_pool; ++i) {
                if (deleted_columns_in_enumeration_pool[i]) ++del_size;
            }
            if (del_size == 0) return;
            node->valid_size = size_pool - del_size;
            if (!if_force) {
                int num_row;
                SAFE_SOLVER(node->solver.getNumRow(&num_row))
                auto ratio_1 = static_cast<double>(node->valid_size) / size_pool;
                auto ratio_2 = static_cast<double>(node->countActiveCuts(duals)) / num_row;
                if (ratio_1 * ratio_2 > LeftThresholdRCFixing4EnumerationPool) {
                    std::cout << "columns pending deletion= " << del_size << " remain= " << size_pool << " col ratio= "
                            << ratio_1 << " row ratio= " << ratio_2 << std::endl;
                    return;
                }
                std::cout << "delete columns= " << del_size << " remain= " << node->valid_size << std::endl;
            }
        } else {
            out_node = node2;
            for (int i = 0; i < size_pool; ++i) {
                if (deleted_columns_in_enumeration_pool[i]) ++del_size;
            }
            node2->valid_size = size_pool - del_size;
        }

        int colIndex = 0;
        std::vector<int> new_col_map(size_pool, -1);
        for (int i = 0; i < size_pool; ++i) {
            if (!deleted_columns_in_enumeration_pool[i]) {
                new_col_map[i] = colIndex++;
            }
        }

        int num_row;
        SAFE_SOLVER(node->solver.getNumRow(&num_row))
        std::vector<int> oper_cstr_index;
        node->findNonActiveCuts(duals, oper_cstr_index);
        SAFE_SOLVER(node->solver.reoptimize(SOLVER_DUAL_SIMPLEX));
        int new_size_pool = size_pool - del_size;
        SAFE_SOLVER(node->solver.getNumRow(&num_row))
        sparseRowMatrixXd tmpMatrix(num_row, new_size_pool);

        int rowIndex = 0;
        std::vector<Eigen::Triplet<double> > triplets;
        triplets.reserve(
            static_cast<size_t>(static_cast<double>(num_row) * static_cast<double>(new_size_pool) *
                                NodeLPDensityEstimation));

        for (auto &it: node->matrix_in_enumeration) {
            for (int i = 0; i < it.rows(); ++i) {
                if (oper_cstr_index[rowIndex] != INVALID_ROW_INDEX) {
                    int j = oper_cstr_index[rowIndex];
                    if (j >= num_row)
                        THROW_RUNTIME_ERROR(
                        "row index out of range!" + std::to_string(j) + " " + std::to_string(num_row));

                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_it(it, i); inner_it; ++
                         inner_it) {
                        if (new_col_map[inner_it.col()] != -1) {
                            triplets.emplace_back(j, new_col_map[inner_it.col()], inner_it.value());
                        }
                    }
                }
                ++rowIndex;
            }
        }

        SAFE_EIGEN(tmpMatrix.setFromTriplets(triplets.begin(), triplets.end());)

        out_node->matrix_in_enumeration.clear();
        out_node->matrix_in_enumeration.push_back(std::move(tmpMatrix));

        colIndex = 0;
        auto size_enumeration_col_pool = static_cast<int>(out_node->index_columns_in_enumeration_column_pool.size());
        for (int i = 0; i < size_enumeration_col_pool; ++i) {
            if (!out_node->deleted_columns_in_enumeration_pool[i]) {
                out_node->cost_for_columns_in_enumeration_column_pool[colIndex] =
                        out_node->cost_for_columns_in_enumeration_column_pool[i];
                out_node->index_columns_in_enumeration_column_pool[colIndex++] =
                        out_node->index_columns_in_enumeration_column_pool[i];
            }
        }

        out_node->cost_for_columns_in_enumeration_column_pool.conservativeResize(colIndex);
        out_node->index_columns_in_enumeration_column_pool.conservativeResize(colIndex);
        out_node->deleted_columns_in_enumeration_pool.assign(colIndex, false);


        out_node->createBasicMatrix();
    }
}
