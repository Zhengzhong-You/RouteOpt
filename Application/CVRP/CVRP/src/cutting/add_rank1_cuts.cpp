//
// Created by You, Zhengzhong on 1/10/24.
//


#include "cvrp.hpp"
#include "get_rank1_matrix.hpp"

using namespace std;


void CVRP::addLimitedMemoryR1CsNodeBased(BbNode *node,
                                         const std::vector<R1c> &full_cuts) {
    vector<int> idx;
    giveMemInNode(node, full_cuts, idx);
    addR1CAtOnce(node, idx);
}


void CVRP::giveMemInNode(BbNode *node,
                         const std::vector<R1c> &full_cuts,
                         std::vector<int> &idx) const {
    idx.clear();
    int numRow = num_row;
    unordered_map<int, int> map_rowIndex_cutIndex;
    for (int i = 0; i < node->r1cs.size(); ++i) {
        map_rowIndex_cutIndex[node->r1cs[i].idx_r1c] = i;
    }
    for (const auto &full_cut: full_cuts) {
        int cut_index = full_cut.idx_r1c;

        if (cut_index == INITIAL_IDX_R1C) {
            R1c r1c = full_cut;
            r1c.idx_r1c = numRow++;
            Rank1CutsSeparator::getMapRhs(r1c.rhs, static_cast<int>(r1c.info_r1c.first.size()), r1c.info_r1c.second);
            node->r1cs.emplace_back(r1c);
            idx.emplace_back(static_cast<int>(node->r1cs.size()) - 1);
        } else {
            cut_index = map_rowIndex_cutIndex.at(cut_index);
            auto &r1c = node->r1cs[cut_index];
            r1c.arc_mem = full_cut.arc_mem;
            idx.emplace_back(cut_index);
        }
    }
}

void CVRP::addR1CAtOnce(BbNode *node,
                        const std::vector<int> &idx) {
    getLimitedR1CPre(node, idx);

    sparseRowMatrixXd mat;
    getLimitedR1CCoeffs(node->getCols(), mat);

    vector<int> solver_ind(num_col);
    vector<double> solver_val(num_col);
    int old_num_row = num_row;
    for (int i = 0; i < idx.size(); ++i) {
        int numnz = 0;
        int row_index = node->r1cs[idx[i]].idx_r1c;
        double rhs = node->r1cs[idx[i]].rhs;
        double coeff = mat.coeff(i, 0);
        if (abs(coeff - rhs) > TOLERANCE) {
            if (abs(coeff) > TOLERANCE) {
                solver_ind[numnz] = 0;
                solver_val[numnz] = rhs;
                ++numnz;
                sparseRowMatrixXd::InnerIterator it(mat, i);
                ++it;
                for (; it; ++it) {
                    solver_ind[numnz] = (int) it.col();
                    solver_val[numnz] = it.value();
                    ++numnz;
                }
            } else {
                solver_ind[numnz] = 0;
                solver_val[numnz] = rhs;
                ++numnz;
                goto HERE;
            }
        } else {
        HERE:
            for (sparseRowMatrixXd::InnerIterator it(mat, i); it; ++it) {
                solver_ind[numnz] = (int) it.col();
                solver_val[numnz] = it.value();
                ++numnz;
            }
        }
        if (row_index >= old_num_row) {
            safe_solver(node->getSolver().addConstraint(numnz,
                solver_ind.data(),
                solver_val.data(),
                SOLVER_LESS_EQUAL,
                rhs,
                nullptr))
            int vind = 0;
            safe_solver(node->getSolver().changeCoeffs(1, &row_index, &vind, &rhs))
            ++num_row;
        } else {
            vector<int> solver_ind2(numnz);
            fill(solver_ind2.begin(), solver_ind2.end(), row_index);
            safe_solver(node->getSolver().XchangeCoeffs(numnz, solver_ind2.data(), solver_ind.data(), solver_val.data()))
        }
    }
}

void CVRP::addLimitedMemoryR1Cs(BbNode *node,
                                std::vector<R1c> &full_cuts) {
    vector<int> idx;
    writeIntoNode(node, full_cuts, idx);
    addR1CAtOnce(node, idx);
}

void CVRP::writeIntoNode(BbNode *node, vector<R1c> &full_cuts, vector<int> &idx) const {
    idx.clear();
    int numRow = num_row;
    unordered_map<int, int> map_rowIndex_cutIndex;
    for (int i = 0; i < node->r1cs.size(); ++i) {
        map_rowIndex_cutIndex[node->r1cs[i].idx_r1c] = i;
    }
    for (auto &cut: full_cuts) {
        int cut_index = cut.idx_r1c;
        if (cut_index == INITIAL_IDX_R1C) {
            cut.idx_r1c = numRow++;
            Rank1CutsSeparator::getMapRhs(cut.rhs, static_cast<int>(cut.info_r1c.first.size()), cut.info_r1c.second);
            node->r1cs.emplace_back(cut);
            idx.emplace_back(static_cast<int>(node->r1cs.size()) - 1);
        } else {
            cut_index = map_rowIndex_cutIndex.at(cut_index);
            auto &r1c = node->r1cs[cut_index];
            r1c.arc_mem = cut.arc_mem;
            idx.emplace_back(cut_index);
        }
    }
}

void CVRP::addR1CAtOnceInEnum(BbNode *node,
                              const std::vector<R1c> &new_cuts) {
    size_t num_nz, numnzP;
    vector<int> ai_col(num_col);
    vector<size_t> solver_beg(2);
    vector<int> solver_ind(num_col);
    vector<double> solver_val(num_col);
    for (auto &c: new_cuts) {
        auto &cut = c.info_r1c;
        num_nz = 0;
        memset(ai_col.data(), 0, sizeof(int) * num_col);
        vector<int> multi;
        int denominator, rhs;
        Rank1CutsSeparator::getMapPlanInfo(multi, denominator, rhs, static_cast<int>(cut.first.size()), cut.second);
        int count = 0;
        for (auto i: cut.first) {
            numnzP = num_col;
            safe_solver(node->getSolver().XgetConstraints(&numnzP,
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
        if (rhs > 0.1) {
            solver_ind[num_nz] = 0;
            solver_val[num_nz++] = rhs;
        }
        for (int i = 1; i < num_col; ++i) {
            int tmp = int(ai_col[i] / denominator);
            if (tmp) {
                solver_ind[num_nz] = i;
                solver_val[num_nz++] = tmp;
            }
        }

        safe_solver(node->getSolver().addConstraint((int)num_nz,
            solver_ind.data(),
            solver_val.data(),
            SOLVER_LESS_EQUAL,
            rhs,
            nullptr))
        R1c r1c;
        r1c.info_r1c = cut;
        r1c.idx_r1c = num_row++;
        r1c.rhs = rhs;
        node->r1cs.emplace_back(r1c);
    }
}




