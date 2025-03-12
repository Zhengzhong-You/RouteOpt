/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include "node.hpp"

namespace RouteOpt::Application::CVRP {
    void BbNode::deleteCuts(std::vector<int> &deleted_cstrs, std::vector<int> &local_cstr_index) {
        std::sort(deleted_cstrs.begin(), deleted_cstrs.end());
        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        local_cstr_index.resize(num_row);
        std::iota(local_cstr_index.begin(), local_cstr_index.end(), 0);
        if (deleted_cstrs.empty()) return;

        for (auto &i: deleted_cstrs) local_cstr_index[i] = INVALID_ROW_INDEX;
        auto stop_sign = deleted_cstrs.end() - 1;
        int delta = 0;
        for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
            ++delta;
            for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
        }
        ++delta;
        for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;

        SAFE_SOLVER(solver.delConstraints(static_cast<int>(deleted_cstrs.size()), deleted_cstrs.data()))
        SAFE_SOLVER(solver.updateModel())

        for (auto i = rccs.begin(); i < rccs.end();) {
            if (local_cstr_index[i->idx_rcc] == INVALID_ROW_INDEX) {
                i = rccs.erase(i);
            } else {
                i->idx_rcc = local_cstr_index[i->idx_rcc];
                ++i;
            }
        }

        for (auto i = r1cs.begin(); i < r1cs.end();) {
            if (local_cstr_index[i->idx_r1c] == INVALID_ROW_INDEX) {
                i = r1cs.erase(i);
            } else {
                i->idx_r1c = local_cstr_index[i->idx_r1c];
                ++i;
            }
        }

        for (auto i = brcs.begin(); i < brcs.end(); ++i) {
            if (i->idx_brc == INVALID_BRC_INDEX) continue;
            i->idx_brc = local_cstr_index[i->idx_brc];
        }
    }

    int BbNode::countActiveCuts(const std::vector<double> &optimal_dual) const {
        if (if_terminate)
            THROW_RUNTIME_ERROR("cannot count active cuts in a terminated node");
        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        if (optimal_dual.empty()) return num_row;
        if (optimal_dual.size() != num_row)
            THROW_RUNTIME_ERROR("optimal dual std::vector size is not equal to num_row");

        int cnt = 0;

        for (auto &rcc: rccs) {
            if (rcc.if_keep) continue;
            if (int idx = rcc.idx_rcc; std::abs(optimal_dual[idx]) < DUAL_TOLERANCE) {
                ++cnt;
            }
        }

        for (auto &r1c: r1cs) {
            if (int idx = r1c.idx_r1c; std::abs(optimal_dual[idx]) < DUAL_TOLERANCE) {
                ++cnt;
            }
        }
        return num_row - cnt;
    }

    void BbNode::findNonActiveCuts(std::vector<double> &optional_optimal_dual, std::vector<int> &cstr_index) {
        if (if_terminate) return;
        std::vector<int> nonactive_cuts;
        std::vector<double> current_dual_vec;
        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        nonactive_cuts.reserve(num_row);
        auto &dual_vec = optional_optimal_dual.empty() ? current_dual_vec : optional_optimal_dual;
        if (!optional_optimal_dual.empty() && optional_optimal_dual.size() != num_row)
            THROW_RUNTIME_ERROR("optimal dual std::vector size is not equal to num_row");

        if (optional_optimal_dual.empty()) {
            SAFE_SOLVER(solver.reoptimize())
            current_dual_vec.resize(num_row);
            SAFE_SOLVER(solver.getDual(0, num_row, current_dual_vec.data()))
        }

        for (auto &rcc: rccs) {
            if (rcc.if_keep) continue;
            int idx = rcc.idx_rcc;
            if (std::abs(dual_vec[idx]) < DUAL_TOLERANCE) {
                nonactive_cuts.emplace_back(idx);
            }
        }

        for (auto &r1c: r1cs) {
            int idx = r1c.idx_r1c;
            if (std::abs(dual_vec[idx]) < DUAL_TOLERANCE) {
                nonactive_cuts.emplace_back(idx);
            }
        }

        cstr_index.clear();
        deleteNonactiveCuts(nonactive_cuts, optional_optimal_dual, cstr_index);
    }

    void BbNode::deleteNonactiveCuts(std::vector<int> &nonactive_cuts, std::vector<double> &optimal_dual_vector,
                                     std::vector<int> &cstr_index) {
        // int num_row;
        // SAFE_SOLVER(solver.getNumRow(&num_row))
        deleteCuts(nonactive_cuts, cstr_index);
        if (nonactive_cuts.empty()) return;
        if (!optimal_dual_vector.empty()) {
            int ccnt = 0;
            for (int i = 0; i < optimal_dual_vector.size(); ++i) {
                if (cstr_index[i] == INVALID_ROW_INDEX) continue;
                optimal_dual_vector[ccnt++] = optimal_dual_vector[i];
            }
            optimal_dual_vector.resize(ccnt);
        }
    }

    bool BbNode::validateCuts(int old_num, bool if_care_lb_improvement) {
        /**
         * since this mode only delete cuts by slack and only new cuts will be tested,
         * therefore, the rcc can be safely deleted
         */
        SAFE_SOLVER(solver.reoptimize())
        bool if_continue = true;
        if (if_care_lb_improvement) {
            double lp_val;
            SAFE_SOLVER(solver.getObjVal(&lp_val))
            if (lp_val - value < CUTTING_BRANCHING_RATIO * br_value_improved) {
                if_continue = false;
                std::cout << "minor (" << lp_val - value << ") improvement after cutting, stop this cutting iteration"
                        <<
                        std::endl;
                std::cout << SMALL_PHASE_SEPARATION;
            }
        }
        int num_row;
        SAFE_SOLVER(solver.getNumRow(&num_row))
        std::vector<double> slack(num_row);
        SAFE_SOLVER(solver.getSlack(0, num_row, slack.data()))
        std::vector<int> deleted_cstrs;

        if (!if_continue) {
            deleted_cstrs.resize(num_row - old_num);
            std::iota(deleted_cstrs.begin(), deleted_cstrs.end(), old_num);
        } else {
            for (int i = old_num; i < num_row; ++i) {
                if (std::abs(slack[i]) > TOLERANCE) {
                    deleted_cstrs.emplace_back(i);
                }
            }
        }

        if (deleted_cstrs.empty()) goto QUIT; {
            std::vector<int> local_cstr_index;
            deleteCuts(deleted_cstrs, local_cstr_index);
        }

    QUIT:
        return if_continue;
    }

    void BbNode::mergeR1Cs() {
        if (if_in_enu_state) return;
        if (r1cs.empty()) return;
        std::set<std::pair<std::vector<int>, int> > r1c_info_set;
        std::vector<int> cind;
        for (auto it = r1cs.rbegin(); it != r1cs.rend(); ++it) {
            if (r1c_info_set.find({it->info_r1c}) != r1c_info_set.end()) {
                cind.emplace_back(it->idx_r1c);
            } else {
                r1c_info_set.emplace(it->info_r1c);
            }
        }
        if (cind.empty()) return;
        std::vector<int> local_cstr_index;
        deleteCuts(cind, local_cstr_index);
    }
}
