/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include "node.hpp"

namespace RouteOpt::Application::CVRP {
    void BbNode::rmLPCols(const std::vector<int> &col_idx) {
        if (col_idx.empty()) return;
        std::vector<int> sort_col_idx = col_idx;
        if (!is_sorted(col_idx.begin(), col_idx.end())) {
            std::sort(sort_col_idx.begin(), sort_col_idx.end());
        }

        int num_col;
        SAFE_SOLVER(solver.getNumCol(&num_col))
        int delta = 0;
        auto stop_sign = sort_col_idx.end() - 1;
        for (auto i = sort_col_idx.begin(); i < stop_sign; ++i) {
            ++delta;
            for (int j = *i + 1; j < *(i + 1); ++j) cols[j - delta] = cols[j];
        }
        ++delta;
        for (int j = *stop_sign + 1; j < num_col; ++j) cols[j - delta] = cols[j];

        SAFE_SOLVER(solver.delVars(sort_col_idx.size(), sort_col_idx.data()))
        SAFE_SOLVER(solver.updateModel())
        SAFE_SOLVER(solver.getNumCol(&num_col))
        cols.resize(num_col);
    }


    void BbNode::cleanIndexColForNode() {
        std::vector<int> col_idx;
        int beg = 1;
        SAFE_SOLVER(solver.reoptimize(SOLVER_PRIMAL_SIMPLEX))
        int num_col;
        SAFE_SOLVER(solver.getNumCol(&num_col));
        std::vector<double> rc(num_col);
        SAFE_SOLVER(solver.getRC(0, num_col, rc.data()))
        std::vector<double> rc_copy(rc.data() + beg, rc.data() + num_col);
        auto n = static_cast<int>((num_col - beg) * COL_KEEP_FRAC);
        nth_element(rc_copy.begin(), rc_copy.begin() + n, rc_copy.end());
        double threshold = std::max(rc_copy[n], TOLERANCE);
        for (int i = beg; i < num_col; ++i) if (rc[i] > threshold) col_idx.emplace_back(i);
        rmLPCols(col_idx);
        SAFE_SOLVER(solver.reoptimize(SOLVER_PRIMAL_SIMPLEX))
    }

    void BbNode::cleanIndexColForNode(double ub, std::vector<double> &optimal_dual_vector) {
        std::vector<int> col_idx;
        int beg = 1; //rc fixing can safely delete the columns
        SAFE_SOLVER(solver.updateModel())
        int num_col;
        SAFE_SOLVER(solver.getNumCol(&num_col))
        auto opt_gap = calculateOptimalGap(ub);
        size_t numnzP;
        SAFE_SOLVER(solver.XgetVars(&numnzP, nullptr, nullptr, nullptr, 0, num_col))
        std::vector<size_t> vbeg(num_col + 1);
        std::vector<int> vind(numnzP);
        std::vector<double> vval(numnzP);
        SAFE_SOLVER(solver.XgetVars(&numnzP, vbeg.data(), vind.data(), vval.data(), 0, num_col))
        std::vector<Eigen::Triplet<double> > triplet(numnzP);
        vbeg[num_col] = numnzP;
        numnzP = 0;
        for (int i = 0; i < num_col; ++i) {
            for (auto j = vbeg[i]; j < vbeg[i + 1]; ++j) {
                triplet[numnzP++] = Eigen::Triplet<double>(vind[j], i, vval[j]);
            }
        }
        auto num_row = static_cast<int>(optimal_dual_vector.size());
        sparseColMatrixXd A(num_row, num_col);
        A.setFromTriplets(triplet.begin(), triplet.end());
        Eigen::Map<Eigen::RowVectorXd> pi(optimal_dual_vector.data(), num_row);
        std::vector<double> cost(num_col);
        SAFE_SOLVER(solver.getObj(0, num_col, cost.data()))
        Eigen::Map<Eigen::RowVectorXd> c(cost.data(), num_col);
        RowVectorXd rc;
        SAFE_EIGEN(rc = c - pi * A;)
        for (int i = beg; i < num_col; ++i) if (rc[i] > opt_gap) col_idx.emplace_back(i);
        rmLPCols(col_idx);
        SAFE_SOLVER(solver.reoptimize(SOLVER_PRIMAL_SIMPLEX))
    }

    void BbNode::deleteColumnByNGMemory(int start, const std::vector<routeOptLong> &ng_mem4_vertex, bool if_full_mem) {
        std::vector<int> col_idx;

        for (int i = start; i < cols.size(); ++i) {
            routeOptLong local_pi = 0;
            for (auto current_node: getCols()[i].col_seq) {
                if (local_pi[current_node]) {
                    col_idx.emplace_back(i);
                    break;
                }
                if (!if_full_mem)local_pi = local_pi & ng_mem4_vertex[current_node];
                local_pi.set(current_node);
            }
        }
        rmLPCols(col_idx);
        SAFE_SOLVER(solver.reoptimize())
        double lp_val;
        SAFE_SOLVER(solver.getObjVal(&lp_val))
        if (if_full_mem) {
            std::cout << "delete non-ele routes lp val= " << lp_val << std::endl;
        } else {
            std::cout << "delete ng-infeasible routes lp val= " << lp_val << std::endl;
        }
    }

    void BbNode::cleanColumnsCapInfeasible(const std::vector<double> &demand, double cap) {
        /**
         * only focus on the lp, the enumeration pool is considered elsewhere
         */

        std::vector<int> col_idx;

        for (int i = 1; i < cols.size(); ++i) {
            double res = 0.;
            for (auto current_node: cols[i].col_seq) {
                res += demand[current_node];
            }
            if (res > cap + TOLERANCE) {
                col_idx.emplace_back(i);
            }
        }
        rmLPCols(col_idx);
        SAFE_SOLVER(solver.reoptimize(SOLVER_PRIMAL_SIMPLEX))
        double lp_val;
        SAFE_SOLVER(solver.getObjVal(&lp_val))
        std::cout << "after delete " << col_idx.size() << " res-cap infeasible columns, lp val= " << lp_val <<
                std::endl;
    }
}
