/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp.hpp"
#include "node.hpp"


namespace RouteOpt::Application::CVRP {
    void BbNode::solveMIP(double ub) {
        int num_col;
        SAFE_SOLVER(solver.getNumCol(&num_col))
        std::vector<char> xtype(num_col, SOLVER_BINARY);
        SAFE_SOLVER(solver.setEnvCutoff(ub + round_up_tolerance))
        SAFE_SOLVER(solver.setVTypeArray(0, num_col, xtype.data()))

        auto eps = TimeSetter::measure([&]() {
            SAFE_SOLVER(solver.mipOptimize())
        });
        printTimeMessage("mip", eps);

        SAFE_SOLVER(solver.setEnvCutoff(SOLVER_INFINITY))

        if_terminate = true;
    }

    void CVRPSolver::terminateByMIP(BbNode *node) {
        BbNode::regenerateEnumMat(node, nullptr, true, optimal_dual_vector);
        auto col_size = node->getIndexColPool().size();
        auto &deleted_columns_in_enumeration_pool = node->refDeletedColumnsInEnumerationColumnPool();
        std::vector<int> added_cols(col_size);
        int cnt = 0;
        for (int i = 0; i < col_size; ++i) {
            if (deleted_columns_in_enumeration_pool[i])continue;
            added_cols[cnt++] = i;
        }
        added_cols.resize(cnt);
        add_column_controller.addColumnsByInspection(added_cols);

        int num_col;
        int num_row;
        size_t nz;
        SAFE_SOLVER(node->refSolver().getNumCol(&num_col))
        SAFE_SOLVER(node->refSolver().getNumRow(&num_row))
        SAFE_SOLVER(node->refSolver().XgetConstraints(&nz, nullptr, nullptr, nullptr, 0, num_row))

        std::cout << "mip: row= " << num_row << " col= " << num_col << " density= "
                << static_cast<double>(nz) / num_row / num_col * 100 << "%" << std::endl;

        node->solveMIP(ub);


        int status;
        SAFE_SOLVER(node->refSolver().getStatus(&status))
        if (status == SOLVER_OPTIMAL) {
            double lp_val;
            SAFE_SOLVER(node->refSolver().getObjVal(&lp_val))
            std::vector<double> X(num_col);
            SAFE_SOLVER(node->refSolver().getX(0, num_col, X.data()))
            bool if_integer, if_feasible;
            updateIntegerSolution(lp_val, X, node->getCols(), if_integer, if_feasible);
            if (!if_integer || !if_feasible) {
                THROW_RUNTIME_ERROR("mip solution is not integer or feasible!");
            }
        } else {
            std::cout << "mip status= " << status << std::endl;
        }
    }
}
