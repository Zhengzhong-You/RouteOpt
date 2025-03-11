/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "node.hpp"

#include "node_macro.hpp"


namespace RouteOpt::Application::CVRP {
    int BbNode::dim{};
    int BbNode::node_idx_counter{};
    sparseRowMatrixXd BbNode::row_basic_matrix{};

    void BbNode::buildModel(int num_vehicle, int dim, Solver *cvrp_solver, BbNode *node) {
        BbNode::dim = dim;
        node->if_root_node = true;
        std::vector<int> solver_beg, solver_ind;
        std::vector<double> solver_val, solver_obj;

        node->cols.emplace_back();
        solver_beg.emplace_back(0);
        auto &col = node->cols.back();
        for (int i = 1; i < dim; ++i) {
            col.col_seq.emplace_back(i);
            solver_ind.emplace_back(i - 1);
            solver_val.emplace_back(1);
        }
        solver_ind.emplace_back(dim - 1);
        solver_val.emplace_back(num_vehicle);
        solver_beg.emplace_back(static_cast<int>(solver_ind.size()));
        col.forward_concatenate_pos = static_cast<int>(col.col_seq.size()) - 1;

        solver_obj.emplace_back(OBJ_ARTIFICIAL);
        std::vector<double> rhs(dim, 1);
        std::vector<char> sense(dim, SOLVER_EQUAL);

        rhs[dim - 1] = num_vehicle;
        sense[dim - 1] = SOLVER_GREATER_EQUAL;

        node->solver.getEnv(cvrp_solver);

        SAFE_SOLVER(node->solver.newModel(MODEL_NAME, 0, nullptr, nullptr, nullptr, nullptr, nullptr))
        SAFE_SOLVER(node->solver.addConstraints(
                dim,
                0,
                nullptr,
                nullptr,
                nullptr,
                sense.data(),
                rhs.data(),
                nullptr)
        )
        SAFE_SOLVER(node->solver.addVars(1,
            static_cast<int>(solver_ind.size()),
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            solver_obj.data(),
            nullptr,
            nullptr,
            nullptr,
            nullptr))
        SAFE_SOLVER(node->solver.updateModel())
    }


    void BbNode::createBasicMatrix() {
        auto size_enumeration_col_pool = static_cast<int>(index_columns_in_enumeration_column_pool.size());
        if (size_enumeration_col_pool > 0) {
            basic_matrix = (matrix_in_enumeration.front()).block(
                0, 0, dim - 1, size_enumeration_col_pool);
        } else {
            basic_matrix = Eigen::SparseMatrix<double>(dim - 1, 0);
        }
    }


    BbNode::~BbNode() {
        if (all_forward_buckets) {
            for (int i = 0; i < dim; ++i) {
                delete[]all_forward_buckets[i];
            }
            delete[]all_forward_buckets;
        }

        if (all_backward_buckets) {
            for (int i = 0; i < dim; ++i) {
                delete[]all_backward_buckets[i];
            }
            delete[]all_backward_buckets;
        }
        solver.freeModel();
    }
}
