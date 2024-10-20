//
// Created by You, Zhengzhong on 12/9/23.
//

#include <branching.hpp>

#include "cvrp.hpp"

#ifdef DELUXING_APPLIED
using namespace std;
using namespace std::chrono;
using namespace Eigen;

void CVRP::applyRCF(BbNode *node, int round, bool if_verbose) {
    cout << "WARNING: RCF can only be applied by GRB model so far!" << endl;
    cstr_index.resize(num_row);
    iota(cstr_index.begin(), cstr_index.end(), 0);
    regenerateEnumMat(node, nullptr, true);
    Solver local_solver{};
    local_solver.getSolver(&node->getSolver());
    size_t numnz = 0;
    int ccnt = 0;
    auto &mat = node->matrix_in_enumeration.front();
    auto &cost = node->cost_for_columns_in_enumeration_column_pool;
    int pool_size = node->size_enumeration_col_pool;
    int old_num_col = num_col;

    SparseMatrix<double, ColMajor> col_mat = mat;

    vector<size_t> solver_beg(pool_size + 1);
    vector<int> solver_ind(col_mat.nonZeros());
    vector<double> solver_val(col_mat.nonZeros());
    vector<double> solver_obj(pool_size);
    for (int i = 0; i < node->size_enumeration_col_pool; ++i) {
        solver_beg[ccnt] = numnz;
        solver_obj[ccnt++] = cost[i];
        for (SparseMatrix<double, ColMajor>::InnerIterator it(col_mat, i); it; ++it) {
            solver_ind[numnz] = (int) it.row();
            solver_val[numnz++] = it.value();
        }
    }
    solver_beg[ccnt] = numnz;
    safe_solver(local_solver.XaddVars(ccnt,
        numnz,
        solver_beg.data(),
        solver_ind.data(),
        solver_val.data(),
        solver_obj.data(),
        nullptr,
        nullptr,
        nullptr,
        nullptr))
    safe_solver(local_solver.updateModel())
    vector<int> idx;
    idx.reserve(num_col + ccnt);

    int beta1 = 1000, beta2 = 1000;
    double time_limit = 100000;
    deLuxing(local_solver.model,
             BaseBranching::ub + round_up_tolerance,
             // real_dim,
             round,
             beta1,
             beta2,
             idx,
             time_limit,
             MIP_GAP_TOLERANCE,
             if_verbose);
    int new_col;
    safe_solver(local_solver.getNumCol(&new_col))
    if (new_col == 0) {
        idx.resize(num_col + ccnt);
        iota(idx.begin(), idx.end(), 0);
    }
    safe_solver(local_solver.freeModel())

    auto if_del = new bool[num_col]();
    int cnt = 0;
    for (auto i: idx) {
        if (i >= num_col) break;
        if_del[i] = true;
        ++cnt;
    }

    vector<int> delete_cols;
    for (int i = 1; i < num_col; ++i) {
        if (if_del[i]) delete_cols.emplace_back(i);
    }
    rmLPCols(node, delete_cols);
    delete[]if_del;
    for (int i = cnt; i < idx.size(); ++i) {
        node->deleted_columns_in_enumeration_pool[idx[i] - old_num_col] = true;
    }
}
#endif
