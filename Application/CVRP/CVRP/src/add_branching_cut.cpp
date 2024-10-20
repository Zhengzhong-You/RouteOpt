#include "cvrp.hpp"
#if SOLUTION_TYPE == 1
#include "best_bound_first_branching.hpp"
#elif SOLUTION_TYPE == 2
#include "depth_first_branching.hpp"
#endif

using namespace std;
using namespace chrono;

void CVRP::addBranchCutToUnsolved(
    BbNode *const node,
    const pair<int, int> &info
) {
    if (node->getIfTerminated()) return;

    int c = num_row;
    Brc bf;
    bf.edge = info;
    bf.idx_br_c = c;
    ++node->getTreeLevel();
    vector<int> solver_ind;
    vector<double> solver_val;
    getNewConstraintCoefficientByEdge(node, info, solver_ind, solver_val);

    bf.br_dir = true;
    auto node2 = new BbNode(node, idx_node + 2, bf, num_buckets_per_vertex);
    safe_solver(addBranchConstraint(solver_ind, solver_val, SOLVER_EQUAL, 1, nullptr, node2->solver))

    bf.br_dir = false;
    bf.idx_br_c = -1;
    node->getBrCs().emplace_back(bf);

    auto it = std::find(solver_ind.begin(), solver_ind.end(), 0);
    if (it != solver_ind.end()) solver_ind.erase(it);

    rmLPCols(node, solver_ind);
    deleteArcByFalseBranchConstraint(node->all_forward_buckets, info);

    symmetry_prohibit_call(deleteArcByFalseBranchConstraint(node->all_backward_buckets, info))

    node->index = ++idx_node;
    ++idx_node;

#if SOLUTION_TYPE == 1
    BestBoundFirstBranching::bbt.push(node2); //solve true first
    BestBoundFirstBranching::bbt.push(node);
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::addNodeIn(DepthFirstBranching::bbt, node2);
  DepthFirstBranching::addNodeIn(DepthFirstBranching::bbt, node);
#endif
}


void CVRP::addBranchConstraint2ColPoolInEnumByColMap(BbNode *node,
                                                     const pair<int, int> &edge) const {
    auto &mat = node->matrix_in_enumeration;
    int size = node->size_enumeration_col_pool;
    if (!size) return;
    auto &mat_last = mat.back();
    mat_last.setZero();
    vector<Eigen::Triplet<double> > triplets(size);
    sparseRowMatrixXd tmp(1, size);
    int cnt = 0;
    int ai = edge.first, aj = edge.second;
    auto &mat0 = mat.front();
    if (ai) {
        tmp = mat0.row(ai - 1) + mat0.row(aj - 1);
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
            if (it.value() > 1.5) {
                for (auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;; ++j) {
                    int current_node = col_pool4_pricing[j];
                    if (!current_node) break;
                    if (current_node == ai) {
                        if (col_pool4_pricing[j + 1] == aj || col_pool4_pricing[j - 1] == aj)
                            triplets[cnt++] = {0, (int) it.col(), 1};
                    }
                }
            }
        }
    } else {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat0, aj - 1); it; ++it) {
            if (it.value() > 0.5) {
                auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;
                if (col_pool4_pricing[j] == aj) {
                    triplets[cnt++] = {0, (int) it.col(), 1};
                    continue;
                }
                for (;; ++j) {
                    int current_node = col_pool4_pricing[j];
                    if (!current_node) break;
                }
                if (col_pool4_pricing[j - 1] == aj) {
                    triplets[cnt++] = {0, (int) it.col(), 1};
                }
            }
        }
    }
    safe_eigen(mat_last.setFromTriplets(triplets.begin(), triplets.end());)
}
