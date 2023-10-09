//
// Created by You, Zhengzhong on 6/17/23.
//

#include "CVRP.hpp"
using namespace std;

#ifdef NominalBranchingInEnu

void CVRP::NominalBranching(BBNODE *node, bool &if_suc) {
//  //get constraints name
//  vector<char *> ConstrName(NumRow);
//  //find NominalBranching
//  for (auto &i : node->NBrCs) {
//    if (i.IdxNBrC != -2) {
//      GRBgetstrattrarray(node->solver.model, GRB_STR_ATTR_CONSTRNAME, i.IdxNBrC, 1, ConstrName.data());
//      if (strstr(ConstrName[0], "NominalBranching") == nullptr) {
//        cout << "error in NominalBranching" << endl;
//        exit(1);
//      } else {
//        cout << "i= " << i.IdxNBrC << " name= " << ConstrName[0] << endl;
//      }
//    }
//  }

  if_suc = false;
  if (node->if_terminated) return;

  BranchingColSet.clear();
  auto &mat = node->MatInEnu;
  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
  Eigen::Map<RowVectorXd> pi(Pi, NumRow);
  RowVectorXd rc = node->Cost4ColsInEnuColPool - pi * mat;
  OptGap = calculateOptGap(node);
  double half_gap = OptGap / 2;
  for (int i = 0; i < node->SizeEnuColPool; ++i) {
    if (node->Deleted_ColsInEnuPool[i]) continue;
    if (rc[i] > half_gap) {
      BranchingColSet.emplace_back(i);
    }
  }
  cout << "NominalBranching" << endl;
  if (BranchingColSet.empty()) {
    cout << "BranchingColSet is empty" << endl;
    return;
  }
  if_suc = true;

  NBrC nbf{-2};//do not change without thinking
  ++node->TreeLevel;
  auto node2 = new BBNODE(node, NumCol, BranchingColSet, IdxNode + 2, nbf);

  int cnt = 0;
  int keep = 1;
  safe_solver(node->solver.SOLVERgetRC(0, NumCol, RC))
  solver_ind[cnt++] = 0;//keep the first column
  for (int i = keep; i < NumCol; ++i) {
    if (RC[i] > half_gap) {
      solver_ind[cnt++] = i;
    } else {
      node2->IdxCols[keep++] = node2->IdxCols[i];
    }
  }

  if (cnt >= 2) {
    safe_solver(node2->solver.SOLVERdelvars(cnt - 1, solver_ind + 1))
  }
  fill_n(solver_val, cnt, 1);
  safe_solver(node->solver.SOLVERaddconstr(cnt, solver_ind, solver_val, SOLVER_EQUAL, 1, nullptr))

  cleanColsInPool(node, node2);
  subBBT2.push(node2);

  nbf.IdxNBrC = NumRow;
  ++NumRow;
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  node->NBrCs.emplace_back(nbf);
  node->Idx = ++IdxNode;
  ++IdxNode;

  mat.conservativeResize(NumRow, Eigen::NoChange);
  mat.row(NumRow - 1).setZero();
  for (auto &i : BranchingColSet) {
    mat(NumRow - 1, i) = 1;
  }
}
#endif
