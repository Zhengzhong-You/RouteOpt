//
// Created by You, Zhengzhong on 6/20/23.
//
#include "CVRP.hpp"

using namespace std;
using namespace chrono;

/**
 * SB_MemSave:
 * 1. solve the nodes first, solve the node with larger heuristic value, and then put them into the stack
 * TODO:
 * all the other things will not be implemented here!
 * this mode is only for the use of RCF!
 */

#ifdef BranchFashion_MemSaving

void CVRP::solveModel() {

  cout << "\n<solveModel>\n\n";

  lateProcessing();
  GloBeg = std::chrono::high_resolution_clock::now();

  buildModel();

#ifdef writeEnumerationTrees
  enumeration_col_idx.reserve(size_t(2e7));
#endif

  BBNODE *node = BBT.top();
  BBT.pop();
  tackleUnsolvedNode(node);
  if (node) {
    BBT.push(node);
    supp_set.emplace(node->Val);
    supp_queue.push(node->Val);
  }

  SB_MemSave(BBT);

  LB_transformed = UB;
  cout << "All remaining nodes pruned! Optimality has been proven!  lb=ub= " << UB << endl;
  cout << BIG_PHASE_SEPARATION;

  GlobalGap = (UB - LB_transformed) / UB;
  Solver.SOLVERfreeenv();

  printOptIntSol();

#ifdef writeEnumerationTrees
  writeEnuCols();
#endif
}

void CVRP::terminateNodeMemSave(BBNODE *&root_node) {
  If_in_Enu_State = true;
  int root_index = root_node->Idx;

  cout << "To terminate the current node, " <<
       "chg node selection strategy by exploring the children nodes of the node first!\n";

  --NumExploredNodes;
  for (auto &r1c : root_node->R1Cs) {
    r1c.Mem.clear();
  }
  for (auto &r1c : root_node->R1Cs_multi) {
    r1c.Mem.clear();
  }

  tackleUnsolvedNode(root_node);
  if (root_node) {
    subBBT.push(root_node);
    supp_set.emplace(root_node->Val);
    supp_queue.push(root_node->Val);
  }

  SB_MemSave(subBBT);

  root_node = nullptr;;
  If_in_Enu_State = false;
}

bool CVRP::checkTimeFail() {
  GloEnd = high_resolution_clock::now();
  GloEps = duration<double>(GloEnd - GloBeg).count();
  if (GloEps > GlobalTimeLimit) {
    cout << SMALL_PHASE_SEPARATION;
    cout << "Ins= " << FileName << "  cap= " << Cap
         << "  shut down due to reaching lm. gt.= " << GlobalTimeLimit << endl;
    cout << "lb= " << LB << "  ub= " << UB << "  gap(ub-lb/ub)= "
         << (UB - LB) / (UB) * 100 << "%  gt= "
         << GloEps << " nd= "
         << NumExploredNodes << "  br= " << NumBr << endl;
    return true;
  }
  return false;
}

void CVRP::printInfo(BBNODE *node) {
  cout << BIG_PHASE_SEPARATION;
  cout << "nd_ind= " << node->Idx << "  nd_col= " << NumCol << "  nd_val= " << node->Val << "  nd_dep= "
       << node->TreeLevel << "  et= " << GloEps << "  lb= " << LB << "  ub= "
       << UB << "  nd_rmn= " << BBT.size() << "  sub_rmn= " << subBBT.size() << endl;
  int count = 0;  // 设置一个计数器
  for (auto &brc : node->BrCs) {
    cout << (brc.BrDir ? "true" : "false") << "(" << brc.Edge.first << "," << brc.Edge.second << ")";
    count++;
    if (count % 8 == 0) cout << endl; else cout << " ";
  }
  cout << endl;
}

void CVRP::resetEnv(BBNODE *node) {
  safe_solver(node->solver.SOLVERupdatemodel())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  if (!If_in_Enu_State) {
    convertVertex2R1CsInOneLP(node);
  }
}

void CVRP::SB_MemSave(StackTree &tree) {
  /**
   * take one -> SB -> cutting -> push into the stack
   */
  while (!tree.empty()) {

    if (checkTimeFail()) goto DELETE;

    BBNODE *node = tree.top();
    tree.pop();

    supp_set.erase(node->Val);
    if (abs(LB - supp_queue.top()) > TOLERANCE)updateLB(supp_queue.top());
    if (supp_set.find(supp_queue.top()) == supp_set.end()) supp_queue.pop();

    resetEnv(node);

    printInfo(node);

    recordOptCol(node, true);//the old data can be the other node, so force rewrite!
    if (If_in_Enu_State) {
      getInfoEdge(node, true);
      writeMap_Edge_ColIdx_in_Enu(node);
    }
    pair<int, int> info;
    do_SB(node, info);
    BBNODE *node2{};
    addBrCut(node, node2, info);
    ++NumBr;

    tackleUnsolvedNode(node);
    tackleUnsolvedNode(node2);

    if (node && node2) {
      if (node->Val < node2->Val) swap(node, node2);
    }

    if (node) {
      supp_set.emplace(node->Val);
      supp_queue.push(node->Val);
      tree.push(node);
    }
    if (node2) {
      supp_set.emplace(node2->Val);
      supp_queue.push(node2->Val);
      tree.push(node2);
    }
  }
  DELETE:
  while (!tree.empty()) {
    BBNODE *extra_node = tree.top();
    tree.pop();
    delete extra_node;
  }
}

void CVRP::addBrCut(
    BBNODE *node,
    BBNODE *&node2,
    const pair<int, int> &info
) {
  int numnz;
  BrC bf;
  bf.Edge = info;
  bf.BrDir = true;
  ++node->TreeLevel;
  BidirLinkList *lnode{}, *rnode{};
  if (!If_in_Enu_State) {
    node->Ptr->giveBirth(lnode, rnode);
    getNewCstrCoeffByEdge(node, info, solver_ind, solver_val, numnz);
    bf.IdxBrC = NumRow;
    node2 = new BBNODE(node, rnode, NumCol, IdxNode + 2, bf, NumBucketsPerVertex);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_GREATER_EQUAL, 1, nullptr, node2->solver))
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
  } else {
    bf.IdxBrC = -1;
    CstrIndex.resize(NumRow);
    iota(CstrIndex.begin(), CstrIndex.end(), 0);
    node2 = new BBNODE(node, NumCol, IdxNode + 2, bf);
    reviseEnuColInfoByBrC(node, node2, bf);
  }

  bf.BrDir = false;
  node->BrCs.emplace_back(bf);

  if (!If_in_Enu_State) {
    node->NumParentCols = NumCol;
    deleteArcByBrC_false(node, info);
    node->Ptr = lnode;
    ++NumRow;
    safe_Hyperparameter(checkCST_LIMIT())
  } else {
    CstrIndex.resize(NumRow);
    iota(CstrIndex.begin(), CstrIndex.end(), 0);
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))//cannot be deleted, recover the env is very important!
    reviseEnuColInfoByBrC(node, node, bf);
  }
  node->Idx = ++IdxNode;
  ++IdxNode;
}

void CVRP::tackleUnsolvedNode(BBNODE *&node) {
  ++NumExploredNodes;
  int numParentCols;

  resetEnv(node);

  if (!If_in_Enu_State) {
    cout << "NumBucketsPerVertex= " << NumBucketsPerVertex << " StepSize= " << StepSize << endl;
  }

  double old_val;
  if (node->Idx) old_val = node->Val;

  bool if_sep_cuts = true;
  if (!If_in_Enu_State) {
    PoolBeg4Pricing = 0;
    LastMaxTimeLabeling = 0;
    ForceNotRollback = true;
    cout << "node->TreeLevel= " << node->TreeLevel << endl;
    solveLPInLabeling(node, true, true, true);
    ForceNotRollback = false;
    if (Rollback == 3) {
      if_sep_cuts = false;
      cout << "labels reach the soft limit!" << endl;
    }
  } else {
    MaxNumEnuColPool = node->SizeEnuColPool;
    solveLPByInspection(node, false, false, true);
  }

  if (node->Idx) {
    auto &brc = node->BrCs.back();
    if (brc.BrDir) {
      RealImprovement_up[brc.Edge].first += node->Val - old_val;
      ++RealImprovement_up[brc.Edge].second;
    } else {
      RealImprovement_down[brc.Edge].first += node->Val - old_val;
      ++RealImprovement_down[brc.Edge].second;
    }
  }

  if (node->if_terminated) {
    delete node;
    node = nullptr;
    goto QUIT;
  } else if (If_in_Enu_State && NumCol + node->SizeEnuColPool <= MaxNumRoute4MIP) {
    terminateByMIP(node);
    if (node->if_terminated) {
      delete node;
      node = nullptr;
      goto QUIT;
    }
  }

  if (!If_in_Enu_State && (!node->Idx || !if_sep_cuts)) {
    final_decision_4_arc_elimination = true;
    eliminateArcs(node);
    final_decision_4_enumeration = true;
    enumerateMIP(node);
    if (!node) goto QUIT;
    cleanIdxCol4Node(node, node->NumParentCols, true);
  }

//  if_sep_cuts = false;
  if (if_sep_cuts) {
    sepHybridCuts(node);
    if (!node) goto QUIT;
    else if (node->if_terminated) {
      delete node;
      node = nullptr;
      goto QUIT;
    } else if (If_in_Enu_State && NumCol + node->SizeEnuColPool <= MaxNumRoute4MIP) {
      terminateByMIP(node);
      if (node->if_terminated) {
        delete node;
        node = nullptr;
        goto QUIT;
      }
    }
  } else {
    findNonactiveCuts(node);
    printCutsInfo(node);
  }

  if (!If_in_Enu_State) {
    /**
      * once mem is written, cannot delete any column before the branching selection and cg!
      */
    numParentCols = node->NumParentCols;

    writeCol2Mem(node);

    if (NumCol > LPCol_FINAL_LIMIT) {
      cleanIdxCol4Node(node, Dim);
      cout << BIG_PHASE_SEPARATION;
      cout << "Run column reduction in memory... ncol= " << NumCol << " are left!" << endl;
      cout << BIG_PHASE_SEPARATION;
      //include constructMap itself
      node->Ptr->becomeParent(node, this);
    } else {
      constructMap(node, numParentCols);
    }
  } else {
#ifdef DELUXING_APPLIED
    cout << "Trying to apply RCF..." << endl;
    applyRCF(node, DELUXING_ROUND, true);
    solveLPByInspection(node, false, false, true);
    if (NumCol + node->SizeEnuColPool <= MaxNumRoute4MIP) {
      terminateByMIP(node);
      if (node->if_terminated) {
        delete node;
        node = nullptr;
        goto QUIT;
      }
    }
#endif
#ifdef writeEnumerationTrees
    writeEnuTree(node);
    delete node;
    node = nullptr;
    goto QUIT;
#endif
    regenerateEnuMat(node, nullptr, true);
  }
  QUIT:
  return;
}

#endif
