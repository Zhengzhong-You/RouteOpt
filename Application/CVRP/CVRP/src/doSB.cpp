//
// Created by Zhengzhong You on 1/31/22.
//
#include "CVRP.hpp"

using namespace std;
using namespace chrono;

bool CVRP::solveSBModel() {

  cout << "\n<solveSBModel>\n\n";

  //initialization
  lateProcessing();
  GloBeg = std::chrono::high_resolution_clock::now();

  buildModel();

#ifdef VERBOSE
  GloEnd = chrono::high_resolution_clock::now();
  GloEps = duration<double>(GloEnd - GloBeg).count();
  cout << "Build model...  create art_var= " << Dim << "  main_cstr= " << RealDim << "  lb= " << LB
       << "  ub= " << UB << "  et= " << GloEps << endl;
  cout << BIG_PHASE_SEPARATION;
#endif

#ifdef if_draw_BBT_graph
  string message;
#endif

#ifdef writeEnumerationTrees
  enumeration_col_idx.reserve(size_t(2e7));
#endif

  while (!BBT.empty()) {
    GloEnd = std::chrono::high_resolution_clock::now();
    GloEps = duration<double>(GloEnd - GloBeg).count();

    BBNODE *node = BBT.top();
    BBT.pop();
    LB = node->Val;
    LB_transformed = ceil_transformed_number_related(LB - TOLERANCE);
    ++NumExploredNodes;

    if (LB_transformed >= UB) {
#ifdef VERBOSE
      LB_transformed = UB;
      cout << "Pruned! Optimality has been proven!  lb=ub= " << UB << endl;
      cout << SMALL_PHASE_SEPARATION;
#endif
#ifdef if_draw_BBT_graph
      message = "Pruned";
      constructBBNodeName(node, message);
#endif
      delete node;
      goto DELETE;
    }

    if (GloEps > GlobalTimeLimit) {
      cout << SMALL_PHASE_SEPARATION;
      cout << "Ins= " << FileName << "  cap= " << Cap
           << "  shut down due to reaching lm. gt.= " << GlobalTimeLimit << endl;
      cout << "lb= " << LB << "  ub= " << UB << "  gap(ub-lb/ub)= "
           << (UB - LB) / (UB) * 100 << "%  gt= "
           << GloEps << " nd= "
           << NumExploredNodes << "  br= " << NumBr << endl;
#ifdef if_draw_BBT_graph
      message = "TimeLimit";
      constructBBNodeName(node, message);
#endif
      delete node;
      goto DELETE;
    }

    //更新关于这个node的宏观lp信息
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))

    cout << "NumBucketsPerVertex= " << NumBucketsPerVertex << " StepSize= " << StepSize << endl;

    convertVertex2R1CsInOneLP(node);
#ifdef  VERBOSE
    cout << BIG_PHASE_SEPARATION;
    cout << "nd_ind= " << node->Idx << "  nd_col= " << NumCol << "  nd_val= " << node->Val << "  nd_dep= "
         << node->TreeLevel << "  et= " << GloEps << "  lb= " << LB << "  ub= "
         << UB << "  nd_rmn= " << BBT.size() << endl;
    for (auto &brc : node->BrCs) {
      cout << (brc.BrDir ? "true" : "false") << "(" << brc.Edge.first << "," << brc.Edge.second << ")" << endl;
    }
#endif

    double old_val;
    if (node->Idx) {
      old_val = node->Val;
    }
#if defined(TrainingDataTreeLevel)
#ifdef TrainingDataTreeLevel
    if (node->TreeLevel > TrainingDataTreeLevel) {
      delete node;
      continue;
    }
#endif
#endif

#ifdef find_missed_solution
    vector<int> data;
    bool if_suc;
    findWhySolutionDisappear(node, this, data, if_suc, true);
    if (!if_suc) {
      cout << "This is not the optimal node!" << endl;
      delete node;
      continue;
    }
#endif

#ifdef useM_dynamically
    auto beg_c = std::chrono::high_resolution_clock::now();
#endif

    //reset PoolBeg4Pricing
    PoolBeg4Pricing = 0;
    LastMaxTimeLabeling = 0;
    ForceNotRollback = true;
    solveLPInLabeling(node, true, true, true);
    ForceNotRollback = false;
    bool if_sep_cuts = true;
//    if_sep_cuts = false;

#if defined(STD_BRANCHING) && !defined(openCutsAndEnumerationAtEachEndNode) || defined(CutAtNonRootNotAllowed)
    if (node->Idx) if_sep_cuts = false;
#endif

#ifdef check_lower_bound
    if (node->Idx) {
      if (node->Val < old_val - TOLERANCE) {
        cerr << "lower bound is not wrong!" << endl;
        cerr << "node->Val= " << node->Val << " old_val= " << old_val << endl;
        safe_Hyperparameter(1);
      }
    }
#endif

#ifdef MachineLearning
    //need to be two place!
    if (!node->Idx) {
      ml.RootValue = node->Val;
    }
#endif

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

    if (Rollback == 3) {
      if_sep_cuts = false;
      cout << "labels reach the soft limit!" << endl;
    }

    if (node->if_terminated) {
#ifdef if_draw_BBT_graph
      constructBBNodeName(node, "Pruned");
#endif
      delete node;
      continue;
    }

    if (!node->Idx) {
      final_decision_4_arc_elimination = true;
      eliminateArcs(node);
      final_decision_4_enumeration = true;
      enumerateMIP(node);
      if (!node) continue;
      cleanIdxCol4Node(node, node->NumParentCols, true);
    }

//    if (node->Idx == 0) if_sep_cuts = false;

//    if_sep_cuts = false;
    if (if_sep_cuts) {
//      double tmp_std = CONFIG::CutsTailOff;
//      if (!node->Idx) CONFIG::CutsTailOff = 0.020;
      sepHybridCuts(node);
//      CONFIG::CutsTailOff = tmp_std;
      if (!node) continue;
      else if (node->if_terminated) {
#ifdef if_draw_BBT_graph
        constructBBNodeName(node, "Pruned");
#endif
        delete node;
        continue;
      }
    } else {
#ifndef Setting_NoBranching_In_Enumeration
      findNonactiveCuts(node);
      convertVertex2R1CsInOneLP(node);
      printCutsInfo(node);
#endif
    }
#ifdef MachineLearning
    if (!node->Idx) {
      ml.RootValue = node->Val;
    }
#endif

#ifdef openCutsAndEnumerationAtEachEndNode
    /**
     * can do since in sepCuts the cg is called
     */
    final_decision_4_arc_elimination = true;
    eliminateArcs(node);
    final_decision_4_enumeration = true;
    enumerateMIP(node);
    if (!node) continue;
#endif

#ifdef Use_heuristic_UB
    findBetterUsingHeurEnumeration(node);
    if (!node) continue;
#endif

//    cout << "exit!" << endl;
//    exit(0);

/**
 * once mem is written, cannot delete any column before the branching selection and cg!
 */
    cout << BIG_PHASE_SEPARATION;

    int numParentCols = node->NumParentCols;

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
#ifdef BranchInDefNotAllowed
    cout << "BranchInDefNotAllowed" << endl;
    exit(1);
#endif

#ifdef useM_dynamically
    if (node->Idx) {
      auto end_c = std::chrono::high_resolution_clock::now();
      auto eps = duration<double>(end_c - beg_c).count();
      node->updateState(eps, node->c, node->BrCs.size() - 1);
      double new_r;
      node->calculateR_star(node->Val - old_val, new_r, this);
      node->updateState(new_r, node->geo_r_star, node->BrCs.size() - 1);
      cout << "node->c= " << node->c << " node->geo_r_star= " << node->geo_r_star << endl;
    }
#endif

    cout << "Begin branching...\n";
    //update the solution accordingly
    recordOptCol(node);

    pair<int, int> info;
#ifdef TestIfPhase0IsWorthy
    do_SB_openTest(node, info);
    goto ADD;
#endif

#ifdef indicate_search_tree
    vector<pair<int, int>> infos = {
        {31, 169},
        {3, 179},
        {75, 193},
        {154, 189},
        {0, 124},
        {11, 168},
        {11, 179},
        {4, 112}
    };
    info = infos[node->TreeLevel];
#else
    do_SB(node, info);
#endif
    ADD:
    addBrCut2Unsolved(node, info);

    ++NumBr;
  }
#ifdef VERBOSE
  LB_transformed = UB;
  cout << "All remaining nodes pruned! Optimality has been proven!  lb=ub= " << UB << endl;
  cout << BIG_PHASE_SEPARATION;
#endif

  DELETE:
  while (!BBT.empty()) {
    BBNODE *extra_node = BBT.top();
    BBT.pop();
    ++NumExploredNodes;
#ifdef if_draw_BBT_graph
    constructBBNodeName(extra_node, message);
#endif
    delete extra_node;
  }
  Solver.SOLVERfreeenv();

  GlobalGap = (UB - LB_transformed) / UB;

#ifdef PrintOptimalSolution
  printOptIntSol();
#endif

#ifdef TestIfPhase0IsWorthy
  Experiment::printPhase0TopNAcc();
#endif

#ifdef writeEnumerationTrees
  writeEnuCols();
#endif
  return false;
}

bool CVRP::addBrCut2Unsolved(
    BBNODE *const node,
    const pair<int, int> &info
) {
  if (node->if_terminated) return false;
  int numnz;
  int c = NumRow;//old_num_cut
  BrC bf;
  bf.Edge = info;
  bf.IdxBrC = c;
  BidirLinkList *lnode{}, *rnode{};
  node->Ptr->giveBirth(lnode, rnode);
  ++node->TreeLevel;
  getNewCstrCoeffByEdge(node, info, solver_ind, solver_val, numnz);

  //2
  bf.BrDir = true;
  auto node2 = new BBNODE(node, rnode, NumCol, IdxNode + 2, bf, NumBucketsPerVertex);
  safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_GREATER_EQUAL, 1, nullptr, node2->solver))
  safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))

  bf.BrDir = false;
  node->BrCs.emplace_back(bf);
#ifdef if_draw_BBT_graph
  //must be after node->BrCs.emplace_back(bf);
  constructBBNodeName(node);
#endif
  //if not use the edge, we delete the left edges
  deleteArcByBrC_false(node, info);
  node->Ptr = lnode;
  node->Idx = ++IdxNode;
  node->NumParentCols = NumCol;
  ++IdxNode;
  ++NumRow;
  safe_Hyperparameter(checkCST_LIMIT())

  //push to BBT
  BBT.push(node2);//solve true first
  BBT.push(node);

  return false;
}

//1 and normal
