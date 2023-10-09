//
// Created by Zhengzhong You on 3/27/23.
//


#include <utility>
#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::CGTesting(BBNODE *node, bool if_exact_CG, bool if_record_product_val,
                     bool if_record_improvement, bool if_force_complete) {
  if (Branch_pair.size() == 1) {
    cout << "CGTesting: Branch_pair.size() == 1, return!" << endl;
#ifdef useM_dynamically
    node->obj_change[Branch_pair[0]] = {0, 0, k};
#endif
    return;
  }
  auto beg = high_resolution_clock::now();
  int BeforeNumRow = NumRow;
  auto org_val = node->Val;
  int cnt = 0;
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);
  int col_start = NumCol;
  auto pseudo_down = &HeuristicImprovement_down;
  auto pseudo_up = &HeuristicImprovement_up;
  if (if_exact_CG) {
    if (!If_in_Enu_State) {
      ForceNotRollback = true;
      if_force_not_regenerate_bucket_graph = true;
    }
    pseudo_down = &RealImprovement_down;
    pseudo_up = &RealImprovement_up;
  }
  if (If_in_Enu_State) {
    sparseRowMatrixXd mat(1, node->SizeEnuColPool);
    node->MatInEnu.push_back(std::move(mat));
  }

  fill(solver_ind2, solver_ind2 + NumCol, BeforeNumRow);
  vector<double> solver_val3(NumCol);
  vector<int> solver_ind1(NumCol);
  iota(solver_ind1.begin(), solver_ind1.end(), 0);
  bool if_changed = false;

  ++NumRow;
  safe_Hyperparameter(checkCST_LIMIT())

#ifdef useM_dynamically
  unordered_map<pair<int, int>, pair<double, double>, PairHasher> obj_change;
#endif
  double max_product = 0;
  for (auto &edge : Branch_pair) {
    int ai = edge.first;
    int aj = edge.second;
    if (!if_force_complete) if (max_product > branch_LP[edge]) continue;
    cout << MID_PHASE_SEPARATION;
    cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
    int numnz;
    double temp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    if (!if_changed) {
      safe_solver(addBranchconstr(numnz,
                                  solver_ind,
                                  solver_val,
                                  SOLVER_LESS_EQUAL,
                                  0,
                                  nullptr,
                                  node->solver))
      if_changed = true;
    } else {
      chgBranchconstr(solver_val3.data(),
                      solver_ind2,
                      solver_ind1.data(),
                      numnz,
                      solver_ind,
                      solver_val,
                      SOLVER_LESS_EQUAL,
                      0,
                      node->solver);
    }
    safe_solver(node->solver.SOLVERupdatemodel())
    node->BrCs.back().Edge = {ai, aj};
    node->BrCs.back().BrDir = false;
    if (If_in_Enu_State) {
      addBrC2ColPoolInEnuBycolmap(node, edge);
      solveLPByInspection(node, true, !if_exact_CG, false);
    } else {
      PoolBeg4Pricing = 0;
      solveLPInLabeling(node, true, if_exact_CG, false);
    }
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif1 = calculateDif(temp_val, org_val);
    cout << SMALL_PHASE_SEPARATION;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    node->BrCs.back().BrDir = true;
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    if (If_in_Enu_State) {
      solveLPByInspection(node, true, !if_exact_CG, false);
    } else {
      PoolBeg4Pricing = 0;
      solveLPInLabeling(node, true, if_exact_CG, false);
    }
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val);
    auto product = dif1 * dif2;
    cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
         << left << product << endl;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))

    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
#ifdef useM_dynamically
    obj_change[edge] = {dif1, dif2};
#endif
    Branch_Val[cnt++] = {edge, product};
    if (product > max_product) {
      max_product = product;
    }
    if (if_record_improvement) {
      (*pseudo_down)[edge].first += dif1;
      ++(*pseudo_down)[edge].second;
      (*pseudo_up)[edge].first += dif2;
      ++(*pseudo_up)[edge].second;
    }
  }

  Branch_Val.resize(cnt);
  safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
  safe_solver(node->solver.SOLVERupdatemodel())
  --NumRow;

  if (if_exact_CG && !If_in_Enu_State) {
    ForceNotRollback = false;
    if_force_not_regenerate_bucket_graph = false;
  }
  if (If_in_Enu_State) {
    node->MatInEnu.pop_back();
  }
  safe_solver(node->solver.SOLVERreoptimize())

  if (if_record_product_val) {
    sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
      return a.second > b.second;
    });
    Branch_Pair_Val = Branch_Val;
    Branch_pair = {Branch_Val[0].first};
  } else {
    auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
      return a.second < b.second;
    });
    Branch_pair = {it->first};
  }
  node->BrCs.pop_back();
#ifdef useM_dynamically
  node->obj_change[Branch_pair[0]] = {obj_change[Branch_pair[0]].first, obj_change[Branch_pair[0]].second, k};
#endif
  /**
   * if_terminated, Val. these two cannot be deleted! since when exact cg could determine if this node is no promising to be solved
   */
  node->if_terminated = false;
  node->Val = org_val;
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "CGTesting spent= " << eps << "s" << endl;
}

#ifdef MASTER_VALVE_ML

void CVRP::useModelInPhase2(BBNODE *node, int num) {
  if (Branch_pair.size() == 1) {
    cout << "useModelInPhase2: only one edge, no need to use model" << endl;
    return;
  }
  auto beg = high_resolution_clock::now();
  getTrainingDataInPhase2(node);
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  transform(Branch_pair.begin(),
            Branch_pair.end(),
            Branch_Val.begin(),
            [&](const auto &edge) {
              return std::make_pair(edge, 0.0);
            });
  ml.predict(Branch_Val, 2);

#ifdef checkSimilar
  int cnt=0;
    for (int i = 1; i < checkSimilar; ++i) {
      if (Branch_Val[i].second/Branch_Val[0].second > 0.8) {
        ++cnt;
      }
    }
    num = cnt+1;
#endif

  Branch_pair.resize(num);
  transform(Branch_Val.begin(),
            Branch_Val.begin() + num,
            Branch_pair.begin(),
            [&](const auto &val) {
              return val.first;
            });
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();

#ifdef useM_dynamically
  tested_Branch_pair.resize(Branch_Val.size());
  transform(Branch_Val.begin(),
            Branch_Val.end(),
            tested_Branch_pair.begin(),
            [&](const auto &val) {
              return val.first;
            });
  cout << "eps= " << eps << endl;
  double average_t = eps / 2 / Branch_Val.size();
  cout << "average time= " << average_t << endl;
  node->updateState(average_t, node->t4oneLP, node->BrCs.size());
  cout << "node->t4oneLP= " << node->t4oneLP << endl;
#endif
  cout << "useModelInPhase2 spent= " << eps << "s" << endl;
}

void CVRP::getTrainingDataInPhase2(BBNODE *node) {
  if (Branch_pair.size() == 1) {
    cout << "getTrainingDataInPhase2: only one edge, no need to get this data!" << endl;
    return;
  }
  int BeforeNumRow = NumRow;
  auto org_val = node->Val;

  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);

  fill(solver_ind2, solver_ind2 + NumCol, BeforeNumRow);
  vector<double> solver_val3(NumCol);
  vector<int> solver_ind1(NumCol);
  iota(solver_ind1.begin(), solver_ind1.end(), 0);
  bool if_changed = false;

  branch_LP.clear();
  ++NumRow;
  safe_Hyperparameter(checkCST_LIMIT())
  for (auto &edge : Branch_pair) {
    int numnz;
    double tmp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    if (!if_changed) {
      safe_solver(addBranchconstr(numnz,
                                  solver_ind,
                                  solver_val,
                                  SOLVER_LESS_EQUAL,
                                  0,
                                  nullptr,
                                  node->solver))
      if_changed = true;
    } else {
      chgBranchconstr(solver_val3.data(),
                      solver_ind2,
                      solver_ind1.data(),
                      numnz,
                      solver_ind,
                      solver_val,
                      SOLVER_LESS_EQUAL,
                      0,
                      node->solver);
    }
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto dif1 = calculateDif(tmp_val, org_val);
    ml.collectResolvingFeatures(this, node, edge, BeforeNumRow, tmp_val, org_val, numnz, false);
#ifdef explore_more_stage2_feature
    veryLightHeuristicCG4ML(node, org_val, BeforeNumRow, edge, false, col_start, 1);
#endif
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto dif2 = calculateDif(tmp_val, org_val);
    ml.collectResolvingFeatures(this, node, edge, BeforeNumRow, tmp_val, org_val, numnz, true);
#ifdef explore_more_stage2_feature
    veryLightHeuristicCG4ML(node, org_val, BeforeNumRow, edge, true, col_start, 1);
#endif
    cout << "( " << edge.first << " , " << edge.second << " ): " << dif1 * dif2 << endl;
    branch_LP[edge] = dif1 * dif2;
    LPTestingImprovement_down[edge].first += dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += dif2;
    ++LPTestingImprovement_up[edge].second;
  }
  --NumRow;
  safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
  safe_solver(node->solver.SOLVERreoptimize())
  node->BrCs.pop_back();
  cout << "getTrainingDataInPhase2!" << endl;
}

void CVRP::veryLightHeuristicCG4ML(BBNODE *node,
                                   double org_val,
                                   int BeforeNumRow,
                                   const pair<int, int> &edge,
                                   bool dir,
                                   int col_start,
                                   int iter) {
  if (If_in_Enu_State) throw runtime_error("veryLightHeuristicCG4ML: If_in_Enu_State");
  node->BrCs.back().Edge = edge;
  node->BrCs.back().BrDir = dir;
  PoolBeg4Pricing = 0;
  runCertainNumberOfHeurCGs(node, iter);
  int add = NumCol - col_start;
  double min_rc, mean_rc;
  if (add == 0) {
    min_rc = 0;
    mean_rc = 0;
  } else {
    min_rc = get<2>(NegativeRCLabelTuple[0]);
    mean_rc = accumulate(NegativeRCLabelTuple.begin(), NegativeRCLabelTuple.begin() + add, 0.0,
                         [](double sum, const auto &p) { return sum + get<2>(p); })
        / (int) NegativeRCLabelTuple.size();
  }
  double tmp_val;
  safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
  ml.collectResolvingFeatures_runCG(this, node, edge, min_rc, mean_rc, add, BeforeNumRow, tmp_val, org_val);
  node->if_Int = false;
  node->if_terminated = false;
  safe_solver(node->solver.SOLVERdelvars(add, const_for_branching + col_start))
  safe_solver(node->solver.SOLVERupdatemodel())
  NumCol = col_start;
}

#endif

void CVRP::addBrC2ColPoolInEnuBycolmap(BBNODE *node,
                                       const pair<int, int> &edge) {
  auto &mat = node->MatInEnu;
  int size = node->SizeEnuColPool;
  if (!size) return;
  auto &colmap = node->map_col_pool;
  auto &mat_last = mat.back();
  mat_last.setZero();
  for (auto i : colmap[edge]) {
    mat_last.insert(0, i) = 1;
  }
}

void CVRP::runCertainNumberOfHeurCGs(BBNODE *node, int iteration) {
  int iter = 0;
  double b4_node_val;
  safe_solver(node->solver.SOLVERgetObjVal(&b4_node_val))
  int env_method;
  bool if_changed = false;
  optimizeLP4OneIter(node, b4_node_val);
  safe_solver(node->solver.SOLVERgetenvMethod(&env_method))
  if (env_method != SOLVER_PRIMAL_SIMPLEX) {
    safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_PRIMAL_SIMPLEX))
    if_changed = true;
  }
  while (++iter <= iteration) {
    if (generateColsByLighterHeur(node) == 0) break;
    optimizeLP4OneIter(node, b4_node_val);
    b4_node_val = LPVal;
  }
  if (if_changed) {
    safe_solver(node->solver.SOLVERsetenvMethod(env_method))
  }
}

