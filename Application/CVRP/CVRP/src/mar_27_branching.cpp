#include <utility>

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::doSB_default(BBNODE *node,
                        pair<int, int> &info) {
  if (node->if_terminated) return;
  //not using exact labeling control
  //step 0: clean the branch_pair
  Branch_pair.clear();
  doSB_Stage0(node);
  doSB_Stage1(node);
  doSB_Stage2(node);
#ifdef OpenHeurCGInUseModel
  doSB_Stage3(node);
#endif
#ifdef Stage1
  generateMap_Edge_SB_rankByExactCG(node);
  tellAccuracy(node->TreeLevel, 0);f
#endif
  info = Branch_pair[0];

  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++BranchChoice[info];
  ++BranchTimes;
}

void CVRP::doSB_Stage0(BBNODE *node) {
  //using old info to select one candidate
  getInfoEdge(node, true);
  vector<tuple<int, int, double>> fracEdges;
#ifdef MachineLearning
  ml.all_fractional_edges.clear();
#endif
#if !defined(MachineLearning) || defined(openZeroPhase)
  vector<tuple<int, int, double>> OldBranch;
  int num_all_frac_edge = 0;
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeVal[i] < 1 - TOLERANCE) {
#ifdef MachineLearning
      ml.all_fractional_edges.emplace_back(node->EdgeTail[i], node->EdgeHead[i], node->EdgeVal[i]);
#endif
      ++num_all_frac_edge;
      int ai = node->EdgeTail[i];
      int aj = node->EdgeHead[i];
      double frac_down = node->EdgeVal[i];
      double frac_up = 1 - frac_down;
      pair<int, int> edge = {ai, aj};
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_up_ptr;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_down_ptr;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator up_iter;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator down_iter;
      int loop = 1;
      //test if in real br
      HERE:
      switch (loop) {
        case 1: {
          improvement_up_ptr = &RealImprovement_up;
          improvement_down_ptr = &RealImprovement_down;
          break;
        }
        case 2: {
#if         !defined(MachineLearning) || defined(OpenHeurCGInUseModel)
          improvement_up_ptr = &HeuristicImprovement_up;
          improvement_down_ptr = &HeuristicImprovement_down;
#endif
          break;
        }
        case 3: {
          improvement_up_ptr = &LPTestingImprovement_up;
          improvement_down_ptr = &LPTestingImprovement_down;
          break;
        }
        default: {
          goto END;
        }
      }
      up_iter = improvement_up_ptr->find(edge);
      down_iter = improvement_down_ptr->find(edge);
      if (up_iter != improvement_up_ptr->end() && down_iter != improvement_down_ptr->end()) {
        double up = up_iter->second.first / up_iter->second.second;
        double down = down_iter->second.first / down_iter->second.second;
        OldBranch.emplace_back(ai, aj, up * frac_up * down * frac_down);
        continue;
      } else {
        ++loop;
        goto HERE;
      }
      END:
      fracEdges.emplace_back(ai, aj, abs(node->EdgeVal[i] - 0.5));
    }
  }
  //sort
  std::sort(OldBranch.begin(), OldBranch.end(), [](const auto &a, const auto &b) {
    return get<2>(a) > get<2>(b);
  });
  //reverse
  std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
    return get<2>(a) < get<2>(b);
  });
#ifdef MachineLearning
  int all_branch_phase1 = min(CONFIG::ML_NumBrCandiInPhase1, num_all_frac_edge);
  int sudo_cap = ceil(CONFIG::ML_Fracsudo * all_branch_phase1);
  sudo_cap = max(sudo_cap, all_branch_phase1 - (int) fracEdges.size());
#else
  int all_branch_phase1 = min(CONFIG::BranPhase1, num_all_frac_edge);
  int sudo_cap = ceil(CONFIG::Frac4sudoCostBranPhase0 * all_branch_phase1);
#endif
  auto num_sudo_edge = min(int(OldBranch.size()), sudo_cap);
  Branch_pair.resize(all_branch_phase1);
  transform(OldBranch.begin(), OldBranch.begin() + num_sudo_edge, Branch_pair.begin(), [](const auto &a) {
    return make_pair(get<0>(a), get<1>(a));
  });
#if defined(Stage1) || defined(GenerateTrainingData_2) || defined(UseModel)
  Branch_pair_from_pseudo.resize(num_sudo_edge);
  transform(OldBranch.begin(), OldBranch.begin() + num_sudo_edge, Branch_pair_from_pseudo.begin(), [](const auto &a) {
    return make_pair(get<0>(a), get<1>(a));
  });
  cout << "pseudo: " << num_sudo_edge << endl;
#endif
  auto candi_cnt = min(all_branch_phase1 - num_sudo_edge, int(fracEdges.size()));
  transform(fracEdges.begin(),
            fracEdges.begin() + candi_cnt,
            Branch_pair.begin() + num_sudo_edge,
            [](const auto &a) {
              return make_pair(get<0>(a), get<1>(a));
            });
  candi_cnt += num_sudo_edge;
  Branch_pair.resize(candi_cnt);
#if defined(Stage1) || defined(GenerateTrainingData_2) || defined(UseModel)
  int frac_candidates = candi_cnt - num_sudo_edge;
  Branch_pair_from_fractional.resize(frac_candidates);
  transform(fracEdges.begin(),
            fracEdges.begin() + frac_candidates,
            Branch_pair_from_fractional.begin(),
            [](const auto &a) {
              return make_pair(get<0>(a), get<1>(a));
            });
  cout << "frac: " << frac_candidates << endl;
#endif
#else
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeVal[i] < 1 - TOLERANCE) {
      fracEdges.emplace_back(node->EdgeTail[i], node->EdgeHead[i], abs(node->EdgeVal[i] - 0.5));
      ml.all_fractional_edges.emplace_back(node->EdgeTail[i], node->EdgeHead[i], node->EdgeVal[i]);
    }
  }
  std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
    return get<2>(a) < get<2>(b);
  });
  Branch_pair.resize(int(fracEdges.size()));
  transform(fracEdges.begin(),
            fracEdges.end(),
            Branch_pair.begin(),
            [](const auto &a) {
              return make_pair(get<0>(a), get<1>(a));
            });
#endif
  cout << "StageZero Branch_pair size: " << Branch_pair.size() << endl;
}

void CVRP::doSB_Stage1(BBNODE *node) {
  auto begin = high_resolution_clock::now();
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  double org_val = node->Val;
  int BeforeNumRow = NumRow;
  int cnt = 0;

#ifdef MachineLearning
  ml.collectEdgeRelatedFeatures(this, node, org_val);
  for (auto &edge : Branch_pair) {
    int numnz;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    ml.collectVariableRelatedFeatures(this, node, edge, BeforeNumRow, numnz, org_val);
  }
  ml.collectInteractFeatures(this);//we do not change here for now, since this feature could be abandoned in the future!
//  transform(Branch_pair.begin(), Branch_pair.end(), Branch_Val.begin(), [](const auto &edge) {
//    return make_pair(edge, 0.0);
//  });
  unordered_set<pair<int, int>, PairHasher> record;
  record.reserve(Branch_pair.size());
#else
  for (auto &edge : Branch_pair) {
    int ai = edge.first;
    int aj = edge.second;
    int numnz;
    double tmp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    ++NumRow;
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto dif1 = calculateDif(tmp_val, org_val);
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto dif2 = calculateDif(tmp_val, org_val);
    auto product = dif1 * dif2;
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    --NumRow;
    Branch_Val[cnt++] = {edge, product};
    LPTestingImprovement_down[edge].first += dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += dif2;
    ++LPTestingImprovement_up[edge].second;
  }
  safe_solver(node->solver.SOLVERreoptimize())
#endif

#if defined( GenerateTrainingData_2 ) || defined( UseModel)
  int all_size = ml.giveTestingNumCandidates(node->TreeLevel);
  auto ratio_pseudo =
      ((double) Branch_pair_from_pseudo.size()
          / (double) (Branch_pair_from_pseudo.size() + Branch_pair_from_fractional.size()));
  int pseudo_size = (int) (all_size * ratio_pseudo);
  int frac_size = all_size - pseudo_size;
  pseudo_size = min(pseudo_size, int(Branch_pair_from_pseudo.size()));
  frac_size = min(frac_size, int(Branch_pair_from_fractional.size()));
//evaluate the model from pseudo cost
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_pseudo(Branch_pair_from_pseudo.size());
  transform(Branch_pair_from_pseudo.begin(),
            Branch_pair_from_pseudo.end(),
            Branch_Val_model_pseudo.begin(),
            [](const auto &edge) {
              return make_pair(edge, 0.0);
            });
  ml.predict(Branch_Val_model_pseudo, 1);
  sort(Branch_Val_model_pseudo.begin(), Branch_Val_model_pseudo.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  Branch_pair.resize(pseudo_size + frac_size);
  transform(Branch_Val_model_pseudo.begin(),
            Branch_Val_model_pseudo.begin() + pseudo_size,
            Branch_pair.begin(),
            [](const auto &a) {
              return a.first;
            });
  //evaluate the model from fractional
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_fractional(Branch_pair_from_fractional.size());
  transform(Branch_pair_from_fractional.begin(),
            Branch_pair_from_fractional.end(),
            Branch_Val_model_fractional.begin(),
            [](const auto &edge) {
              return make_pair(edge, 0.0);
            });
  ml.predict(Branch_Val_model_fractional, 2);
  sort(Branch_Val_model_fractional.begin(), Branch_Val_model_fractional.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  transform(Branch_Val_model_fractional.begin(),
            Branch_Val_model_fractional.begin() + frac_size,
            Branch_pair.begin() + pseudo_size,
            [](const auto &a) {
              return a.first;
            });
#endif

//#if  defined(Stage1)
//  //  int tmp_size = node->Idx ? CONFIG::UseModelBranNormalStage4ModelNormal : CONFIG::UseModelBranRootStageModelNormal;
//  //  int size1 = min(tmp_size, int(Branch_Val.size()));
//  //  ml.predict(Branch_Val, 1);
//  //  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
//  //    return a.second < b.second;
//  //  });
//  //  unordered_set<pair<int, int>, PairHasher> record;
//  //  record.reserve(Branch_Val.size());
//  //  // insert to record by size1 using transform
//  //  transform(Branch_Val.begin(), Branch_Val.begin() + size1, inserter(record, record.begin()), [](const auto &a) {
//  //    return a.first;
//  //  });
//  //  tmp_size = node->Idx ? CONFIG::UseModelBranNormalStage4ModelPseudo : CONFIG::UseModelBranRootStageModelPseudo;
//  //  int size2 = min(tmp_size, int(Branch_Val.size()));
//  //  ml.predict(Branch_Val, 2);
//  //  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
//  //    return a.second < b.second;
//  //  });
//  //  transform(Branch_Val.begin(), Branch_Val.begin() + size2, inserter(record, record.begin()), [](const auto &a) {
//  //    return a.first;
//  //  });
//  //  int size = int(record.size());
//  //  Branch_Val.resize(size);
//  //  //transform the record to Branch_Val
//  //  transform(record.begin(), record.end(), Branch_Val.begin(), [](const auto &a) {
//  //    return make_pair(a, 0.0);
//  //  });
//  //evaluate the model from pseudo cost
//    std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_pseudo(Branch_pair_from_pseudo.size());
//    transform(Branch_pair_from_pseudo.begin(),
//              Branch_pair_from_pseudo.end(),
//              Branch_Val_model_pseudo.begin(),
//              [](const auto &edge) {
//                return make_pair(edge, 0.0);
//              });
//    ml.predict(Branch_Val_model_pseudo, 1);
//    sort(Branch_Val_model_pseudo.begin(), Branch_Val_model_pseudo.end(), [](const auto &a, const auto &b) {
//      return a.second < b.second;
//    });
//#ifdef Stage1
//    int size1 = min(NumCandidates, int(Branch_Val_model_pseudo.size()));
//#else
//    int tmp_size = node->Idx ? CONFIG::UseModelBranNormalStage4ModelNormal : CONFIG::UseModelBranRootStageModelNormal;
//    int size1 = min(tmp_size, int(Branch_Val_model_pseudo.size())); ?
//#endif
//    Branching_selected_by_model1.resize(size1);
//    transform(Branch_Val_model_pseudo.begin(),
//              Branch_Val_model_pseudo.begin() + size1,
//              Branching_selected_by_model1.begin(),
//              [](const auto &a) {
//                return a.first;
//              });
//    //evaluate the model from fractional
//    std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_fractional(Branch_pair_from_fractional.size());
//    transform(Branch_pair_from_fractional.begin(),
//              Branch_pair_from_fractional.end(),
//              Branch_Val_model_fractional.begin(),
//              [](const auto &edge) {
//                return make_pair(edge, 0.0);
//              });
//    ml.predict(Branch_Val_model_fractional, 2);
//    sort(Branch_Val_model_fractional.begin(), Branch_Val_model_fractional.end(), [](const auto &a, const auto &b) {
//      return a.second < b.second;
//    });
//#ifdef Stage1
//    int size2 = min(NumCandidates, int(Branch_pair_from_fractional.size()));
//#else
//    tmp_size = node->Idx ? CONFIG::UseModelBranNormalStage4ModelPseudo : CONFIG::UseModelBranRootStageModelPseudo;
//    size2 = min(tmp_size, int(Branch_Val_model_fractional.size()));
//#endif
//    Branching_selected_by_model2.resize(size2);
//    transform(Branch_Val_model_fractional.begin(),
//              Branch_Val_model_fractional.begin() + size2,
//              Branching_selected_by_model2.begin(),
//              [](const auto &a) {
//                return a.first;
//              });
//#elif defined( GenerateTrainingData_1 )
//  int size = int(Branch_Val.size());
//#else
//  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
//    return a.second > b.second;
//  });
//  int size = min(CONFIG::BranPhase2, int(Branch_Val.size()));
//#endif
//  Branch_pair.resize(size);
//  transform(Branch_Val.begin(),
//            Branch_Val.begin() + size,
//            Branch_pair.begin(),
//            [](const auto &a) {
//              return a.first;
//            });
//#ifdef AccuracyTest
//  tellAccuracy(1);
//#endif
#ifndef MachineLearning
  int size = min(CONFIG::BranPhase2, int(Branch_Val.size()));
  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second > b.second;
  });
  Branch_pair.resize(size);
  transform(Branch_Val.begin(),
            Branch_Val.begin() + size,
            Branch_pair.begin(),
            [](const auto &a) {
              return a.first;
            });
#endif
  cout << "StageOne Branch_pair size= " << Branch_pair.size() << endl;
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - begin).count();
  cout << "StageOne spent= " << eps << "s" << endl;
}

#ifdef OpenHeurCGInUseModel
void CVRP::doSB_Stage3(BBNODE *node) {
  auto beg = high_resolution_clock::now();
  int BeforeNumRow = NumRow;
  auto org_val = node->Val;
  int cnt = 0;
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);
  int col_start = NumCol;
  auto store_lpoptindex = node->Idx4LPSolsInColPool;
  int edge_idx = 0;
  for (; edge_idx < Branch_pair.size(); ++edge_idx) {
    auto &edge = Branch_pair[edge_idx];
    int ai = edge.first;
    int aj = edge.second;
    cout << MID_PHASE_SEPARATION;
    cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
    int numnz;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_Hyperparameter(checkCST_LIMIT())
    ++NumRow;
    node->BrCs.back().Edge = {ai, aj};
    node->BrCs.back().BrDir = false;
    PoolBeg4Pricing = 0;
    solveLP(node, true, false, false);
    node->if_Int = false;
    node->if_terminated = false;
    double temp_val;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif1 = calculateDif(temp_val, org_val, true);
    cout << SMALL_PHASE_SEPARATION;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    node->BrCs.back().BrDir = true;
    PoolBeg4Pricing = 0;
    solveLP(node, true, false, false);
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val, true);
    auto product = dif1 * dif2;
    cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
         << left << product << endl;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    --NumRow;
    Branch_Val[cnt++] = {edge, product};
    //write data to HeuristicBranch
    HeuristicImprovement_down[edge].first += dif1;
    ++HeuristicImprovement_down[edge].second;
    HeuristicImprovement_up[edge].first += dif2;
    ++HeuristicImprovement_up[edge].second;
//    auto max_it = max_element(Branch_Val.begin(),
//                              Branch_Val.begin() + cnt,
//                              [](const auto &a, const auto &b) {
//                                return a.second < b.second;
//                              });
//    if (max_it - Branch_Val.begin() == cnt - 1) {
//      cout << "The best edge is the " << cnt << " -th!" << endl;
//      cout << "The best edge is ( " << edge.first << " , " << edge.second << " )\t";
//      cout << "The best product is " << max_it->second << endl;
//    } else {
//      break;
//    }
  }
  safe_solver(node->solver.SOLVERreoptimize())
  auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  Branch_pair = {it->first};
  node->BrCs.pop_back();
  node->Idx4LPSolsInColPool = store_lpoptindex;
  node->Val = org_val;
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "StageThree spent= " << eps << "s" << endl;
}
#endif

#ifdef OnlyUseModel1
void CVRP::doSB_Stage2(BBNODE *node) {
  auto beg = high_resolution_clock::now();
  int BeforeNumRow = NumRow;
  auto org_val = node->Val;
  int cnt = 0;
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);
  int col_start = NumCol;
  auto store_lpoptindex = node->Idx4LPSolsInColPool;
  int edge_idx = 0;
  for (; edge_idx < Branch_pair.size(); ++edge_idx) {
    auto &edge = Branch_pair[edge_idx];
    int ai = edge.first;
    int aj = edge.second;
    cout << MID_PHASE_SEPARATION;
    cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
    int numnz;
    double temp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_Hyperparameter(checkCST_LIMIT())
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto LP_dif1 = calculateDif(temp_val, org_val);
    ++NumRow;
    node->BrCs.back().Edge = {ai, aj};
    node->BrCs.back().BrDir = false;
    PoolBeg4Pricing = 0;
    solveLP(node, true, false, false);
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif1 = calculateDif(temp_val, org_val, true);
    cout << SMALL_PHASE_SEPARATION;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    node->BrCs.back().BrDir = true;
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto LP_dif2 = calculateDif(temp_val, org_val);
    PoolBeg4Pricing = 0;
    solveLP(node, true, false, false);
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val, true);
    auto product = dif1 * dif2;
    cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
         << left << product << endl;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    --NumRow;
    Branch_Val[cnt++] = {edge, product};
    //can see LP's info
    LPTestingImprovement_down[edge].first += LP_dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += LP_dif2;
    ++LPTestingImprovement_up[edge].second;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  cout << "Best= " << it - Branch_Val.begin() << endl;
  Branch_pair = {it->first};
  node->BrCs.pop_back();
  node->Idx4LPSolsInColPool = store_lpoptindex;
  node->Val = org_val;

  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "StageTwo spent= " << eps << "s" << endl;
}
#elif defined(Model1_LP_Heuristic)
void CVRP::doSB_Stage2(BBNODE *node) {
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  double org_val = node->Val;
  int BeforeNumRow = NumRow;
  int cnt = 0;
  //resolve LP and obtain the sub_selection
  for (auto &edge : Branch_pair) {
    int ai = edge.first;
    int aj = edge.second;
    int numnz;
    double tmp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    ++NumRow;
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto dif1 = calculateDif(tmp_val, org_val);
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto dif2 = calculateDif(tmp_val, org_val);
    auto product = dif1 * dif2;
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    --NumRow;
    Branch_Val[cnt++] = {edge, product};
    LPTestingImprovement_down[edge].first += dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += dif2;
    ++LPTestingImprovement_up[edge].second;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second > b.second;
  });
  int size = min(CONFIG::BranPhase2, int(Branch_Val.size()));
  Branch_pair.resize(size);
  transform(Branch_Val.begin(),
            Branch_Val.begin() + size,
            Branch_pair.begin(),
            [](const auto &a) {
              return a.first;
            });
  //use heuristic CG to get the exact solution
  Branch_Val.resize(size);
  cnt = 0;
  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);
  int col_start = NumCol;
  auto store_lpoptindex = node->Idx4LPSolsInColPool;
  int edge_idx = 0;
  for (; edge_idx < Branch_pair.size(); ++edge_idx) {
    auto &edge = Branch_pair[edge_idx];
    int ai = edge.first;
    int aj = edge.second;
    cout << MID_PHASE_SEPARATION;
    cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
    int numnz;
    double temp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_Hyperparameter(checkCST_LIMIT())
    ++NumRow;
    node->BrCs.back().Edge = {ai, aj};
    node->BrCs.back().BrDir = false;
    PoolBeg4Pricing = 0;
    solveLP(node, true, false, false);
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif1 = calculateDif(temp_val, org_val, true);
    cout << SMALL_PHASE_SEPARATION;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    node->BrCs.back().BrDir = true;
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    PoolBeg4Pricing = 0;
    solveLP(node, true, false, false);
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val, true);
    auto product = dif1 * dif2;
    cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
         << left << product << endl;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    --NumRow;
    Branch_Val[cnt++].second = product;
    //write data to HeuristicBranch
//    HeuristicImprovement_down[edge].first += dif1;
//    ++HeuristicImprovement_down[edge].second;
//    HeuristicImprovement_up[edge].first += dif2;
//    ++HeuristicImprovement_up[edge].second;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  Branch_pair = {it->first};
}
#else
void CVRP::doSB_Stage2(BBNODE *node) {
  auto beg = high_resolution_clock::now();
  int BeforeNumRow = NumRow;
  auto org_val = node->Val;
#ifdef Stage1
  // run lp testing for initialization of pseudo cost
  std::vector<std::pair<int, int >> Branch_;
  int size1, size2;
  if (Branch_pair_from_fractional.size() > Branch_pair_from_pseudo.size()) {
    size1 = min(NumCandidates_2, int(Branch_pair_from_pseudo.size()));
    Branch_.assign(Branch_pair_from_pseudo.begin(), Branch_pair_from_pseudo.begin() + size1);
    size2 = min(2 * NumCandidates_2 - size1, int(Branch_pair_from_fractional.size()));
    Branch_.insert(Branch_.end(), Branch_pair_from_fractional.begin(), Branch_pair_from_fractional.begin() + size2);
    cout << "pseudo= " << size1 << "  frac= " << size2 << endl;
  } else {
    size1 = min(NumCandidates_2, int(Branch_pair_from_fractional.size()));
    Branch_.assign(Branch_pair_from_fractional.begin(), Branch_pair_from_fractional.begin() + size1);
    size2 = min(2 * NumCandidates_2 - size1, int(Branch_pair_from_pseudo.size()));
    Branch_.insert(Branch_.end(), Branch_pair_from_pseudo.begin(), Branch_pair_from_pseudo.begin() + size2);
    cout << "pseudo= " << size2 << "  frac= " << size1 << endl;
  }
  for (auto &edge : Branch_) {
    int ai = edge.first;
    int aj = edge.second;
    int numnz;
    double tmp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    ++NumRow;
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto left_val = tmp_val;
    auto dif1 = calculateDif(tmp_val, org_val);
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto right_val = tmp_val;
    auto dif2 = calculateDif(tmp_val, org_val);
    auto product = dif1 * dif2;
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    --NumRow;
    LPTestingImprovement_down[edge].first += dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += dif2;
    ++LPTestingImprovement_up[edge].second;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  ml.EdgeTmpInfo.clear();
  return;
#endif
  int cnt = 0;
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
#ifdef UseModel
  transform(Branch_pair.begin(),
            Branch_pair.end(),
            Branch_Val.begin(),
            [&](const auto &edge) {
              return std::make_pair(edge, 0.0);
            });
  int edge_idx = 0;
  for (; edge_idx < Branch_Val.size(); ++edge_idx) {
    auto &edge = Branch_Val[edge_idx].first;
    int ai = edge.first;
    int aj = edge.second;
    int numnz;
    double tmp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    ++NumRow;
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto left_val = tmp_val;
    auto dif1 = calculateDif(tmp_val, org_val);
    ml.collectResolvingFeatures(this, node, edge, BeforeNumRow, tmp_val, org_val, numnz, false);
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    auto right_val = tmp_val;
    auto dif2 = calculateDif(tmp_val, org_val);
    ml.collectResolvingFeatures(this, node, edge, BeforeNumRow, tmp_val, org_val, numnz, true);
    auto product = dif1 * dif2;
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    //collectResolvingFeatures
    if (ceil_transformed_number_related(left_val) >= UB || ceil_transformed_number_related(right_val) >= UB) {
      ml.EdgeTmpInfo[edge].ResolvingLPFeatures.emplace_back("unbalance_left_and_right",
                                                            abs(left_val - right_val) / (left_val + right_val));
    } else {
      ml.EdgeTmpInfo[edge].ResolvingLPFeatures.emplace_back("unbalance_left_and_right", 0);
    }
    //
    --NumRow;
    LPTestingImprovement_down[edge].first += dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += dif2;
    ++LPTestingImprovement_up[edge].second;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  ml.predict(Branch_Val, 3);
  auto it = min_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
#ifdef OpenHeurCGInUseModel
  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  int size = min(CONFIG::MaxUseModelBranPhase2, int(Branch_Val.size()));
  Branch_pair.resize(size);
  transform(Branch_Val.begin(),
            Branch_Val.begin() + size,
            Branch_pair.begin(),
            [](const auto &a) {
              return a.first;
            });
#else
  Branch_pair = {it->first};
#endif
#ifdef AccuracyTest
  //sort Branch_Val and give the real rank!
  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  cout << "phase2 real rank: ";
  for (auto &i : Branch_Val) {
    cout << Map_Edge_SB_rank[i.first].second << " ";
  }
  cout << endl;
  tellAccuracy(2);
  for (auto &map_it : Map_Edge_SB_rank) {
    if (map_it.second.second == 1) {
      cout << "choose the best!" << endl;
      Branch_pair = {map_it.first};
      break;
    }
  }
#endif
#else
  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);
  int col_start = NumCol;
  auto store_lpoptindex = node->Idx4LPSolsInColPool;
  int edge_idx = 0;
  for (; edge_idx < Branch_pair.size(); ++edge_idx) {
    auto &edge = Branch_pair[edge_idx];
    int ai = edge.first;
    int aj = edge.second;
    cout << MID_PHASE_SEPARATION;
    cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
    int numnz;
    double temp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_Hyperparameter(checkCST_LIMIT())
    ++NumRow;
#ifdef GenerateTrainingData_2
    double tmp_val;
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    double left_val = tmp_val;
    ml.collectResolvingFeatures(this, node, edge, BeforeNumRow, tmp_val, org_val, numnz, false);
#endif
    node->BrCs.back().Edge = {ai, aj};
    node->BrCs.back().BrDir = false;
    PoolBeg4Pricing = 0;
#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2) || defined(Stage1)
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto LP_dif1 = calculateDif(temp_val, org_val);
    ForceNotRollback = true;
    if_force_not_regenerate_bucket_graph = true;
    solveLP(node, true, true, false);
    ForceNotRollback = false;
    if_force_not_regenerate_bucket_graph = false;
#else
    solveLPInLabeling(node, true, false, false);
#endif
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
#ifdef GenerateTrainingData_2
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
    double right_val = tmp_val;
    ml.collectResolvingFeatures(this, node, edge, BeforeNumRow, tmp_val, org_val, numnz, true);
    if (ceil_transformed_number_related(left_val) >= UB || ceil_transformed_number_related(right_val) >= UB) {
      ml.EdgeTmpInfo[edge].ResolvingLPFeatures.emplace_back("unbalance_left_and_right",
                                                            abs(left_val - right_val) / (left_val + right_val));
    } else {
      ml.EdgeTmpInfo[edge].ResolvingLPFeatures.emplace_back("unbalance_left_and_right", 0);
    }
#endif
    PoolBeg4Pricing = 0;
#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2) || defined(Stage1)
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto LP_dif2 = calculateDif(temp_val, org_val);
    ForceNotRollback = true;
    if_force_not_regenerate_bucket_graph = true;
    solveLP(node, true, true, false);
    ForceNotRollback = false;
    if_force_not_regenerate_bucket_graph = false;
#else
    solveLPInLabeling(node, true, false, false);
#endif
    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val);
    auto product = dif1 * dif2;
    cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
         << left << product << endl;
#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2)
    ml.EdgeTmpInfo[edge].SB_scores = product;
#endif
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    --NumRow;
    Branch_Val[cnt++] = {edge, product};
    //write data to HeuristicBranch
#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2) || defined(Stage1)
    LPTestingImprovement_down[edge].first += LP_dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += LP_dif2;
    ++LPTestingImprovement_up[edge].second;
#else
    HeuristicImprovement_down[edge].first += dif1;
    ++HeuristicImprovement_down[edge].second;
    HeuristicImprovement_up[edge].first += dif2;
    ++HeuristicImprovement_up[edge].second;
#endif
  }
  safe_solver(node->solver.SOLVERreoptimize())
  auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  Branch_pair = {it->first};
  node->BrCs.pop_back();
  node->Idx4LPSolsInColPool = store_lpoptindex;
  node->Val = org_val;
#endif
#ifdef GenerateTrainingData_1
  ml.writeTrainingLPFile();
#elif defined(GenerateTrainingData_2)
  ml.writeTrainingExactFile();
#endif
#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2)
  ml.EdgeTmpInfo.clear();
  ++ml.QID;
#endif
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "StageTwo spent= " << eps << "s" << endl;
}
#endif

void CVRP::doSB_enumeration(BBNODE *node, pair<int, int> &info) {
  if (node->if_terminated) return;
  ++BranchTimes;

  Branch_pair.clear();
  auto tmp_phase1 = CONFIG::BranPhase1;
  CONFIG::BranPhase1 = CONFIG::BranPhase2;
  doSB_Stage0(node);
  CONFIG::BranPhase1 = tmp_phase1;

  //stage 1 in enu
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  double org_val = node->Val;
  int BeforeNumRow = NumRow;
  int cnt = 0;
  for (auto &edge : Branch_pair) {
    int ai = edge.first;
    int aj = edge.second;
    int numnz;
    double temp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif1 = calculateDif(temp_val, org_val);
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val);
    auto product = dif1 * dif2;
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    Branch_Val[cnt++] = {edge, product};
    LPTestingImprovement_down[edge].first += dif1;
    ++LPTestingImprovement_down[edge].second;
    LPTestingImprovement_up[edge].first += dif2;
    ++LPTestingImprovement_up[edge].second;
  }

  safe_solver(node->solver.SOLVERreoptimize())
  auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
    return a.second < b.second;
  });

  info = it->first;
  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++BranchChoice[info];
}

//we don't collect data from enumeration branch thus we don't need to write data to file
//enumeration branch write another version of this function



#ifdef AccuracyTest
void CVRP::generateMap_Edge_SB_rankByExactCG(BBNODE *node) {
  Map_Edge_SB_rank.clear();
  getInfoEdge(node, true);
  int BeforeNumRow = NumRow;
  auto org_val = node->Val;
  BrC bf;
  bf.IdxBrC = NumRow;
  node->BrCs.emplace_back(bf);
  int col_start = NumCol;

  auto store_lpoptindex = node->Idx4LPSolsInColPool;
  int edge_idx = 0;
  unordered_set<pair<int, int>, PairHasher> eval_branching_candidates;
  eval_branching_candidates.reserve(Branching_selected_by_model1.size() + Branching_selected_by_model2.size());
  eval_branching_candidates.insert(Branching_selected_by_model1.begin(), Branching_selected_by_model1.end());
  eval_branching_candidates.insert(Branching_selected_by_model2.begin(), Branching_selected_by_model2.end());
  auto beg = high_resolution_clock::now();
  for (auto &edge : eval_branching_candidates) {
    int ai = edge.first;
    int aj = edge.second;
    cout << MID_PHASE_SEPARATION;
    cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
    int numnz;
    double temp_val;
    getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
    safe_solver(addBranchconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_Hyperparameter(checkCST_LIMIT())
    ++NumRow;
    node->BrCs.back().Edge = {ai, aj};
    node->BrCs.back().BrDir = false;
    PoolBeg4Pricing = 0;

    ForceNotRollback = true;
    if_force_not_regenerate_bucket_graph = true;
    solveLP(node, true, true, false);
    ForceNotRollback = false;
    if_force_not_regenerate_bucket_graph = false;

    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif1 = calculateDif(temp_val, org_val, true);
    cout << SMALL_PHASE_SEPARATION;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    node->BrCs.back().BrDir = true;
    PoolBeg4Pricing = 0;
    ForceNotRollback = true;
    if_force_not_regenerate_bucket_graph = true;
    solveLP(node, true, true, false);
    ForceNotRollback = false;
    if_force_not_regenerate_bucket_graph = false;

    node->if_Int = false;
    node->if_terminated = false;
    safe_solver(node->solver.SOLVERgetObjVal(&temp_val))
    auto dif2 = calculateDif(temp_val, org_val, true);
    auto product = dif1 * dif2;
    cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
         << left << product << endl;
    safe_solver(node->solver.SOLVERdelvars(NumCol - col_start, const_for_branching + col_start))
    safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
    --NumRow;
    Map_Edge_SB_rank[edge] = {product, 0};
  }
  vector<pair<pair<int, int>, pair<double, int>>> vec_map(Map_Edge_SB_rank.begin(), Map_Edge_SB_rank.end());
  sort(vec_map.begin(), vec_map.end(), [](const auto &a, const auto &b) {
    return a.second.first > b.second.first;
  });
  int cnt = 1;
  max_SB_scores = vec_map[0].second.first;
  double std_sb_score = max_SB_scores * AccuracyTest_alpha;
  for (auto &tmp : vec_map) {
    if (tmp.second.first < std_sb_score) ++cnt;
    Map_Edge_SB_rank[tmp.first].second = cnt;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  node->BrCs.pop_back();
  node->Idx4LPSolsInColPool = store_lpoptindex;
  node->Val = org_val;
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "generate Map= " << eps << "s" << endl;
  Branch_pair = {vec_map[0].first};
};

void CVRP::tellAccuracy(int node_dep, int phase) {
  //find the minimum rank from Map_Edge_SB_rank
#ifdef Stage1
// 2 models, check candidates selected by model 1 and model 2, print their rank and ratio, more detail, better.
// later information can be checked by python
  write2Files(node_dep);
#else
  int min_rank = INT_MAX;
  double ratio;
  vector<int> rank_vec(Branch_pair.size());
  int cnt = 0;
  for (auto &edge : Branch_pair) {
    auto rank = Map_Edge_SB_rank[edge].second;
    rank_vec[cnt++] = rank;
    if (rank < min_rank) {
      min_rank = rank;
      ratio = Map_Edge_SB_rank[edge].first / max_SB_scores;
    }
  }
  if (phase == 1) {
    rank_phase1 = min_rank;
    average_rank_phase1.first += min_rank;
    average_rank_phase1.second++;
    average_ratio_phase1.first += ratio;
    average_ratio_phase1.second++;
  } else if (phase == 2) {
    rank_phase2 = min_rank;
    average_rank_phase2.first += min_rank;
    average_rank_phase2.second++;
    average_ratio_phase2.first += ratio;
    average_ratio_phase2.second++;
  };
  cout << "Best rank of phase" << phase << "= " << min_rank << ", ratio= " << ratio << endl;
  cout << "ordered rank and corresponding alpha: " << endl;
  for (auto &edge : Branch_pair) {
    cout << "(" << Map_Edge_SB_rank[edge].second << "," << Map_Edge_SB_rank[edge].first / max_SB_scores << ")" << " ";
  }
  cout << endl;
#endif
}

void CVRP::write2Files(int node_dep) {
  // if AccuracyTest_outputDir does not exist, create it
  if (!std::filesystem::exists(AccuracyTest_outputDir)) {
    std::filesystem::create_directory(AccuracyTest_outputDir);
  }
  ofstream fout;
  string out_path = AccuracyTest_outputDir + "/" + FileName + ".txt";
  // 2 models, check candidates selected by model 1 and model 2, print their rank and ratio, more detail, better.
  // later information can be checked by python
  fout.open(out_path, ios::app);
  fout << "node_dep= " << node_dep << endl;
  fout << "max_SB_scores= " << max_SB_scores << endl;
  fout << "model_1= " << endl;
  for (auto &edge : Branching_selected_by_model1) {
    fout << "(" << Map_Edge_SB_rank[edge].second << "," << Map_Edge_SB_rank[edge].first << ")" << " ";
  }
  fout << endl;
  fout << "model_2= " << endl;
  for (auto &edge : Branching_selected_by_model2) {
    fout << "(" << Map_Edge_SB_rank[edge].second << "," << Map_Edge_SB_rank[edge].first << ")" << " ";
  }
  fout << endl;
}
#endif