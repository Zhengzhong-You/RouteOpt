//
// Created by Zhengzhong You on 5/31/22.
//

#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;

void CVRP::sepRCCs(BBNODE *&node) {
  int numnz;

  CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

  CMGR_CreateCMgr(&MyCutsCMP, 100);
  CMGR_CreateCMgr(&MyOldCutsCMP, 100);

  int oldNum = NumRow;
  while (true) {
    getInfoEdge(node, false);
    if (node->if_Int) {
      break;
    }

    //Because the separation functor considers n+1 as the depot, we do the corresponding converse
    for (int i = 1; i <= node->NumEdges; ++i) {
      if (!node->EdgeTail[i]) {
        node->EdgeTail[i] = Dim;
      } else break;
    }

    CAPSEP_SeparateCapCuts(RealDim, Demand, Cap, node->NumEdges, node->EdgeTail,
                           node->EdgeHead, node->EdgeVal, MyOldCutsCMP,
                           MAXNOOFCUTS, TOLERANCE, TOLERANCE,
                           &if_IntNFeasible, &MaxVio, MyCutsCMP);

    //Because the separation functor considers n+1 as the depot, we do the corresponding converse
    for (int i = 1; i <= node->NumEdges; ++i) {
      if (node->EdgeTail[i] == Dim) {
        node->EdgeTail[i] = 0;
      } else break;
    }

    if (!MyCutsCMP->Size) break;
    //Read cuts and add them to lp
    int cnt = 0;
    for (int i = 0; i < MyCutsCMP->Size; ++i) {
      RCC rcc;
      auto &tmp_customerInfo = rcc.InfoRCCCustomer;
      for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
        tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
      }

      if (tmp_customerInfo.size() <= Dim / 2) {
        rcc.FormRCC = true;
        rcc.RHS = MyCutsCMP->CPL[i]->RHS;
      } else {
        rcc.FormRCC = false;
        rcc.RHS = RealDim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
        auto &tmp_NoCustomerInfo = rcc.InfoRCCOutsideCustomer;
        vector<bool> tmp(Dim, false);
        for (int j : tmp_customerInfo) {
          tmp[j] = true;
        }
        for (int j = 1; j < Dim; ++j) {
          if (!tmp[j]) {
            tmp_NoCustomerInfo.emplace_back(j);
          }
        }
      }

      if (std::find(node->RCCs.begin(), node->RCCs.end(), rcc) != node->RCCs.end()) {
        continue;
      }

      ++cnt;
      rcc.IdxRCC = NumRow;
      node->RCCs.emplace_back(rcc);

      getCoeffRCC(node, rcc, solver_ind, solver_val, numnz);
      safe_solver(
          node->solver.SOLVERaddconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, rcc.RHS, nullptr))
      safe_solver(node->solver.SOLVERupdatemodel())
      safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
      safe_Hyperparameter(checkCST_LIMIT())
    }

    if (!cnt) break;

    for (int i = 0; i < MyCutsCMP->Size; ++i) {
      CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
    }

#ifdef VERBOSE
    cout << "Cuts added...  rcc= " << cnt << endl;
#endif
    MyCutsCMP->Size = 0;

    solveLPInLabeling(node);

    if (node->if_terminated) {
      delete node;
      node = nullptr;
      goto QUIT;
    } else {
      eliminateArcs(node);
      enumerateMIP(node);
      if (!node) goto QUIT;
      cleanIdxCol4Node(node, node->NumParentCols, true);
    }
    LB = node->Val;
    LB_transformed = ceil_transformed_number_related(LB - TOLERANCE);
    if (ceil_transformed_number_related(node->Val - TOLERANCE) + TOLERANCE >= UB || node->if_Int) {
      node->if_terminated = true;
      cout << TERMINATED_MESSAGE_SEP_RCC;
      break;
    }

#ifdef VERBOSE
    cout << "gap= " << (double(UB - LB) / UB > TOLERANCE ? double(UB - LB) / UB * 100 : 0) << "%\n";
    cout << SMALL_PHASE_SEPARATION;
#endif
  }

  QUIT:
  if (node) {
//    deleteNewAddedNonActiveCutsBySlack(node, oldNum);
    findNonactiveCuts(node);
    convertVertex2R1CsInOneLP(node);
    solveLPInLabeling(node);
  }
  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  CMGR_FreeMemCMgr(&MyCutsCMP);
}

void CVRP::generateRCCs(BBNODE *node) {
  int numnz, cnt = 0;

  CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
  CMGR_CreateCMgr(&MyCutsCMP, 100);
  CMGR_CreateCMgr(&MyOldCutsCMP, 100);

  getInfoEdge(node, false);
  //Because the separation functor considers n+1 as the depot, we do the corresponding converse
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (!node->EdgeTail[i]) {
      node->EdgeTail[i] = Dim;
    } else break;
  }

  CAPSEP_SeparateCapCuts(RealDim, Demand, Cap, node->NumEdges, node->EdgeTail,
                         node->EdgeHead, node->EdgeVal, MyOldCutsCMP,
                         MAXNOOFCUTS, TOLERANCE, TOLERANCE,
                         &if_IntNFeasible, &MaxVio, MyCutsCMP);
  //Because the separation functor considers n+1 as the depot, we do the corresponding converse
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeTail[i] == Dim) {
      node->EdgeTail[i] = 0;
    } else break;
  }

  if (!MyCutsCMP->Size) goto QUIT;
  cnt = 0;
  //Read cuts and add them to lp
  for (int i = 0; i < MyCutsCMP->Size; ++i) {
    RCC rcc;
    auto &tmp_customerInfo = rcc.InfoRCCCustomer;
    for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
      tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
    }

    if (tmp_customerInfo.size() <= Dim / 2) {
      rcc.FormRCC = true;
      rcc.RHS = MyCutsCMP->CPL[i]->RHS;
    } else {
      rcc.FormRCC = false;
      rcc.RHS = RealDim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
      auto &tmp_NoCustomerInfo = rcc.InfoRCCOutsideCustomer;
      vector<bool> tmp(Dim, false);
      for (int j : tmp_customerInfo) tmp[j] = true;
      for (int j = 1; j < Dim; ++j) {
        if (!tmp[j]) {
          tmp_NoCustomerInfo.emplace_back(j);
        }
      }
    }

    if (std::find(node->RCCs.begin(), node->RCCs.end(), rcc) != node->RCCs.end()) {
      continue;
    }

    ++cnt;

    rcc.IdxRCC = NumRow;
    node->RCCs.emplace_back(rcc);

    getCoeffRCC(node, rcc, solver_ind, solver_val, numnz);
    safe_solver(node->solver.SOLVERaddconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, rcc.RHS, nullptr))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
    safe_Hyperparameter(checkCST_LIMIT())
  }
  QUIT:
  for (int i = 0; i < MyCutsCMP->Size; ++i) {
    CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
  }
  MyCutsCMP->Size = 0;

  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  CMGR_FreeMemCMgr(&MyCutsCMP);
}

///---------------------------------------------------------------------------------------------------------------------

void CVRP::sepHybridCuts(BBNODE *&node) {
  //For CVRP, we consider RCC as well as other robust cuts and SR3
  int cnt_tail_off = 0;
  double standard;
  int count = 0;
  int RCCnum, R1Cnum;
  bool if_opt = false;
  std::chrono::high_resolution_clock::time_point t1, t2;
  double duration;
//  int aggressive_level;// right now, the policy is already the best one!
  //1. add RCC 2. if validation of rcc < 0.1, it shows tail-off

  if (If_in_Enu_State) goto HYBRID;
  if (!node->Idx) {
#ifdef SOLVER_VRPTW
    augmentNGRound(node);
#endif

    cout << "Begin RCC separation...\n";

    sepRCCs(node);

    cout << MID_PHASE_SEPARATION;
    cout << "RCC tail-off!\nSwitch on Hybrid...\n";

    if (!node || node->if_terminated) {
      if_opt = true;
      goto QUIT;
    }
    safe_Hyperparameter(Rollback)
    cnt_tail_off = -CONFIG::ExtraCuts4Root;

    setTailOffStd_N_RollBackStd();
  } else {
#ifdef VERBOSE
    cout << "Begin Hybrid separation...\n";
#endif
  }

  HYBRID:
//  aggressive_level = 0;
  while (true) {

    bool if_soft_limit = false;
    //record old cuts info
    //because this is extremely unusual, we seek the way of minimum damage
//    auto old_rcc = node->RCCs;
//    auto old_r1c = node->R1Cs;
//    auto old_r1c_multi = node->R1Cs_multi;
//    auto old_brc = node->BrCs;
//    vector<size_t> old_col_idx(node->IdxCols, node->IdxCols + NumCol);

    ++count;
    cout << BIG_PHASE_SEPARATION;
    cout << "Hybrid separation round " << count << endl;

    double prior_nodeVal = node->Val;

    int oldNum = NumRow;

    auto beg = chrono::high_resolution_clock::now();

    if (If_in_Enu_State) {
      generateRCCsInEnu(node);
    } else {
      generateRCCs(node);
    }

//#ifndef NewWay2AddCuts
//    generateAllR1Cs(node, false);
//#else
////    if (cnt_tail_off > 0 && aggressive_level < 2) {
////      cnt_tail_off = 0;
////      ++aggressive_level;
////    }
////    NewgenerateAllR1Cs(node, aggressive_level);
//#endif

    NewgenerateAllR1Cs(node, 2);

#ifdef checkR1C_Mul_State
    checkR1C_mul_look(node, oldNum);
#endif

    int cuts_sum = NumRow - oldNum;
    if (!cuts_sum) {
      cout << "sep fails due to zero cuts" << endl;
      if (if_fill_mem) {
        cout << "fill mem! resolve the cg!" << endl;
        if_fill_mem = false;
        cnt_tail_off = numeric_limits<int>::max();
        goto SOLVE_CG;
      }
      goto QUIT;
    }
    cout << "cuts_sum: " << cuts_sum << endl;
    t1 = chrono::high_resolution_clock::now();
    safe_solver(node->solver.SOLVERreoptimize())
    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration<double>(t2 - t1).count();
    cout << "re-optimize time: " << duration << endl;

//    deleteNonActiveCutsByDual(node);
//    deleteNewAddedNonActiveCutsBySlack(node, oldNum, true);
//    deleteNewAddedActiveCutsByDual_N_Mem(node, oldNum);
//    deleteNonActiveCutsBySlack(node, false, true);

    deleteNonActiveCutsSafely(node, oldNum);

#ifdef checkR1C_Mul_State
    checkR1C_mul_look(node, oldNum);
#endif

    printCutsInfo(node);

    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration<double>(t2 - t1).count();
    cout << "Hybrid separation time: " << duration << endl;

#ifdef checkR1CTotally
    checkR1Ctotally(node);
#endif

    SOLVE_CG:
    if (If_in_Enu_State) {
      solveLPByInspection(node, false, false, true);
      if (!node->Idx) updateLB(node->Val);
      if (node->if_terminated)goto QUIT;
//      findNonactiveCuts(node);///forbidden here!
      if (node->SizeEnuColPool + NumCol <= MaxNumRoute4MIP) goto QUIT;
    } else {
      convertVertex2R1CsInOneLP(node);
      solveLPInLabeling(node);
      if (!node->Idx) updateLB(node->Val);
      if (Rollback == 1) {
        cout << BIG_PHASE_SEPARATION;
        rollbackEaseWay(node, oldNum);
//      rollback2PreState(node, old_rcc, old_r1c, old_r1c_multi, old_brc, old_col_idx);
        cout << BIG_PHASE_SEPARATION;
        goto QUIT;
      } else if (Rollback == 3) {
        if_soft_limit = true;
      }

      if (node->if_terminated) {
        delete node;
        node = nullptr;
        goto QUIT;
      }
      eliminateArcs(node);
      enumerateMIP(node);
      if (!node) goto QUIT;
      cleanIdxCol4Node(node, node->NumParentCols, true);
      findNonactiveCuts(node);
    }
    cout << "after deleting those inactivate cuts! We left with cols: " << NumCol << " rows: " << NumRow << endl;
#ifdef checkR1C_Mul_State
    checkR1C_mul_look(node, 0);
#endif

    standard = calculateGapImprovement(node->Val, prior_nodeVal);
    cout << "local gap= " << (double(UB - node->Val) / UB > TOLERANCE ? double(UB - node->Val) / UB * 100 : 0) << endl;
    cout << "gap improved= " << standard << endl;
    cout << SMALL_PHASE_SEPARATION;
    if (UB <= LB_transformed) {
      if_opt = true;
      goto QUIT;
    }
    if (if_soft_limit) {
      cout << "labels reach the soft limit!" << endl;
      goto QUIT;
    }
    if (standard < TOLERANCE) goto QUIT;
    else if (standard < CONFIG::CutsTailOff[0]) {
      ++cnt_tail_off;
      cout << "cnt_tail_off= " << cnt_tail_off << endl;
      if (cnt_tail_off >= CONFIG::CutsTolerance)goto QUIT;
    }
  }

  QUIT:
  if (!If_in_Enu_State) {
    if (node && !if_opt) {
      convertVertex2R1CsInOneLP(node);
      cout << "we clean the col pool!" << endl;
      cleanIdxCol4Node(node, node->NumParentCols);
#ifdef openCutsAndEnumerationAtEachEndNode
      cout << "we resolve the LP and we force not generate new bucket graph!" << endl;
    if_force_not_regenerate_bucket_graph = true;
    solveLPInLabeling(node);
    if_force_not_regenerate_bucket_graph = false;
#endif
      cout << "Hybrid tali-off! " << "Stop cut-separation!" << " Last column reduction... ncol= " << NumCol
           << endl;
      cout << BIG_PHASE_SEPARATION;
    } else {
#ifdef VERBOSE
      cout << "Terminate the node!\n" << "Stop cut-separation!" << endl;
      cout << SMALL_PHASE_SEPARATION;
#endif
    }
  } else {
    cout << "Hybrid tali-off! " << "Stop cut-separation!" << " Last column reduction... ncol= " << NumCol
         << endl;
    cout << BIG_PHASE_SEPARATION;
  }
  if (node) {
    recordOptCol(node);//very essential
    if (!node->Idx) updateLB(node->Val);
  }
#ifdef check_if_cuts_are_legal
  if (node) checkIfCutsLegal(node);
#endif
}

void CVRP::generateR1C3s(const vector<vector<int>> &ele_routes,
                         const vector<vector<int>> &non_ele_routes,
                         const vector<double> &frac_ele_routes,
                         const vector<double> &frac_non_ele_routes,
                         vector<pair<vector<int>, double>> &cut_info) const {

  int num_cut = 0;
  int ai, aj, ak;
  int num_ele_routes = (int) ele_routes.size(), num_non_ele_routes = (int) non_ele_routes.size();
  int cnt;
  double frac;
  yzzLong tmp_long, PI_ai, PI_aj, PI_ai_aj;
  vector<int> aux(3);

  vector<vector<int>> vertex2non_ele_routes(Dim);
  for (int i = 0; i < Dim; ++i) vertex2non_ele_routes.reserve(num_non_ele_routes);

  unordered_map<yzzLong, tuple<int, int, int, double>> three_combinations;
  three_combinations.reserve(int(pow_self(MaxLengthEleRoute, 3) / 6));
  unordered_map<yzzLong, tuple<int, int, int, double >>::iterator
      iter3;
  unordered_map<yzzLong, tuple<int, int, double >> two_combinations;
  two_combinations.reserve(int(pow_self(MaxLengthEleRoute, 2) / 2));
  unordered_map<yzzLong, tuple<int, int, double >>::iterator
      iter2;
  unordered_map<yzzLong, tuple<int, int, int, double>> all_combinations;
  all_combinations.reserve(int(pow_self(MaxLengthEleRoute, 3) / 6));
  unordered_map<yzzLong, tuple<int, int, int, double >>::iterator
      all_iter;

  for (int i = 0; i < non_ele_routes.size(); ++i) {
    for (auto j : non_ele_routes[i]) {
      vertex2non_ele_routes[j].emplace_back(i);
    }
  }

  //find all possible combinations in ele-routes
  for (int i = 0; i < num_ele_routes; ++i) {
    auto &seq = ele_routes[i];
    cnt = (int) ele_routes[i].size();
    frac = frac_ele_routes[i];
    //write data into two combinations
    for (int j = 0; j < cnt; ++j) {
      for (int k = j + 1; k < cnt; ++k) {
        tmp_long = 0;
        tmp_long.set(seq[j]);
        tmp_long.set(seq[k]);
        iter2 = two_combinations.find(tmp_long);
        if (iter2 == two_combinations.end()) {
          if (seq[j] < seq[k])two_combinations.emplace(tmp_long, make_tuple(seq[j], seq[k], frac));
          else two_combinations.emplace(tmp_long, make_tuple(seq[k], seq[j], frac));
        } else get<2>(iter2->second) += frac;
      }
    }
    //write data into three combinations
    for (int j = 0; j < cnt; ++j) {
      for (int k = j + 1; k < cnt; ++k) {
        for (int l = k + 1; l < cnt; ++l) {
          tmp_long = 0;
          tmp_long.set(seq[j]);
          tmp_long.set(seq[k]);
          tmp_long.set(seq[l]);
          iter3 = three_combinations.find(tmp_long);
          if (iter3 == three_combinations.end()) {
            //sort
            aux[0] = seq[j];
            aux[1] = seq[k];
            aux[2] = seq[l];
            std::stable_sort(aux.begin(), aux.end());
            three_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
          } else get<3>(iter3->second) += frac;
        }
      }
    }
  }
  //e.g. find 1,2 and 1,3. we use 1,2 + 1,3 + 2,3 -2* 1,2,3 to determine the violation (of course if 2,3 and 1,2,3 exist)
  for (auto &two_comb : two_combinations) {
    ai = get<0>(two_comb.second);
    aj = get<1>(two_comb.second);
    PI_ai_aj = 0;
    PI_ai_aj.set(ai);
    PI_ai_aj.set(aj);
    PI_ai = 0;
    PI_aj = 0;
    PI_ai.set(ai);
    PI_aj.set(aj);
    //find pair
    for (int j = 1; j < ai; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      //find if counted already
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      //find 1,3
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 2,3
      tmp_long = PI_aj;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 1,2,3
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      iter3 = three_combinations.find(tmp_long);
      if (iter3 != three_combinations.end()) frac -= 2 * get<3>(iter3->second);
      aux[0] = ai;
      aux[1] = aj;
      aux[2] = j;
      std::stable_sort(aux.begin(), aux.end());
      all_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
    }
    for (int j = ai + 1; j < aj; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      //find if counted already
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      //find 1,3
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 2,3
      tmp_long = PI_aj;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 1,2,3
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      iter3 = three_combinations.find(tmp_long);
      if (iter3 != three_combinations.end()) frac -= 2 * get<3>(iter3->second);
      //if frac
      aux[0] = ai;
      aux[1] = aj;
      aux[2] = j;
      std::stable_sort(aux.begin(), aux.end());
      all_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
    }
    for (int j = aj + 1; j < Dim; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      //find if counted already
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      //find 1,3
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 2,3
      tmp_long = PI_aj;
      tmp_long.set(j) = 1;
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 1,2,3
      tmp_long = PI_ai_aj;
      tmp_long.set(j) = 1;
      iter3 = three_combinations.find(tmp_long);
      if (iter3 != three_combinations.end()) frac -= 2 * get<3>(iter3->second);
      aux[0] = ai;
      aux[1] = aj;
      aux[2] = j;
      std::stable_sort(aux.begin(), aux.end());
      all_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
    }
  }
  //traverse the combinations to calculate the final violations (check if there is a need for optimization)

  vector<double> tmp_com_info_non_ele(num_non_ele_routes);
  vector<tuple<int, int, int, double>> violated_combs;
  violated_combs.reserve(all_combinations.size());
  for (auto &all_comb : all_combinations) {
    auto &comb_frac = get<3>(all_comb.second);
    memset(tmp_com_info_non_ele.data(), 0, sizeof(double) * num_non_ele_routes);
    ai = get<0>(all_comb.second);
    aj = get<1>(all_comb.second);
    ak = get<2>(all_comb.second);
    for (auto item : vertex2non_ele_routes[ai])++tmp_com_info_non_ele[item];
    for (auto item : vertex2non_ele_routes[aj])++tmp_com_info_non_ele[item];
    for (auto item : vertex2non_ele_routes[ak])++tmp_com_info_non_ele[item];
    for (int i = 0; i < tmp_com_info_non_ele.size(); ++i)
      comb_frac += floor(tmp_com_info_non_ele[i] / 2 + TOLERANCE) * frac_non_ele_routes[i];
    if (comb_frac - 1 > TOLERANCE)violated_combs.emplace_back(ai, aj, ak, comb_frac);
  }

  std::stable_sort(violated_combs.begin(),
                   violated_combs.end(),
                   [](const std::tuple<int, int, int, double> &a, const std::tuple<int, int, int, double> &b) {
                     return get<3>(a) > get<3>(b);
                   });
  int limit = min((int) violated_combs.size(), CONFIG::MaxNumR1C3PerRound);
  transform(violated_combs.begin(), violated_combs.begin() + limit, back_inserter(cut_info),
            [](const std::tuple<int, int, int, double> &a) {
              return make_pair(vector<int>{get<0>(a), get<1>(a), get<2>(a)}, get<3>(a) - 1);
            });
}

void CVRP::generateR1C1s(const vector<vector<int>> &non_ele_routes,
                         const vector<double> &frac_none_ele_routes,
                         vector<pair<vector<int>, double>> &cut_info) const {
//  return;

  yzzLong tmp_long;
  set<int> re_visited_vertices;
  unordered_map<int, set<int>> vertex2non_ele_routes;
  vertex2non_ele_routes.reserve(Dim);

  for (int r = 0; r < non_ele_routes.size(); ++r) {
    tmp_long = 0;
    for (auto i : non_ele_routes[r]) {
      if (tmp_long[i]) {
        re_visited_vertices.emplace(i);
        vertex2non_ele_routes[i].emplace(r);
      } else tmp_long.set(i);
    }
  }

  vector<pair<int, double >> vio_mem_vec_re_vertices(re_visited_vertices.size());//v & vio
  int cnt = 0;
  for (auto i : re_visited_vertices) {
    double vio = 0;
    for (auto r : vertex2non_ele_routes[i]) {
      double times = 0;
      for (auto node : non_ele_routes[r]) {
        if (node == i)++times;
      }
      vio += floor(times / 2 + TOLERANCE) * frac_none_ele_routes[r];
    }
    vio_mem_vec_re_vertices[cnt++] = {i, vio};
  }

  std::stable_sort(vio_mem_vec_re_vertices.begin(),
                   vio_mem_vec_re_vertices.end(),
                   [](const pair<int, double> &a, const pair<int, double> &b) {
                     return a.second > b.second;
                   });

  int limit = min((int) vio_mem_vec_re_vertices.size(), CONFIG::MaxNumR1CPerRound);
  transform(vio_mem_vec_re_vertices.begin(), vio_mem_vec_re_vertices.begin() + limit, back_inserter(cut_info),
            [](const pair<int, double> &a) {
              return make_pair(vector<int>{a.first}, a.second);
            });
}

void CVRP::generateHighDimR1Cs(const vector<vector<int>> &routes,
                               const vector<double> &frac_routes,
                               const vector<int> &cut_type,
                               vector<pair<vector<int>, double>> &cut_info) {
// we generate high dimensional rank-1 cuts by solving a MIP
// solver
  int num_routes = (int) routes.size();
  SOLVER solver{};
  solver.SOLVERgetenv(&Solver);
  int numrow = num_routes + 1;
  int numcol = num_routes + RealDim;
  //the first node->Idx4LPSolsInColPool.size() cols are fractional solutions, and the next RealDim cols are vertices
  vector<double> obj(numcol, 0);
  for (int i = 0; i < num_routes; ++i) obj[i] = frac_routes[i];
  vector<char> xtype(numcol, SOLVER_BINARY);
  std::fill(xtype.begin(), xtype.begin() + num_routes, SOLVER_INTEGER);

  safe_solver(solver.SOLVERnewmodel("HighDimensionalR1CsMIP",
                                    numcol,
                                    obj.data(),
                                    nullptr,
                                    nullptr,
                                    xtype.data(),
                                    nullptr))
  safe_solver(solver.SOLVERsetModelSense(SOLVER_MAX_SENSE))
  // add variables type
  vector<char> sense(numrow, SOLVER_LESS_EQUAL);
  sense[numrow - 1] = SOLVER_EQUAL;
  vector<double> rhs(numrow, 0);
  //add constraints
  size_t numnz = 0;
  int ccnt = 0;
  vector<double> dim(Dim, 0);
  int num_first_batch_col = num_routes - 1;
  for (auto &r : routes) {
    solver_beg[ccnt] = numnz;
    solver_ind[numnz] = ccnt;
    solver_val[numnz++] = 1;
    ccnt++;
    for (auto i : r) ++dim[i];
    for (int j = 1; j < Dim; ++j) {
      if (dim[j] != 0) {
        solver_ind[numnz] = num_first_batch_col + j;
        solver_val[numnz++] = -dim[j] / 2;
      }
    }
    memset(dim.data(), 0, Dim * sizeof(double));
  }
  solver_beg[numrow - 1] = numnz;
  for (int i = num_first_batch_col + 1; i < numcol; ++i) {
    solver_ind[numnz] = i;
    solver_val[numnz++] = 1;
  }
  solver_beg[numrow] = numnz;
  safe_solver(solver.SOLVERXaddconstrs(numrow,
                                       numnz,
                                       solver_beg,
                                       solver_ind,
                                       solver_val,
                                       sense.data(),
                                       rhs.data(),
                                       nullptr))
  safe_solver(solver.SOLVERupdatemodel())
  safe_solver(solver.SOLVERsetenvTimeLimit(MAXTIMELIMT4R1CMIP))
  for (auto i : cut_type) {
    // we only need to change the coefficient of the last row
    double new_rhs = i;
    safe_solver(solver.SOLVERsetRhs(numrow - 1, 1, &new_rhs))
    safe_solver(solver.SOLVERoptimize())
    int SolCount;
    double PoolObjVal;
    safe_solver(solver.SOLVERgetSolCount(&SolCount))
    for (int j = 0; j < SolCount; ++j) {
      safe_solver(solver.SOLVERgetObjFromPool(j, &PoolObjVal))
      if (PoolObjVal - int(i / 2) < TOLERANCE) break;
      safe_solver(solver.SOLVERgetSolFromPool(j, num_first_batch_col + 1, RealDim, X))
      vector<int> cut(i);
      int tmp_cnt = 0;
      for (int k = 0; k < RealDim; ++k) {
        if (X[k] > 0.5) {
          cut[tmp_cnt++] = k + 1;
        }
      }
      cut_info.emplace_back(cut, PoolObjVal - int(i / 2));
    }
    cout << "generate R1C" << i << "= " << SolCount << endl;
  }
  safe_solver(solver.SOLVERsetenvTimeLimit(MAXTIMELIMT4MIP))
  safe_solver(solver.SOLVERfreemodel())
}

void CVRP::addR1C(BBNODE *node, const vector<int> &cut, const set<int> &mem, int cut_index) {
  size_t num_nz = 0;
  yzzLong v_comb = 0;
  for (auto i : cut) v_comb.set(i);
  int rhs = int(cut.size() / 2);
  if (rhs > 0.1) {// is not 0
    solver_ind[num_nz] = 0;
    solver_val[num_nz++] = rhs;
  }
  for (int i = 1; i < node->NumParentCols; ++i) {
    int coeff = 0;
    int times = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int current_node = ColPool4Mem[j];
      if (!current_node) break;
      if (v_comb[current_node]) {
        ++times;
        if (times == 2) {
          ++coeff;
          times = 0;
        }
      } else if (mem.find(current_node) == mem.end()) {
        times = 0;
      }
    }
    if (coeff) {
      solver_val[num_nz] = coeff;
      solver_ind[num_nz++] = i;
    }
  }
  for (int i = node->NumParentCols; i < NumCol; ++i) {
    int coeff = 0;
    int times = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int current_node = ColPool4Pricing[j];
      if (!current_node) break;
      if (v_comb[current_node]) {
        ++times;
        if (times == 2) {
          ++coeff;
          times = 0;
        }
      } else if (mem.find(current_node) == mem.end()) {
        times = 0;
      }
    }
    if (coeff) {
      solver_val[num_nz] = coeff;
      solver_ind[num_nz++] = i;
    }
  }

  if (cut_index == numeric_limits<int>::max()) {
    safe_solver(node->solver.SOLVERaddconstr((int) num_nz,
                                             solver_ind,
                                             solver_val,
                                             SOLVER_LESS_EQUAL,
                                             rhs,
                                             nullptr))
    R1C r1c3;
    r1c3.InfoR1C = cut;
    r1c3.Mem.assign(mem.begin(), mem.end());
    r1c3.IdxR1C = NumRow++;
    r1c3.RHS = rhs;
    node->R1Cs.emplace_back(r1c3);
  } else {
    fill(solver_ind2, solver_ind2 + num_nz, node->R1Cs[cut_index].IdxR1C);
    node->R1Cs[cut_index].Mem.assign(mem.begin(), mem.end());
    safe_solver(node->solver.SOLVERXchgcoeffs(num_nz, solver_ind2, solver_ind, solver_val))
  }
}

//void CVRP::findMem4R1Cs(int mode, yzzLong v_comb,
//                        set<int> &mem,
//                        vector<vector<int>> &routes,
//                        vector<double> &frac_routes) {
//  int times = 0;
//  set<int> aux;
//  if (mode == 1) {
//    //this is default mode
//    for (auto &route : routes) {
//      times = 0;
//      aux.clear();
//      for (auto v : route) {
//        if (v_comb[v]) {
//          ++times;
//          if (times == 2) {
//            for (auto v1 : aux) mem.emplace(v1);
//            times = 0;
//            aux.clear();
//          }
//        } else if (times) {
//          aux.emplace(v);
//        }
//      }
//    }
//  } else if (mode == 2) {
//    vector<pair<int, double>> idx_frac(frac_routes.size());
//    for (int i = 0; i < frac_routes.size(); ++i) {
//      idx_frac[i] = {i, frac_routes[i]};
//    }
//    std::sort(idx_frac.begin(), idx_frac.end(), [](const auto &a, const auto &b) {
//      return a.second > b.second;
//    });
//    for (auto &i : idx_frac) {
//      int idx = i.first;
//      vector<set<int>> tmp_segment_route;
//      times = 0;
//      aux.clear();
//      for (auto v : routes[idx]) {
//        if (v_comb[v]) {
//          ++times;
//          if (times == 2) {
//            tmp_segment_route.emplace_back(aux);
//            times = 1;//change to 1! Notice!
//            aux.clear();
//          }
//        } else if (times) {
//          aux.emplace(v);
//        }
//      }
//      if (!tmp_segment_route.empty()) {
//        vector<vector<int>> data;
//        combinationUtilOuter(data, (int) tmp_segment_route.size());
//        int min_size = MaxInt;
//        int min_idx;
//        for (int c = 0; c < data.size(); ++c) {
//          int num_new_added = 0;
//          for (auto &j : data[c]) {
//            for (auto v : tmp_segment_route[j]) {
//              if (mem.find(v) == mem.end()) ++num_new_added;
//            }
//          }
//          if (num_new_added < min_size) {
//            min_size = num_new_added;
//            min_idx = c;
//          }
//        }
//        for (auto &j : data[min_idx]) {
//          for (auto v : tmp_segment_route[j]) {
//            mem.emplace(v);
//          }
//        }
//      }
//    }
//  } else if (mode == 3) {
//    auto beg = std::chrono::high_resolution_clock::now();
//    findMem4R1CsMode3ByEnumeration_N_MIP(v_comb, mem, routes, );
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> diff = end - beg;
//    if (diff.count() > 0.01) {
//      cout << "findMem4R1CsMode3ByEnumeration_N_MIP takes " << diff.count() << " seconds" << endl;
//    }
//  }
//}

void CVRP::findMem4R1CsMode3ByEnumeration_N_MIP(yzzLong v_comb,
                                                std::set<int> &mem,
                                                const std::vector<std::vector<int>> &routes,
                                                bool &if_suc) {
  if_suc = true;
  int times;
  set<int> aux;
  vector<vector<vector<int>>> vec_data;
  vector<vector<set<int>>> vec_segment_route;
  yzzLong mem_long = 0;
  for (auto &i : routes) {
    vector<set<int>> tmp_segment_route;
    times = 0;
    aux.clear();
    for (auto v : i) {
      if (v_comb[v]) {
        ++times;
        if (times == 2) {
          tmp_segment_route.emplace_back(aux);
          times = 1;//change to 1! Notice!
          aux.clear();
        }
      } else if (times) {
        aux.emplace(v);
      }
    }
    if (!tmp_segment_route.empty()) {
      if (tmp_segment_route.size() % 2) {//odd number of segments //choice could only be one
        for (int j = 0; j < tmp_segment_route.size(); j += 2) {
          for (auto v : tmp_segment_route[j]) {
            mem_long.set(v);
          }
        }
      } else {//even number of segments //we add the common used vertices in the mem
        vector<int> arr(tmp_segment_route.size());
        int n = (int) arr.size();
        int r = (int) (floor((n + 1) / 2) + TOLERANCE);
        vector<vector<int>> plans;
        vector<int> tmp(r);
        iota(arr.begin(), arr.end(), 0);
        combinationUtil(arr, tmp, plans, 0, n - 1, 0, r);
        vec_data.emplace_back(plans);
        vec_segment_route.emplace_back(tmp_segment_route);
        auto bs = yzzLong{}.set();
        for (auto &p : plans) {
          yzzLong tmp_bs = 0;
          for (auto j : p) {
            for (auto v : tmp_segment_route[j]) {
              tmp_bs.set(v);
            }
          }
          bs &= tmp_bs;
        }
        if (bs.any()) {
          for (int v = 0; v < Dim; ++v) {
            if (bs[v]) {
              mem_long.set(v);
            }
          }
        }
      }
    }
  }

  int cnt = 1;
  //more reduction if it is zero then we know what to be added!
  for (int i = 0; i < vec_data.size();) {
    bool if_clear = false;
    for (auto &j : vec_data[i]) {
      bool if_all_satis = true;
      for (auto k : j) {
        for (auto l : vec_segment_route[i][k]) {
          if (!mem_long[l]) {
            if_all_satis = false;
            goto outside;
          }
        }
      }
      outside:
      if (if_all_satis) {//clear this vec
        vec_data.erase(vec_data.begin() + i);
        vec_segment_route.erase(vec_segment_route.begin() + i);
        if_clear = true;
        break;
      }
    }
    if (!if_clear) {
      cnt *= (int) vec_data[i].size();
      ++i;
    }
  }

  if (cnt != 1) {
    //we delete the repeated ones
    for (auto &r : vec_segment_route) {
      for (auto &s : r) {
        for (auto i = s.begin(); i != s.end();) {
          if (mem_long[*i]) {
            i = s.erase(i);
          } else ++i;
        }
      }
    }
  }

  if (mem_long.count() > Rank1MemSizeLimit) {
    if_suc = false;
    return;
  }

  for (int i = 1; i < Dim; ++i) {
    if (mem_long[i]) mem.emplace(i);
  }

  if (cnt == 1) {
    return;
  } else if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
    //use enumeration
    vector<int> tmp;
    set<int> new_mem;
    int record_min = MaxInt;
    combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
    mem = new_mem;
  } else {
    //use MIP
    getMemByMIP(vec_data, vec_segment_route, mem, if_suc);
  }
}

void CVRP::combinationUtil(const std::vector<int> &arr,
                           std::vector<int> &tmp,
                           std::vector<std::vector<int>> &data,
                           int start,
                           int end,
                           int index,
                           int r) {
  if (index == r) {
    data.emplace_back(tmp);
    return;
  }

  for (int i = start; end - i >= 2 * (r - index - 1); i++) {
    tmp[index] = arr[i];
    combinationUtil(arr, tmp, data, i + 2,
                    end, index + 1, r);
  }
}

void CVRP::combinationUtilOuter(vector<vector<int>> &data, int n) {
  if (recorded_combinations4R1C.find(n) != recorded_combinations4R1C.end()) {
    data = recorded_combinations4R1C[n];
    return;
  }
  vector<int> arr(n);
  iota(arr.begin(), arr.end(), 0);
  int r = int((n + 1) / 2);
  vector<int> tmp(r);
  combinationUtil(arr, tmp, data, 0, n - 1, 0, r);
  recorded_combinations4R1C[n] = data;
}

void CVRP::combinations(const vector<vector<vector<int>>> &array,
                        const vector<vector<set<int>>> &vec_segment,
                        int i,
                        const vector<int> &accum,
                        const set<int> &mem,
                        int &record_min,
                        set<int> &new_mem) {
  if (i == array.size()) {
    int num = 0;
    auto tmp_mem = mem;
    for (int j = 0; j < array.size(); ++j) {
      for (auto k : array[j][accum[j]]) {
        for (auto l : vec_segment[j][k]) {
          if (tmp_mem.find(l) == tmp_mem.end()) {
            tmp_mem.emplace(l);
            num++;
          }
        }
      }
    }
    if (num < record_min) {
      record_min = num;
      new_mem = tmp_mem;
    }
  } else {
    for (int j = 0; j < array[i].size(); ++j) {
      vector<int> tmp(accum);
      tmp.emplace_back(j);
      combinations(array, vec_segment, i + 1, tmp, mem, record_min, new_mem);
    }
  }
}

void CVRP::getMemByMIP(const vector<vector<vector<int>>> &array,
                       const vector<vector<set<int>>> &vec_segment,
                       set<int> &mem, bool &if_suc) {
//  cout << "we don't use MIP to solve the problem!" << endl;
//  return;
  unordered_map<int, int> new_mem_map;
  unordered_map<int, int> re_new_mem_map;
  new_mem_map.reserve(Dim);

  int idx = 0;
  for (auto &r : vec_segment) {
    for (auto &s : r) {
      for (auto i : s) {
        if (new_mem_map.find(i) == new_mem_map.end()) {
          new_mem_map.emplace(i, idx);
          re_new_mem_map.emplace(idx, i);
          ++idx;
        }
      }
    }
  }

  int num_r_p = 0;
  for (auto &r : array) {
    num_r_p += (int) r.size();
  }

  int half_num_row = (int) new_mem_map.size();
  int last_half = (int) array.size();
  int num_row = half_num_row + last_half;
  vector<char> sense(num_row, SOLVER_LESS_EQUAL);
  fill(sense.begin() + half_num_row, sense.end(), SOLVER_EQUAL);
  vector<double> rhs(num_row, 0);
  fill(rhs.begin() + half_num_row, rhs.end(), 1);

  SOLVER solver{};
  solver.SOLVERgetenv(&Solver);//need load environment
  cout << "we build model= getMemByMIP_.lp to get Mem!" << endl;
  //then we create a new model.
  safe_solver(solver.SOLVERnewmodel("getMemByMIP_.lp", 0, nullptr, nullptr, nullptr, nullptr, nullptr))
  safe_solver(solver.SOLVERaddconstrs(num_row, 0, nullptr, nullptr, nullptr, sense.data(), rhs.data(), nullptr))

  size_t nzcnt = 0;
  int ccnt = 0;
  int last_num_idx = half_num_row;
  for (int i = 0; i < array.size(); ++i, ++last_num_idx) {
    for (auto &p : array[i]) {
      //new variable
      solver_beg[ccnt++] = nzcnt;
      yzzLong tmp = 0;
      for (auto n : p) {
        for (auto j : vec_segment[i][n]) {
          if (tmp[j]) continue;
          tmp.set(j);
          solver_ind[nzcnt++] = new_mem_map[j];
        }
      }
      solver_ind[nzcnt++] = last_num_idx;
    }
  }
  memset(solver_obj, 0, sizeof(double) * ccnt);
  fill_n(solver_val, nzcnt, 1);
  //y variables
  fill_n(solver_obj + ccnt, half_num_row, 1);
  iota(solver_beg + ccnt, solver_beg + ccnt + half_num_row, nzcnt);
  ccnt += half_num_row;
  iota(solver_ind + nzcnt, solver_ind + nzcnt + half_num_row, 0);
  fill_n(solver_val + nzcnt, half_num_row, -num_r_p);
  nzcnt += half_num_row;
  vector<char> xtype(ccnt, SOLVER_BINARY);
  safe_solver(solver.SOLVERXaddvars(ccnt,
                                    nzcnt,
                                    solver_beg,
                                    solver_ind,
                                    solver_val,
                                    solver_obj,
                                    nullptr,
                                    nullptr,
                                    xtype.data(),
                                    nullptr))
// set time limit
  safe_solver(solver.SOLVERsetenvTimeLimit(TimeLimit4MIPFindMem))
  safe_solver(solver.SOLVERoptimize())
  safe_solver(solver.SOLVERsetenvTimeLimit(MAXTIMELIMT4MIP))
  int status, left;
  safe_solver(solver.SOLVERgetStatus(&status))
  if (status == SOLVER_TIME_LIMIT) {
    cout << "time limit for getMemByMIP" << endl;
    if_suc = false;
    goto here;
  }
  left = ccnt - num_r_p;
  safe_solver(solver.SOLVERgetX(num_r_p, left, X))
  for (int i = 0; i < left; ++i) {
    if (X[i] > 0.5) {
      mem.emplace(re_new_mem_map[i]);
    }
  }
  here:
  safe_solver(solver.SOLVERfreemodel())
}

struct other_ {
  int beg{};
  int tor{};
  vector<int> left_c;
  vector<int> mem_c;
  other_(int beg, int tor, vector<int> left_c, vector<int> mem_c) :
      beg(beg), tor(tor), left_c(std::move(left_c)), mem_c(std::move(mem_c)) {}
  other_() = default;
};

void CVRP::findR1C_multi(const vector<vector<int>> &routes,
                         const vector<double> &frac_routes,
                         const vector<vector<int>> &vbasis,
                         vector<pair<vector<int>, double>> &cut_info,
                         vector<tuple<vector<int>, int, double>> &multi_cut_info) {
//  throw runtime_error(
//      "could possibly have error about ~new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));\n"
//      "          new_cij.emplace_back(choice_swap_i_j.second);~");
  vector<vector<int>> v_r_map(Dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      v_r_map[i].emplace_back(r);
    }
  }

  int route_size = (int) routes.size();
  vector<vector<vector<int>>> c(Dim, vector<vector<int>>(route_size));
  vector<vector<vector<int>>> w_no_c(Dim, vector<vector<int>>(route_size));
  auto if_vis = new bool[route_size];
  unordered_map<int, vector<pair<int, int>>> c_map;//size, c_index
  for (int i = 1; i < Dim; ++i) {
    memset(if_vis, 0, sizeof(bool) * route_size);
    yzzLong wc = 0;// c within i
    for (auto &j : v_r_map[i]) {
      if_vis[j] = true;
      for (auto v : routes[j]) {
        if (Rank1SepHeurMem4Vertex[i][v]) wc.set(v);
      }
    }
    for (int j = 0; j < routes.size(); ++j) {
      if (if_vis[j]) continue;
      yzzLong tmp_c = 0;
      for (auto v : routes[j]) {
        if (wc[v]) {
          tmp_c.set(v);//r no i but c within i
        }
      }
      tmp_c.set(i);
      int c_size = (int) tmp_c.count();// basic c set
      if (c_size < 3 || c_size > CONFIG::MaxHeurInitialCsetSize4RowRank1C) continue;
      c_map[c_size].emplace_back(i, j);
      c[i][j].resize(c_size);
      w_no_c[i][j].reserve(Dim);
      c_size = 0;
      for (int k = 1; k < Dim; ++k) {
        if (tmp_c[k]) c[i][j][c_size++] = k;
        else if (wc[k]) w_no_c[i][j].emplace_back(k);//candidates set
      }
    }
  }
  delete[] if_vis;

  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> generated_rank1_multi_pool;
  vector<rank1_multi_label> rank1_multi_label_pool(Initial_rank1_multi_label_pool_size);
  unordered_map<int, int> num_operations;
  int label_cnt = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  double total_time = 0;
  for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
    for (int len = 4; len <= CONFIG::MaxRowRank1; ++len) {
      for (auto &p : c_map[len]) {
        int i = p.first, j = p.second;// i and non-visited routes idx: j
        vector<int> new_cij;
        double vio;
        bool if_succeed;
        //calculate the best vio while using the multiplier plan and the set c[i][j] and the corresponding multiplier
        auto beg10 = std::chrono::high_resolution_clock::now();
        findCombs4rankCutSet2getBestMultiplier(plan_idx,
                                               c[i][j],
                                               new_cij,
                                               i,
                                               vio,
                                               v_r_map,
                                               frac_routes,
                                               if_succeed);
        auto end10 = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double>(end10 - beg10).count();
        if (!if_succeed) continue;//if negative, then no such multiplier
        //revise the add... function to get the start position as the parameter(the above set is sorted again!)
        vector<pair<int, double>> move_vio(4);
        move_vio[0] = {0, vio};
        double new_vio = vio;
        int choice_add_j;
        int choice_remove_j;
        pair<int, int> choice_swap_i_j;
        addInRank1HeurSep(plan_idx, new_cij, w_no_c[i][j], new_vio, choice_add_j, v_r_map, frac_routes);
        move_vio[1] = {1, new_vio};
        new_vio = vio;
        removeInRank1HeurSep(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
        move_vio[2] = {2, new_vio};
        new_vio = vio;
        swapInRank1HeurSep(plan_idx, new_cij, w_no_c[i][j], new_vio, choice_swap_i_j, v_r_map, frac_routes);
        move_vio[3] = {3, new_vio};
        stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                         const pair<int, double> &b) {
          return a.second > b.second;
        });
        int best_move = move_vio[0].first;
        ++num_operations[best_move];
        if (best_move == 0) {
          if (move_vio[0].second > TOLERANCE)
            generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij, plan_idx, move_vio[0].second);
        } else if (best_move == 1) {
          new_cij.emplace_back(choice_add_j);
          auto new_w_no_c = w_no_c[i][j];
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_add_j));
          rank1_multi_label_pool[label_cnt++] = rank1_multi_label{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  'a'};
          if (label_cnt == rank1_multi_label_pool.size()) {
            rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
          }
        } else if (best_move == 2) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
          rank1_multi_label_pool[label_cnt++] = rank1_multi_label{new_cij, w_no_c[i][j], plan_idx, move_vio[0].second,
                                                                  'r'};
          if (label_cnt == rank1_multi_label_pool.size()) {
            rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
          }
        } else if (best_move == 3) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
          new_cij.emplace_back(choice_swap_i_j.second);
          auto new_w_no_c = w_no_c[i][j];
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_swap_i_j.second));
          rank1_multi_label_pool[label_cnt++] = rank1_multi_label{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  's'};
          if (label_cnt == rank1_multi_label_pool.size()) {
            rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
          }
        }
      }
    }
  }
  cout << "total time: " << total_time << endl;
  auto end = std::chrono::high_resolution_clock::now();
  cout << "time for generating rank1 multi label pool: " << std::chrono::duration<double>(end - beg).count() << endl;
  //extend the rank1_multi_label_pool
  int num_extend = 0;
  beg = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < label_cnt;) {
    auto &label = rank1_multi_label_pool[i];
    auto &vio = label.vio;
    auto &new_cij = label.c;
    auto &w_no_cij = label.w_no_c;
    auto &plan_idx = label.plan_idx;
    if (label.search_dir == 'a') {
      vector<pair<int, double>> move_vio(3);
      move_vio[0] = {0, vio};
      double new_vio = vio;
      int choice_add_j;
      pair<int, int> choice_swap_i_j;
      addInRank1HeurSep(plan_idx, new_cij, w_no_cij, new_vio, choice_add_j, v_r_map, frac_routes);
      move_vio[1] = {1, new_vio};
      new_vio = vio;
      swapInRank1HeurSep(plan_idx, new_cij, w_no_cij, new_vio, choice_swap_i_j, v_r_map, frac_routes);
      move_vio[2] = {3, new_vio};
      stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                       const pair<int, double> &b) {
        return a.second > b.second;
      });

      int best_move = move_vio[0].first;
      if (best_move != 0) {
        num_extend++;
      }
      if (best_move == 0) {
        if (move_vio[0].second > TOLERANCE) {
          generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij,
                                                                        plan_idx,
                                                                        move_vio[0].second);
        }
        ++i;
      } else if (best_move == 1) {
        new_cij.emplace_back(choice_add_j);
        w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_add_j));
        vio = move_vio[0].second;
      } else if (best_move == 3) {
        new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
        new_cij.emplace_back(choice_swap_i_j.second);
        w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_swap_i_j.second));
        vio = move_vio[0].second;
      }
    } else if (label.search_dir == 'r') {
      vector<pair<int, double>> move_vio(3);
      move_vio[0] = {0, vio};
      double new_vio = vio;
      int choice_remove_j;
      pair<int, int> choice_swap_i_j;
      removeInRank1HeurSep(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
      move_vio[1] = {2, new_vio};
      new_vio = vio;
      swapInRank1HeurSep(plan_idx, new_cij, w_no_cij, new_vio, choice_swap_i_j, v_r_map, frac_routes);
      move_vio[2] = {3, new_vio};
      stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                       const pair<int, double> &b) {
        return a.second > b.second;
      });

      int best_move = move_vio[0].first;
      if (best_move != 0) {
        num_extend++;
      }
      if (best_move == 0) {
        if (move_vio[0].second > TOLERANCE) {
          generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij,
                                                                        plan_idx,
                                                                        move_vio[0].second);
        }
        ++i;
      } else if (best_move == 2) {
        new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
        vio = move_vio[0].second;
      } else if (best_move == 3) {
        new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
        new_cij.emplace_back(choice_swap_i_j.second);
        w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_swap_i_j.second));
        vio = move_vio[0].second;
      }
    } else if (label.search_dir == 's') {
      vector<pair<int, double>> move_vio(3);
      move_vio[0] = {0, vio};
      double new_vio = vio;
      int choice_add_j;
      int choice_remove_j;
      addInRank1HeurSep(plan_idx, new_cij, w_no_cij, new_vio, choice_add_j, v_r_map, frac_routes);
      move_vio[1] = {1, new_vio};
      new_vio = vio;
      removeInRank1HeurSep(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
      move_vio[2] = {2, new_vio};
      stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                       const pair<int, double> &b) {
        return a.second > b.second;
      });

      int best_move = move_vio[0].first;
      if (best_move != 0) {
        num_extend++;
      }
      if (best_move == 0) {
        if (move_vio[0].second > TOLERANCE) {
          generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij,
                                                                        plan_idx,
                                                                        move_vio[0].second);
        }
        ++i;
      } else if (best_move == 1) {
        new_cij.emplace_back(choice_add_j);
        w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_add_j));
        vio = move_vio[0].second;
      } else if (best_move == 2) {
        new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
        vio = move_vio[0].second;
      }
    }
  }
  end = std::chrono::high_resolution_clock::now();
  cout << "num_extend: " << num_extend << endl;
  for (auto &i : num_operations) {
    switch (i.first) {
      case 0:cout << "self ";
        break;
      case 1:cout << "add ";
        break;
      case 2:cout << "remove ";
        break;
      case 3:cout << "swap ";
        break;
    }
    cout << i.second << endl;
  }
  cout << "rank1 multi pool generation time: " << std::chrono::duration<double>(end - beg).count() << endl;
  //sort the generated pool
  for (auto i = generated_rank1_multi_pool.begin(); i != generated_rank1_multi_pool.end();) {
    if (i->first > CONFIG::MaxRowRank1) {
      i = generated_rank1_multi_pool.erase(i);
    } else {
      ++i;
    }
  }

  for (auto &i : generated_rank1_multi_pool) {
    stable_sort(i.second.begin(), i.second.end(),
                [](auto &a, auto &b) {
                  return get<2>(a) > get<2>(b);
                });
    double vio_std = get<2>(i.second[0]) * CONFIG::CutVioFactor;
    //if cut's vio is less than vio_std, then cut is deleted
    i.second.erase(remove_if(i.second.begin(), i.second.end(),
                             [vio_std](auto &a) {
                               return get<2>(a) < vio_std;
                             }), i.second.end());
  }
  //get v_r_map related to vbasis
  v_r_map.clear();
  v_r_map.resize(Dim);
  for (int r = 0; r < vbasis.size(); ++r) {
    for (auto &i : vbasis[r]) {
      v_r_map[i].emplace_back(r);
    }
  }
  map<int, vector<vector<int>>> test_group;
  for (auto &i : generated_rank1_multi_pool) {
    auto &test = test_group[i.first];
    double b4_vio = get<2>(i.second[0]);
    vector<pair<int, int>> tmp;//pair of <idx ,plan_idx>
    for (int j = 1; j < i.second.size(); ++j) {
      double now_vio = get<2>(i.second[j]);
      int plan_idx = get<1>(i.second[j - 1]);//old_plan
      if (abs(now_vio - b4_vio) > TOLERANCE) {
        b4_vio = now_vio;
        if (!tmp.empty()) {
          tmp.emplace_back(j - 1, plan_idx);
          sort(tmp.begin(), tmp.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
            return a.second < b.second;
          }//no need stable
          );
          vector<int> tmp_idx(tmp.size());
          transform(tmp.begin(), tmp.end(), tmp_idx.begin(), [](const pair<int, int> &a) {
            return a.first;
          });
          test.emplace_back(std::move(tmp_idx));
          tmp.clear();
        }
      } else {
        tmp.emplace_back(j - 1, plan_idx);
      }
    }
  }
//  for (auto &i : generated_rank1_multi_pool) {
//    auto &test = test_group[i.first];
//    double b4_vio = get<2>(i.second[0]);
//    vector<int> tmp;//pair of <idx ,plan_idx>
//    for (int j = 1; j < i.second.size(); ++j) {
//      double now_vio = get<2>(i.second[j]);
//      int plan_idx = get<1>(i.second[j]);
//      if (abs(now_vio - b4_vio) > TOLERANCE) {
//        b4_vio = now_vio;
//        if (!tmp.empty()) {
//          tmp.emplace_back(j - 1);
//          vector<int> tmp_idx = tmp;
//          test.emplace_back(std::move(tmp_idx));
//          tmp.clear();
//        }
//      } else {
//        tmp.emplace_back(j - 1);
//      }
//    }
//  }
//we check if the two plans about the basis routes are the same
  vector<int> num_vis(vbasis.size());
  for (auto &i : test_group) {
    auto if_del = new bool[generated_rank1_multi_pool[i.first].size()]();
    for (auto &j : i.second) {
      vector<vector<int>> tmp;
      for (auto &k : j) {
        auto &cut = generated_rank1_multi_pool[i.first][k];
        memset(num_vis.data(), 0, sizeof(int) * num_vis.size());
        const auto &plan = get<0>(map_rank1_multiplier[i.first][get<1>(cut)]);
        auto deno = get<1>(map_rank1_multiplier[i.first][get<1>(cut)]);
        int cnt = 0;
        for (auto &l : get<0>(cut)) {
          for (auto &m : v_r_map[l]) {
            num_vis[m] += plan[cnt];
          }
          ++cnt;
        }
        transform(num_vis.begin(), num_vis.end(), num_vis.begin(),
                  [deno](int a) { return a / deno; });
        if (std::find(tmp.begin(), tmp.end(), num_vis) == tmp.end()) {
          tmp.emplace_back(num_vis);
        } else if_del[k] = true;
      }
    }
    int len = 0, keep = 0;
    for (int j = 0; j < generated_rank1_multi_pool[i.first].size(); ++j) {
      if (!if_del[j]) {
//      if (true) {
        auto &plan_idx = get<1>(generated_rank1_multi_pool[i.first][j]);
        if (!plan_idx) {
          //construct the basic plan
          cut_info.emplace_back(get<0>(generated_rank1_multi_pool[i.first][j]),
                                get<2>(generated_rank1_multi_pool[i.first][j]));
        } else {
          multi_cut_info.emplace_back(get<0>(generated_rank1_multi_pool[i.first][j]),
                                      plan_idx,
                                      get<2>(generated_rank1_multi_pool[i.first][j]));
        }
        ++len;
        if (len == CONFIG::MaxNumR1CPerRound) break;
      }
    }
    delete[] if_del;
  }

#ifdef DETAILED_RANK1_PRINT_INFO
  //now we find the memory and violation
  for (auto &i : generated_rank1_multi_pool) {
    cout << "----------\nrow= " << i.first << endl;
    vector<int> plans(7, 0);
    for (auto &j : i.second) {
      plans[get<1>(j)]++;
    }
    double sum_vio =
        accumulate(i.second.begin(), i.second.end(), 0.0, [](double a, const tuple<vector<int>, int, double> &b) {
          return a + get<2>(b);
        });
    sum_vio /= (int) i.second.size();
    cout << "aver_vio= " << sum_vio << endl;
    double max_vio = get<2>(*max_element(i.second.begin(), i.second.end(), [](const tuple<vector<int>, int, double> &a,
                                                                              const tuple<vector<int>,
                                                                                          int,
                                                                                          double> &b) {
      return get<2>(a) < get<2>(b);
    }));
    cout << "max_vio= " << max_vio << endl;
    cout << "plans:" << endl;
    int cnt = 0;
    for (auto &j : plans) {
      cout << "idx= " << cnt++ << " num= " << j << endl;
    }
  }
#endif
}

//void CVRP::generateAllR1Cs(BBNODE *node, bool if_use_MIP) {
//  //old cuts
//  unordered_map<yzzLong, std::pair<std::set<int>, int>> R1CsPool;//mem & cut_index
//  for (int i = 0; i < node->R1Cs.size(); ++i) {
//    yzzLong tmp_long = 0;
//    for (auto j : node->R1Cs[i].InfoR1C) tmp_long.set(j);
//    R1CsPool[tmp_long] = {set<int>(node->R1Cs[i].Mem.begin(), node->R1Cs[i].Mem.end()), i};
//  }
//  unordered_map<pair<vector<int>, int>, std::pair<std::set<int>, int>, Rank1_multi_pair_Hasher>
//      R1C_multi_Pool;//mem & cut_index
//  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
//    R1C_multi_Pool[node->R1Cs_multi[i].InfoR1C] =
//        {set<int>(node->R1Cs_multi[i].Mem.begin(), node->R1Cs_multi[i].Mem.end()), i};
//  }
//
//  vector<pair<vector<int>, double>> cut_info;
//  vector<tuple<vector<int>, int, double>> multi_cut_info;
//  vector<vector<int>> routes;
//  vector<double> frac_routes;
//  vector<int> route;
//  route.reserve(Dim);
//  for (int i = 0; i < node->NumParentColsInLPSols; ++i) {
//    route.clear();
//    for (auto j = node->Idx4LPSolsInColPool[i].first + 1;; ++j) {
//      int curr_node = ColPool4Mem[j];
//      if (!curr_node)break;
//      route.emplace_back(curr_node);
//    }
//    routes.emplace_back(route);
//    frac_routes.emplace_back(node->Idx4LPSolsInColPool[i].second);
//  }
//  for (int i = node->NumParentColsInLPSols; i < node->Idx4LPSolsInColPool.size(); ++i) {
//    route.clear();
//    for (auto j = node->Idx4LPSolsInColPool[i].first + 1;; ++j) {
//      int curr_node = ColPool4Pricing[j];
//      if (!curr_node)break;
//      route.emplace_back(curr_node);
//    }
//    routes.emplace_back(route);
//    frac_routes.emplace_back(node->Idx4LPSolsInColPool[i].second);
//  }
//
//  vector<vector<int>> ele_routes;//contain no depot
//  vector<vector<int>> non_ele_routes;
//  vector<double> frac_ele_routes;
//  vector<double> frac_non_ele_routes;
//
//  for (int r = 0; r < routes.size(); ++r) {
//    yzzLong tmp_long = 0;
//    bool if_ele = true;
//    for (auto i : routes[r]) {
//      if (tmp_long[i]) {
//        if_ele = false;
//        goto outside;
//      }
//      tmp_long.set(i);
//    }
//    outside:
//    if (if_ele) {
//      ele_routes.emplace_back(routes[r]);
//      frac_ele_routes.emplace_back(frac_routes[r]);
//    } else {
//      non_ele_routes.emplace_back(routes[r]);
//      frac_non_ele_routes.emplace_back(frac_routes[r]);
//    }
//  }
//
//  auto beg2 = chrono::high_resolution_clock::now();
//
//  generateR1C1s(non_ele_routes, frac_non_ele_routes, cut_info);
//
//  auto end2 = chrono::high_resolution_clock::now();
//  auto duration2 = chrono::duration<double>(end2 - beg2).count();
//  cout << "R1C1s time: " << duration2 << endl;
//
//  generateR1C3s(ele_routes, non_ele_routes, frac_ele_routes, frac_non_ele_routes, cut_info);
//
//  auto end3 = chrono::high_resolution_clock::now();
//  auto duration3 = chrono::duration<double>(end3 - end2).count();
//  cout << "R1C3s time: " << duration3 << endl;
//
//  if (if_use_MIP) {
//    vector<int> in_cut_type;
//    for (int i = 5; i <= CONFIG::MaxRowRank1; i += 2) {
//      in_cut_type.emplace_back(i);
//    }
//    generateHighDimR1Cs(routes, frac_routes, in_cut_type, cut_info);
//  } else {
//    vector<vector<int>> VBasis;
//    for (int i = 0; i < NumCol; ++i) {
//      int v_idx = node->VBasis[i];
//      if (v_idx) continue;
//      int *col_pool;
//      if (v_idx < node->NumParentCols) {
//        col_pool = ColPool4Mem;
//      } else col_pool = ColPool4Pricing;
//      route.clear();
//      for (auto j = node->IdxCols[v_idx] + 1;; ++j) {
//        int curr_node = col_pool[j];
//        if (!curr_node)break;
//        route.emplace_back(curr_node);
//      }
//      VBasis.emplace_back(route);
//    }
//    auto beg = std::chrono::high_resolution_clock::now();
//
//    findR1C_multi(routes, frac_routes, VBasis, cut_info, multi_cut_info);
//
//    auto end = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration<double>(end - beg).count();
//    cout << "findR1C_multi time: " << duration << endl;
//  }
//
//  cout << "cut_info size: " << cut_info.size() << endl;
//  cout << "multi_cut_info size: " << multi_cut_info.size() << endl;
//
//  auto beg = std::chrono::high_resolution_clock::now();
//
//  vector<tuple<int, set<int>, double, int, int>> cut_info_set;//cnt & mem & vio & index & used_mem
//  int cnt = 0;
//  for (auto &cut : cut_info) {
//    yzzLong v_comb = 0;
//    for (auto i : cut.first) {
//      v_comb.set(i);
//    }
//    set<int> mem = {};
//    int old_mem_size = 0;
//    int index = numeric_limits<int>::max();
//    if (R1CsPool.find(v_comb) != R1CsPool.end()) {
//      mem = R1CsPool[v_comb].first;
//      index = R1CsPool[v_comb].second;
//      old_mem_size = (int) mem.size();
//    }
//    bool if_suc;
//    //get memory for the new cuts
//    findMem4R1CsMode3ByEnumeration_N_MIP(v_comb, mem, routes, if_suc);
//    //add cuts into lp
//    if (if_suc) {
//      cut_info_set.emplace_back(cnt, mem, cut.second, index, mem.size() - old_mem_size + cut.first.size());
//    }
//    ++cnt;
//  }
//
//  addSelectedR1C_N_multiCuts<vector<pair<vector<int>, double>>, false>(node, cut_info_set, cut_info);
//
//  cout << "cut_info_set size: " << cut_info_set.size() << endl;
//  cut_info_set.clear();
//  cnt = 0;
//  for (auto &cut_pair : multi_cut_info) {
//    auto cut = make_pair(get<0>(cut_pair), get<1>(cut_pair));
//    set<int> mem = {};
//    int old_mem_size = 0;
//    int index = numeric_limits<int>::max();
//    if (R1C_multi_Pool.find(cut) != R1C_multi_Pool.end()) {
//      mem = R1C_multi_Pool[cut].first;
//      index = R1C_multi_Pool[cut].second;
//      old_mem_size = (int) mem.size();
//    }
//    bool if_suc;
//    //got memory for the new cuts
//    findMem4Rank1_multi(routes, cut, mem, if_suc);
//    //add cuts into lp
//    if (if_suc) {
//      cut_info_set.emplace_back(cnt, mem, get<2>(cut_pair), index, mem.size() - old_mem_size + cut.first.size());
//    }
//    ++cnt;
//  }
//  cout << "cut_info_set size: " << cut_info_set.size() << endl;
//
//  addSelectedR1C_N_multiCuts<vector<tuple<vector<int>, int, double>>, false>(node, cut_info_set, multi_cut_info);
//
//  auto end = std::chrono::high_resolution_clock::now();
//  auto duration = std::chrono::duration<double>(end - beg).count();
//  cout << "add cuts time: " << duration << endl;
//
//  safe_Hyperparameter(checkMaxNum_R1Cs((int) node->R1Cs.size()))
//  safe_Hyperparameter((int) node->R1Cs_multi.size() > MaxNum_R1C_multi)
//  safe_Hyperparameter(NumRow > CST_LIMIT)
//}

void CVRP::addR1C_multi(BBNODE *node, const pair<vector<int>, int> &cut, const set<int> &mem, int cut_index) {
  size_t num_nz = 0;
  yzzLong v_comb = 0;
  int size = (int) cut.first.size();
  const auto &plan = get<0>(map_rank1_multiplier[size][cut.second]);
  auto denominator = get<1>(map_rank1_multiplier[size][cut.second]);
  auto rhs = get<2>(map_rank1_multiplier[size][cut.second]);
  if (rhs > 0.1) {// is not 0
    solver_ind[num_nz] = 0;
    solver_val[num_nz++] = rhs;
  }

  unordered_map<int, int> v_2_mul;
  for (int i = 0; i < size; ++i) {
    v_2_mul[cut.first[i]] = plan[i];
  }

  for (int i = 1; i < node->NumParentCols; ++i) {
    int coeff = 0;
    int times = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int current_node = ColPool4Mem[j];
      if (!current_node) break;
      if (v_2_mul.find(current_node) != v_2_mul.end()) {
        times += v_2_mul[current_node];
        if (times >= denominator) {
          ++coeff;
          times -= denominator;
        }
      } else if (mem.find(current_node) == mem.end()) {
        times = 0;
      }
    }
    if (coeff) {
      solver_val[num_nz] = coeff;
      solver_ind[num_nz++] = i;
    }
  }
  for (int i = node->NumParentCols; i < NumCol; ++i) {
    int coeff = 0;
    int times = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int current_node = ColPool4Pricing[j];
      if (!current_node) break;
      if (v_2_mul.find(current_node) != v_2_mul.end()) {
        times += v_2_mul[current_node];
        if (times >= denominator) {
          ++coeff;
          times -= denominator;
        }
      } else if (mem.find(current_node) == mem.end()) {
        times = 0;
      }
    }
    if (coeff) {
      solver_val[num_nz] = coeff;
      solver_ind[num_nz++] = i;
    }
  }

  if (cut_index == numeric_limits<int>::max()) {
    safe_solver(node->solver.SOLVERaddconstr((int) num_nz,
                                             solver_ind,
                                             solver_val,
                                             SOLVER_LESS_EQUAL,
                                             rhs,
                                             nullptr))
    R1C_multi r1c;
    r1c.InfoR1C = cut;
    r1c.Mem.assign(mem.begin(), mem.end());
    r1c.IdxR1C = NumRow++;
    r1c.RHS = rhs;
    node->R1Cs_multi.emplace_back(r1c);
  } else {
    fill(solver_ind2, solver_ind2 + num_nz, node->R1Cs_multi[cut_index].IdxR1C);
    node->R1Cs_multi[cut_index].Mem.assign(mem.begin(), mem.end());
    safe_solver(node->solver.SOLVERXchgcoeffs(num_nz, solver_ind2, solver_ind, solver_val))
  }
}

void CVRP::findMem4Rank1_multi(const vector<vector<int>> &routes,
                               const std::vector<std::unordered_map<int, int>> &v_r_map,
                               const pair<vector<int>, int> &cut_pair,
                               set<int> &mem,
                               bool &if_suc) {
  //find vec_data and vec_segment_route
  if_suc = true;
  auto &cut = cut_pair.first;
  auto plan_idx = cut_pair.second;
  int size = (int) cut.size();
  const auto &multi = get<0>(map_rank1_multiplier[size][plan_idx]);
  auto denominator = get<1>(map_rank1_multiplier[size][plan_idx]);
  vector<vector<vector<int>>> vec_data;
  vector<vector<set<int>>> vec_segment_route;
  unordered_map<int, int> map_cut_mul;
  vector<int> num_vis_times(routes.size(), 0);
  for (int i = 0; i < cut.size(); ++i) {
    map_cut_mul[cut[i]] = multi[i];
    for (auto &pr : v_r_map[cut[i]]) {
      num_vis_times[pr.first] += multi[i] * pr.second;
    }
  }
  transform(num_vis_times.begin(), num_vis_times.end(), num_vis_times.begin(),
            [denominator](int x) { return int(x / denominator); });
  yzzLong mem_long = 0;
  int num = 0;
  for (auto &i : routes) {
    if (num_vis_times[num++] == 0) continue;
    vector<vector<int>> data;
    vector<int> vis;
    vector<set<int>> segment_route;
    set<int> tmp_seg;
    for (auto &j : i) {
      if (map_cut_mul.find(j) != map_cut_mul.end()) {
        vis.emplace_back(map_cut_mul[j]);
        segment_route.emplace_back(tmp_seg);
        tmp_seg.clear();
      } else {
        tmp_seg.insert(j);
      }
    }
    if (!segment_route.empty())
      segment_route.erase(segment_route.begin());//remove the first one
    findPlan4Rank1_multi(vis, denominator, mem_long, segment_route, data);
    if (!data.empty()) {
      vec_data.emplace_back(data);
      vec_segment_route.emplace_back(segment_route);
    }
  }

  //deep clean
  size_t cnt = 1;
  for (int i = 0; i < vec_data.size();) {
    bool if_clear = false;
    for (auto &j : vec_data[i]) {
      bool if_all_satis = true;
      for (auto k : j) {
        for (auto l : vec_segment_route[i][k]) {
          if (!mem_long[l]) {
            if_all_satis = false;
            goto outside;
          }
        }
      }
      outside:
      if (if_all_satis) {//clear this vec
        vec_data.erase(vec_data.begin() + i);
        vec_segment_route.erase(vec_segment_route.begin() + i);
        if_clear = true;
        break;
      }
    }
    if (!if_clear) {
      cnt *= (int) vec_data[i].size();
      ++i;
    }
  }

  if (cnt != 1) {
    //we delete the repeated ones
    for (auto &r : vec_segment_route) {
      for (auto &s : r) {
        for (auto i = s.begin(); i != s.end();) {
          if (mem_long[*i]) {
            i = s.erase(i);
          } else ++i;
        }
      }
    }
  }

  if (mem_long.count() > Rank1MemSizeLimit) {
    if_suc = false;
    return;
  }

  for (int i = 1; i < Dim; ++i) {
    if (mem_long[i]) {
      mem.emplace(i);
    }
  }

  cnt = 1;
  for (auto &i : vec_data) {// in this way we don't have to worry about overflow!
    cnt *= i.size();
    if (cnt >= FIND_MEM_USE_ENUMERATION_OR_MIP) break;
  }

  if (cnt == 1) {
    return;
  } else if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
    //use enumeration
    vector<int> tmp;
    set<int> new_mem;
    int record_min = MaxInt;
    combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
    mem = new_mem;
  } else {
    //TODO: perfect the use! right now we just never use it!
//    throw runtime_error("we don't use MIP to find mem now!");
    //use MIP
    getMemByMIP(vec_data, vec_segment_route, mem, if_suc);
  }
}

void CVRP::construct_mem(BBNODE *node, const vector<vector<int>> &routes,
                         const std::vector<std::unordered_map<int, int>> &v_r_map,
                         const vector<pair<vector<int>, int>> &cuts,
                         vector<tuple<vector<int>, int, int, set<int>>> &full_cuts) {//cut self, plan idx, cut idx, mem
  unordered_map<pair<vector<int>, int>, std::pair<std::set<int>, int>, Rank1_multi_pair_Hasher>
      R1C_multi_Pool;//mem & cut_index
  for (int i = 0; i < node->R1Cs.size(); ++i) {
    R1C_multi_Pool[{node->R1Cs[i].InfoR1C, 0}] =
        {set < int > (node->R1Cs[i].Mem.begin(), node->R1Cs[i].Mem.end()), i};
  }
  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
    R1C_multi_Pool[node->R1Cs_multi[i].InfoR1C] =
        {set < int > (node->R1Cs_multi[i].Mem.begin(), node->R1Cs_multi[i].Mem.end()), i};
  }
  full_cuts.reserve(cuts.size());
  int num_r1c = MaxNum_R1Cs - node->R1Cs.size() - 1;
  int num_mul_r1c = MaxNum_R1C_multi - node->R1Cs_multi.size() - 1;
  for (auto &cut : cuts) {
    set<int> mem = {};
    int index = numeric_limits<int>::max();
    auto it_find = R1C_multi_Pool.find(cut);
    bool if_suc;
    if (it_find != R1C_multi_Pool.end()) {
      mem = it_find->second.first;
      index = it_find->second.second;
    }
    if (!cut.second) {
      if (--num_r1c < 0) continue;
      yzzLong v_comb = 0;
      for (auto i : cut.first) v_comb.set(i);
      findMem4R1CsMode3ByEnumeration_N_MIP(v_comb, mem, routes, if_suc);
    } else {
      if (--num_mul_r1c < 0) continue;
      findMem4Rank1_multi(routes, v_r_map, cut, mem, if_suc);
    }
    if (if_suc) {
      full_cuts.emplace_back(cut.first, cut.second, index, mem);
    }
  }
}

void CVRP::NewgenerateAllR1Cs(BBNODE *node, int aggressive_level) {
//  for (int i=0;i<NumCol;++i){
//	safe_solver(node->solver.SOLVERsetColUpper(i, 1))
//  }
  //old cuts
  cut_record.clear();
#ifdef debugVio
  vio_map.clear();
#endif
  vector<pair<vector<int>, int>> cuts;
  vector<vector<int>> routes;
  vector<double> frac_routes;
  vector<int> route;
  route.reserve(Dim);
  //contain no depot
  vector<vector<int>> ele_routes;
  vector<vector<int>> non_ele_routes;
  vector<double> frac_ele_routes;
  vector<double> frac_non_ele_routes;

  if (If_in_Enu_State) {
    for (auto &i : node->Idx4LPSolsInColPool) {
      route.clear();
      for (auto j = i.first + 1;; ++j) {
        int curr_node = ColPool4Pricing[j];
        if (!curr_node)break;
        route.emplace_back(curr_node);
      }
      routes.emplace_back(route);
      frac_routes.emplace_back(i.second);
    }
    //sort
    std::vector<size_t> indices(frac_routes.size());
    iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&frac_routes](size_t i1, size_t i2) {
      return frac_routes[i1] > frac_routes[i2];
    });
    std::vector<std::vector<int>> temp_routes = routes;
    std::vector<double> temp_frac_routes = frac_routes;
    transform(indices.begin(), indices.end(), routes.begin(), [&temp_routes](size_t i) {
      return temp_routes[i];
    });
    transform(indices.begin(), indices.end(), frac_routes.begin(), [&temp_frac_routes](size_t i) {
      return temp_frac_routes[i];
    });
  } else {
    for (int i = 0; i < node->NumParentColsInLPSols; ++i) {
      route.clear();
      if (abs(node->Idx4LPSolsInColPool[i].second - 1) < TOLERANCE)continue;
      for (auto j = node->Idx4LPSolsInColPool[i].first + 1;; ++j) {
        int curr_node = ColPool4Mem[j];
        if (!curr_node)break;
        route.emplace_back(curr_node);
      }
      routes.emplace_back(route);
      frac_routes.emplace_back(node->Idx4LPSolsInColPool[i].second);
    }
    for (int i = node->NumParentColsInLPSols; i < node->Idx4LPSolsInColPool.size(); ++i) {
      route.clear();
      if (abs(node->Idx4LPSolsInColPool[i].second - 1) < TOLERANCE)continue;
      for (auto j = node->Idx4LPSolsInColPool[i].first + 1;; ++j) {
        int curr_node = ColPool4Pricing[j];
        if (!curr_node)break;
        route.emplace_back(curr_node);
      }
      routes.emplace_back(route);
      frac_routes.emplace_back(node->Idx4LPSolsInColPool[i].second);
    }
    //sort
    std::vector<size_t> indices(frac_routes.size());
    iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&frac_routes](size_t i1, size_t i2) {
      return frac_routes[i1] > frac_routes[i2];
    });
    std::vector<std::vector<int>> temp_routes = routes;
    std::vector<double> temp_frac_routes = frac_routes;
    transform(indices.begin(), indices.end(), routes.begin(), [&temp_routes](size_t i) {
      return temp_routes[i];
    });
    transform(indices.begin(), indices.end(), frac_routes.begin(), [&temp_frac_routes](size_t i) {
      return temp_frac_routes[i];
    });

    for (int r = 0; r < routes.size(); ++r) {
      yzzLong tmp_long = 0;
      bool if_ele = true;
      for (auto i : routes[r]) {
        if (tmp_long[i]) {
          if_ele = false;
          goto outside;
        }
        tmp_long.set(i);
      }
      outside:
      if (if_ele) {
        ele_routes.emplace_back(routes[r]);
        frac_ele_routes.emplace_back(frac_routes[r]);
      } else {
        non_ele_routes.emplace_back(routes[r]);
        frac_non_ele_routes.emplace_back(frac_routes[r]);
      }
    }
  }

  if (!If_in_Enu_State) {
#ifndef debugVio
    fill_mem_first(node, routes, frac_routes, cuts);
#else
    cout << "fill_mem_first has been closed!" << endl;
#endif
  }

  NewgenerateR1C1s(non_ele_routes, frac_non_ele_routes, cuts);

  NewgenerateR1C3s(routes, frac_routes, cuts);// unnecessarily enumeration!//can be supplied by findR1C_multi_aggressive

  auto beg = std::chrono::high_resolution_clock::now();
  std::vector<std::unordered_map<int, int>> v_r_map;

  switch (aggressive_level) {
    case 0:newfindR1C_multi(routes, frac_routes, v_r_map, cuts);
      break;
    case 1:findR1C_multi_aggressive(routes, frac_routes, v_r_map, cuts);
      break;
    case 2:crazySearch(routes, frac_routes, v_r_map, cuts);
      break;
    default: throw std::runtime_error("aggressive_level error");
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double>(end - beg).count();

  cout << "findR1C_multi time: " << duration << endl;

  cout << "cuts.size() = " << cuts.size() << endl;

  if (If_in_Enu_State) {
    addR1C_atOnceInEnu(node, cuts);
  } else {
    vector<tuple<vector<int>, int, int, set<int>>> full_cuts;// cut, plan, index, mem
    construct_mem(node, routes, v_r_map, cuts, full_cuts);
    addR1C_atOnce(node, full_cuts);
  }

  cout << "now num_row=  " << NumRow << endl;
  safe_Hyperparameter(NumRow > CST_LIMIT)
}

void CVRP::findPlan4Rank1_multi(const vector<int> &vis, int denominator, yzzLong &mem,
                                vector<set<int>> &segment, vector<vector<int>> &plan) {
  int sum = accumulate(vis.begin(), vis.end(), 0);
  int mod = sum % denominator;
  vector<int> key = vis;
  key.emplace_back(mod);
  auto &other2 = rank1_multi_mem_plan_map[key];
  if (other2.empty()) {
    list<other_> other;//we can not use
    other.emplace_back(0, mod, vis, vector<int>{});
    for (auto i = other.begin(); i != other.end(); ++i) {
      auto &o = *i;
      int cnt = 0;
      int tor = o.tor;
      int beg = o.beg;
      for (int j = 0; j < o.left_c.size(); ++j) {
        cnt += o.left_c[j];
        if (cnt >= denominator) {
          cnt -= denominator;
        }
        if (cnt) {
          if (cnt <= tor && o.left_c.begin() + j + 1 < o.left_c.end()) {
            other.emplace_back(beg + j + 1,
                               tor - cnt,
                               vector<int>(o.left_c.begin() + j + 1, o.left_c.end()),
                               o.mem_c);
          }
          int rem = beg + j;
          if (rem != vis.size() - 1)
            o.mem_c.emplace_back(beg + j);//can not change the sequence!
        }
      }
    }
    other2.resize(other.size());
    transform(other.begin(), other.end(), other2.begin(), [](const other_ &o) {
      return o.mem_c;
    });
  }

  //clean segment
  for (int i = 1; i < Dim; ++i) {
    if (mem[i]) {
      for (auto &j : segment) {
        j.erase(i);
      }
    }
  }
  vector<set<int>> num_vis(Dim);
  for (int i = 0; i < other2.size(); ++i) {
    for (auto j : other2[i]) {
      for (auto k : segment[j]) {
        num_vis[k].emplace(i);
      }
    }
  }
  for (int i = 1; i < Dim; ++i) {
    if (num_vis[i].size() == other2.size()) {
      mem.set(i);
      for (auto &j : segment) {
        j.erase(i);
      }
    }
  }

  vector<pair<yzzLong, vector<int>>> mem_other;
  for (const auto &i : other2) {
    yzzLong p_mem = 0;
    for (auto j : i) {
      for (auto k : segment[j]) {
        p_mem.set(k);
      }
    }
    if (p_mem.none()) {
      plan.clear();
      goto QUIT;
    }
    for (auto &j : mem_other) {
      if (((p_mem & j.first) ^ p_mem).none()) {
        j = {p_mem, i};
        goto NEXT;
      }
    }
    mem_other.emplace_back(p_mem, i);
    NEXT:;
  }
  plan.resize(mem_other.size());
  transform(mem_other.begin(), mem_other.end(), plan.begin(),
            [](const auto &i) { return i.second; });
  QUIT:
  return;
}

void CVRP::findCombs4rankCutSet2getBestMultiplier(int plan_idx,
                                                  const vector<int> &cset,
                                                  vector<int> &new_cset,
                                                  int c,
                                                  double &new_vio,
                                                  const vector<vector<int>> &v_r_map,
                                                  const vector<double> &frac_routes,
                                                  bool &if_succeed) {
  //basic_info: <multiplier, num> this num is deducted by one since c is excluded
  //and new cset[0] must be c
  //we sort the basic_info by the num
  if_succeed = false;
  auto &plan = map_rank1_multiplier[(int) cset.size()][plan_idx];
  if (!get<1>(plan)) {
    return;
  }
  if_succeed = true;
  unordered_map<int, int> basic_info_map;
  for (auto i : get<0>(plan)) {
    ++basic_info_map[i];
  }
  int denominator = get<1>(plan);
  vector<pair<int, int>> basic_info(basic_info_map.begin(), basic_info_map.end());
  if (basic_info.size() == 1) {
    //multi can only be 1
    vector<double> num_vis(frac_routes.size(), 0);
    int mul = basic_info[0].first;
    for (auto i : cset) {
      for (auto j : v_r_map[i]) {
        num_vis[j] += mul;
      }
    }
    transform(num_vis.begin(),
              num_vis.end(),
              frac_routes.begin(),
              num_vis.begin(),
              [denominator](auto &a, auto &b) {
                return int(a / denominator + TOLERANCE) * b;
              });
    new_vio = accumulate(num_vis.begin(), num_vis.end(), double(-get<2>(plan)));
    new_cset = cset;
    return;
  }
  auto sorted_info = basic_info;
  stable_sort(sorted_info.begin(), sorted_info.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
    return a.second > b.second;
  });
  vector<int> sorted_num_info(sorted_info.size());
  transform(sorted_info.begin(), sorted_info.end(), sorted_num_info.begin(), [](const pair<int, int> &a) {
    return a.second;
  });
  vector<vector<vector<int>>> map;
  pure_map(map, sorted_num_info);
  int n = accumulate(sorted_num_info.begin(), sorted_num_info.end(), 0, [](int a, const auto &b) { return a + b; });
  vector<int> tmp(n);
  vector<double> frac_(frac_routes.size());
  vector<int> num_of_vis(frac_routes.size());
  double best_vio = 0;
  vector<int> best_plan;
  if (sorted_num_info.size() == 2) {
    for (auto &i : map[0]) {
      fill(tmp.begin(), tmp.end(), 1);
      for (auto j : i) tmp[j] = 0;
      memset(num_of_vis.data(), 0, sizeof(int) * num_of_vis.size());
      for (int j = 0; j < cset.size(); ++j) {
        int mul = sorted_info[tmp[j]].first;
        for (auto k : v_r_map[cset[j]]) {
          num_of_vis[k] += mul;
        }
      }
      transform(num_of_vis.begin(),
                num_of_vis.end(),
                frac_routes.begin(),
                frac_.begin(),
                [denominator](auto &a, auto &b) {
                  return int(a / denominator) * b;
                });
      double vio = accumulate(frac_.begin(), frac_.end(), double(-get<2>(plan)));
      if (vio > best_vio) {
        best_vio = vio;
        best_plan = tmp;
      }
    }
  } else if (sorted_num_info.size() == 3) {
    vector<int> supp(n - sorted_num_info[0]);
    for (auto &i : map[0]) {
      fill(tmp.begin(), tmp.end(), 2);
      for (auto j : i) tmp[j] = 0;
      int cnt = 0;
      for (int j = 0; j < n; ++j) if (tmp[j] == 2) supp[cnt++] = j;
      for (auto &j : map[1]) {
        auto tmp2 = tmp;
        for (auto k : j) tmp2[supp[k]] = 1;
        memset(num_of_vis.data(), 0, sizeof(int) * num_of_vis.size());
        for (int l = 0; l < cset.size(); ++l) {
          int mul = sorted_info[tmp2[l]].first;
          for (auto k : v_r_map[cset[l]]) {
            num_of_vis[k] += mul;
          }
        }
        transform(num_of_vis.begin(),
                  num_of_vis.end(),
                  frac_routes.begin(),
                  frac_.begin(),
                  [denominator](auto &a, auto &b) {
                    return int(a / denominator) * b;
                  });
        double vio = accumulate(frac_.begin(), frac_.end(), double(-get<2>(plan)));
        if (vio > best_vio) {
          best_vio = vio;
          best_plan = tmp2;//here is tmp2
        }
      }
    }
  }
  new_vio = best_vio;
  new_cset.resize(cset.size());
  int cnt = 0;
  vector<pair<int, int>> set_pair(cset.size());
  for (int i = 0; i < best_plan.size(); ++i) {
    set_pair[cnt++] = {cset[i], sorted_info[best_plan[i]].first};
  }
  std::stable_sort(set_pair.begin(), set_pair.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
    return a.second > b.second;
  });
  transform(set_pair.begin(), set_pair.end(), new_cset.begin(), [](const pair<int, int> &a) {
    return a.first;
  });
}

void CVRP::addInRank1HeurSep(int plan_idx,
                             const vector<int> &c,
                             const vector<int> &w_no_c,
                             double &new_vio,
                             int &choice_j,
                             const vector<vector<int>> &v_r_map,
                             const vector<double> &frac_routes) {
  int new_c_size = (int) c.size() + 1;
  const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
  if (new_c_size > CONFIG::MaxRowRank1 || !get<1>(plan)) {
    new_vio = -Dim;
    return;
  }

  choice_j = 0;

  int route_size = (int) frac_routes.size();
  int denominator = get<1>(plan);
  vector<int> num_use_route(route_size, 0);
  vector<double> frac_(route_size);
  for (int i = 0; i < c.size(); ++i) {
    for (auto j : v_r_map[c[i]]) {
      num_use_route[j] += get<0>(plan)[i];
    }
  }

  int add_idx = (int) c.size();
  for (auto j : w_no_c) {
    auto tmp = num_use_route;
    for (auto r : v_r_map[j]) {
      ++tmp[r];//coeff must be 1
    }
    transform(tmp.begin(), tmp.end(), frac_routes.begin(), frac_.begin(), [denominator](auto &a, auto &b) {
      return int(a / denominator) * b;
    });
    double af_vio = accumulate(frac_.begin(), frac_.end(), double(-get<2>(plan)));
    if (af_vio > new_vio + TOLERANCE) {
      new_vio = af_vio;
      choice_j = j;
    }
  }

  new_vio -= TOLERANCE;//penalty

  if (!choice_j) {
    new_vio = -Dim;
  }
}

void CVRP::removeInRank1HeurSep(int plan_idx,
                                const vector<int> &c,
                                double &new_vio,
                                int &choice_j,
                                const vector<vector<int>> &v_r_map,
                                const vector<double> &frac_routes) {
  int new_c_size = (int) c.size() - 1;
  const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
  if (new_c_size < 3 || !get<1>(plan)) {
    new_vio = -Dim;
    return;
  }
  choice_j = 0;
  int route_size = (int) frac_routes.size();
  int denominator = get<1>(plan);
  vector<int> num_use_route(route_size, 0);
  vector<double> frac_(route_size);
  int tie = 0;
  for (; tie < c.size(); ++tie) {
    if (get<0>(plan)[tie] == 1) break;// we only delete the node with coefficient 1
    for (auto j : v_r_map[c[tie]]) {
      num_use_route[j] += get<0>(plan)[tie];
    }
  }
  vector<int> movable_c(c.begin() + tie, c.end());
  for (auto j : movable_c) {
    auto tmp = num_use_route;
    for (auto r : v_r_map[j]) {
      --tmp[r];
    }
    transform(tmp.begin(), tmp.end(), frac_routes.begin(), frac_.begin(), [denominator](auto &a, auto &b) {
      return int(a / denominator) * b;
    });
    double af_vio = accumulate(frac_.begin(), frac_.end(), double(-get<2>(plan)));
    if (af_vio > new_vio + TOLERANCE) {
      new_vio = af_vio;
      choice_j = j;
    }
  }

  new_vio += TOLERANCE;//penalty

  if (!choice_j) {
    new_vio = -Dim;
  }
}

void CVRP::swapInRank1HeurSep(int plan_idx,
                              const vector<int> &c,
                              const vector<int> &w_no_c,
                              double &new_vio,
                              pair<int, int> &choice_i_j,
                              const vector<vector<int>> &v_r_map,
                              const vector<double> &frac_routes) {
  const auto &plan = map_rank1_multiplier[(int) c.size()][plan_idx];
  int new_c_size = (int) c.size();
  if ((new_c_size < 3 || new_c_size > CONFIG::MaxRowRank1) || !get<1>(plan)) {
    new_vio = -Dim;
    return;
  }

  choice_i_j.second = 0;
  int route_size = (int) frac_routes.size();
  int denominator = get<1>(plan);
  vector<int> num_use_route(route_size, 0);
  vector<double> frac_(route_size);
  int tie = 0;
  for (; tie < c.size(); ++tie) {
    if (get<0>(plan)[tie] == 1) break;// we only delete the node with coefficient 1
    for (auto j : v_r_map[c[tie]]) {
      num_use_route[j] += get<0>(plan)[tie];
    }
  }

  vector<int> w_no_c2;
  for (int i = tie; i < c.size(); ++i) {
    for (auto j : v_r_map[c[i]]) {
      ++num_use_route[j];
    }//fix these nodes
    w_no_c2.emplace_back(c[i]);
  }

  for (auto old_i : w_no_c2) {
    for (auto j : w_no_c) {
      auto tmp = num_use_route;
      for (auto r : v_r_map[old_i]) {
        --tmp[r];
      }
      for (auto r : v_r_map[j]) {
        ++tmp[r];
      }
      transform(tmp.begin(), tmp.end(), frac_routes.begin(), frac_.begin(), [denominator](auto &a, auto &b) {
        return int(a / denominator) * b;
      });
      double af_vio = accumulate(frac_.begin(), frac_.end(), double(-get<2>(plan)));

      if (af_vio > new_vio + TOLERANCE) {
        new_vio = af_vio;
        choice_i_j = {old_i, j};
      }
    }
  }
  if (!choice_i_j.second) {
    new_vio = -Dim;
  }
}

void CVRP::combinationUtil_addOne(const std::vector<int> &arr,
                                  std::vector<int> &tmp,
                                  std::vector<std::vector<int>> &data,
                                  int start,
                                  int end,
                                  int index,
                                  int r) {
  if (index == r) {
    data.emplace_back(tmp);
    return;
  }

  for (int i = start; end - i >= (r - index - 1); i++) {
    tmp[index] = arr[i];
    combinationUtil_addOne(arr, tmp, data, i + 1,
                           end, index + 1, r);
  }
}

void CVRP::comb_all(int n, int r, vector<vector<int>> &data) {
  vector<int> arr(n);
  iota(arr.begin(), arr.end(), 0);
  vector<int> tmp(r);
  combinationUtil_addOne(arr, tmp, data, 0, (int) arr.size() - 1, 0, r);
}

void CVRP::pure_map(vector<vector<vector<int>>> &map, const vector<int> &sorted_info) {
  int n, mul_idx;
  bool if_find = true;
  if (pureMap.find(sorted_info) != pureMap.end()) {
    map = pureMap[sorted_info];
    goto QUIT;
  } else {
    if_find = false;
  }

  if (sorted_info.size() == 1) {
    map.emplace_back(1, vector<int>(sorted_info[0], 0));
    goto QUIT;
  }

  n = accumulate(sorted_info.begin(), sorted_info.end(), 0, [](int a, const auto &b) { return a + b; });
  mul_idx = 0;
  map.resize(sorted_info.size());
  GET_LEFT:
  comb_all(n, sorted_info[mul_idx], map[mul_idx]);
  n -= sorted_info[mul_idx];
  ++mul_idx;
  if (mul_idx != sorted_info.size()) {
    goto GET_LEFT;
  }
  QUIT:
  if (!if_find) pureMap[sorted_info] = map;
}

void CVRP::NewgenerateR1C1s(const vector<vector<int>> &non_ele_routes,
                            const vector<double> &frac_none_ele_routes,
                            std::vector<std::pair<std::vector<int>, int>> &cuts) {
  if (non_ele_routes.empty())
    return;
//  return;
  yzzLong tmp_long;
  set<int> re_visited_vertices;
  unordered_map<int, set<int>> vertex2non_ele_routes;
  vertex2non_ele_routes.reserve(Dim);

  for (int r = 0; r < non_ele_routes.size(); ++r) {
    tmp_long = 0;
    for (auto i : non_ele_routes[r]) {
      if (tmp_long[i]) {
        re_visited_vertices.emplace(i);
        vertex2non_ele_routes[i].emplace(r);
      } else tmp_long.set(i);
    }
  }

  vector<pair<int, double >> vio_mem_vec_re_vertices(re_visited_vertices.size());//v & vio
  int cnt = 0;
  for (auto i : re_visited_vertices) {
    double vio = 0;
    for (auto r : vertex2non_ele_routes[i]) {
      double times = 0;
      for (auto node : non_ele_routes[r]) {
        if (node == i)++times;
      }
      vio += floor(times / 2 + TOLERANCE) * frac_none_ele_routes[r];
    }
    vio_mem_vec_re_vertices[cnt++] = {i, vio};
  }

  sort(vio_mem_vec_re_vertices.begin(),
       vio_mem_vec_re_vertices.end(),
       [](const pair<int, double> &a, const pair<int, double> &b) {
         return a.second > b.second;
       });

  vector<pair<vector<int>, int>> tmp_cut(vio_mem_vec_re_vertices.size());
  transform(vio_mem_vec_re_vertices.begin(), vio_mem_vec_re_vertices.end(), tmp_cut.begin(),
            [](const pair<int, double> &a) {
              return make_pair(vector<int>{a.first}, 0);
            });
#ifdef debugVio
  for (auto &cut : vio_mem_vec_re_vertices) {
    yzzLong tmp = 0;
    tmp.set(cut.first);
    if (vio_map.find(tmp) == vio_map.end()) {
      vio_map[tmp].resize(7, -numeric_limits<double>::max());
    }
    vio_map[tmp][0] = cut.second;
  }
#endif
  selectCuts(tmp_cut, cuts, CONFIG::MaxNumR1CPerRound);
}

void CVRP::NewgenerateR1C3s(const vector<vector<int>> &routes,
                            const vector<double> &frac_routes,
                            std::vector<std::pair<std::vector<int>, int>> &cuts) {
  /**
   * we only search for cuts that could potentially make the cuts have a small memory
   * for the routes, we can construct a residual graph
   * the cuts set {i,j,k}. for i, (i,j) or (i,k) must exist in the graph
   */

  int allNumRoutes = (int) routes.size();
  vector<unordered_map<int, int>> visited_times_by_routes(Dim);
  vector<unordered_set<int>> i_connections(Dim);//j,k that connects i
  int routeNumber = 0;
  for (const auto &ele_route : routes) {
    for (int j = 0; j < ele_route.size() - 1; ++j) {
      i_connections[ele_route[j]].emplace(ele_route[j + 1]);
      i_connections[ele_route[j + 1]].emplace(ele_route[j]);
      ++visited_times_by_routes[ele_route[j]][routeNumber];
    }
    ++visited_times_by_routes[ele_route.back()][routeNumber];
    ++routeNumber;
  }
  vector<vector<int>> i_connections_vec(Dim);
  transform(i_connections.begin(), i_connections.end(), i_connections_vec.begin(),
            [](const unordered_set<int> &ele) -> vector<int> {
              return {ele.begin(), ele.end()};
            });
  for (auto &i : i_connections_vec) {
    sort(i.begin(), i.end());
  }
  vector<int> v_tmp;
  vector<int> v_tmp2(allNumRoutes);//find the number of times visited the cut
  vector<pair<vector<int>, double>> cuts_pool;
  cuts_pool.reserve(2048);
  for (int i = 1; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      v_tmp.clear();
      if (i_connections[i].find(j) != i_connections[i].end()) {
        //if i and j are connected
        set_union(i_connections_vec[i].begin(), i_connections_vec[i].end(),
                  i_connections_vec[j].begin(), i_connections_vec[j].end(),
                  back_inserter(v_tmp));
      } else {
        set_intersection(i_connections_vec[i].begin(), i_connections_vec[i].end(),
                         i_connections_vec[j].begin(), i_connections_vec[j].end(),
                         back_inserter(v_tmp));
      }
      for (int k : v_tmp) {
        if (k <= i || k <= j) continue;
        //calculate the violation of i,j,k
        double vio = -1;
        memset(v_tmp2.data(), 0, sizeof(int) * allNumRoutes);
        for (auto &ele : visited_times_by_routes[i]) {
          v_tmp2[ele.first] += ele.second;
        }
        for (auto &ele : visited_times_by_routes[j]) {
          v_tmp2[ele.first] += ele.second;
        }
        for (auto &ele : visited_times_by_routes[k]) {
          v_tmp2[ele.first] += ele.second;
        }
        for (int l = 0; l < allNumRoutes; ++l) {
          vio += int(v_tmp2[l] / 2) * frac_routes[l];
        }
        if (vio > CutVioTolerance) {
          cuts_pool.emplace_back(vector<int>{i, j, k}, vio);
        }
      }
    }
  }

  sort(cuts_pool.begin(), cuts_pool.end(),
       [](const pair<vector<int>, double> &a, const pair<vector<int>, double> &b) -> bool {
         return a.second > b.second;
       });
  vector<pair<vector<int>, int>> tmp_cut(cuts_pool.size());
  transform(cuts_pool.begin(), cuts_pool.end(), tmp_cut.begin(),
            [](const pair<vector<int>, double> &ele) -> pair<vector<int>, int> {
              return {ele.first, 0};//the plan is 0
            });
#ifdef debugVio
  for (auto &cut : cuts_pool) {
    yzzLong tmp = 0;
    for (int i : cut.first) {
      tmp.set(i);
    }
    if (vio_map.find(tmp) == vio_map.end()) {
      vio_map[tmp].resize(7, -numeric_limits<double>::max());
    }
    vio_map[tmp][0] = cut.second;
  }
#endif
  selectCuts(tmp_cut, cuts, CONFIG::MaxNumR1C3PerRound);
}

void CVRP::construct_v_r_vec(const vector<vector<int>> &routes,
                             vector<vector<int>> &v_r_vec,
                             vector<unordered_map<int, int>> &v_r_map) const {
  int NumRoutes = (int) routes.size();
  v_r_vec.resize(Dim, vector<int>(NumRoutes, 0));
  v_r_map.resize(Dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      ++v_r_map[i][r];
    }
  }
  for (int i = 1; i < Dim; ++i) {
    for (auto &ele : v_r_map[i]) {
      v_r_vec[i][ele.first] = ele.second;
    }
  }
}

template<bool if_symmetry>
void CVRP::construct_v_resourceGap(const vector<vector<int>> &routes, vector<double> &v_resourceGap) {
  vector<pair<double, int>> supp_v_resourceGap(Dim);
  for (auto &r : routes) {
    int past_node = 0;
    int cur_node;
    double res = 0;
    for (int it : r) {
      cur_node = it;
      increaseMainResourceConsumption(res, res, past_node, cur_node);
      if (res > MeetPointResourceInBiDir) break;
      supp_v_resourceGap[cur_node].first += MeetPointResourceInBiDir - res;
      ++supp_v_resourceGap[cur_node].second;
      past_node = cur_node;
    }
    past_node = 0;
    res = if_symmetry ? 0 : MaxMainResource;
    for (auto it = r.rbegin(); it != r.rend(); ++it) {
      cur_node = *it;
      if_symmetry ? increaseMainResourceConsumption(res, res, past_node, cur_node)
                  : decreaseMainResourceConsumption(res, res, past_node, cur_node);
      if ((if_symmetry ? MeetPointResourceInBiDir - res : res - MeetPointResourceInBiDir) < 0) break;
      supp_v_resourceGap[cur_node].first += if_symmetry ? MeetPointResourceInBiDir - res
                                                        : res - MeetPointResourceInBiDir;
      ++supp_v_resourceGap[cur_node].second;
      past_node = cur_node;
    }
  }
  v_resourceGap.resize(Dim);
  transform(supp_v_resourceGap.begin(), supp_v_resourceGap.end(), v_resourceGap.begin(),
            [](const pair<double, int> &p) -> double { return p.first / p.second; });
}

void CVRP::construct_seed(const vector<vector<int>> &routes,
                          const vector<vector<int>> &v_r_vec,
                          const vector<double> &v_resourceGap,
                          vector<vector<int>> &seed) {
  unordered_set < yzzLong > seed_pool;
  seed_pool.reserve(1024);
  int route_size = (int) routes.size();
  vector<vector<vector<int>>> c(Dim, vector<vector<int>>(route_size));
  unordered_map<int, vector<pair<int, int>>> c_map;//size, c_index
  vector<int> no_route;
  no_route.reserve(route_size);
  for (int i = 1; i < Dim; ++i) {
    yzzLong wc = 0;// c within i
    no_route.clear();
    for (int r = 0; r < route_size; ++r) {
      if (!v_r_vec[i][r]) no_route.emplace_back(r);
      else {
        for (auto v : routes[r]) {
          if (Rank1SepHeurMem4Vertex[i][v]) wc.set(v);
        }
      }
    }
    for (auto &r : no_route) {
      vector<pair<int, double>> tmp_seed;
      yzzLong yzz_seed = 0;
      for (auto v : routes[r]) {
        if (wc[v] && !yzz_seed.test(v)) {
          tmp_seed.emplace_back(v, v_resourceGap[v]);
          yzz_seed.set(v);
        }
      }
      tmp_seed.emplace_back(i, v_resourceGap[i]);
      yzz_seed.set(i);
      if (tmp_seed.size() < 3 || tmp_seed.size() > CONFIG::MaxRowRank1)
        continue;// in this case, we use do not allow extension!
      if (seed_pool.find(yzz_seed) != seed_pool.end()) continue;
      sort(tmp_seed.begin(), tmp_seed.end(), [](const pair<int, double> &a, const pair<int, double> &b) -> bool {
        return a.second < b.second;
      });
      seed_pool.emplace(yzz_seed);
      //insert tmp_seed.first to seed
      seed.emplace_back();
      seed.back().resize(tmp_seed.size());
      transform(tmp_seed.begin(), tmp_seed.end(), seed.back().begin(), [](const pair<int, double> &p) -> int {
        return p.first;
      });
    }
  }
//  cout << "这里有测试！" << endl;
//  for (auto &cut : seed) {
//    unordered_set<int> a;
//    for (auto &v : cut) {
//      if (a.find(v) != a.end()) {
//        throw std::runtime_error("seed has duplicate vertex!");
//      }
//      a.emplace(v);
//    }
//  }
}

void CVRP::construct_cuts_by_resourceGap_only(const vector<unordered_map<int, int>> &v_r_map,
                                              const vector<vector<int>> &seed,
                                              const vector<double> &frac_routes,
                                              vector<pair<vector<int>, int>> &cuts) {
  vector<double> num_times_vis_routes(frac_routes.size());
  vector<tuple<vector<int>, int, double>> cuts_pool;
  int std_size = (CONFIG::MaxRowRank1 - 3) * CONFIG::MaxNumR1CPerRound;
  cuts_pool.reserve(std_size);
  //max_vio= (n-1)/n, n is the denominator, the gap is less, is the better.
  double vio_threshold = -1;
  int route_size = (int) frac_routes.size();
  for (auto &cut : seed) {
    int cut_size = (int) cut.size();
    double max_vio = vio_threshold;
    int max_idx = -1;
    for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
      auto &plan = map_rank1_multiplier[cut_size][plan_idx];
      if (!get<1>(plan)) continue;
      memset(num_times_vis_routes.data(), 0, sizeof(double) * route_size);
      auto &coeff = get<0>(plan);
      int denominator = get<1>(plan);
//      double self_cmp = double(denominator - 1) / denominator;
      double rhs = get<2>(plan);//must be double
      for (int i = 0; i < cut_size; ++i) {
        int c = cut[i];
        for (auto &p : v_r_map[c]) {
          num_times_vis_routes[p.first] += p.second * coeff[i];
        }
      }
      transform(num_times_vis_routes.begin(),
                num_times_vis_routes.end(),
                frac_routes.begin(),
                num_times_vis_routes.begin(),
                [denominator](double d, double f) -> double { return int(d / denominator) * f; });
      double vio =
          accumulate(num_times_vis_routes.begin(),
                     num_times_vis_routes.end(),
                     -rhs);
      if (vio > CutVioTolerance) {
//        vio -= self_cmp;
//now we only compare the vio!
        if (vio > max_vio) {
          max_vio = vio;
          max_idx = plan_idx;
        }
      }
    };
    if (max_idx == -1) continue;
    cuts_pool.emplace_back(cut, max_idx, max_vio);
    if (cuts_pool.size() >= std_size) {
      nth_element(cuts_pool.begin(), cuts_pool.begin() + std_size - 1, cuts_pool.end(),
                  [](const tuple<vector<int>, int, double> &a,
                     const tuple<vector<int>, int, double> &b) -> bool {
                    return get<2>(a) > get<2>(b);
                  });
      vio_threshold = get<2>(cuts_pool[std_size - 1]);
    }
  }
  sort(cuts_pool.begin(), cuts_pool.end(), [](const tuple<vector<int>, int, double> &a,
                                              const tuple<vector<int>, int, double> &b) -> bool {
    return get<2>(a) > get<2>(b);
  });
  if (cuts_pool.size() > std_size) cuts_pool.resize(std_size);
  cuts.resize(cuts.size() + cuts_pool.size());
  transform(cuts_pool.begin(),
            cuts_pool.end(),
            cuts.end() - (int) cuts_pool.size(),
            [](const tuple<vector<int>, int, double> &t) -> pair<vector<int>, int> {
              return make_pair(get<0>(t), get<1>(t));
            });
//  cout << "cut= " << endl;
//  for (auto &cut : cuts_pool) {
//    for (auto &c : get<0>(cut)) {
//      cout << c << " ";
//    }
//    cout << get<2>(cut) << endl;
//  }
//  cout << "cuts_pool= " << cuts_pool.size() << endl;
}

void CVRP::construct_cuts(const std::vector<std::unordered_map<int, int>> &v_r_map,
                          const vector<vector<int>> &seed,
                          const vector<double> &frac_routes,
                          vector<pair<vector<int>, int>> &cuts) {
//  time1 = time2 = time3 = time4 = time5 = 0;
  vector<tuple<vector<int>, int, double>> cuts_pool;
  int std_size = (CONFIG::MaxRowRank1 - 3) * CONFIG::MaxNumR1CPerRound;
  cuts_pool.reserve(std_size);
  //max_vio= (n-1)/n, n is the denominator, the gap is less, is the better.
  double vio;
  int route_size = (int) frac_routes.size();
//  double total_time = 0;
  int number = 0;
  for (auto &cut : seed) {
    int cut_size = (int) cut.size();
    double max_vio = TOLERANCE;
    int max_idx = -1;
    vector<int> best_cut;
    for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
      auto &plan = map_rank1_multiplier[cut_size][plan_idx];
      if (!get<1>(plan)) continue;
      vector<int> new_cut;
//      auto beg = chrono::high_resolution_clock::now();
      heuristicFindBestPermuatation4Oneplan(v_r_map, frac_routes, cut, plan_idx, new_cut, vio);
//      auto end = chrono::high_resolution_clock::now();
//      total_time += chrono::duration<double>(end - beg).count();
      if (vio > CutVioTolerance) {
        if (vio > max_vio) {
          max_vio = vio;
          max_idx = plan_idx;
          best_cut = new_cut;
        }
      }
    }
    if (max_idx == -1) continue;
    ++number;
    cuts_pool.emplace_back(best_cut, max_idx, max_vio);
  }
  sort(cuts_pool.begin(), cuts_pool.end(), [](const tuple<vector<int>, int, double> &a,
                                              const tuple<vector<int>, int, double> &b) -> bool {
    return get<2>(a) > get<2>(b);
  });
  vector<pair<vector<int>, int>> tmp_cuts(cuts_pool.size());
  transform(cuts_pool.begin(),
            cuts_pool.end(),
            tmp_cuts.begin(),
            [](const tuple<vector<int>, int, double> &t) -> pair<vector<int>, int> {
              return make_pair(get<0>(t), get<1>(t));
            });
  selectCuts(tmp_cuts, cuts, std_size);
//  cout << "max_vio_time= " << total_time << endl;
//  cout << "time1= " << time1 << ", time2= " << time2 << ", time3= " << time3 << ", time4= " << time4 << ", time5= "
//       << time5 << endl;
}

void CVRP::heuristicFindBestPermuatation4Oneplan(
    const std::vector<std::unordered_map<int, int>> &v_r_map,
    const vector<double> &frac_routes,
    const vector<int> &cut, int plan_idx,
    vector<int> &new_cut,
    double &vio) {
  bool is_sorted = std::is_sorted(frac_routes.begin(), frac_routes.end(), greater<>());
  if (!is_sorted) throw runtime_error("frac_routes is not sorted!");
  //get plan information
  int cut_size = (int) cut.size();
  auto &plan = map_rank1_multiplier[cut_size][plan_idx];
  if (!get<1>(plan)) throw runtime_error("plan is not valid!");
  auto &coeff = get<0>(plan);
  int denominator = get<1>(plan);
  double rhs = get<2>(plan);//must be double
  unordered_map<int, int> map_coeff_kind_num;
  map_coeff_kind_num.reserve(cut_size);
  for (int i = 0; i < cut_size; ++i) ++map_coeff_kind_num[coeff[i]];
  vector<double> num_times_vis_routes(frac_routes.size(), 0);
//  auto beg = chrono::high_resolution_clock::now();
  if (map_coeff_kind_num.size() == 1) {
    new_cut = cut;
    for (int i = 0; i < cut_size; ++i) {
      int c = cut[i];
      for (auto &pr : v_r_map[c]) {
        num_times_vis_routes[pr.first] += pr.second * coeff[i];
      }
    }
    transform(num_times_vis_routes.begin(),
              num_times_vis_routes.end(),
              frac_routes.begin(),
              num_times_vis_routes.begin(),
              [denominator](double d, double f) -> double { return int(d / denominator) * f; });
    vio = accumulate(num_times_vis_routes.begin(),
                     num_times_vis_routes.end(),
                     -rhs);
//    auto end = chrono::high_resolution_clock::now();
//    time1 += chrono::duration<double>(end - beg).count();
    return;
  }
//  beg = chrono::high_resolution_clock::now();
  // the order is important! If the best route even smaller than rhs, then we return the vio as 0.
  vector<vector<pair<int, int>>> tmp_vec_r_v(frac_routes.size());
  for (auto c : cut) {
    for (auto &pr : v_r_map[c]) {
      tmp_vec_r_v[pr.first].emplace_back(c, pr.second);
    }
  }
  vector<pair<int, vector<pair<int, int>>>> vec_r_v;
  vec_r_v.reserve(frac_routes.size());
  for (int r = 0; r < frac_routes.size(); ++r) {
    if (tmp_vec_r_v[r].size() < 2) continue;
    vec_r_v.emplace_back(r, std::move(tmp_vec_r_v[r]));
  }

  //notice! here must use stable_sort!
  for (auto &pr : vec_r_v) {
    auto &vec_v = pr.second;
    bool if_sorted = std::is_sorted(vec_v.begin(), vec_v.end(),
                                    [](const pair<int, int> &a,
                                       const pair<int, int> &b) -> bool {
                                      return a.second > b.second;
                                    });
    if (!if_sorted) {
      stable_sort(vec_v.begin(), vec_v.end(),
                  [](const pair<int, int> &a,
                     const pair<int, int> &b) -> bool {
                    return a.second > b.second;
                  });
    }
  }
//  auto end = chrono::high_resolution_clock::now();
//  time2 += chrono::duration<double>(end - beg).count();
  new_cut.clear();
  unordered_map<int, int> new_cut_coeff;
  auto cp_coeff = coeff;
  for (auto &pr : vec_r_v) {
    auto &vec_v = pr.second;
    //in order and reverse the order about the unknown vertex!
    int cnt = 0, cnt2 = 0;
    int ccnt = 0;
    vector<pair<int, int>> v_left;
    for (auto &i : vec_v) {
      int c = i.first;
//      cout << "c= " << c << " times= " << i.second << endl;
//      cout << "cp_coeff= ";
//      for (auto &j : cp_coeff) cout << j << " ";
//      cout << endl;
      auto if_find = new_cut_coeff.find(c);
      if (if_find != new_cut_coeff.end()) {
        int times = i.second * new_cut_coeff[c];
        cnt += times;
        cnt2 += times;
      } else {
        cnt += i.second * cp_coeff[ccnt];
        cnt2 += i.second * cp_coeff[cp_coeff.size() - 1 - ccnt];
        v_left.emplace_back(c, cp_coeff[ccnt]);
        ++ccnt;
      }
    }
//    cout << "----------------" << endl;
    cnt /= denominator;
    cnt2 /= denominator;
    if (cnt != cnt2) {
//      cout << "cnt= " << cnt << " cnt2= " << cnt2 << endl;
      //we know the partial sequence
      for (auto &i : v_left) {
        new_cut.emplace_back(i.first);
        --map_coeff_kind_num[i.second];
        new_cut_coeff[i.first] = i.second;
        cp_coeff.erase(find(cp_coeff.begin(), cp_coeff.end(), i.second));
      }
      for (auto it = map_coeff_kind_num.begin(); it != map_coeff_kind_num.end();) {
        if (it->second == 0) {
          it = map_coeff_kind_num.erase(it);
        } else {
          ++it;
        }
      }
//      cout << "map_coeff_kind_num.size()= " << map_coeff_kind_num.size() << endl;
      if (map_coeff_kind_num.size() <= 1) {
        for (auto &c : cut) {
          if (std::find(new_cut.begin(), new_cut.end(), c) == new_cut.end())
            new_cut.emplace_back(c);
        }
      }
//      for (auto &c : new_cut) {
//        cout << c << " ";
//        cout << "coeff= " << new_cut_coeff[c] << endl;
//      }
//      cout << endl;
    }
    if (new_cut.size() == cut.size()) goto QUIT;
  }
  //keep the old sequence
  new_cut = cut;
  QUIT:
//  beg = chrono::high_resolution_clock::now();
  for (int i = 0; i < cut_size; ++i) {
    int c = new_cut[i];
    for (auto &pr : v_r_map[c]) {
      num_times_vis_routes[pr.first] += pr.second * coeff[i];
    }
  }
  transform(num_times_vis_routes.begin(),
            num_times_vis_routes.end(),
            frac_routes.begin(),
            num_times_vis_routes.begin(),
            [denominator](double d, double f) -> double { return int(d / denominator) * f; });
  vio = accumulate(num_times_vis_routes.begin(),
                   num_times_vis_routes.end(),
                   -rhs);
//  end = chrono::high_resolution_clock::now();
//  time5 += chrono::duration<double>(end - beg).count();
}

void CVRP::newfindR1C_multi(const std::vector<std::vector<int>> &routes,
                            const std::vector<double> &frac_routes,
                            std::vector<std::unordered_map<int, int>> &v_r_map,
                            std::vector<std::pair<std::vector<int>, int>> &cuts) {
  /**
   * (1) for each i, find the r passed i and the times, and the routes don't pass i
   * and the average resource in the bucket, the closet to the meet point, the better
   * (2) search strategy is quite straightforward,
   * find j's in the same route as a set, and for each non-vis route find intersection,
   * combing i, the set is initialized as the basic set C
   * (3) the sequence could be decided by rule (1)
   * for every multiplier, we search for the best violation, we only collect the seed if the violation is above std
   * (4) the std is decided by the number of cuts, std= min(cut_pool), if cut_pool is not full, then it is zero
   * (5) when all the seeds have been collected, we tried to use dominance rule to avoid the permutation,
   * (6) in this version, we do not consider to do the extension!
   * (7) if for all the routes, the coefficient is the same, we abandon the corresponding cut
   */

  /**
   * vector<vector<int>> v_r_vec;//v_r_vec[i] is the number of times the routes passed i
   * vector<double> v_resourceGap;//v_resourceGap[i] is the resource gap of i, Mean(abs(MeetPoint - aver_res_now))
   */

  vector<vector<int>> v_r_vec;
  vector<double> v_resourceGap;
  construct_v_r_vec(routes, v_r_vec, v_r_map);
#ifdef SYMMETRY_PROHIBIT
  construct_v_resourceGap<false>(routes, v_resourceGap);
#else
  construct_v_resourceGap<true>(routes, v_resourceGap);
#endif
  vector<vector<int>> seed;
  construct_seed(routes, v_r_vec, v_resourceGap, seed);
//  cout << "seed size: " << seed.size() << endl;
//  for (auto &s : seed) {
//    for (auto v : s) {
//      cout << v << " ";
//    }
//    cout << endl;
//  }
//  cout << "pass 167 first route:" << endl;
//  int cnt = 0;
//  for (auto &i : v_r_vec[167]) {
//    if (i != 0) {
//      for (auto j : routes[cnt]) {
//        cout << j << " ";
//      }
//      break;
//    }
//    ++cnt;
//  }
//  cout << "233: " << v_resourceGap[233] << endl;
//  cout << "167: " << v_resourceGap[167] << endl;
//  cout << "11: " << v_resourceGap[11] << endl;
//  cout << "find all 233 routes: " << endl;
//  cnt = 0;
//  for (auto &i : v_r_vec[233]) {
//    if (i != 0) {
//      for (auto j : routes[cnt]) {
//        cout << j << " ";
//      }
//      cout << endl;
//    }
//    ++cnt;
//  }
//  cout << "find all 11 routes: " << endl;
//  cnt = 0;
//  for (auto &i : v_r_vec[11]) {
//    if (i != 0) {
//      for (auto j : routes[cnt]) {
//        cout << j << " ";
//      }
//      cout << endl;
//    }
//    ++cnt;
//  }
//  cout << endl;
  construct_cuts(v_r_map, seed, frac_routes, cuts);
}

void CVRP::addR1C_atOnce(BBNODE *node,
                         const vector<tuple<vector<int>, int, int, set<int>>> &full_cuts
) {
  //we only traverse the sequence for once!
  vector<vector<int>> coeffs(full_cuts.size(), vector<int>(NumCol, 0));
  vector<int> states(full_cuts.size());
  vector<vector<pair<int, int>>> v_cuts(Dim);//(idx, increase)
  vector<vector<bool>> v_Mem(Dim, vector<bool>(full_cuts.size(), false));
  vector<int> denominator(full_cuts.size());
  vector<int> rhs(full_cuts.size());
  for (int i = 0; i < full_cuts.size(); ++i) {
    auto &cut = full_cuts[i];
    const auto &big_plan = map_rank1_multiplier[(int) get<0>(cut).size()][get<1>(cut)];
    const auto &plan = get<0>(big_plan);
    denominator[i] = get<1>(big_plan);
    rhs[i] = get<2>(big_plan);
    int cnt = 0;
    for (auto &c : get<0>(cut)) {
      v_cuts[c].emplace_back(i, plan[cnt]);
      v_Mem[c][i] = true;
      ++cnt;
    }
    for (auto &m : get<3>(cut)) {
      v_Mem[m][i] = true;
    }
    coeffs[i][0] = rhs[i];
  }

  vector<vector<int>> v_lostMem(Dim);
  for (int i = 1; i < Dim; ++i) {
    for (int j = 0; j < full_cuts.size(); ++j) {
      if (!v_Mem[i][j]) {
        v_lostMem[i].emplace_back(j);
      }
    }
  }

  int beg = 1, end = node->NumParentCols, if_rep = false;
  int *col_pool = ColPool4Mem;
  here:
  for (int i = beg; i < end; ++i) {
    memset(states.data(), 0, sizeof(int) * full_cuts.size());
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int current_node = col_pool[j];
      if (!current_node) break;
      for (auto &p : v_cuts[current_node]) {
        states[p.first] += p.second;
        if (states[p.first] >= denominator[p.first]) {
          ++coeffs[p.first][i];
          states[p.first] -= denominator[p.first];
        }
      }
      for (auto &p : v_lostMem[current_node]) states[p] = 0;
    }
  }

  if (!if_rep) {
    if_rep = true;
    beg = end;
    end = NumCol;
    col_pool = ColPool4Pricing;
    goto here;
  }

#ifdef debugVio
  vector<double> frac(NumCol, 0);
  unordered_map<size_t, double> tmp1;
  for (auto &i : node->Idx4LPSolsInColPool) {
    tmp1[i.first] = i.second;
  }
  for (int i = 0; i < NumCol; ++i) {
    auto idx = node->IdxCols[i];
    if (tmp1.find(idx) != tmp1.end()) {
      frac[i] = tmp1[idx];
    }
  }
  vector<double> tmp_sum(NumCol);
#endif
  for (int i = 0; i < full_cuts.size(); ++i) {
    int numnz = 0;
    for (int j = 0; j < NumCol; ++j) {
      if (coeffs[i][j] != 0) {
        solver_ind[numnz] = j;
        solver_val[numnz] = coeffs[i][j];
        ++numnz;
      }
    }
#ifdef debugVio
    transform(coeffs[i].begin(), coeffs[i].end(), frac.begin(), tmp_sum.begin(), multiplies<>());
    double sum = accumulate(tmp_sum.begin(), tmp_sum.end(), 0.0);
    auto vio = sum - rhs[i];
    auto &cut = get<0>(full_cuts[i]);
    auto &plan_idx = get<1>(full_cuts[i]);
    yzzLong tmp = 0;
    for (auto &j : cut) {
      tmp.set(j);
    }
    if (abs(vio - vio_map[tmp][plan_idx]) > 0.01) {
      cout << "sum: " << sum << " rhs: " << rhs[i] << " vio: " << vio << " vio_map: " << vio_map[tmp][plan_idx] << endl;
//      throw runtime_error("cut_vio3 error");
      cout << "cut: ";
      for (auto &j : cut) {
        cout << j << " ";
      }
      cout << "plan_idx: " << plan_idx << endl;
      if (vio > vio_map[tmp][plan_idx] + 0.01) {
        cout << "new error!" << endl;
//        throw runtime_error("new error!");
      }
    }
#endif
    int cut_index = get<2>(full_cuts[i]);
    if (cut_index == numeric_limits<int>::max()) {
      safe_solver(node->solver.SOLVERaddconstr(numnz,
                                               solver_ind,
                                               solver_val,
                                               SOLVER_LESS_EQUAL,
                                               rhs[i],
                                               nullptr))
      if (get<1>(full_cuts[i]) == 0) {
        R1C r1c3;
        r1c3.InfoR1C = get<0>(full_cuts[i]);
        r1c3.Mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
        r1c3.IdxR1C = NumRow++;
        r1c3.RHS = rhs[i];
        node->R1Cs.emplace_back(r1c3);
      } else {
        R1C_multi r1c;
        r1c.InfoR1C = make_pair(get<0>(full_cuts[i]), get<1>(full_cuts[i]));
        r1c.Mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
        r1c.IdxR1C = NumRow++;
        r1c.RHS = rhs[i];
        node->R1Cs_multi.emplace_back(r1c);
      }
    } else {
      if (get<1>(full_cuts[i]) == 0) {
        fill(solver_ind2, solver_ind2 + numnz, node->R1Cs[cut_index].IdxR1C);
        node->R1Cs[cut_index].Mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
      } else {
        fill(solver_ind2, solver_ind2 + numnz, node->R1Cs_multi[cut_index].IdxR1C);
        node->R1Cs_multi[cut_index].Mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
      }
      safe_solver(node->solver.SOLVERXchgcoeffs(numnz, solver_ind2, solver_ind, solver_val))
    }
  }
}

void CVRP::generateR1C3s_aggressive(const vector<vector<int>> &ele_routes,
                                    const vector<vector<int>> &non_ele_routes,
                                    const vector<double> &frac_ele_routes,
                                    const vector<double> &frac_non_ele_routes,
                                    vector<pair<vector<int>, int>> &cuts) const {
  auto beg = chrono::high_resolution_clock::now();
  int num_cut = 0;
  int ai, aj, ak;
  int num_ele_routes = (int) ele_routes.size(), num_non_ele_routes = (int) non_ele_routes.size();
  int cnt;
  double frac;
  yzzLong tmp_long, PI_ai, PI_aj, PI_ai_aj;
  vector<int> aux(3);

  vector<vector<int>> vertex2non_ele_routes(Dim);
  for (int i = 0; i < Dim; ++i) vertex2non_ele_routes.reserve(num_non_ele_routes);

  unordered_map<yzzLong, tuple<int, int, int, double>> three_combinations;
  three_combinations.reserve(int(pow_self(MaxLengthEleRoute, 3) / 6));
  unordered_map<yzzLong, tuple<int, int, int, double >>::iterator
      iter3;
  unordered_map<yzzLong, tuple<int, int, double >> two_combinations;
  two_combinations.reserve(int(pow_self(MaxLengthEleRoute, 2) / 2));
  unordered_map<yzzLong, tuple<int, int, double >>::iterator
      iter2;
  unordered_map<yzzLong, tuple<int, int, int, double>> all_combinations;
  all_combinations.reserve(int(pow_self(MaxLengthEleRoute, 3) / 6));
  unordered_map<yzzLong, tuple<int, int, int, double >>::iterator
      all_iter;

  for (int i = 0; i < non_ele_routes.size(); ++i) {
    for (auto j : non_ele_routes[i]) {
      vertex2non_ele_routes[j].emplace_back(i);
    }
  }

  //find all possible combinations in ele-routes
  for (int i = 0; i < num_ele_routes; ++i) {
    auto &seq = ele_routes[i];
    cnt = (int) ele_routes[i].size();
    frac = frac_ele_routes[i];
    //write data into two combinations
    for (int j = 0; j < cnt; ++j) {
      for (int k = j + 1; k < cnt; ++k) {
        tmp_long = 0;
        tmp_long.set(seq[j]);
        tmp_long.set(seq[k]);
        iter2 = two_combinations.find(tmp_long);
        if (iter2 == two_combinations.end()) {
          if (seq[j] < seq[k])two_combinations.emplace(tmp_long, make_tuple(seq[j], seq[k], frac));
          else two_combinations.emplace(tmp_long, make_tuple(seq[k], seq[j], frac));
        } else get<2>(iter2->second) += frac;
      }
    }
    //write data into three combinations
    for (int j = 0; j < cnt; ++j) {
      for (int k = j + 1; k < cnt; ++k) {
        for (int l = k + 1; l < cnt; ++l) {
          tmp_long = 0;
          tmp_long.set(seq[j]);
          tmp_long.set(seq[k]);
          tmp_long.set(seq[l]);
          iter3 = three_combinations.find(tmp_long);
          if (iter3 == three_combinations.end()) {
            //sort
            aux[0] = seq[j];
            aux[1] = seq[k];
            aux[2] = seq[l];
            std::stable_sort(aux.begin(), aux.end());
            three_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
          } else get<3>(iter3->second) += frac;
        }
      }
    }
  }
  //e.g. find 1,2 and 1,3. we use 1,2 + 1,3 + 2,3 -2* 1,2,3 to determine the violation (of course if 2,3 and 1,2,3 exist)
  for (auto &two_comb : two_combinations) {
    ai = get<0>(two_comb.second);
    aj = get<1>(two_comb.second);
    PI_ai_aj = 0;
    PI_ai_aj.set(ai);
    PI_ai_aj.set(aj);
    PI_ai = 0;
    PI_aj = 0;
    PI_ai.set(ai);
    PI_aj.set(aj);
    //find pair
    for (int j = 1; j < ai; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      //find if counted already
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      //find 1,3
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 2,3
      tmp_long = PI_aj;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 1,2,3
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      iter3 = three_combinations.find(tmp_long);
      if (iter3 != three_combinations.end()) frac -= 2 * get<3>(iter3->second);
      aux[0] = ai;
      aux[1] = aj;
      aux[2] = j;
      std::stable_sort(aux.begin(), aux.end());
      all_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
    }
    for (int j = ai + 1; j < aj; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      //find if counted already
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      //find 1,3
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 2,3
      tmp_long = PI_aj;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 1,2,3
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      iter3 = three_combinations.find(tmp_long);
      if (iter3 != three_combinations.end()) frac -= 2 * get<3>(iter3->second);
      //if frac
      aux[0] = ai;
      aux[1] = aj;
      aux[2] = j;
      std::stable_sort(aux.begin(), aux.end());
      all_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
    }
    for (int j = aj + 1; j < Dim; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      //find if counted already
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      //find 1,3
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 2,3
      tmp_long = PI_aj;
      tmp_long.set(j) = 1;
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      //find 1,2,3
      tmp_long = PI_ai_aj;
      tmp_long.set(j) = 1;
      iter3 = three_combinations.find(tmp_long);
      if (iter3 != three_combinations.end()) frac -= 2 * get<3>(iter3->second);
      aux[0] = ai;
      aux[1] = aj;
      aux[2] = j;
      std::stable_sort(aux.begin(), aux.end());
      all_combinations.emplace(tmp_long, make_tuple(aux[0], aux[1], aux[2], frac));
    }
  }
  //traverse the combinations to calculate the final violations (check if there is a need for optimization)

  vector<double> tmp_com_info_non_ele(num_non_ele_routes);
  vector<tuple<int, int, int, double>> violated_combs;
  violated_combs.reserve(all_combinations.size());
  for (auto &all_comb : all_combinations) {
    auto &comb_frac = get<3>(all_comb.second);
    memset(tmp_com_info_non_ele.data(), 0, sizeof(double) * num_non_ele_routes);
    ai = get<0>(all_comb.second);
    aj = get<1>(all_comb.second);
    ak = get<2>(all_comb.second);
    for (auto item : vertex2non_ele_routes[ai])++tmp_com_info_non_ele[item];
    for (auto item : vertex2non_ele_routes[aj])++tmp_com_info_non_ele[item];
    for (auto item : vertex2non_ele_routes[ak])++tmp_com_info_non_ele[item];
    for (int i = 0; i < tmp_com_info_non_ele.size(); ++i)
      comb_frac += floor(tmp_com_info_non_ele[i] / 2 + TOLERANCE) * frac_non_ele_routes[i];
    if (comb_frac - 1 > TOLERANCE)violated_combs.emplace_back(ai, aj, ak, comb_frac);
  }

  std::stable_sort(violated_combs.begin(),
                   violated_combs.end(),
                   [](const std::tuple<int, int, int, double> &a, const std::tuple<int, int, int, double> &b) {
                     return get<3>(a) > get<3>(b);
                   });
  int limit = min((int) violated_combs.size(), CONFIG::MaxNumR1C3PerRound);
  cuts.resize(cuts.size() + limit);
  transform(violated_combs.begin(), violated_combs.begin() + limit, cuts.end() - limit,
            [](const std::tuple<int, int, int, double> &a) {
              return make_pair(vector<int>{get<0>(a), get<1>(a), get<2>(a)}, 0);
            });
  auto end = chrono::high_resolution_clock::now();
  cout << "R1C3 time: " << chrono::duration<double>(end - beg).count() << endl;
}

void CVRP::checkR1Ctotally(BBNODE *node) {
  cout << "WARNING: checkR1Ctotally" << endl;
  int num_rows = int(node->R1Cs.size() + node->R1Cs_multi.size());
  vector<vector<int>> coeffs(num_rows, vector<int>(NumCol, 0));
  vector<int> states(num_rows);
  vector<vector<pair<int, int>>> v_cuts(Dim);//(idx, increase)
  vector<vector<bool>> v_Mem(Dim, vector<bool>(num_rows, false));
  vector<int> denominator(num_rows);
  vector<int> rhs(num_rows);
  vector<tuple<vector<int>, int, int, vector<int>>> full_cuts;// cut, plan, idx, mem
  for (auto &r1c : node->R1Cs) {
    full_cuts.emplace_back(r1c.InfoR1C, 0, r1c.IdxR1C, r1c.Mem);
  }
  for (auto &r1c : node->R1Cs_multi) {
    full_cuts.emplace_back(r1c.InfoR1C.first, r1c.InfoR1C.second, r1c.IdxR1C, r1c.Mem);
  }
  for (int i = 0; i < num_rows; ++i) {
    auto &cut = full_cuts[i];
    const auto &big_plan = map_rank1_multiplier[(int) get<0>(cut).size()][get<1>(cut)];
    const auto &plan = get<0>(big_plan);
    denominator[i] = get<1>(big_plan);
    rhs[i] = get<2>(big_plan);
    int cnt = 0;
    for (auto &c : get<0>(cut)) {
      v_cuts[c].emplace_back(i, plan[cnt]);
      v_Mem[c][i] = true;
      ++cnt;
    }
    for (auto &m : get<3>(cut)) {
      v_Mem[m][i] = true;
    }
    coeffs[i][0] = rhs[i];
  }

  vector<vector<int>> v_lostMem(Dim);
  for (int i = 1; i < Dim; ++i) {
    for (int j = 0; j < full_cuts.size(); ++j) {
      if (!v_Mem[i][j]) {
        v_lostMem[i].emplace_back(j);
      }
    }
  }

  int beg = 1, end = node->NumParentCols, if_rep = false;
  int *col_pool = ColPool4Mem;
  here:
  for (int i = beg; i < end; ++i) {
    memset(states.data(), 0, sizeof(int) * full_cuts.size());
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int current_node = col_pool[j];
      if (!current_node) break;
      for (auto &p : v_cuts[current_node]) {
        states[p.first] += p.second;
        if (states[p.first] >= denominator[p.first]) {
          ++coeffs[p.first][i];
          states[p.first] -= denominator[p.first];
        }
      }
      for (auto &p : v_lostMem[current_node]) states[p] = 0;
    }
  }

  if (!if_rep) {
    if_rep = true;
    beg = end;
    end = NumCol;
    col_pool = ColPool4Pricing;
    goto here;
  }

  for (int i = 0; i < full_cuts.size(); ++i) {
    int numnz;
    int idx = get<2>(full_cuts[i]);
    //take the coefficient of this row
    safe_solver(node->solver.SOLVERgetconstrs(&numnz, solver_ind2, solver_ind, solver_val, idx, 1))
    vector<int> coes(NumCol, 0);
    for (int j = 0; j < numnz; ++j) {
      coes[solver_ind[j]] = int(solver_val[j] + 1e-6);
    }
    if (coes != coeffs[i]) {
      cout << "-------------------------------" << endl;
      cout << "wrong in " << idx << endl;
      int j = 0;
      for (; j < NumCol; ++j) {
        if (coes[j] != coeffs[i][j]) {
          cout << j << " lp= " << coes[j] << " real= " << coeffs[i][j] << endl;
          break;
        }
      }
      cout << "assume j is in pricing!" << endl;
      for (auto k = node->IdxCols[j] + 1;; ++k) {
        int current_node = ColPool4Pricing[k];
        if (!current_node) break;
        cout << current_node << " ";
      }
      cout << endl;
      if (i < node->R1Cs.size()) {
        cout << "R1C" << endl;
      } else {
        cout << "R1C_multi" << endl;
      }
      cout << "cuts: ";
      for (auto &c : get<0>(full_cuts[i]))cout << c << " ";
      cout << endl;
      cout << "plan: " << get<1>(full_cuts[i]) << endl;
      cout << "mem: ";
      for (auto &m : get<3>(full_cuts[i]))cout << m << " ";
      cout << endl;
      exit(0);
    }
  }
}

//void CVRP::addSelectedR1C_N_multiCuts(BBNODE *node,
//                                      vector<tuple<int, set<int>, double, int>> &cut_info_set,
//                                      T &cut_info) {
//  double vio_std;
//  vector<pair<int, double>> list(cut_info_set.size());
//  if constexpr (is_same<T, vector<pair<vector<int>, double>>>::value) {
//    vio_std = 0.5;
//    //calculate the vio and the mem size prior: (vio_std-vio)*(cut+mem_size), less is better
//    //and sort the cut_info_set by the prior
//    for (int i = 0; i < cut_info_set.size(); ++i) {
//      list[i] = {i, (vio_std - get<2>(cut_info_set[i])) *
//          double((get<1>(cut_info_set[i]).size() + cut_info[get<0>(cut_info_set[i])].first.size()))};
//    }
//    sort(list.begin(), list.end(), [](auto &a, auto &b) { return a.second < b.second; });
//    //add the cut from top 150
//    int num_r1c1 = 0, num_r1c3 = 0;
//    for (auto &pr : list) {
//      auto &cut = cut_info_set[pr.first];
//      if (cut_info[get<0>(cut)].first.size() == 1 && ++num_r1c1 > CONFIG::MaxNumR1CPerRoundSecondSelect) {
//        continue;
//      }
//      if (cut_info[get<0>(cut)].first.size() == 3 && ++num_r1c3 > CONFIG::MaxNumR1C3PerRoundSecondSelect) {
//        continue;
//      }
//      cout << "add cut: " << get<1>(cut).size() << " " << get<2>(cut) << " " << pr.second << endl;
//      addR1C(node, cut_info[get<0>(cut)].first, get<1>(cut), get<3>(cut));
//    }
//  } else {
//    vector<double> vio_vec(cut_info_set.size());
//    int cnt = 0;
//    for (auto &i : cut_info_set) {
//      auto &cut = cut_info[get<0>(i)];
//      auto deno = get<1>(map_rank1_multiplier[(int) get<0>(cut).size()][get<1>(cut)]);
//      vio_vec[cnt++] = 1. - 1. / deno;
//    }
//    for (int i = 0; i < cut_info_set.size(); ++i) {
//      list[i] = {i, (vio_vec[i] - get<2>(cut_info_set[i])) *
//          double((get<1>(cut_info_set[i]).size() + get<0>(cut_info[get<0>(cut_info_set[i])]).size()))};
//    }
//    sort(list.begin(), list.end(), [](auto &a, auto &b) { return a.second < b.second; });
//    //add the cut from top 150
//    unordered_map<int, int> num_r1c;
//    for (auto &pr : list) {
//      auto &cut = cut_info_set[pr.first];
//      if (++num_r1c[get<0>(cut_info[get<0>(cut)]).size()] > CONFIG::MaxNumR1CPerRoundSecondSelect) {
//        continue;
//      }
//      cout << "add cut: " << get<1>(cut).size() << " " << get<2>(cut) << endl;
//      addR1C_multi(node,
//                   make_pair(get<0>(cut_info[get<0>(cut)]), get<1>(cut_info[get<0>(cut)])),
//                   get<1>(cut),
//                   get<3>(cut));
//    }
//  }
//}


//void CVRP::addSelectedR1C_N_multiCuts(BBNODE *node,
//                                      vector<tuple<int, set<int>, double, int, int>> &cut_info_set,
//    //cut_idx, mem, vio, lp_index(if max, then new cut!)
//                                      T &cut_info) {
////sort by vio first and then by sie of mem, when sort by vio, the tolerance is 1e-6
//  sort(cut_info_set.begin(), cut_info_set.end(), [](const tuple<int, set<int>, double, int, int> &a,
//                                                    const tuple<int, set<int>, double, int, int> &b) {
//    if (abs(get<2>(a) - get<2>(b)) < 1e-6) {
//      return get<4>(a) < get<4>(b);
//    } else return get<2>(a) > get<2>(b);
//  });
//  //for each mem size, record the max vio
//  unordered_map<int, double> max_vio;
//  for (auto &cut : cut_info_set) {
//    int mem_size = get<4>(cut);
//    if (max_vio.find(mem_size) == max_vio.end()) {
//      max_vio[mem_size] = get<2>(cut);
//    } else {
//      max_vio[mem_size] = max(max_vio[mem_size], get<2>(cut));
//    }
//  }
//  //subtract by 1e-6 for each element
//  transform(cut_info_set.begin(), cut_info_set.end(), cut_info_set.begin(),
//            [](tuple<int, set<int>, double, int, int> &a) {
//              get<2>(a) -= 1e-6;
//              return a;
//            });
//  //we check if the cut's mem is greater than corresponding max vio,
//  //it should be greater than the max vio, when faced with less mem size
//  for (auto &cut : cut_info_set) {
//    int mem_size = get<4>(cut);
//    bool if_should_add = true;
//    for (int i = mem_size - 1; i >= 0; --i) {
//      if (max_vio.find(i) != max_vio.end() && get<2>(cut) < max_vio[i]) {
//        if_should_add = false;
//        break;
//      }
//    }
//    if (if_should_add) {
//      //if cutinfo is vector<pair<vector<int>, double>>, we use addR1C
//      //else we use addR1C_multi
//      if constexpr (is_same<T, vector<pair<vector<int>, double>>>::value) {
//        addR1C(node, cut_info[get<0>(cut)].first, get<1>(cut), get<3>(cut));
//      } else {
//        addR1C_multi(node,
//                     make_pair(get<0>(cut_info[get<0>(cut)]), get<1>(cut_info[get<0>(cut)])),
//                     get<1>(cut),
//                     get<3>(cut));
//      }
//    }
//  }
//}





void CVRP::construct_v_r_vec_aggressive(const vector<vector<int>> &routes,
                                        vector<vector<int>> &v_r_map,
                                        std::vector<std::unordered_map<int, int>> &v_r_map2,
                                        vector<vector<pair<vector<int>, vector<int>>>> &c_n_w_no_c,
                                        unordered_map<int, vector<pair<int, int>>> &c_map) {
  v_r_map.resize(Dim);
  v_r_map2.resize(Dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      v_r_map[i].emplace_back(r);
      ++v_r_map2[i][r];
    }
  }
  int route_size = (int) routes.size();
  c_n_w_no_c.resize(Dim);
  for (int i = 1; i < Dim; ++i) {
    if (v_r_map[i].empty()) continue;
    yzzLong wc = 0;// c within i
    for (auto &j : v_r_map2[i]) {
      for (auto v : routes[j.first]) {
        if (Rank1SepHeurMem4Vertex[i][v]) {
          wc.set(v);
        }
      }
    }
    unordered_map<yzzLong, int> seed_map;
    seed_map.reserve(routes.size());
    unordered_set<int> bad_c_size;
//    unordered_set<int> bad_c_size2;
    vector<pair<vector<int>, vector<int>>> sub_c_n_w_no_c;
    sub_c_n_w_no_c.reserve(16);

    for (int r = 0; r < routes.size(); ++r) {
      if (v_r_map2[i].find(r) != v_r_map2[i].end()) continue;
      yzzLong tmp_c = 0;
      for (auto v : routes[r]) {
        if (wc[v]) {
          tmp_c.set(v);//r no i but c within i
        }
      }
      tmp_c.set(i);
      int c_size = (int) tmp_c.count();// basic c set
      if (c_size < 3 || c_size > CONFIG::MaxHeurInitialCsetSize4RowRank1C) continue;
      if (bad_c_size.find(c_size) != bad_c_size.end()) continue;
//      if (bad_c_size2.find(c_size) != bad_c_size2.end()) continue;
      bad_c_size.insert(map_rank1_multiplier_dominance[c_size].begin(), map_rank1_multiplier_dominance[c_size].end());
//      bad_c_size2.insert(c_size);
      if (seed_map.find(tmp_c) == seed_map.end()) {
        seed_map[tmp_c] = (int) sub_c_n_w_no_c.size();
        yzzLong w_no_c = wc ^ tmp_c;
        sub_c_n_w_no_c.emplace_back(vector<int>(c_size), vector<int>((int) w_no_c.count()));
        int cnt = 0, cnt2 = 0;
        for (int k = 1; k < Dim; ++k) {
          if (wc.test(k)) {
            if (tmp_c.test(k)) sub_c_n_w_no_c.back().first[cnt++] = k;
            else sub_c_n_w_no_c.back().second[cnt2++] = k;
          }
        }
      }
    }
    if (sub_c_n_w_no_c.empty()) continue;
    auto &cset = c_n_w_no_c[i];
    cset.reserve(sub_c_n_w_no_c.size());
    for (auto &j : sub_c_n_w_no_c) {
      if (bad_c_size.find((int) j.first.size()) != bad_c_size.end()) continue;
      cset.emplace_back(j);
      c_map[(int) j.first.size()].emplace_back(i, (int) cset.size() - 1);
    }
//    for (auto &j : cset) {
//      cout << "cset: ";
//      for (auto &k : j.first) {
//        cout << k << " ";
//      }
//      cout << "| ";
//      for (auto &k : j.second) {
//        cout << k << " ";
//      }
//      cout << endl;
//    }
//    cout << "-------------------------------" << endl;
  }
}

void CVRP::construct_operations_aggressive(
    rank1_multi_label &label,
    const vector<vector<int>> &v_r_map,
    const vector<double> &frac_routes,
    int &i,
    std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool
) {
  auto &vio = label.vio;
  auto &new_cij = label.c;
  auto &w_no_cij = label.w_no_c;
  auto &plan_idx = label.plan_idx;

  auto dir = label.search_dir;
  int best_move;
  double best_move_vio;
  int choice_add_j, choice_remove_j;
  pair<int, int> choice_swap_i_j;
  vector<pair<int, double>> move_vio(4);
  for (int j = 1; j < 4; ++j) {
    move_vio[j] = {j, -numeric_limits<double>::max()};
  }
  move_vio[0] = {0, vio};
  double new_vio;
  if (dir == 'a' || dir == 's') {
    new_vio = vio;
    addInRank1HeurSep(plan_idx, new_cij, w_no_cij, new_vio, choice_add_j, v_r_map, frac_routes);
    move_vio[1].second = new_vio;
  }
  if (dir == 'r' || dir == 's') {
    new_vio = vio;
    removeInRank1HeurSep(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
    move_vio[2].second = new_vio;
  }
  if (dir == 'a' || dir == 'r') {
    new_vio = vio;
    swapInRank1HeurSep(plan_idx, new_cij, w_no_cij, new_vio, choice_swap_i_j, v_r_map, frac_routes);
    move_vio[3].second = new_vio;
  }
  //need stable sort
  stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                   const pair<int, double> &b) {
    return a.second > b.second;
  });
  best_move = move_vio[0].first;
  best_move_vio = move_vio[0].second;

  if (best_move == 0) {
    if (best_move_vio > TOLERANCE && new_cij.size() <= CONFIG::MaxRowRank1 && new_cij.size() > 2) {
      generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij,
                                                                    plan_idx,
                                                                    best_move_vio);
    }
    ++i;
  } else if (best_move == 1) {
    new_cij.emplace_back(choice_add_j);
    w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_add_j));
  } else if (best_move == 2) {
    new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
  } else if (best_move == 3) {
//    if (std::find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first) == new_cij.end()) {
//      cout << "choice_swap_i_j.first= " << choice_swap_i_j.first << " choice_swap_i_j.second= "
//           << choice_swap_i_j.second << endl;
//      cout << "new_cij: ";
//      for (auto &j : new_cij) {
//        cout << j << " ";
//      }
//      cout << endl;
//      cout << "w_no_cij: ";
//      for (auto &j : w_no_cij) {
//        cout << j << " ";
//      }
//      cout << endl;
//      throw std::runtime_error("choice_swap_i_j.first not in new_cij");
//    }
    new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
    new_cij.emplace_back(choice_swap_i_j.second);
    w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_swap_i_j.second));
  }
  vio = best_move_vio;
}

void CVRP::construct_seed_aggressive(
    const vector<vector<int>> &v_r_map,
    const vector<double> &frac_routes,
    std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool,
    vector<rank1_multi_label> &rank1_multi_label_pool,
    std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map,
    std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> &c_n_w_no_c,
    int &num_label
) {
  rank1_multi_label_pool.resize(Initial_rank1_multi_label_pool_size);
  num_label = 0;
  for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
    for (int len = 3; len <= CONFIG::MaxHeurInitialCsetSize4RowRank1C; ++len) {
      for (auto &p : c_map[len]) {
        int i = p.first, j = p.second;// i and non-visited idx j
        vector<int> new_cij;
        double vio;
        bool if_succeed;
        //calculate the best vio while using the multiplier plan and the set c[i][j] and the corresponding multiplier
        findCombs4rankCutSet2getBestMultiplier(plan_idx,
                                               c_n_w_no_c[i][j].first,
                                               new_cij,
                                               i,
                                               vio,
                                               v_r_map,
                                               frac_routes,
                                               if_succeed);
        if (!if_succeed) continue;//if negative, then no such multiplier
        //revise the add... function to get the start position as the parameter(the above set is sorted again!)
        vector<pair<int, double>> move_vio(4);
        move_vio[0] = {0, vio};
        double new_vio = vio;
        int choice_add_j;
        int choice_remove_j;
        pair<int, int> choice_swap_i_j;
        addInRank1HeurSep(plan_idx, new_cij, c_n_w_no_c[i][j].second, new_vio, choice_add_j, v_r_map, frac_routes);
        move_vio[1] = {1, new_vio};
        new_vio = vio;
        removeInRank1HeurSep(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
        move_vio[2] = {2, new_vio};
        new_vio = vio;
        swapInRank1HeurSep(plan_idx, new_cij, c_n_w_no_c[i][j].second, new_vio, choice_swap_i_j, v_r_map, frac_routes);
        move_vio[3] = {3, new_vio};
        //need stable sort
        stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                         const pair<int, double> &b) {
          return a.second > b.second;
        });
        int best_move = move_vio[0].first;
        if (best_move == 0) {
          if (move_vio[0].second > TOLERANCE && (new_cij.size() <= CONFIG::MaxRowRank1 && new_cij.size() > 2))
            // 3 cut is allowed to be added!
            generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij, plan_idx, move_vio[0].second);
          continue;
        } else if (best_move == 1) {
          new_cij.emplace_back(choice_add_j);
          auto new_w_no_c = c_n_w_no_c[i][j].second;
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_add_j));
          rank1_multi_label_pool[num_label++] = rank1_multi_label{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  'a'};
        } else if (best_move == 2) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
          rank1_multi_label_pool[num_label++] =
              rank1_multi_label{new_cij, c_n_w_no_c[i][j].second, plan_idx, move_vio[0].second,
                                'r'};
        } else if (best_move == 3) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
          new_cij.emplace_back(choice_swap_i_j.second);
          auto new_w_no_c = c_n_w_no_c[i][j].second;
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_swap_i_j.second));
          rank1_multi_label_pool[num_label++] = rank1_multi_label{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  's'};
        }
        if (num_label == rank1_multi_label_pool.size()) {
          rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
        }
      }
    }
  }
}

void CVRP::findR1C_multi_aggressive(const vector<vector<int>> &routes,
                                    const vector<double> &frac_routes,
                                    std::vector<std::unordered_map<int, int>> &v_r_map,
                                    vector<pair<vector<int>, int>> &cuts) {
  vector<vector<int>> v_r_map_special;
  std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> c_n_w_no_c;
  std::unordered_map<int, std::vector<std::pair<int, int>>> c_map;
  construct_v_r_vec_aggressive(routes, v_r_map_special, v_r_map, c_n_w_no_c, c_map);

  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> generated_rank1_multi_pool;
  vector<rank1_multi_label> rank1_multi_label_pool;
  int num_label;
  construct_seed_aggressive(v_r_map_special,
                            frac_routes,
                            generated_rank1_multi_pool,
                            rank1_multi_label_pool,
                            c_map,
                            c_n_w_no_c,
                            num_label);
  for (int i = 0; i < num_label;) {
    construct_operations_aggressive(rank1_multi_label_pool[i],
                                    v_r_map_special,
                                    frac_routes,
                                    i,
                                    generated_rank1_multi_pool);
  }
  construct_cuts_aggressive(generated_rank1_multi_pool, cuts);
}

void CVRP::checkR1C_mul_look(BBNODE *node, int start) {
  unordered_map<int, int> tmp_map;
  unordered_map<yzzLong, pair<int, int>> tmp_map2;
  bool if_print_all = false;
  int idx = 0;
  for (auto &r1c : node->R1Cs) {
    if (r1c.IdxR1C >= start) {
      if (if_print_all) {
        for (auto c : r1c.InfoR1C) {
          cout << c << " ";
        }
        cout << endl;
      }
      ++tmp_map[0];
    }
    yzzLong tmp_yzz = 0;
    for (auto c : r1c.InfoR1C) {
      tmp_yzz.set(c);
    }
    if (tmp_map2.find(tmp_yzz) != tmp_map2.end() && tmp_map2[tmp_yzz].first == 0) {
      auto &w1 = node->R1Cs[tmp_map2[tmp_yzz].second].InfoR1C;
      for (auto &i : w1) {
        cout << i << " ";
      }
      cout << endl;
      for (auto c : r1c.InfoR1C) {
        cout << c << " ";
      }
      cout << endl;
//      throw std::runtime_error("same r1c");
    }
    tmp_map2[tmp_yzz] = {0, idx};
    ++idx;
  }
  cout << "-----------------" << endl;
  idx = 0;
  for (auto &r1c : node->R1Cs_multi) {
    if (r1c.IdxR1C >= start) {
      if (if_print_all) {
        for (auto c : r1c.InfoR1C.first) {
          cout << c << " ";
        }
        cout << " p: " << r1c.InfoR1C.second << endl;
      }
      ++tmp_map[r1c.InfoR1C.second];
    }
    yzzLong tmp_yzz = 0;
    for (auto c : r1c.InfoR1C.first) {
      tmp_yzz.set(c);
    }
    if (tmp_map2.find(tmp_yzz) != tmp_map2.end() && tmp_map2[tmp_yzz].first == r1c.InfoR1C.second) {
      auto &w1 = node->R1Cs_multi[tmp_map2[tmp_yzz].second].InfoR1C.first;
      for (auto &i : w1) {
        cout << i << " ";
      }
      cout << " p: " << node->R1Cs_multi[tmp_map2[tmp_yzz].second].InfoR1C.second << endl;
      for (auto c : r1c.InfoR1C.first) {
        cout << c << " ";
      }
      cout << " p: " << r1c.InfoR1C.second << endl;
//      throw std::runtime_error("same r1c_multi");
    }
    tmp_map2[tmp_yzz] = {r1c.InfoR1C.second, idx};
    ++idx;
  }

  vector<pair<int, int>> tmp_vec(tmp_map.begin(), tmp_map.end());
  sort(tmp_vec.begin(), tmp_vec.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
    return a.first < b.first;
  });
  for (auto &i : tmp_vec) {
    cout << i.first << " " << i.second << endl;
  }
}

void CVRP::construct_cuts_aggressive(
    std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool,
    vector<pair<vector<int>, int>> &cuts) {

  for (auto &i : generated_rank1_multi_pool) {
    sort(i.second.begin(), i.second.end(),
         [](auto &a, auto &b) {
           return get<2>(a) > get<2>(b);
         });
  }

//  unordered_map<int, double> map_plan_vio;
  for (auto &i : generated_rank1_multi_pool) {
    unordered_set < yzzLong > cut_set;
    unordered_set<int> p_set;
    double vio_std = get<2>(i.second[0]) * CONFIG::CutVioFactor;
//    cout << i.first << " : " << get<2>(i.second[0]) << " size= " << i.second.size() << endl;
    vector<pair<vector<int>, int>> tmp_cuts;
    for (auto &j : i.second) {
      if (get<2>(j) < vio_std) break;
      yzzLong key = 0;
      for (auto c : get<0>(j)) {
        key.set(c);
      }
      if (cut_set.find(key) != cut_set.end() && p_set.find(get<1>(j)) != p_set.end()) continue;
//      if (map_plan_vio[get<1>(j)] < get<2>(j)) {
//        map_plan_vio[get<1>(j)] = get<2>(j);
//      }
      tmp_cuts.emplace_back(get<0>(j), get<1>(j));
      cut_set.insert(key);
      p_set.insert(get<1>(j));
#ifdef debugVio
      if (vio_map.find(key) == vio_map.end()) vio_map[key].resize(7, -numeric_limits<double>::max());
      vio_map[key][get<1>(j)] = get<2>(j);
#endif
    }
    selectCuts(tmp_cuts, cuts, CONFIG::MaxNumR1CPerRound);
  }
//  cout << "-------------------" << endl;
//
//  vector<pair<int, double>> map_vec(map_plan_vio.begin(), map_plan_vio.end());
//  sort(map_vec.begin(), map_vec.end(), [](auto &a, auto &b) {
//    return a.first < b.first;
//  });
//  for (auto &i : map_vec) {
//    cout << i.first << " : " << i.second << endl;
//  }
}

/**
 * May 18....
 */

/**
 * 先按照 seed 找到 集合 为 4 的情况！如果效果好的话，再看看能不能多加一点！
 */

void CVRP::construct_seed_mode3(const vector<vector<int>> &routes,
                                const vector<vector<int>> &v_r_vec,
                                const vector<double> &v_resourceGap,
                                vector<vector<int>> &seed) {
  throw runtime_error("do not even consider! not implemented");
  //this should only be called for once! after that, the rest should be stored!
  for (int i = 1; i < Dim; ++i) {
    yzzLong tmp = 0;
    tmp = Rank1SepHeurMem4Vertex[i];
    for (int j = i + 1; j < Dim; ++j) {
      if (!tmp.test(j)) continue;
      tmp &= Rank1SepHeurMem4Vertex[j];
      for (int k = j + 1; k < Dim; ++k) {
        if (!tmp.test(k)) continue;
        tmp &= Rank1SepHeurMem4Vertex[k];
        for (int l = k + 1; l < Dim; ++l) {
          if (!tmp.test(l)) continue;
          tmp &= Rank1SepHeurMem4Vertex[l];
          for (int p = l + 1; p < Dim; ++p) {
            if (!tmp.test(p)) continue;
            if (!Rank1SepHeurMem4Vertex[p].test(i) ||
                !Rank1SepHeurMem4Vertex[p].test(j) ||
                !Rank1SepHeurMem4Vertex[p].test(k) ||
                !Rank1SepHeurMem4Vertex[p].test(l))
              continue;
            seed.emplace_back(vector<int>{i, j, k, l, p});
          }
        }
      }
    }
  }
  cout << "seed.size()= " << seed.size() << endl;
//  for (auto &i : seed) {
//    for (auto j : i) {
//      cout << j << " ";
//    }
//    cout << endl;
//  }
//  cout << "size= " << seed.size() << endl;
}

void CVRP::NewfindR1C_multi_mode3(const vector<vector<int>> &routes,
                                  const vector<double> &frac_routes,
                                  std::vector<std::unordered_map<int, int>> &v_r_map,
                                  vector<pair<vector<int>, int>> &cuts) {
  throw runtime_error("do not even consider! not implemented");
  vector<vector<int>> v_r_vec;
  vector<double> v_resourceGap;
  construct_v_r_vec(routes, v_r_vec, v_r_map);
#ifdef SYMMETRY_PROHIBIT
  construct_v_resourceGap<false>(routes, v_resourceGap);
#else
  construct_v_resourceGap<true>(routes, v_resourceGap);
#endif
  vector<vector<int>> seed;
  construct_seed_mode3(routes, v_r_vec, v_resourceGap, seed);
  construct_cuts(v_r_map, seed, frac_routes, cuts);
}

void CVRP::fill_mem_first(BBNODE *node,
                          const std::vector<std::vector<int>> &routes,
                          const std::vector<double> &frac_routes,
                          std::vector<std::pair<std::vector<int>, int>> &cuts) {
  if (!cuts.empty()) throw runtime_error("cuts should be empty at first!");
  ResetCutMem.clear();
  /**
   * there's a v_map, could use the v_map outside, but could easily lead to bad results!
   */
  std::vector<std::unordered_map<int, int>> v_r_map(Dim);
  int r_idx = 0;
  for (auto &i : routes) {
    for (auto j : i) {
      ++v_r_map[j][r_idx];
    }
    ++r_idx;
  }
  int num_add_mem = 0;
  vector<double> vis_times(routes.size());//need to be double
  int idx = 0;
  for (auto &r1c : node->R1Cs) {
    memset(&vis_times[0], 0, sizeof(double) * routes.size());
    for (auto &v : r1c.InfoR1C) {
      for (auto &r : v_r_map[v]) {
        vis_times[r.first] += r.second;
      }
    }
    transform(vis_times.begin(), vis_times.end(), frac_routes.begin(), vis_times.begin(), [](auto &a, auto &b) {
      return int(a / 2 + TOLERANCE) * b;
    });
    auto vio = accumulate(vis_times.begin(), vis_times.end(), double(-int((int) r1c.InfoR1C.size() / 2 + TOLERANCE)));
    if (vio > TOLERANCE) {
      cuts.emplace_back(r1c.InfoR1C, 0);
      yzzLong tmp = 0;
      for (auto &v : r1c.InfoR1C) {
        tmp.set(v);
      }
      cut_record[tmp].insert(0);
      ++num_add_mem;
//#ifdef debugVio
//      if (vio_map.find(tmp) == vio_map.end()) vio_map[tmp].resize(7, -numeric_limits<double>::max());
//      if (vio_map[tmp][0] != -numeric_limits<double>::max())continue;
//      vio_map[tmp][0] = vio;
//#endif
      tmp = 0;
      for (auto &v : r1c.Mem) {
        tmp.set(v);
      }
      ResetCutMem.emplace_back(true, idx, tmp);
    }
    ++idx;
  }

  idx = 0;
  for (auto &r1c : node->R1Cs_multi) {
    memset(&vis_times[0], 0, sizeof(double) * routes.size());
    const auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
    const auto &coeff = get<0>(plan);
    auto deno = get<1>(plan);
    auto rhs = get<2>(plan);
    int cnt = 0;
    for (auto &v : r1c.InfoR1C.first) {
      for (auto &r : v_r_map[v]) {
        vis_times[r.first] += r.second * coeff[cnt];
      }
      ++cnt;
    }
    transform(vis_times.begin(), vis_times.end(), frac_routes.begin(), vis_times.begin(), [deno](auto &a, auto &b) {
      return int(a / deno + TOLERANCE) * b;
    });
    auto vio = accumulate(vis_times.begin(), vis_times.end(), -double(rhs));
    if (vio > TOLERANCE) {
      cuts.emplace_back(r1c.InfoR1C);
      yzzLong tmp = 0;
      for (auto &v : r1c.InfoR1C.first) {
        tmp.set(v);
      }
      cut_record[tmp].insert(r1c.InfoR1C.second);
      ++num_add_mem;
//#ifdef debugVio
//      if (vio_map.find(tmp) == vio_map.end()) vio_map[tmp].resize(7, 0);
//      if (vio_map[tmp][0] != -numeric_limits<double>::max())continue;
//      vio_map[tmp][r1c.InfoR1C.second] = vio;
//#endif
      tmp = 0;
      for (auto &v : r1c.Mem) {
        tmp.set(v);
      }
      ResetCutMem.emplace_back(false, idx, tmp);
    }
    ++idx;
  }
  cout << "num_add_mem= " << num_add_mem << endl;
  if (num_add_mem) if_fill_mem = true;
  else if_fill_mem = false;
}

void CVRP::selectCuts(const vector<pair<vector<int>, int>> &tmp_cuts,
                      vector<pair<vector<int>, int>> &cuts,
                      int numCuts) {
  /**
   * this function select valid cuts from all cuts
   */
  numCuts = min(numCuts, (int) tmp_cuts.size());

  for (auto &cut : tmp_cuts) {
    yzzLong tmp = 0;
    for (auto i : cut.first) {
      tmp.set(i);
    }

    if (cut_record[tmp].find(cut.second) != cut_record[tmp].end()) continue;
    cut_record[tmp].insert(cut.second);

    /**
     * make the cuts exists one and only one (very important for avoiding repeated cuts)
     */
    int size = (int) cut.first.size();
    const auto &coeff = get<0>(map_rank1_multiplier[size][cut.second]);
    vector<vector<int>> tmp_cut(coeff[0] + 1);
    for (int i = 0; i < size; ++i) {
      tmp_cut[coeff[i]].emplace_back(cut.first[i]);
    }
    for (int i = 0; i <= coeff[0]; ++i) {
      if (tmp_cut[i].size() > 1) sort(tmp_cut[i].begin(), tmp_cut[i].end());
    }
    vector<int> new_cut(size);
    int index = 0;
    for (int i = coeff[0]; i >= 0; --i) {
      for (auto j : tmp_cut[i]) {
        new_cut[index++] = j;
      }
    }
    cuts.emplace_back(new_cut, cut.second);
    numCuts--;
    if (numCuts == 0) break;
  }
}



