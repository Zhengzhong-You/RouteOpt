
#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;

void CVRP::separateRCCs(BbNode *&node) {
  int numnz;

  CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

  CMGR_CreateCMgr(&MyCutsCMP, 100);
  CMGR_CreateCMgr(&MyOldCutsCMP, 100);
  int numRe=2*num_col;
  vector<int> solver_ind(numRe);
  vector<double> solver_val(numRe);
  int oldNum = num_row;
  while (true) {
    getEdgeInfo(node, false);
    if (node->is_integer) {
      break;
    }

    for (int i = 1; i <= node->num_edges; ++i) {
      if (!node->edge_tail[i]) {
        node->edge_tail[i] = dim;
      } else break;
    }

    CAPSEP_SeparateCapCuts(real_dim, demand, cap, node->num_edges, node->edge_tail.data(),
                           node->edge_head.data(), node->edge_value.data(), MyOldCutsCMP,
                           MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                           &if_int_n_feasible, &max_vio, MyCutsCMP);

    for (int i = 1; i <= node->num_edges; ++i) {
      if (node->edge_tail[i] == dim) {
        node->edge_tail[i] = 0;
      } else break;
    }

    if (!MyCutsCMP->Size) break;
    int cnt = 0;
    for (int i = 0; i < MyCutsCMP->Size; ++i) {
      Rcc rcc;
      auto &tmp_customerInfo = rcc.info_rcc_customer;
      for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
        tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
      }

      if (tmp_customerInfo.size() <= dim / 2) {
        rcc.form_rcc = true;
        rcc.rhs = MyCutsCMP->CPL[i]->RHS;
      } else {
        rcc.form_rcc = false;
        rcc.rhs = real_dim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
        auto &tmp_NoCustomerInfo = rcc.info_rcc_outside_customer;
        vector<bool> tmp(dim, false);
        for (int j : tmp_customerInfo) {
          tmp[j] = true;
        }
        for (int j = 1; j < dim; ++j) {
          if (!tmp[j]) {
            tmp_NoCustomerInfo.emplace_back(j);
          }
        }
      }

      if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) {
        continue;
      }

      ++cnt;
      rcc.idx_rcc = num_row;
      node->rccs.emplace_back(rcc);
	  if (num_col> numRe) {
		  numRe = 2*num_col;
		  solver_ind.resize(numRe);
		  solver_val.resize(numRe);
	  }
      getCoefficientRCC(node, rcc, solver_ind.data(), solver_val.data(), numnz);
	  safe_solver(
          node->solver.addConstraint(numnz, solver_ind.data(), solver_val.data(), SOLVER_LESS_EQUAL, rcc.rhs, nullptr))
      safe_solver(node->solver.updateModel())
      safe_solver(node->solver.getNumRow(&num_row))
      safe_Hyperparameter(checkCSTLimit())
    }

    if (!cnt) break;

    for (int i = 0; i < MyCutsCMP->Size; ++i) {
      CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
    }

#ifdef VERBOSE_MODE
    cout << "Cuts added...  rcc= " << cnt << endl;
#endif
    MyCutsCMP->Size = 0;

    solveLPInLabeling(node);

    if (node->is_terminated) {
      delete node;
      node = nullptr;
      goto QUIT;
    } else {
      eliminateArcs(node);
      enumerateMIP(node);
      if (!node) goto QUIT;
      cleanIndexColForNode(node, node->num_parent_cols, true);
    }
    lb = node->value;
    lb_transformed = ceilTransformedNumberRelated(lb - TOLERANCE);
    if (ceilTransformedNumberRelated(node->value - TOLERANCE) + TOLERANCE >= ub || node->is_integer) {
      node->is_terminated = true;
      cout << TERMINATED_MESSAGE_SEP_RCC;
      break;
    }

#ifdef VERBOSE_MODE
    cout << "gap= " << (double(ub - lb) / ub > TOLERANCE ? double(ub - lb) / ub * 100 : 0) << "%\n";
    cout << SMALL_PHASE_SEPARATION;
#endif
  }

  QUIT:
  if (node) {
    findNonActiveCuts(node);
    convertVertexToR1CsInOneLP(node);
    solveLPInLabeling(node);
  }
  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  CMGR_FreeMemCMgr(&MyCutsCMP);
}

void CVRP::generateRCCs(BbNode *node) {
  int numnz, cnt = 0;

  CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
  CMGR_CreateCMgr(&MyCutsCMP, 100);
  CMGR_CreateCMgr(&MyOldCutsCMP, 100);
  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
#ifdef SOLVER_VRPTW
  if(if_force_keep_rcc) popArcGraph(node);
#endif
  getEdgeInfo(node, false);
  for (int i = 1; i <= node->num_edges; ++i) {
    if (!node->edge_tail[i]) {
      node->edge_tail[i] = dim;
    } else break;
  }

  CAPSEP_SeparateCapCuts(real_dim, demand, cap, node->num_edges, node->edge_tail.data(),
                         node->edge_head.data(), node->edge_value.data(), MyOldCutsCMP,
                         MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                         &if_int_n_feasible, &max_vio, MyCutsCMP);

  for (int i = 1; i <= node->num_edges; ++i) {
    if (node->edge_tail[i] == dim) {
      node->edge_tail[i] = 0;
    } else break;
  }

  if (!MyCutsCMP->Size) goto QUIT;
  cnt = 0;
  for (int i = 0; i < MyCutsCMP->Size; ++i) {
    Rcc rcc;
    auto &tmp_customerInfo = rcc.info_rcc_customer;
    for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
      tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
    }

    if (tmp_customerInfo.size() <= dim / 2) {
      rcc.form_rcc = true;
      rcc.rhs = MyCutsCMP->CPL[i]->RHS;
    } else {
      rcc.form_rcc = false;
      rcc.rhs = real_dim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
      auto &tmp_NoCustomerInfo = rcc.info_rcc_outside_customer;
      vector<bool> tmp(dim, false);
      for (int j : tmp_customerInfo) tmp[j] = true;
      for (int j = 1; j < dim; ++j) {
        if (!tmp[j]) {
          tmp_NoCustomerInfo.emplace_back(j);
        }
      }
    }

    if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) {
      continue;
    }

    ++cnt;

    rcc.idx_rcc = num_row;
#ifdef SOLVER_VRPTW
	if (if_force_keep_rcc) rcc.if_keep=true;
#endif
    node->rccs.emplace_back(rcc);

    getCoefficientRCC(node, rcc, solver_ind.data(), solver_val.data(), numnz);
    safe_solver(node->solver.addConstraint(numnz, solver_ind.data(), solver_val.data(), SOLVER_LESS_EQUAL, rcc.rhs, nullptr))
    safe_solver(node->solver.updateModel())
    safe_solver(node->solver.getNumRow(&num_row))
    safe_Hyperparameter(checkCSTLimit())
  }
  QUIT:
  for (int i = 0; i < MyCutsCMP->Size; ++i) {
    CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
  }
  MyCutsCMP->Size = 0;

  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  CMGR_FreeMemCMgr(&MyCutsCMP);
}


void CVRP::separateHybridCuts(BbNode *&node) {
  int cnt_tail_off = 0;
  double standard;
  int count = 0;
  int RCCnum, R1Cnum;
  bool if_opt = false;
  std::chrono::high_resolution_clock::time_point t1, t2;
  double duration;

  if (if_in_enu_state) goto HYBRID;
  if (!node->index) {
#ifdef SOLVER_VRPTW
    augmentNonGuillotineRound(node);
#endif

    cout << "Begin Rcc separation...\n";

    separateRCCs(node);

    cout << MID_PHASE_SEPARATION;
    cout << "Rcc tail-off!\nSwitch on Hybrid...\n";

    if (!node || node->is_terminated) {
      if_opt = true;
      goto QUIT;
    }
    safe_Hyperparameter(rollback)

    setTailOffStandardAndRollBackStandard();
  } else {
#ifdef VERBOSE_MODE
    cout << "Begin Hybrid separation...\n";
#endif
  }

  HYBRID:
  while (true) {

    bool if_soft_limit = false;
    ++count;
    cout << BIG_PHASE_SEPARATION;
    cout << "Hybrid separation round " << count << endl;

    double prior_nodeVal = node->value;

    int oldNum = num_row;

    auto beg = chrono::high_resolution_clock::now();

    if (if_in_enu_state) {
      generateRCCsInEnum(node);
    } else {
      generateRCCs(node);
    }

    newGenerateAllR1Cs(node, 2);

    int cuts_sum = num_row - oldNum;
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
    safe_solver(node->solver.reoptimize())
    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration<double>(t2 - t1).count();
    cout << "re-optimize time: " << duration << endl;

    deleteNonActiveCutsSafely(node, oldNum);

    printCutsInformation(node);

    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration<double>(t2 - t1).count();
    cout << "Hybrid separation time: " << duration << endl;

    SOLVE_CG:
    if (if_in_enu_state) {
      solveLPByInspection(node, false, false, true);
      if (!node->index) updateLowerBound(node->value);
      if (node->is_terminated)goto QUIT;
      if (node->size_enumeration_col_pool + num_col <= max_num_route4_mip) goto QUIT;
    } else {
      convertVertexToR1CsInOneLP(node);
      solveLPInLabeling(node);
      if (!node->index) updateLowerBound(node->value);
      if (rollback == 1) {
        cout << BIG_PHASE_SEPARATION;
        rollbackEasyWay(node, oldNum);
        cout << BIG_PHASE_SEPARATION;
        goto QUIT;
      } else if (rollback == 3) {
        if_soft_limit = true;
      }

      if (node->is_terminated) {
        delete node;
        node = nullptr;
        goto QUIT;
      }
      eliminateArcs(node);
      enumerateMIP(node);
      if (!node) goto QUIT;
      cleanIndexColForNode(node, node->num_parent_cols, true);
      findNonActiveCuts(node);
    }
    cout << "after deleting those inactivate cuts! We left with cols: " << num_col << " rows: " << num_row << endl;

    standard = calculateGapImprovement(node->value, prior_nodeVal);
    cout << "local gap= " << (double(ub - node->value) / ub > TOLERANCE ? double(ub - node->value) / ub * 100 : 0) << endl;
    cout << "gap improved= " << standard << endl;
    cout << SMALL_PHASE_SEPARATION;
    if (ub <= lb_transformed) {
      if_opt = true;
      goto QUIT;
    }
    if (if_soft_limit) {
      cout << "labels reach the soft limit!" << endl;
      goto QUIT;
    }
    if (standard < TOLERANCE) goto QUIT;
    else if (standard < Config::CutsTailOff[0]) {
      ++cnt_tail_off;
      cout << "cnt_tail_off= " << cnt_tail_off << endl;
      if (cnt_tail_off >= Config::CutsTolerance)goto QUIT;
    }
  }

  QUIT:
  if (!if_in_enu_state) {
    if (node && !if_opt) {
      convertVertexToR1CsInOneLP(node);
      cout << "we clean the col pool!" << endl;
      cleanIndexColForNode(node, node->num_parent_cols);
      cout << "Hybrid tali-off! " << "Stop cut-separation!" << " Last column reduction... ncol= " << num_col
           << endl;
      cout << BIG_PHASE_SEPARATION;
    } else {
      cout << "Terminate the node!\n" << "Stop cut-separation!" << endl;
      cout << SMALL_PHASE_SEPARATION;
    }
  } else {
    cout << "Hybrid tali-off! " << "Stop cut-separation!" << " Last column reduction... ncol= " << num_col
         << endl;
    cout << BIG_PHASE_SEPARATION;
  }
  if (node) {
    recordOptimalColumn(node);//very essential
    if (!node->index) updateLowerBound(node->value);
  }
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
  int MaxLengthEleRoute= aver_route_length*2;
  yzzLong tmp_long, PI_ai, PI_aj, PI_ai_aj;
  vector<int> aux(3);

  vector<vector<int>> vertex2non_ele_routes(dim);
  for (int i = 0; i < dim; ++i) vertex2non_ele_routes.reserve(num_non_ele_routes);

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

  for (int i = 0; i < num_ele_routes; ++i) {
    auto &seq = ele_routes[i];
    cnt = (int) ele_routes[i].size();
    frac = frac_ele_routes[i];
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
    for (int j = 0; j < cnt; ++j) {
      for (int k = j + 1; k < cnt; ++k) {
        for (int l = k + 1; l < cnt; ++l) {
          tmp_long = 0;
          tmp_long.set(seq[j]);
          tmp_long.set(seq[k]);
          tmp_long.set(seq[l]);
          iter3 = three_combinations.find(tmp_long);
          if (iter3 == three_combinations.end()) {
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
    for (int j = 1; j < ai; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      tmp_long = PI_aj;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
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
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      tmp_long = PI_aj;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
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
    for (int j = aj + 1; j < dim; ++j) {
      tmp_long = PI_ai_aj;
      tmp_long.set(j);
      all_iter = all_combinations.find(tmp_long);
      if (all_iter != all_combinations.end()) continue;
      frac = get<2>(two_comb.second);
      tmp_long = PI_ai;
      tmp_long.set(j);
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
      tmp_long = PI_aj;
      tmp_long.set(j) = 1;
      iter2 = two_combinations.find(tmp_long);
      if (iter2 != two_combinations.end()) frac += get<2>(iter2->second);
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
  int limit = min((int) violated_combs.size(), Config::MaxNumR1C3PerRound);
  transform(violated_combs.begin(), violated_combs.begin() + limit, back_inserter(cut_info),
            [](const std::tuple<int, int, int, double> &a) {
              return make_pair(vector<int>{get<0>(a), get<1>(a), get<2>(a)}, get<3>(a) - 1);
            });
}

void CVRP::generateR1C1s(const vector<vector<int>> &non_ele_routes,
                         const vector<double> &frac_none_ele_routes,
                         vector<pair<vector<int>, double>> &cut_info) const {

  yzzLong tmp_long;
  set<int> re_visited_vertices;
  unordered_map<int, set<int>> vertex2non_ele_routes;
  vertex2non_ele_routes.reserve(dim);

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

  int limit = min((int) vio_mem_vec_re_vertices.size(), Config::MaxNumR1CPerRound);
  transform(vio_mem_vec_re_vertices.begin(), vio_mem_vec_re_vertices.begin() + limit, back_inserter(cut_info),
            [](const pair<int, double> &a) {
              return make_pair(vector<int>{a.first}, a.second);
            });
}

void CVRP::findMemoryForR1CsMode3ByEnumerationAndMIP(yzzLong v_comb,
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
          for (int v = 0; v < dim; ++v) {
            if (bs[v]) {
              mem_long.set(v);
            }
          }
        }
      }
    }
  }

  int cnt = 1;
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

  if (mem_long.count() > rank1_mem_size_limit) {
    if_suc = false;
    return;
  }

  for (int i = 1; i < dim; ++i) {
    if (mem_long[i]) mem.emplace(i);
  }

  if (cnt == 1) {
    return;
  } else if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
    vector<int> tmp;
    set<int> new_mem;
    int record_min = MAX_INT;
    combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
    mem = new_mem;
  } else {
    getMemoryByMIP(vec_data, vec_segment_route, mem, if_suc);
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
  if (recorded_combinations4_r1c.find(n) != recorded_combinations4_r1c.end()) {
    data = recorded_combinations4_r1c[n];
    return;
  }
  vector<int> arr(n);
  iota(arr.begin(), arr.end(), 0);
  int r = int((n + 1) / 2);
  vector<int> tmp(r);
  combinationUtil(arr, tmp, data, 0, n - 1, 0, r);
  recorded_combinations4_r1c[n] = data;
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

void CVRP::getMemoryByMIP(const vector<vector<vector<int>>> &array,
                       const vector<vector<set<int>>> &vec_segment,
                       set<int> &mem, bool &if_suc) {
  unordered_map<int, int> new_mem_map;
  unordered_map<int, int> re_new_mem_map;
  new_mem_map.reserve(dim);

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
  int local_num_row = half_num_row + last_half;
  vector<char> sense(local_num_row, SOLVER_LESS_EQUAL);
  fill(sense.begin() + half_num_row, sense.end(), SOLVER_EQUAL);
  vector<double> rhs(local_num_row, 0);
  fill(rhs.begin() + half_num_row, rhs.end(), 1);

  Solver local_solver{};
  local_solver.getEnv(&solver);//need load environment
  cout << "we build model= getMemByMIP_.lp to get mem!" << endl;
  safe_solver(local_solver.newModel("getMemByMIP_.lp", 0, nullptr, nullptr, nullptr, nullptr, nullptr))
  safe_solver(local_solver.addConstraints(local_num_row, 0, nullptr, nullptr, nullptr, sense.data(), rhs.data(), nullptr))

  int last_num_idx = half_num_row;
  vector<size_t> solver_beg;
  vector<int> solver_ind;
  vector<double> solver_val, solver_obj;
  int numRe=10000;
  solver_beg.reserve(numRe);
  solver_ind.reserve(numRe);
  solver_val.reserve(numRe);
  solver_obj.reserve(numRe);
  for (int i = 0; i < array.size(); ++i, ++last_num_idx) {
    for (auto &p : array[i]) {
	  solver_beg.emplace_back((int) solver_ind.size());
      yzzLong tmp = 0;
      for (auto n : p) {
        for (auto j : vec_segment[i][n]) {
          if (tmp[j]) continue;
          tmp.set(j);
		  solver_ind.emplace_back(new_mem_map[j]);
        }
      }
	  solver_ind.emplace_back(last_num_idx);
    }
  }
  int old_ccnt=(int) solver_beg.size();
  size_t old_nzcnt= solver_ind.size();
  solver_obj.assign(old_ccnt, 0);
  solver_val.assign(old_nzcnt, 1);
  solver_obj.resize(old_ccnt+ half_num_row, 1);
  solver_beg.resize(old_ccnt+ half_num_row+1);
  iota(solver_beg.begin()+old_ccnt, solver_beg.end(), old_nzcnt);
  old_ccnt += half_num_row;
  solver_ind.resize(old_nzcnt+ half_num_row);
  iota(solver_ind.begin()+ old_nzcnt, solver_ind.end(), 0);
  old_nzcnt+= half_num_row;
  solver_val.resize(old_nzcnt, -num_r_p);
  vector<char> xtype(old_ccnt, SOLVER_BINARY);
  safe_solver(local_solver.XaddVars(old_ccnt,
                                    old_nzcnt,
                                    solver_beg.data(),
                                    solver_ind.data(),
                                    solver_val.data(),
                                    solver_obj.data(),
                                    nullptr,
                                    nullptr,
                                    xtype.data(),
                                    nullptr))
  safe_solver(local_solver.setEnvTimeLimit(TIME_LIMIT_FOR_MIP_FIND_MEM))
  safe_solver(local_solver.optimize())
  safe_solver(local_solver.setEnvTimeLimit(MAX_TIME_LIMIT_FOR_MIP))
  int status, left= old_ccnt-num_r_p;
  safe_solver(local_solver.getStatus(&status))
  vector<double> X(left);
  if (status == SOLVER_TIME_LIMIT) {
    cout << "time limit for getMemoryByMIP" << endl;
    if_suc = false;
    goto here;
  }
  safe_solver(local_solver.getX(num_r_p, left, X.data()))
  for (int i = 0; i < left; ++i) {
    if (X[i] > 0.5) {
      mem.emplace(re_new_mem_map[i]);
    }
  }
  here:
  safe_solver(local_solver.freeModel())
}

struct other_ {
  int beg{};
  int tor{};
  vector<int> left_c{};
  vector<int> mem_c{};
  other_(int beg, int tor, vector<int> left_c, vector<int> mem_c) :
      beg(beg), tor(tor), left_c(std::move(left_c)), mem_c(std::move(mem_c)) {}
  other_() = default;
};

void CVRP::findR1CMulti(const vector<vector<int>> &routes,
                         const vector<double> &frac_routes,
                         const vector<vector<int>> &vbasis,
                         vector<pair<vector<int>, double>> &cut_info,
                         vector<tuple<vector<int>, int, double>> &multi_cut_info) {
  vector<vector<int>> v_r_map(dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      v_r_map[i].emplace_back(r);
    }
  }

  int route_size = (int) routes.size();
  vector<vector<vector<int>>> c(dim, vector<vector<int>>(route_size));
  vector<vector<vector<int>>> w_no_c(dim, vector<vector<int>>(route_size));
  auto if_vis = new bool[route_size];
  unordered_map<int, vector<pair<int, int>>> c_map;//size, c_index
  for (int i = 1; i < dim; ++i) {
    memset(if_vis, 0, sizeof(bool) * route_size);
    yzzLong wc = 0;// c within i
    for (auto &j : v_r_map[i]) {
      if_vis[j] = true;
      for (auto v : routes[j]) {
        if (rank1_sep_heur_mem4_vertex[i][v]) wc.set(v);
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
      if (c_size < 3 || c_size > Config::MaxHeurInitialCsetSize4RowRank1C) continue;
      c_map[c_size].emplace_back(i, j);
      c[i][j].resize(c_size);
      w_no_c[i][j].reserve(dim);
      c_size = 0;
      for (int k = 1; k < dim; ++k) {
        if (tmp_c[k]) c[i][j][c_size++] = k;
        else if (wc[k]) w_no_c[i][j].emplace_back(k);//candidates set
      }
    }
  }
  delete[] if_vis;

  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> generated_rank1_multi_pool;
  vector<Rank1MultiLabel> rank1_multi_label_pool(INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE);
  unordered_map<int, int> num_operations;
  int label_cnt = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  double total_time = 0;
  for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
    for (int len = 4; len <= Config::MaxRowRank1; ++len) {
      for (auto &p : c_map[len]) {
        int i = p.first, j = p.second;// i and non-visited routes idx: j
        vector<int> new_cij;
        double vio;
        bool if_succeed;
        auto beg10 = std::chrono::high_resolution_clock::now();
        findCombinationsForRankCutSetToGetBestMultiplier(plan_idx,
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
        vector<pair<int, double>> move_vio(4);
        move_vio[0] = {0, vio};
        double new_vio = vio;
        int choice_add_j;
        int choice_remove_j;
        pair<int, int> choice_swap_i_j;
        addInRank1HeuristicSeparate(plan_idx, new_cij, w_no_c[i][j], new_vio, choice_add_j, v_r_map, frac_routes);
        move_vio[1] = {1, new_vio};
        new_vio = vio;
        removeInRank1HeuristicSeparate(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
        move_vio[2] = {2, new_vio};
        new_vio = vio;
        swapInRank1HeuristicSeparate(plan_idx, new_cij, w_no_c[i][j], new_vio, choice_swap_i_j, v_r_map, frac_routes);
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
          rank1_multi_label_pool[label_cnt++] = Rank1MultiLabel{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  'a'};
          if (label_cnt == rank1_multi_label_pool.size()) {
            rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
          }
        } else if (best_move == 2) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
          rank1_multi_label_pool[label_cnt++] = Rank1MultiLabel{new_cij, w_no_c[i][j], plan_idx, move_vio[0].second,
                                                                  'r'};
          if (label_cnt == rank1_multi_label_pool.size()) {
            rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
          }
        } else if (best_move == 3) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
          new_cij.emplace_back(choice_swap_i_j.second);
          auto new_w_no_c = w_no_c[i][j];
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_swap_i_j.second));
          rank1_multi_label_pool[label_cnt++] = Rank1MultiLabel{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
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
      addInRank1HeuristicSeparate(plan_idx, new_cij, w_no_cij, new_vio, choice_add_j, v_r_map, frac_routes);
      move_vio[1] = {1, new_vio};
      new_vio = vio;
      swapInRank1HeuristicSeparate(plan_idx, new_cij, w_no_cij, new_vio, choice_swap_i_j, v_r_map, frac_routes);
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
      removeInRank1HeuristicSeparate(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
      move_vio[1] = {2, new_vio};
      new_vio = vio;
      swapInRank1HeuristicSeparate(plan_idx, new_cij, w_no_cij, new_vio, choice_swap_i_j, v_r_map, frac_routes);
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
      addInRank1HeuristicSeparate(plan_idx, new_cij, w_no_cij, new_vio, choice_add_j, v_r_map, frac_routes);
      move_vio[1] = {1, new_vio};
      new_vio = vio;
      removeInRank1HeuristicSeparate(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
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
  for (auto i = generated_rank1_multi_pool.begin(); i != generated_rank1_multi_pool.end();) {
    if (i->first > Config::MaxRowRank1) {
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
    double vio_std = get<2>(i.second[0]) * Config::CutVioFactor;
    i.second.erase(remove_if(i.second.begin(), i.second.end(),
                             [vio_std](auto &a) {
                               return get<2>(a) < vio_std;
                             }), i.second.end());
  }
  v_r_map.clear();
  v_r_map.resize(dim);
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
        auto &plan_idx = get<1>(generated_rank1_multi_pool[i.first][j]);
        if (!plan_idx) {
          cut_info.emplace_back(get<0>(generated_rank1_multi_pool[i.first][j]),
                                get<2>(generated_rank1_multi_pool[i.first][j]));
        } else {
          multi_cut_info.emplace_back(get<0>(generated_rank1_multi_pool[i.first][j]),
                                      plan_idx,
                                      get<2>(generated_rank1_multi_pool[i.first][j]));
        }
        ++len;
        if (len == Config::MaxNumR1CPerRound) break;
      }
    }
    delete[] if_del;
  }
}


void CVRP::findMemoryForRank1Multi(const vector<vector<int>> &routes,
                               const std::vector<std::unordered_map<int, int>> &v_r_map,
                               const pair<vector<int>, int> &cut_pair,
                               set<int> &mem,
                               bool &if_suc) {
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
    findPlanForRank1Multi(vis, denominator, mem_long, segment_route, data);
    if (!data.empty()) {
      vec_data.emplace_back(data);
      vec_segment_route.emplace_back(segment_route);
    }
  }

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

  if (mem_long.count() > rank1_mem_size_limit) {
    if_suc = false;
    return;
  }

  for (int i = 1; i < dim; ++i) {
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
    vector<int> tmp;
    set<int> new_mem;
    int record_min = MAX_INT;
    combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
    mem = new_mem;
  } else {
    getMemoryByMIP(vec_data, vec_segment_route, mem, if_suc);
  }
}

void CVRP::constructMemory(BbNode *node, const vector<vector<int>> &routes,
                         const std::vector<std::unordered_map<int, int>> &v_r_map,
                         const vector<pair<vector<int>, int>> &cuts,
                         vector<tuple<vector<int>, int, int, set<int>>> &full_cuts) {//cut self, plan idx, cut idx, mem
  unordered_map<pair<vector<int>, int>, std::pair<std::set<int>, int>, Rank1MultiPairHasher>
      R1C_multi_Pool;//mem & cut_index
  for (int i = 0; i < node->r1cs.size(); ++i) {
    R1C_multi_Pool[{node->r1cs[i].info_r1c, 0}] =
        {set < int > (node->r1cs[i].mem.begin(), node->r1cs[i].mem.end()), i};
  }
  for (int i = 0; i < node->r1cs_multi.size(); ++i) {
    R1C_multi_Pool[node->r1cs_multi[i].info_r1c] =
        {set < int > (node->r1cs_multi[i].mem.begin(), node->r1cs_multi[i].mem.end()), i};
  }
  full_cuts.reserve(cuts.size());
  int num_r1c = MAX_NUM_R1CS - node->r1cs.size() - 1;
  int num_mul_r1c = MAX_NUM_R1C_MULTI - node->r1cs_multi.size() - 1;
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
      findMemoryForR1CsMode3ByEnumerationAndMIP(v_comb, mem, routes, if_suc);
    } else {
      if (--num_mul_r1c < 0) continue;
      findMemoryForRank1Multi(routes, v_r_map, cut, mem, if_suc);
    }
    if (if_suc) {
      full_cuts.emplace_back(cut.first, cut.second, index, mem);
    }
  }
}

void CVRP::newGenerateAllR1Cs(BbNode *node, int aggressive_level) {

  cut_record.clear();
  vector<pair<vector<int>, int>> cuts;
  vector<vector<int>> routes;
  vector<double> frac_routes;
  vector<int> route;
  route.reserve(dim);
  vector<vector<int>> ele_routes;
  vector<vector<int>> non_ele_routes;
  vector<double> frac_ele_routes;
  vector<double> frac_non_ele_routes;

  if (if_in_enu_state) {
    for (auto &i : node->index_for_lp_solutions_in_column_pool) {
      route.clear();
      for (auto j = i.first + 1;; ++j) {
        int curr_node = col_pool4_pricing[j];
        if (!curr_node)break;
        route.emplace_back(curr_node);
      }
      routes.emplace_back(route);
      frac_routes.emplace_back(i.second);
    }
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
    for (int i = 0; i < node->num_parent_cols_in_lp_solutions; ++i) {
      route.clear();
      if (abs(node->index_for_lp_solutions_in_column_pool[i].second - 1) < TOLERANCE)continue;
      for (auto j = node->index_for_lp_solutions_in_column_pool[i].first + 1;; ++j) {
        int curr_node = col_pool4_mem[j];
        if (!curr_node)break;
        route.emplace_back(curr_node);
      }
      routes.emplace_back(route);
      frac_routes.emplace_back(node->index_for_lp_solutions_in_column_pool[i].second);
    }
    for (int i = node->num_parent_cols_in_lp_solutions; i < node->index_for_lp_solutions_in_column_pool.size(); ++i) {
      route.clear();
      if (abs(node->index_for_lp_solutions_in_column_pool[i].second - 1) < TOLERANCE)continue;
      for (auto j = node->index_for_lp_solutions_in_column_pool[i].first + 1;; ++j) {
        int curr_node = col_pool4_pricing[j];
        if (!curr_node)break;
        route.emplace_back(curr_node);
      }
      routes.emplace_back(route);
      frac_routes.emplace_back(node->index_for_lp_solutions_in_column_pool[i].second);
    }
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

  if (!if_in_enu_state) {
	fillMemoryFirst(node, routes, frac_routes, cuts);
  }

  newGenerateR1C1s(non_ele_routes, frac_non_ele_routes, cuts);

  newGenerateR1C3s(routes, frac_routes, cuts);// unnecessarily enumeration!//can be supplied by findR1CMultiAggressive

  auto beg = std::chrono::high_resolution_clock::now();
  std::vector<std::unordered_map<int, int>> v_r_map;

  switch (aggressive_level) {
    case 0:newFindR1CMulti(routes, frac_routes, v_r_map, cuts);
      break;
    case 1:findR1CMultiAggressive(routes, frac_routes, v_r_map, cuts);
      break;
    case 2:searchCrazy(routes, frac_routes, v_r_map, cuts);
      break;
    default: throw std::runtime_error("aggressive_level error");
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double>(end - beg).count();

  cout << "findR1CMulti time: " << duration << endl;

  cout << "cuts.size() = " << cuts.size() << endl;

  if (if_in_enu_state) {
    addR1CAtOnceInEnum(node, cuts);
  } else {
    vector<tuple<vector<int>, int, int, set<int>>> full_cuts;// cut, plan, index, mem
    constructMemory(node, routes, v_r_map, cuts, full_cuts);
    addR1CAtOnce(node, full_cuts);
  }

  cout << "now num_row=  " << num_row << endl;
  safe_Hyperparameter(num_row > CST_LIMIT)
}

void CVRP::findPlanForRank1Multi(const vector<int> &vis, int denominator, yzzLong &mem,
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

  for (int i = 1; i < dim; ++i) {
    if (mem[i]) {
      for (auto &j : segment) {
        j.erase(i);
      }
    }
  }
  vector<set<int>> num_vis(dim);
  for (int i = 0; i < other2.size(); ++i) {
    for (auto j : other2[i]) {
      for (auto k : segment[j]) {
        num_vis[k].emplace(i);
      }
    }
  }
  for (int i = 1; i < dim; ++i) {
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

void CVRP::findCombinationsForRankCutSetToGetBestMultiplier(int plan_idx,
                                                  const vector<int> &cset,
                                                  vector<int> &new_cset,
                                                  int c,
                                                  double &new_vio,
                                                  const vector<vector<int>> &v_r_map,
                                                  const vector<double> &frac_routes,
                                                  bool &if_succeed) {
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
  generatePureMap(map, sorted_num_info);
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

void CVRP::addInRank1HeuristicSeparate(int plan_idx,
                             const vector<int> &c,
                             const vector<int> &w_no_c,
                             double &new_vio,
                             int &choice_j,
                             const vector<vector<int>> &v_r_map,
                             const vector<double> &frac_routes) {
  int new_c_size = (int) c.size() + 1;
  const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
  if (new_c_size > Config::MaxRowRank1 || !get<1>(plan)) {
    new_vio = -dim;
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
    new_vio = -dim;
  }
}

void CVRP::removeInRank1HeuristicSeparate(int plan_idx,
                                const vector<int> &c,
                                double &new_vio,
                                int &choice_j,
                                const vector<vector<int>> &v_r_map,
                                const vector<double> &frac_routes) {
  int new_c_size = (int) c.size() - 1;
  const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
  if (new_c_size < 3 || !get<1>(plan)) {
    new_vio = -dim;
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
    new_vio = -dim;
  }
}

void CVRP::swapInRank1HeuristicSeparate(int plan_idx,
                              const vector<int> &c,
                              const vector<int> &w_no_c,
                              double &new_vio,
                              pair<int, int> &choice_i_j,
                              const vector<vector<int>> &v_r_map,
                              const vector<double> &frac_routes) {
  const auto &plan = map_rank1_multiplier[(int) c.size()][plan_idx];
  int new_c_size = (int) c.size();
  if ((new_c_size < 3 || new_c_size > Config::MaxRowRank1) || !get<1>(plan)) {
    new_vio = -dim;
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
    new_vio = -dim;
  }
}

void CVRP::combinationUtilAddOne(const std::vector<int> &arr,
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
    combinationUtilAddOne(arr, tmp, data, i + 1,
                           end, index + 1, r);
  }
}

void CVRP::combineAll(int n, int r, vector<vector<int>> &data) {
  vector<int> arr(n);
  iota(arr.begin(), arr.end(), 0);
  vector<int> tmp(r);
  combinationUtilAddOne(arr, tmp, data, 0, (int) arr.size() - 1, 0, r);
}

void CVRP::generatePureMap(std::vector<std::vector<std::vector<int>>> &map, const std::vector<int> &sorted_info) {
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
  combineAll(n, sorted_info[mul_idx], map[mul_idx]);
  n -= sorted_info[mul_idx];
  ++mul_idx;
  if (mul_idx != sorted_info.size()) {
    goto GET_LEFT;
  }
  QUIT:
  if (!if_find) pureMap[sorted_info] = map;
}

void CVRP::newGenerateR1C1s(const vector<vector<int>> &non_ele_routes,
                            const vector<double> &frac_none_ele_routes,
                            std::vector<std::pair<std::vector<int>, int>> &cuts) {
  if (non_ele_routes.empty())
    return;
  yzzLong tmp_long;
  set<int> re_visited_vertices;
  unordered_map<int, set<int>> vertex2non_ele_routes;
  vertex2non_ele_routes.reserve(dim);

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
  chooseCuts(tmp_cut, cuts, Config::MaxNumR1CPerRound);
}

void CVRP::newGenerateR1C3s(const vector<vector<int>> &routes,
                            const vector<double> &frac_routes,
                            std::vector<std::pair<std::vector<int>, int>> &cuts) {
  /**
   * we only search for cuts that could potentially make the cuts have a small memory
   * for the routes, we can construct a residual graph
   * the cuts set {i,j,k}. for i, (i,j) or (i,k) must exist in the graph
   */

  int allNumRoutes = (int) routes.size();
  vector<unordered_map<int, int>> visited_times_by_routes(dim);
  vector<unordered_set<int>> i_connections(dim);//j,k that connects i
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
  vector<vector<int>> i_connections_vec(dim);
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
  for (int i = 1; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      v_tmp.clear();
      if (i_connections[i].find(j) != i_connections[i].end()) {
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
        if (vio > CUT_VIO_TOLERANCE) {
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
  chooseCuts(tmp_cut, cuts, Config::MaxNumR1C3PerRound);
}

void CVRP::constructVRVec(const vector<vector<int>> &routes,
                             vector<vector<int>> &v_r_vec,
                             vector<unordered_map<int, int>> &v_r_map) const {
  int NumRoutes = (int) routes.size();
  v_r_vec.resize(dim, vector<int>(NumRoutes, 0));
  v_r_map.resize(dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      ++v_r_map[i][r];
    }
  }
  for (int i = 1; i < dim; ++i) {
    for (auto &ele : v_r_map[i]) {
      v_r_vec[i][ele.first] = ele.second;
    }
  }
}

template<bool if_symmetry>
void CVRP::construct_v_resourceGap(const vector<vector<int>> &routes, vector<double> &v_resourceGap) {
  vector<pair<double, int>> supp_v_resourceGap(dim);
  for (auto &r : routes) {
    int past_node = 0;
    int cur_node;
    double res = 0;
    for (int it : r) {
      cur_node = it;
      increaseMainResourceConsumption(res, res, past_node, cur_node);
      if (res > meet_point_resource_in_bi_dir) break;
      supp_v_resourceGap[cur_node].first += meet_point_resource_in_bi_dir - res;
      ++supp_v_resourceGap[cur_node].second;
      past_node = cur_node;
    }
    past_node = 0;
    res = if_symmetry ? 0 : max_main_resource;
    for (auto it = r.rbegin(); it != r.rend(); ++it) {
      cur_node = *it;
      if_symmetry ? increaseMainResourceConsumption(res, res, past_node, cur_node)
                  : decreaseMainResourceConsumption(res, res, past_node, cur_node);
      if ((if_symmetry ? meet_point_resource_in_bi_dir - res : res - meet_point_resource_in_bi_dir) < 0) break;
      supp_v_resourceGap[cur_node].first += if_symmetry ? meet_point_resource_in_bi_dir - res
                                                        : res - meet_point_resource_in_bi_dir;
      ++supp_v_resourceGap[cur_node].second;
      past_node = cur_node;
    }
  }
  v_resourceGap.resize(dim);
  transform(supp_v_resourceGap.begin(), supp_v_resourceGap.end(), v_resourceGap.begin(),
            [](const pair<double, int> &p) -> double { return p.first / p.second; });
}

void CVRP::constructSeed(const vector<vector<int>> &routes,
                          const vector<vector<int>> &v_r_vec,
                          const vector<double> &v_resourceGap,
                          vector<vector<int>> &seed) {
  unordered_set < yzzLong > seed_pool;
  seed_pool.reserve(1024);
  int route_size = (int) routes.size();
  vector<vector<vector<int>>> c(dim, vector<vector<int>>(route_size));
  unordered_map<int, vector<pair<int, int>>> c_map;//size, c_index
  vector<int> no_route;
  no_route.reserve(route_size);
  for (int i = 1; i < dim; ++i) {
    yzzLong wc = 0;// c within i
    no_route.clear();
    for (int r = 0; r < route_size; ++r) {
      if (!v_r_vec[i][r]) no_route.emplace_back(r);
      else {
        for (auto v : routes[r]) {
          if (rank1_sep_heur_mem4_vertex[i][v]) wc.set(v);
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
      if (tmp_seed.size() < 3 || tmp_seed.size() > Config::MaxRowRank1)
        continue;// in this case, we use do not allow extension!
      if (seed_pool.find(yzz_seed) != seed_pool.end()) continue;
      sort(tmp_seed.begin(), tmp_seed.end(), [](const pair<int, double> &a, const pair<int, double> &b) -> bool {
        return a.second < b.second;
      });
      seed_pool.emplace(yzz_seed);
      seed.emplace_back();
      seed.back().resize(tmp_seed.size());
      transform(tmp_seed.begin(), tmp_seed.end(), seed.back().begin(), [](const pair<int, double> &p) -> int {
        return p.first;
      });
    }
  }
}

void CVRP::constructCutsByResourceGapOnly(const vector<unordered_map<int, int>> &v_r_map,
                                              const vector<vector<int>> &seed,
                                              const vector<double> &frac_routes,
                                              vector<pair<vector<int>, int>> &cuts) {
  vector<double> num_times_vis_routes(frac_routes.size());
  vector<tuple<vector<int>, int, double>> cuts_pool;
  int std_size = (Config::MaxRowRank1 - 3) * Config::MaxNumR1CPerRound;
  cuts_pool.reserve(std_size);
  double vio_threshold = -1;
  int route_size = (int) frac_routes.size();
  for (auto &cut : seed) {
    int cut_size = (int) cut.size();
    double local_max_vio = vio_threshold;
    int max_idx = -1;
    for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
      auto &plan = map_rank1_multiplier[cut_size][plan_idx];
      if (!get<1>(plan)) continue;
      memset(num_times_vis_routes.data(), 0, sizeof(double) * route_size);
      auto &coeff = get<0>(plan);
      int denominator = get<1>(plan);
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
                [denominator](double d, double local_f) -> double { return int(d / denominator) * local_f; });
      double vio =
          accumulate(num_times_vis_routes.begin(),
                     num_times_vis_routes.end(),
                     -rhs);
      if (vio > CUT_VIO_TOLERANCE) {
        if (vio > local_max_vio) {
		  local_max_vio = vio;
          max_idx = plan_idx;
        }
      }
    };
    if (max_idx == -1) continue;
    cuts_pool.emplace_back(cut, max_idx, local_max_vio);
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
}

void CVRP::constructCuts(const std::vector<std::unordered_map<int, int>> &v_r_map,
                          const vector<vector<int>> &seed,
                          const vector<double> &frac_routes,
                          vector<pair<vector<int>, int>> &cuts) {
  vector<tuple<vector<int>, int, double>> cuts_pool;
  int std_size = (Config::MaxRowRank1 - 3) * Config::MaxNumR1CPerRound;
  cuts_pool.reserve(std_size);
  double vio;
  int route_size = (int) frac_routes.size();
  int number = 0;
  for (auto &cut : seed) {
    int cut_size = (int) cut.size();
    double local_max_vio = TOLERANCE;
    int max_idx = -1;
    vector<int> best_cut;
    for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
      auto &plan = map_rank1_multiplier[cut_size][plan_idx];
      if (!get<1>(plan)) continue;
      vector<int> new_cut;
      heuristicFindBestPermutationForOnePlan(v_r_map, frac_routes, cut, plan_idx, new_cut, vio);
      if (vio > CUT_VIO_TOLERANCE) {
        if (vio > local_max_vio) {
		  local_max_vio = vio;
          max_idx = plan_idx;
          best_cut = new_cut;
        }
      }
    }
    if (max_idx == -1) continue;
    ++number;
    cuts_pool.emplace_back(best_cut, max_idx, local_max_vio);
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
  chooseCuts(tmp_cuts, cuts, std_size);
}

void CVRP::heuristicFindBestPermutationForOnePlan(
    const std::vector<std::unordered_map<int, int>> &v_r_map,
    const vector<double> &frac_routes,
    const vector<int> &cut, int plan_idx,
    vector<int> &new_cut,
    double &vio) {
  bool is_sorted = std::is_sorted(frac_routes.begin(), frac_routes.end(), greater<>());
  if (!is_sorted) throw runtime_error("frac_routes is not sorted!");
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
              [denominator](double d, double local_f) -> double { return int(d / denominator) * local_f; });
    vio = accumulate(num_times_vis_routes.begin(),
                     num_times_vis_routes.end(),
                     -rhs);
    return;
  }
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
  new_cut.clear();
  unordered_map<int, int> new_cut_coeff;
  auto cp_coeff = coeff;
  for (auto &pr : vec_r_v) {
    auto &vec_v = pr.second;
    int cnt = 0, cnt2 = 0;
    int ccnt = 0;
    vector<pair<int, int>> v_left;
    for (auto &i : vec_v) {
      int c = i.first;
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
    cnt /= denominator;
    cnt2 /= denominator;
    if (cnt != cnt2) {
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
      if (map_coeff_kind_num.size() <= 1) {
        for (auto &c : cut) {
          if (std::find(new_cut.begin(), new_cut.end(), c) == new_cut.end())
            new_cut.emplace_back(c);
        }
      }
    }
    if (new_cut.size() == cut.size()) goto QUIT;
  }
  new_cut = cut;
  QUIT:
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
            [denominator](double d, double local_f) -> double { return int(d / denominator) * local_f; });
  vio = accumulate(num_times_vis_routes.begin(),
                   num_times_vis_routes.end(),
                   -rhs);
}

void CVRP::newFindR1CMulti(const std::vector<std::vector<int>> &routes,
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
  constructVRVec(routes, v_r_vec, v_r_map);
#ifdef SYMMETRY_PROHIBIT
  construct_v_resourceGap<false>(routes, v_resourceGap);
#else
  construct_v_resourceGap<true>(routes, v_resourceGap);
#endif
  vector<vector<int>> seed;
  constructSeed(routes, v_r_vec, v_resourceGap, seed);
  constructCuts(v_r_map, seed, frac_routes, cuts);
}

void CVRP::addR1CAtOnce(BbNode *node,
                         const vector<tuple<vector<int>, int, int, set<int>>> &full_cuts
) {
  vector<vector<int>> coeffs(full_cuts.size(), vector<int>(num_col, 0));
  vector<int> states(full_cuts.size());
  vector<vector<pair<int, int>>> v_cuts(dim);//(idx, increase)
  vector<vector<bool>> v_Mem(dim, vector<bool>(full_cuts.size(), false));
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

  vector<vector<int>> v_lostMem(dim);
  for (int i = 1; i < dim; ++i) {
    for (int j = 0; j < full_cuts.size(); ++j) {
      if (!v_Mem[i][j]) {
        v_lostMem[i].emplace_back(j);
      }
    }
  }

  int beg = 1, end = node->num_parent_cols, if_rep = false;
  int *col_pool = col_pool4_mem;
  here:
  for (int i = beg; i < end; ++i) {
    memset(states.data(), 0, sizeof(int) * full_cuts.size());
    for (auto j = node->index_columns[i] + 1;; ++j) {
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
    end = num_col;
    col_pool = col_pool4_pricing;
    goto here;
  }
  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
  for (int i = 0; i < full_cuts.size(); ++i) {
    int numnz = 0;
    for (int j = 0; j < num_col; ++j) {
	  if (coeffs[i][j] != 0) {
		solver_ind[numnz] = j;
		solver_val[numnz] = coeffs[i][j];
		++numnz;
	  }
	}
    int cut_index = get<2>(full_cuts[i]);
    if (cut_index == numeric_limits<int>::max()) {
      safe_solver(node->solver.addConstraint(numnz,
                                               solver_ind.data(),
                                               solver_val.data(),
                                               SOLVER_LESS_EQUAL,
                                               rhs[i],
                                               nullptr))
      if (get<1>(full_cuts[i]) == 0) {
        R1c r1c3;
        r1c3.info_r1c = get<0>(full_cuts[i]);
        r1c3.mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
        r1c3.idx_r1c = num_row++;
        r1c3.rhs = rhs[i];
        node->r1cs.emplace_back(r1c3);
      } else {
        R1cMulti r1c;
        r1c.info_r1c = make_pair(get<0>(full_cuts[i]), get<1>(full_cuts[i]));
        r1c.mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
        r1c.idx_r1c = num_row++;
        r1c.rhs = rhs[i];
        node->r1cs_multi.emplace_back(r1c);
      }
    } else {
	  vector<int> solver_ind2(numnz);
      if (get<1>(full_cuts[i]) == 0) {
        fill(solver_ind2.begin(), solver_ind2.end(), node->r1cs[cut_index].idx_r1c);
        node->r1cs[cut_index].mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
      } else {
		fill(solver_ind2.begin(), solver_ind2.end(), node->r1cs_multi[cut_index].idx_r1c);
        node->r1cs_multi[cut_index].mem.assign(get<3>(full_cuts[i]).begin(), get<3>(full_cuts[i]).end());
      }
      safe_solver(node->solver.XchangeCoeffs(numnz, solver_ind2.data(), solver_ind.data(), solver_val.data()))
    }
  }
}


void CVRP::checkR1CTotally(BbNode *node) {
  cout << "WARNING: checkR1CTotally" << endl;
  int num_rows = int(node->r1cs.size() + node->r1cs_multi.size());
  vector<vector<int>> coeffs(num_rows, vector<int>(num_col, 0));
  vector<int> states(num_rows);
  vector<vector<pair<int, int>>> v_cuts(dim);//(idx, increase)
  vector<vector<bool>> v_Mem(dim, vector<bool>(num_rows, false));
  vector<int> denominator(num_rows);
  vector<int> rhs(num_rows);
  vector<tuple<vector<int>, int, int, vector<int>>> full_cuts;// cut, plan, idx, mem
  for (auto &r1c : node->r1cs) {
    full_cuts.emplace_back(r1c.info_r1c, 0, r1c.idx_r1c, r1c.mem);
  }
  for (auto &r1c : node->r1cs_multi) {
    full_cuts.emplace_back(r1c.info_r1c.first, r1c.info_r1c.second, r1c.idx_r1c, r1c.mem);
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

  vector<vector<int>> v_lostMem(dim);
  for (int i = 1; i < dim; ++i) {
    for (int j = 0; j < full_cuts.size(); ++j) {
      if (!v_Mem[i][j]) {
        v_lostMem[i].emplace_back(j);
      }
    }
  }

  int beg = 1, end = node->num_parent_cols, if_rep = false;
  int *col_pool = col_pool4_mem;
  here:
  for (int i = beg; i < end; ++i) {
    memset(states.data(), 0, sizeof(int) * full_cuts.size());
    for (auto j = node->index_columns[i] + 1;; ++j) {
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
    end = num_col;
    col_pool = col_pool4_pricing;
    goto here;
  }


  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
  vector<int> solver_ind2(2);
  for (int i = 0; i < full_cuts.size(); ++i) {
    int numnz;
    int idx = get<2>(full_cuts[i]);
	safe_solver(node->solver.getConstraints(&numnz, solver_ind2.data(), solver_ind.data(), solver_val.data(), idx, 1))
    vector<int> coes(num_col, 0);
    for (int j = 0; j < numnz; ++j) {
      coes[solver_ind[j]] = int(solver_val[j] + 1e-6);
    }
    if (coes != coeffs[i]) {
      cout << "-------------------------------" << endl;
      cout << "wrong in " << idx << endl;
      int j = 0;
      for (; j < num_col; ++j) {
        if (coes[j] != coeffs[i][j]) {
          cout << j << " lp= " << coes[j] << " real= " << coeffs[i][j] << endl;
          break;
        }
      }
      cout << "assume j is in pricing!" << endl;
      for (auto k = node->index_columns[j] + 1;; ++k) {
        int current_node = col_pool4_pricing[k];
        if (!current_node) break;
        cout << current_node << " ";
      }
      cout << endl;
      if (i < node->r1cs.size()) {
        cout << "R1c" << endl;
      } else {
        cout << "R1cMulti" << endl;
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

void CVRP::constructVRVecAggressive(const vector<vector<int>> &routes,
                                        vector<vector<int>> &v_r_map,
                                        std::vector<std::unordered_map<int, int>> &v_r_map2,
                                        vector<vector<pair<vector<int>, vector<int>>>> &c_n_w_no_c,
                                        unordered_map<int, vector<pair<int, int>>> &c_map) {
  v_r_map.resize(dim);
  v_r_map2.resize(dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      v_r_map[i].emplace_back(r);
      ++v_r_map2[i][r];
    }
  }
  int route_size = (int) routes.size();
  c_n_w_no_c.resize(dim);
  for (int i = 1; i < dim; ++i) {
    if (v_r_map[i].empty()) continue;
    yzzLong wc = 0;// c within i
    for (auto &j : v_r_map2[i]) {
      for (auto v : routes[j.first]) {
        if (rank1_sep_heur_mem4_vertex[i][v]) {
          wc.set(v);
        }
      }
    }
    unordered_map<yzzLong, int> seed_map;
    seed_map.reserve(routes.size());
    unordered_set<int> bad_c_size;
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
      if (c_size < 3 || c_size > Config::MaxHeurInitialCsetSize4RowRank1C) continue;
      if (bad_c_size.find(c_size) != bad_c_size.end()) continue;
      bad_c_size.insert(map_rank1_multiplier_dominance[c_size].begin(), map_rank1_multiplier_dominance[c_size].end());
      if (seed_map.find(tmp_c) == seed_map.end()) {
        seed_map[tmp_c] = (int) sub_c_n_w_no_c.size();
        yzzLong w_no_c = wc ^ tmp_c;
        sub_c_n_w_no_c.emplace_back(vector<int>(c_size), vector<int>((int) w_no_c.count()));
        int cnt = 0, cnt2 = 0;
        for (int k = 1; k < dim; ++k) {
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
  }
}

void CVRP::constructOperationsAggressive(
    Rank1MultiLabel &label,
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
    addInRank1HeuristicSeparate(plan_idx, new_cij, w_no_cij, new_vio, choice_add_j, v_r_map, frac_routes);
    move_vio[1].second = new_vio;
  }
  if (dir == 'r' || dir == 's') {
    new_vio = vio;
    removeInRank1HeuristicSeparate(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
    move_vio[2].second = new_vio;
  }
  if (dir == 'a' || dir == 'r') {
    new_vio = vio;
    swapInRank1HeuristicSeparate(plan_idx, new_cij, w_no_cij, new_vio, choice_swap_i_j, v_r_map, frac_routes);
    move_vio[3].second = new_vio;
  }
  stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                   const pair<int, double> &b) {
    return a.second > b.second;
  });
  best_move = move_vio[0].first;
  best_move_vio = move_vio[0].second;

  if (best_move == 0) {
    if (best_move_vio > TOLERANCE && new_cij.size() <= Config::MaxRowRank1 && new_cij.size() > 2) {
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
    new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
    new_cij.emplace_back(choice_swap_i_j.second);
    w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), choice_swap_i_j.second));
  }
  vio = best_move_vio;
}

void CVRP::constructSeedAggressive(
    const vector<vector<int>> &v_r_map,
    const vector<double> &frac_routes,
    std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool,
    vector<Rank1MultiLabel> &rank1_multi_label_pool,
    std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map,
    std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> &c_n_w_no_c,
    int &num_label
) {
  rank1_multi_label_pool.resize(INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE);
  num_label = 0;
  for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
    for (int len = 3; len <= Config::MaxHeurInitialCsetSize4RowRank1C; ++len) {
      for (auto &p : c_map[len]) {
        int i = p.first, j = p.second;// i and non-visited idx j
        vector<int> new_cij;
        double vio;
        bool if_succeed;
        findCombinationsForRankCutSetToGetBestMultiplier(plan_idx,
                                               c_n_w_no_c[i][j].first,
                                               new_cij,
                                               i,
                                               vio,
                                               v_r_map,
                                               frac_routes,
                                               if_succeed);
        if (!if_succeed) continue;//if negative, then no such multiplier
        vector<pair<int, double>> move_vio(4);
        move_vio[0] = {0, vio};
        double new_vio = vio;
        int choice_add_j;
        int choice_remove_j;
        pair<int, int> choice_swap_i_j;
        addInRank1HeuristicSeparate(plan_idx, new_cij, c_n_w_no_c[i][j].second, new_vio, choice_add_j, v_r_map, frac_routes);
        move_vio[1] = {1, new_vio};
        new_vio = vio;
        removeInRank1HeuristicSeparate(plan_idx, new_cij, new_vio, choice_remove_j, v_r_map, frac_routes);
        move_vio[2] = {2, new_vio};
        new_vio = vio;
        swapInRank1HeuristicSeparate(plan_idx, new_cij, c_n_w_no_c[i][j].second, new_vio, choice_swap_i_j, v_r_map, frac_routes);
        move_vio[3] = {3, new_vio};
        stable_sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                                         const pair<int, double> &b) {
          return a.second > b.second;
        });
        int best_move = move_vio[0].first;
        if (best_move == 0) {
          if (move_vio[0].second > TOLERANCE && (new_cij.size() <= Config::MaxRowRank1 && new_cij.size() > 2))
            generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(new_cij, plan_idx, move_vio[0].second);
          continue;
        } else if (best_move == 1) {
          new_cij.emplace_back(choice_add_j);
          auto new_w_no_c = c_n_w_no_c[i][j].second;
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_add_j));
          rank1_multi_label_pool[num_label++] = Rank1MultiLabel{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  'a'};
        } else if (best_move == 2) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_remove_j));
          rank1_multi_label_pool[num_label++] =
              Rank1MultiLabel{new_cij, c_n_w_no_c[i][j].second, plan_idx, move_vio[0].second,
                                'r'};
        } else if (best_move == 3) {
          new_cij.erase(find(new_cij.begin(), new_cij.end(), choice_swap_i_j.first));
          new_cij.emplace_back(choice_swap_i_j.second);
          auto new_w_no_c = c_n_w_no_c[i][j].second;
          new_w_no_c.erase(find(new_w_no_c.begin(), new_w_no_c.end(), choice_swap_i_j.second));
          rank1_multi_label_pool[num_label++] = Rank1MultiLabel{new_cij, new_w_no_c, plan_idx, move_vio[0].second,
                                                                  's'};
        }
        if (num_label == rank1_multi_label_pool.size()) {
          rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
        }
      }
    }
  }
}

void CVRP::findR1CMultiAggressive(const vector<vector<int>> &routes,
                                    const vector<double> &frac_routes,
                                    std::vector<std::unordered_map<int, int>> &v_r_map,
                                    vector<pair<vector<int>, int>> &cuts) {
  vector<vector<int>> v_r_map_special;
  std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> c_n_w_no_c;
  std::unordered_map<int, std::vector<std::pair<int, int>>> c_map;
  constructVRVecAggressive(routes, v_r_map_special, v_r_map, c_n_w_no_c, c_map);

  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> generated_rank1_multi_pool;
  vector<Rank1MultiLabel> rank1_multi_label_pool;
  int num_label;
  constructSeedAggressive(v_r_map_special,
                            frac_routes,
                            generated_rank1_multi_pool,
                            rank1_multi_label_pool,
                            c_map,
                            c_n_w_no_c,
                            num_label);
  for (int i = 0; i < num_label;) {
    constructOperationsAggressive(rank1_multi_label_pool[i],
                                    v_r_map_special,
                                    frac_routes,
                                    i,
                                    generated_rank1_multi_pool);
  }
  constructCutsAggressive(generated_rank1_multi_pool, cuts);
}

void CVRP::checkR1CMultiLook(BbNode *node, int start) {
  unordered_map<int, int> tmp_map;
  unordered_map<yzzLong, pair<int, int>> tmp_map2;
  bool if_print_all = false;
  int idx = 0;
  for (auto &r1c : node->r1cs) {
    if (r1c.idx_r1c >= start) {
      if (if_print_all) {
        for (auto c : r1c.info_r1c) {
          cout << c << " ";
        }
        cout << endl;
      }
      ++tmp_map[0];
    }
    yzzLong tmp_yzz = 0;
    for (auto c : r1c.info_r1c) {
      tmp_yzz.set(c);
    }
    if (tmp_map2.find(tmp_yzz) != tmp_map2.end() && tmp_map2[tmp_yzz].first == 0) {
      auto &w1 = node->r1cs[tmp_map2[tmp_yzz].second].info_r1c;
      for (auto &i : w1) {
        cout << i << " ";
      }
      cout << endl;
      for (auto c : r1c.info_r1c) {
        cout << c << " ";
      }
      cout << endl;
    }
    tmp_map2[tmp_yzz] = {0, idx};
    ++idx;
  }
  cout << "-----------------" << endl;
  idx = 0;
  for (auto &r1c : node->r1cs_multi) {
    if (r1c.idx_r1c >= start) {
      if (if_print_all) {
        for (auto c : r1c.info_r1c.first) {
          cout << c << " ";
        }
        cout << " p: " << r1c.info_r1c.second << endl;
      }
      ++tmp_map[r1c.info_r1c.second];
    }
    yzzLong tmp_yzz = 0;
    for (auto c : r1c.info_r1c.first) {
      tmp_yzz.set(c);
    }
    if (tmp_map2.find(tmp_yzz) != tmp_map2.end() && tmp_map2[tmp_yzz].first == r1c.info_r1c.second) {
      auto &w1 = node->r1cs_multi[tmp_map2[tmp_yzz].second].info_r1c.first;
      for (auto &i : w1) {
        cout << i << " ";
      }
      cout << " p: " << node->r1cs_multi[tmp_map2[tmp_yzz].second].info_r1c.second << endl;
      for (auto c : r1c.info_r1c.first) {
        cout << c << " ";
      }
      cout << " p: " << r1c.info_r1c.second << endl;
    }
    tmp_map2[tmp_yzz] = {r1c.info_r1c.second, idx};
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

void CVRP::constructCutsAggressive(
    std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool,
    vector<pair<vector<int>, int>> &cuts) {

  for (auto &i : generated_rank1_multi_pool) {
    sort(i.second.begin(), i.second.end(),
         [](auto &a, auto &b) {
           return get<2>(a) > get<2>(b);
         });
  }

  for (auto &i : generated_rank1_multi_pool) {
    unordered_set < yzzLong > cut_set;
    unordered_set<int> p_set;
    double vio_std = get<2>(i.second[0]) * Config::CutVioFactor;
    vector<pair<vector<int>, int>> tmp_cuts;
    for (auto &j : i.second) {
      if (get<2>(j) < vio_std) break;
      yzzLong key = 0;
      for (auto c : get<0>(j)) {
        key.set(c);
      }
      if (cut_set.find(key) != cut_set.end() && p_set.find(get<1>(j)) != p_set.end()) continue;
      tmp_cuts.emplace_back(get<0>(j), get<1>(j));
      cut_set.insert(key);
      p_set.insert(get<1>(j));
#ifdef debugVio
      if (vio_map.find(key) == vio_map.end()) vio_map[key].resize(7, -numeric_limits<double>::max());
      vio_map[key][get<1>(j)] = get<2>(j);
#endif
    }
    chooseCuts(tmp_cuts, cuts, Config::MaxNumR1CPerRound);
  }
}

/**
 * May 18....
 */

/**
 * 先按照 seed 找到 集合 为 4 的情况！如果效果好的话，再看看能不能多加一点！
 */

void CVRP::constructSeedMode3(const vector<vector<int>> &routes,
                                const vector<vector<int>> &v_r_vec,
                                const vector<double> &v_resourceGap,
                                vector<vector<int>> &seed) {
  throw runtime_error("do not even consider! not implemented");
  for (int i = 1; i < dim; ++i) {
    yzzLong tmp = 0;
    tmp = rank1_sep_heur_mem4_vertex[i];
    for (int j = i + 1; j < dim; ++j) {
      if (!tmp.test(j)) continue;
      tmp &= rank1_sep_heur_mem4_vertex[j];
      for (int k = j + 1; k < dim; ++k) {
        if (!tmp.test(k)) continue;
        tmp &= rank1_sep_heur_mem4_vertex[k];
        for (int l = k + 1; l < dim; ++l) {
          if (!tmp.test(l)) continue;
          tmp &= rank1_sep_heur_mem4_vertex[l];
          for (int p = l + 1; p < dim; ++p) {
            if (!tmp.test(p)) continue;
            if (!rank1_sep_heur_mem4_vertex[p].test(i) ||
                !rank1_sep_heur_mem4_vertex[p].test(j) ||
                !rank1_sep_heur_mem4_vertex[p].test(k) ||
                !rank1_sep_heur_mem4_vertex[p].test(l))
              continue;
            seed.emplace_back(vector<int>{i, j, k, l, p});
          }
        }
      }
    }
  }
  cout << "seed.size()= " << seed.size() << endl;
}

void CVRP::newFindR1CMultiMode3(const vector<vector<int>> &routes,
                                  const vector<double> &frac_routes,
                                  std::vector<std::unordered_map<int, int>> &v_r_map,
                                  vector<pair<vector<int>, int>> &cuts) {
  throw runtime_error("do not even consider! not implemented");
  vector<vector<int>> v_r_vec;
  vector<double> v_resourceGap;
  constructVRVec(routes, v_r_vec, v_r_map);
#ifdef SYMMETRY_PROHIBIT
  construct_v_resourceGap<false>(routes, v_resourceGap);
#else
  construct_v_resourceGap<true>(routes, v_resourceGap);
#endif
  vector<vector<int>> seed;
  constructSeedMode3(routes, v_r_vec, v_resourceGap, seed);
  constructCuts(v_r_map, seed, frac_routes, cuts);
}

void CVRP::fillMemoryFirst(BbNode *node,
                          const std::vector<std::vector<int>> &routes,
                          const std::vector<double> &frac_routes,
                          std::vector<std::pair<std::vector<int>, int>> &cuts) {
  if (!cuts.empty()) throw runtime_error("cuts should be empty at first!");
  reset_cut_mem.clear();
  /**
   * there's a v_map, could use the v_map outside, but could easily lead to bad results!
   */
  std::vector<std::unordered_map<int, int>> v_r_map(dim);
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
  for (auto &r1c : node->r1cs) {
    memset(&vis_times[0], 0, sizeof(double) * routes.size());
    for (auto &v : r1c.info_r1c) {
      for (auto &r : v_r_map[v]) {
        vis_times[r.first] += r.second;
      }
    }
    transform(vis_times.begin(), vis_times.end(), frac_routes.begin(), vis_times.begin(), [](auto &a, auto &b) {
      return int(a / 2 + TOLERANCE) * b;
    });
    auto vio = accumulate(vis_times.begin(), vis_times.end(), double(-int((int) r1c.info_r1c.size() / 2 + TOLERANCE)));
    if (vio > TOLERANCE) {
      cuts.emplace_back(r1c.info_r1c, 0);
      yzzLong tmp = 0;
      for (auto &v : r1c.info_r1c) {
        tmp.set(v);
      }
      cut_record[tmp].insert(0);
      ++num_add_mem;
      tmp = 0;
      for (auto &v : r1c.mem) {
        tmp.set(v);
      }
      reset_cut_mem.emplace_back(true, idx, tmp);
    }
    ++idx;
  }

  idx = 0;
  for (auto &r1c : node->r1cs_multi) {
    memset(&vis_times[0], 0, sizeof(double) * routes.size());
    const auto &plan = map_rank1_multiplier[(int) r1c.info_r1c.first.size()][r1c.info_r1c.second];
    const auto &coeff = get<0>(plan);
    auto deno = get<1>(plan);
    auto rhs = get<2>(plan);
    int cnt = 0;
    for (auto &v : r1c.info_r1c.first) {
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
      cuts.emplace_back(r1c.info_r1c);
      yzzLong tmp = 0;
      for (auto &v : r1c.info_r1c.first) {
        tmp.set(v);
      }
      cut_record[tmp].insert(r1c.info_r1c.second);
      ++num_add_mem;
      tmp = 0;
      for (auto &v : r1c.mem) {
        tmp.set(v);
      }
      reset_cut_mem.emplace_back(false, idx, tmp);
    }
    ++idx;
  }
  cout << "num_add_mem= " << num_add_mem << endl;
  if (num_add_mem) if_fill_mem = true;
  else if_fill_mem = false;
}

void CVRP::chooseCuts(const vector<pair<vector<int>, int>> &tmp_cuts,
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



