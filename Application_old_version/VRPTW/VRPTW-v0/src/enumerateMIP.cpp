#include <iostream>

#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace chrono;
using namespace Eigen;


void CVRP::solveMIP(BbNode *const node, bool if_inEnu) {
  auto beg = high_resolution_clock::now();

  char *xtype = new char[num_col];
  fill_n(xtype, num_col, SOLVER_BINARY);
  safe_solver(node->solver.setEnvCutoff(ub + round_up_tolerance))
  safe_solver(node->solver.setVTypeArray(0, num_col, xtype))

  if (if_inEnu) {
    if_mip_enumeration_suc = true;
    if (num_col > Config::MinNumRoute4MIP) {
      safe_solver(node->solver.setEnvTimeLimit(Config::MIPInEnumerationTimeLimit))
    }
  }

  safe_solver(node->solver.mipOptimize())

  int status;
  safe_solver(node->solver.getStatus(&status))
  if (status == SOLVER_INFEASIBLE || status == SOLVER_INF_OR_UNBD) {
    cout << "Model is infeasible after pre_solving!" << endl;
    cout << "status: " << status << endl;
    lp_val = 1e100;
  } else {
    safe_solver(node->solver.getObjVal(&lp_val))
  }
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end-beg).count();
  cout << "MIP used " << eps << " seconds!" << endl;

  safe_solver(node->solver.setEnvCutoff(1e100))

  if (if_inEnu) {
    if (status == SOLVER_TIME_LIMIT) {
      if_mip_enumeration_suc = false;
      cout << "MIP time limit reached!" << endl;
      max_num_route4_mip = max(int(Config::MIPInEnumerationTimeLimitChgFactor * num_col), Config::MinNumRoute4MIP);
      cout << "max_num_route4_mip = " << max_num_route4_mip << endl;
    }
    safe_solver(node->solver.setEnvTimeLimit(MAX_TIME_LIMIT_FOR_MIP))
    if (ub > ceilTransformedNumberRelated(lp_val - TOLERANCE) + TOLERANCE) {
	  vector<double> X(num_col);
      safe_solver(node->solver.getX(0, num_col, X.data()))
      ub = ceilTransformedNumberRelated(lp_val - TOLERANCE);
      cout << "solve MIP get " << ub << endl;
      ip_opt_sol.clear();
      for (int i = 0; i < num_col; ++i) {
        if (abs(X[i] - 1) < TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(dim);
          tmp.emplace_back(0);
          for (auto j = node->index_columns[i] + 1;; ++j) {
            if (!col_pool4_pricing[j]) break;
            tmp.emplace_back(col_pool4_pricing[j]);
          }
          tmp.emplace_back(0);
          ip_opt_sol.emplace_back(std::move(tmp));
        }
      }
    }

#ifdef MASTER_VALVE_ML
    updateUpperBoundEdgeSolution();
#endif

    if (status == SOLVER_TIME_LIMIT) {
      fill_n(xtype, num_col, SOLVER_CONTINUOUS);
      safe_solver(node->solver.setVTypeArray(0, num_col, xtype))
      safe_solver(node->solver.reoptimize())
    }
  } else {
    if (ub > ceilTransformedNumberRelated(lp_val - TOLERANCE) + TOLERANCE) {
	  vector<double> X(num_col);
      safe_solver(node->solver.getX(0, num_col, X.data()))
      ub = ceilTransformedNumberRelated(lp_val - TOLERANCE);
      cout << "solve MIP get " << ub << endl;
      ip_opt_sol.clear();
      for (int i = 0; i < node->num_parent_cols; ++i) {
        if (abs(X[i] - 1) < TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(dim);
          tmp.emplace_back(0);
          for (auto j = node->index_columns[i] + 1;; ++j) {
            if (!col_pool4_mem[j]) break;
            tmp.emplace_back(col_pool4_mem[j]);
          }
          tmp.emplace_back(0);
          ip_opt_sol.emplace_back(std::move(tmp));
        }
      }
      for (int i = node->num_parent_cols; i < num_col; ++i) {
        if (abs(X[i] - 1) < TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(dim);
          tmp.emplace_back(0);
          for (auto j = node->index_columns[i] + 1;; ++j) {
            if (!col_pool4_pricing[j]) break;
            tmp.emplace_back(col_pool4_pricing[j]);
          }
          tmp.emplace_back(0);
          ip_opt_sol.emplace_back(std::move(tmp));
        }
      }
    }

#ifdef MASTER_VALVE_ML
    updateUpperBoundEdgeSolution();
#endif

    fill_n(xtype, num_col, SOLVER_CONTINUOUS);
    safe_solver(node->solver.setVTypeArray(0, num_col, xtype))
    safe_solver(node->solver.reoptimize())
  }
  delete[]xtype;
}

void CVRP::enumerateMIP(BbNode *&node) {
#if SETTING==2
  cout<<"In setting II, we do not do enumeration!"<<endl;
  return;
#endif
  determineIfEnumeration(node);
  if (!final_decision_4_enumeration) return;
  final_decision_4_enumeration = false;

  cout << BIG_PHASE_SEPARATION;
  cout << "try EnumerateRoutes..." << endl;

  if (abs(meet_point_resource_in_bi_dir_enu) < TOLERANCE) meet_point_resource_in_bi_dir_enu = meet_point_resource_in_bi_dir;
  bool if_succeed;
  double time_labeling;

  PtrAllR1CS ptrAllR1Cs(node, this);

  auto beg = high_resolution_clock::now();

  rollback = 0;
  if_succeed = enumerateRoutes(node, ptrAllR1Cs);

  auto end = high_resolution_clock::now();
  auto eps = duration_cast<milliseconds>(end - beg);
  time_labeling = (double) eps.count() * 1e-3;

  if (if_succeed) {
    cout << "enumeration time= " << time_labeling << " and succeed!" << endl;
  } else {
#if VERBOSE_MODE==1
	cout << "RollBack= " << rollback << endl;
#endif
    cout << "enumeration time= " << time_labeling << " but failed!" << endl;
  }

  double gap = (ub - node->value) / ub;
  if (if_succeed) {
    if (enumeration_mode) {
      count4_tolerance4_try_enumeration_when_arc_elimination_fails = 0;
    }
    if (gap > gap_tolerance4_arc_elimination_n_enumeration) {
      gap_tolerance4_arc_elimination_n_enumeration = gap;
    }
    last_enumeration_fail_gap /= Config::EnumerationFailFactor;
#ifdef SYMMETRY_PROHIBIT
    double dif = abs(num_forward_labels_in_enu - num_backward_labels_in_enu);
    double over = dif / min(num_forward_labels_in_enu, num_backward_labels_in_enu);
    if (over > Config::NumberOfOverLabelsInMeetPoint) {
      if (num_forward_labels_in_enu > num_backward_labels_in_enu) {
        meet_point_resource_in_bi_dir_enu *= (1 - Config::MeetPointFactor);
      } else {
        meet_point_resource_in_bi_dir_enu *= (1 + Config::MeetPointFactor);
      }
      cout << "meet_point_resource_in_bi_dir_enu= " << meet_point_resource_in_bi_dir_enu << endl;
    }
#endif
    if_enumeration_suc = true;
#ifdef BRANCH_FASHION_MEM_SAVING
    terminateNodeMemSave(node);
#else
    terminateNode(node);
#endif
  } else {
    last_enumeration_fail_gap = gap * Config::EnumerationFailFactor;//set a smaller gap!
    if_enumeration_suc = false;
  }
#if VERBOSE_MODE==1
  cout << "last_enumeration_fail_gap= " << last_enumeration_fail_gap << endl;
  cout << BIG_PHASE_SEPARATION;
#endif
}

bool CVRP::enumerateRoutes(BbNode *const node, const PtrAllR1CS &ptrAllR1Cs) {
  if (enumeration_mode) {//do without arc elimination
	cout<<"in this version, we do not support enumeration without arc elimination!"<<endl;
	return false;
  }
  int Max_labels = Config::MaxNumLabelInEnumeration;
  int Max_routes_all = Config::MaxNumRouteInEnumeration;
  int Max_routes_phase1 = (int) sqrt_self((float) Max_routes_all);
  int num_routes_now = 0;
  auto &cost_m = node->cost_for_columns_in_enumeration_column_pool;
  auto &ptr = node->index_columns_in_enumeration_column_pool;
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;

  int index;
  opt_gap = calculateOptimalGap(node);
  unordered_map<yzzLong, tuple<Label *, Label *, double>> Tags;
  Tags.reserve(Max_routes_all);

  auto copy_Forward_bucket = new vector<Label *> *[dim];
#ifdef SYMMETRY_PROHIBIT
  auto copy_Backward_bucket = new vector<Label *> *[dim];//these memory will be freed inside the forward function
#endif

    for (int i = 0; i < dim; ++i) {
      copy_Forward_bucket[i] = new vector<Label *>[num_buckets_per_vertex];
#ifdef SYMMETRY_PROHIBIT
      copy_Backward_bucket[i] = new vector<Label *>[num_buckets_per_vertex];
#endif
      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        copy_Forward_bucket[i][b].assign(label_array_in_forward_sense[i][b].first.begin(),
                                         label_array_in_forward_sense[i][b].first.begin()
                                             + label_array_in_forward_sense[i][b].second);
#ifdef SYMMETRY_PROHIBIT
        copy_Backward_bucket[i][b].assign(label_array_in_backward_sense[i][b].first.begin(),
                                          label_array_in_backward_sense[i][b].first.begin()
                                              + label_array_in_backward_sense[i][b].second);
#endif
      }
    }

  prior_pool_beg4_pricing = pool_beg4_pricing;
  if (checkPricingPool()) reallocatePricingPool();

  auto beg = high_resolution_clock::now();

  int status =
#ifdef SYMMETRY_PROHIBIT
      enumerateHalfwardRoutes<true, false>(node, r1c_to_pi, r1c_multi_to_pi, Tags,
                                           copy_Backward_bucket, num_routes_now);
#else
      enumerateHalfwardRoutes<true, true>(node, r1c_to_pi, r1c_multi_to_pi, Tags, copy_Forward_bucket, num_routes_now);
#endif

  auto end = high_resolution_clock::now();
  auto eps = duration<double> (end-beg).count();
#if VERBOSE_MODE==1
  cout << "Half Forward time= " << eps << "s" << endl;
#endif
  if (status || rollback) {
    if (status == 1)
      cout << "the number of labels in Forward reached its limit!" << endl;
    return false;
  }

#ifdef SYMMETRY_PROHIBIT
  beg = high_resolution_clock::now();

  status = enumerateHalfwardRoutes<false, false>(node,
                                                 r1c_to_pi,
                                                 r1c_multi_to_pi,
                                                 Tags,
                                                 copy_Forward_bucket,
                                                 num_routes_now);

  end = high_resolution_clock::now();
  eps = duration<double> (end-beg).count();
  #if VERBOSE_MODE==1
  cout << "Half Backward time= " << eps << "s" << endl;
  #endif

  if (status || rollback) {
    if (status == 1)
      cout << "the number of labels in Backward reached its limit!" << endl;
    return false;
  }
#endif


  beg = high_resolution_clock::now();

  status = concatenateRoutesPriorForwardInEnumeration(node, r1c_to_pi, r1c_multi_to_pi, Tags, num_routes_now);

  end = high_resolution_clock::now();
  eps = duration<double> (end-beg).count();
#if VERBOSE_MODE==1
  cout << "Concatenate time= " << eps << "s" << endl;
#endif

  if (status) {
    cout << "the number of routes reached its limit!" << endl;
    return false;
  }

  priceLabeling(node);
  organizeColumnsInMemoryToPricing(node);//reset colPool for pricing! be careful!
  cost_m.resize(num_routes_now);
  ptr.resize(num_routes_now);

  auto PricingWarning = (size_t) (0.9 * (double) mem4_pricing);
  index = 0;
  Label *ki, *li, *p;
  vector<int> seq(dim+1);
  for (auto &tag : Tags) {
    ki = get<0>(tag.second);
    li = get<1>(tag.second);

    cost_m[index] = get<2>(tag.second);
    ptr[index++] = pool_beg4_pricing;

	int cnt = 0;
	p = ki;
	while (p) {
	  seq[cnt++] = p->end_vertex;
	  p = p->p_label;
	}
	for (int k = 0; k < cnt; ++k) {
	  col_pool4_pricing[pool_beg4_pricing++] = seq[cnt - 1 - k];
	}
	if (li) {
	  p = li;
	  while (p) {
		col_pool4_pricing[pool_beg4_pricing++] = p->end_vertex;
		p = p->p_label;
	  }
	} else {
	  col_pool4_pricing[pool_beg4_pricing++] = 0;
	}
#ifdef SOLVER_VRPTW
	double local_cap=0;
	for(size_t i = ptr[index-1]+1; i < pool_beg4_pricing; ++i){
	  local_cap+=demand[col_pool4_pricing[i]];
	}
	if(local_cap> cap){
	  --index;
	 pool_beg4_pricing=ptr[index];
	  continue;
	}
#endif
    if (pool_beg4_pricing >= PricingWarning) {
      cout << SMALL_PHASE_SEPARATION;
      cout << "Warning: the pricing pool is almost full!" << endl;
      cout << "pool_beg4_pricing=" << pool_beg4_pricing << endl;
      cout << "mem4_pricing=" << mem4_pricing << endl;
      cout << "we reallocate the pricing pool!" << endl;
      reallocatePricingPool();
      PricingWarning = (size_t) (0.9 * (double) mem4_pricing);
      cout << "the new mem4_pricing=" << mem4_pricing << endl;
      cout << SMALL_PHASE_SEPARATION;
    }
  }
  num_routes_now=index;
  cost_m.resize(num_routes_now);
  ptr.resize(num_routes_now);
  node->size_enumeration_col_pool = num_routes_now;
  node->valid_size = num_routes_now;
  cout<<"now the number of routes is "<<num_routes_now<<endl;

  cleanColumnsNonElement(node);
#ifdef SOLVER_VRPTW
  for(auto &rcc:node->rccs)rcc.if_keep=false;
  cleanColumnsNonFeasible(node);
#endif
  deleteBranchCutsAndR1C1s(node);
  recoverR1CsInEnum(node);

  generateVertex2IndexColsAndEdge2IndexCols(node);
  return true;
}


int CVRP::concatenateRoutesPriorForwardInEnumeration(BbNode *node,
                                                        const double *r1c_to_pi,
                                                        const double *r1c_multi_to_pi,
                                                        unordered_map<yzzLong, tuple<Label *, Label *, double>> &Tags,
                                                        int &num_routes_now) {
  int status = 0;
  double path_rc, path_cost;
  yzzLong tmp_PI;
#ifdef SYMMETRY_PROHIBIT
  populateRC2TillThisBinNRC2Bin<false>(node);
#else
  populateRC2TillThisBinNRC2Bin<true>(node);//use function here!
#endif
  for (auto &label_list : concatenate_labels_in_forward_cg) {
    int i = label_list.first.first;
    int j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      double tmp_mainResource = pr.second;
      double tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
#ifdef SYMMETRY_PROHIBIT
      int arr_bj = int((tmp_mainResource) / step_size);
      if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap)
        continue;

      if (rc2_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap &&
          (std::find(node->all_backward_buckets[j][arr_bj].bucket_arcs.begin(),
                     node->all_backward_buckets[j][arr_bj].bucket_arcs.end(), i)
              != node->all_backward_buckets[j][arr_bj].bucket_arcs.end())) {
        auto &label_arr = label_array_in_backward_sense[j][arr_bj].first;
        auto &label_valid_num = label_array_in_backward_sense[j][arr_bj].second;
#else
      int arr_bj = int((max_main_resource - tmp_mainResource) / step_size);
      if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
        continue;

      if (rc2_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap
          ) {
        auto &label_arr = label_array_in_forward_sense[j][arr_bj].first;
        auto &label_valid_num = label_array_in_forward_sense[j][arr_bj].second;
#endif
        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
          auto &kj = label_arr[vec_index];
          path_rc = kj->rc + tmp_rc;
          if (path_rc > opt_gap) break;
#ifdef SYMMETRY_PROHIBIT
          if (tmp_mainResource > kj->sum_main_resource) continue;
#else
          if (tmp_mainResource + kj->sum_main_resource > max_main_resource) continue;
#endif
          if ((ki->pi & kj->pi).any()) continue;

          if (ki->num_valid_rank1_cut < kj->num_valid_rank1_cut) {
            for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
              if (kj->rank1_cut_mem[ki->valid_rank1_cut[l]]) {
                path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
                if (path_rc > opt_gap) goto here;
              }
            }
          } else {
            for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
              if (ki->rank1_cut_mem[kj->valid_rank1_cut[l]]) {
                path_rc -= r1c_to_pi[kj->valid_rank1_cut[l]];
                if (path_rc > opt_gap) goto here;
              }
            }
          }

          if (ki->num_valid_rank1_cut_multi < kj->num_valid_rank1_cut_multi) {
            for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
              int tmp_cut = ki->valid_rank1_cut_multi[l];
              if (kj->rank1_cut_mem_multi[tmp_cut] +
                  ki->rank1_cut_mem_multi[tmp_cut]
                  >= r1c_multi_denominator_in_cg[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > opt_gap) goto here;
              }
            }
          } else {
            for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
              int tmp_cut = kj->valid_rank1_cut_multi[l];
              if (ki->rank1_cut_mem_multi[tmp_cut] +
                  kj->rank1_cut_mem_multi[tmp_cut]
                  >= r1c_multi_denominator_in_cg[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > opt_gap) goto here;
              }
            }
          }
          path_cost = ki->cost + cost_mat4_vertex[i][j] + kj->cost;
          tmp_PI = ki->pi | kj->pi;
          if (Tags.find(tmp_PI) == Tags.end()) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
            ++num_routes_now;
            if (num_routes_now > Config::MaxNumRouteInEnumeration) {
              status = 1;
              goto QUIT;
            }
          } else if (get<2>(Tags[tmp_PI]) > path_cost) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
          }
          here:;
        }
      }
#ifdef SYMMETRY_PROHIBIT
      for (++arr_bj; arr_bj < num_buckets_per_vertex; ++arr_bj) {
        if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap)
          break;

        if (rc2_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap ||
            (std::find(node->all_backward_buckets[j][arr_bj].bucket_arcs.begin(),
                       node->all_backward_buckets[j][arr_bj].bucket_arcs.end(), i)
                == node->all_backward_buckets[j][arr_bj].bucket_arcs.end()))
          continue;

        auto &label_arr = label_array_in_backward_sense[j][arr_bj].first;
        auto &label_valid_num = label_array_in_backward_sense[j][arr_bj].second;
#else
      for (--arr_bj; arr_bj >= 0; --arr_bj) {
        if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
          break;

        if (rc2_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap
            )
          continue;

        auto &label_arr = label_array_in_forward_sense[j][arr_bj].first;
        auto &label_valid_num = label_array_in_forward_sense[j][arr_bj].second;
#endif
        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
          auto &kj = label_arr[vec_index];
          path_rc = kj->rc + tmp_rc;
          if (path_rc > opt_gap) break;

          if ((ki->pi & kj->pi).any()) continue;

          if (ki->num_valid_rank1_cut < kj->num_valid_rank1_cut) {
            for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
              if (kj->rank1_cut_mem[ki->valid_rank1_cut[l]]) {
                path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
                if (path_rc > opt_gap) goto here2;
              }
            }
          } else {
            for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
              if (ki->rank1_cut_mem[kj->valid_rank1_cut[l]]) {
                path_rc -= r1c_to_pi[kj->valid_rank1_cut[l]];
                if (path_rc > opt_gap) goto here2;
              }
            }
          }

          if (ki->num_valid_rank1_cut_multi < kj->num_valid_rank1_cut_multi) {
            for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
              int tmp_cut = ki->valid_rank1_cut_multi[l];
              if (kj->rank1_cut_mem_multi[tmp_cut] +
                  ki->rank1_cut_mem_multi[tmp_cut]
                  >= r1c_multi_denominator_in_cg[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > opt_gap) goto here2;
              }
            }
          } else {
            for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
              int tmp_cut = kj->valid_rank1_cut_multi[l];
              if (ki->rank1_cut_mem_multi[tmp_cut] +
                  kj->rank1_cut_mem_multi[tmp_cut]
                  >= r1c_multi_denominator_in_cg[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > opt_gap) goto here2;
              }
            }
          }
          path_cost = ki->cost + cost_mat4_vertex[i][j] + kj->cost;
          tmp_PI = ki->pi | kj->pi;
          if (Tags.find(tmp_PI) == Tags.end()) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
            ++num_routes_now;
            if (num_routes_now > Config::MaxNumRouteInEnumeration) {
              status = 1;
              goto QUIT;
            }
          } else if (get<2>(Tags[tmp_PI]) > path_cost) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
          }
          here2:;
        }
      }
    }
  }
  QUIT:
  cout << "num_routes_now= " << num_routes_now << endl;
  return status;
}

#ifdef SYMMETRY_PROHIBIT
int CVRP::enumerateHalfBackwardRoutes(BbNode *const node,
                                      const double *r1c_to_pi,
                                      const double *r1c_multi_to_pi,
                                      vector<Label *> **copy_bucket) {
  int status = 0;
  int edgemap;
  num_backward_labels_in_enu = 0;
  bool if_keep, if_break;
  double path_rc, path_cost;
  initializeLabels(node, 2, false, {true, 2, false});

  auto beg = high_resolution_clock::now();
  auto end = beg;
  auto b4_end = beg;
  auto af_end = beg;
  double eps;
  double eps2;
  double left_time = Config::HardTimeThresholdInAllEnumeration;

  for (int b = num_buckets_per_vertex - 1; b >= 0; --b) {
    int i = 1;
    STILL_EXIST:
    for (; i < dim; ++i) {
      end = high_resolution_clock::now();
      eps = duration<double>(end - b4_end).count();
      if (eps > left_time) {
        status = 2;
        goto outside;
      }
      auto &valid_num = if_exist_extra_labels_in_backward_sense[i][b].second;
      if (!valid_num) continue;
      auto &label_array = if_exist_extra_labels_in_backward_sense[i][b].first;
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->is_extended) continue;
        ki->is_extended = true;
        for (int j : node->all_backward_buckets[i][b].bucket_arcs) {
          if (ki->pi[j]) continue;
          auto &tmp_mainResource = all_label[idx_glo].sum_main_resource;
          if (!decreaseMainResourceConsumption(ki->sum_main_resource, tmp_mainResource, i, j)) continue;
          if (tmp_mainResource < meet_point_resource_in_bi_dir_enu) {
            continue;
          }
          auto &tmp_rc = all_label[idx_glo].rc;
          tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];//real rc
          if_keep = false;

          int arr_bj = int(tmp_mainResource / step_size);
          if (tmp_rc + rc2_till_this_bin_in_forward_sense[j][arr_bj] > opt_gap) continue;
          if (tmp_rc + rc2_bin_in_forward_sense[j][arr_bj] < opt_gap) {
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if (tmp_mainResource < kkj->sum_main_resource) continue;
              if ((ki->pi & kkj->pi).any()) continue;
              path_rc = tmp_rc + kkj->rc;
              if (path_rc > opt_gap) break;
              if (ki->num_valid_rank1_cut < kkj->num_valid_rank1_cut) {
                for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
                  if (kkj->rank1_cut_mem[ki->valid_rank1_cut[l]])path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut; ++l) {
                  if (ki->rank1_cut_mem[kkj->valid_rank1_cut[l]])path_rc -= r1c_to_pi[kkj->valid_rank1_cut[l]];
                }
              }
              if (ki->num_valid_rank1_cut_multi < kkj->num_valid_rank1_cut_multi) {
                for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = ki->valid_rank1_cut_multi[l];
                  if (kkj->rank1_cut_mem_multi[tmp_cut] +
                      ki->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = kkj->valid_rank1_cut_multi[l];
                  if (ki->rank1_cut_mem_multi[tmp_cut] +
                      kkj->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              }
              if (path_rc < opt_gap) {
                if_keep = true;
                goto outside1;
              }
            }
          }
          for (--arr_bj; arr_bj >= 0; --arr_bj) {
            if (tmp_rc + rc2_till_this_bin_in_forward_sense[j][arr_bj] > opt_gap) break;
            if (tmp_rc + rc2_bin_in_forward_sense[j][arr_bj] > opt_gap) continue;
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if ((ki->pi & kkj->pi).any()) continue;
              path_rc = tmp_rc + kkj->rc;
              if (path_rc > opt_gap) break;
              if (ki->num_valid_rank1_cut < kkj->num_valid_rank1_cut) {
                for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
                  if (kkj->rank1_cut_mem[ki->valid_rank1_cut[l]])path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut; ++l) {
                  if (ki->rank1_cut_mem[kkj->valid_rank1_cut[l]])path_rc -= r1c_to_pi[kkj->valid_rank1_cut[l]];
                }
              }
              if (ki->num_valid_rank1_cut_multi < kkj->num_valid_rank1_cut_multi) {
                for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = ki->valid_rank1_cut_multi[l];
                  if (kkj->rank1_cut_mem_multi[tmp_cut] +
                      ki->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = kkj->valid_rank1_cut_multi[l];
                  if (ki->rank1_cut_mem_multi[tmp_cut] +
                      kkj->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              }

              if (path_rc < opt_gap) {
                if_keep = true;
                goto outside1;
              }
            }
          }
          outside1:
          if (!if_keep) continue;
          int bj = int(tmp_mainResource / step_size);
          auto &labelList_j = label_array_in_backward_sense[j][bj].first;
          auto &valid_num_j = label_array_in_backward_sense[j][bj].second;
          auto &tmp_PI = all_label[idx_glo].pi;
          auto &tmp_Cost = all_label[idx_glo].cost;
          auto &tmp_Rank1CutMem = all_label[idx_glo].rank1_cut_mem;
          auto &tmp_num_valid_rank1_cut = all_label[idx_glo].num_valid_rank1_cut;
          auto &tmp_valid_rank1_cut = all_label[idx_glo].valid_rank1_cut;
          tmp_PI = ki->pi;
          tmp_PI.set(j);
          tmp_Cost = ki->cost + cost_mat4_vertex[i][j];
          tmp_Rank1CutMem = ki->rank1_cut_mem;
          for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_Rank1CutMem[l] = false;
              tmp_rc -= r1c_to_pi[l];
            } else tmp_Rank1CutMem[l] = true;
          }
          tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);

          auto &tmp_Rank1CutMem_multi = all_label[idx_glo].rank1_cut_mem_multi;
          auto &tmp_num_valid_rank1_cut_multi = all_label[idx_glo].num_valid_rank1_cut_multi;
          auto &tmp_valid_rank1_cut_multi = all_label[idx_glo].valid_rank1_cut_multi;
          copy(ki->rank1_cut_mem_multi, ki->rank1_cut_mem_multi + num_valid_r1c_multi_in_cg, tmp_Rank1CutMem_multi);
          for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
            int tmp_cut = get<0>(l);
            tmp_Rank1CutMem_multi[tmp_cut] += get<1>(l);
            if (tmp_Rank1CutMem_multi[tmp_cut] >= get<2>(l)) {
              tmp_rc -= r1c_multi_to_pi[tmp_cut];
              tmp_Rank1CutMem_multi[tmp_cut] -= get<2>(l);
            }
          }
          for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;

          if_break = false;

          for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
            auto &kj = labelList_j[vec_index_j];
#ifdef CAPACITY_AS_MAIN_RESOURCE
            if (abs(kj->sum_main_resource - tmp_mainResource) > TOLERANCE) {
              ++vec_index_j;
              continue;
            }
            if ((kj->pi ^ tmp_PI).none()) {
              if (kj->cost > tmp_Cost) {
                kj->is_extended = true;
                kj = labelList_j[--valid_num_j];
                --num_backward_labels_in_enu;
              } else {
                if_break = true;
                break;
              }
            } else ++vec_index_j;
#else
            if (kj->cost > tmp_Cost) {
              if (kj->sum_main_resource < tmp_mainResource) {
                if ((kj->pi ^ tmp_PI).none()) {
                  kj->is_extended = true;
                  kj = labelList_j[--valid_num_j];
                  --num_backward_labels_in_enu;
                } else ++vec_index_j;
              } else ++vec_index_j;
            } else {
              if (kj->sum_main_resource > tmp_mainResource) {
                if ((kj->pi ^ tmp_PI).none()) {
                  if_break = true;
                  break;
                } else ++vec_index_j;
              } else ++vec_index_j;
            }
#endif
          }
          if (if_break) continue;

          tmp_num_valid_rank1_cut = 0;
          for (auto l : get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
            }
          }

          tmp_num_valid_rank1_cut_multi = 0;
          for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
            if (tmp_Rank1CutMem_multi[l]) {
              tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
            }
          }

          labelList_j[valid_num_j++] = all_label + idx_glo;
          if (valid_num_j == labelList_j.size()) {
            labelList_j.resize(labelList_j.size() * 2);
          }

		  all_label[idx_glo].p_label=ki;
          all_label[idx_glo].end_vertex = j;
          all_label[idx_glo].is_extended = false;
          auto &bucket = if_exist_extra_labels_in_backward_sense[j][bj];
          bucket.first[bucket.second++] = all_label + idx_glo;
          if (bucket.second == bucket.first.size()) {
            bucket.first.resize(bucket.first.size() * 2);
          }

          if (++num_backward_labels_in_enu > Config::MaxNumLabelInEnumeration) {
            status = 3;//all labels limit
            goto outside;
          }
          ++idx_glo;//can be put here, because once go outside, the function will end
          if (idx_glo == label_assign) {
            rollback = 2;
            goto outside;
          }
        }
      }
      valid_num = 0;
    }
    for (i = 1; i < dim; ++i) {
      if (if_exist_extra_labels_in_backward_sense[i][b].second)
        goto STILL_EXIST;
    }
    af_end = high_resolution_clock::now();
    eps2 = duration<double>(af_end - b4_end).count();
    eps = duration<double>(af_end - beg).count();
    left_time = (Config::HardTimeThresholdInAllEnumeration - eps) / (b + 1);
    if (eps2 > left_time) {
      status = 2;
      goto outside;
    }
    b4_end = af_end;
  }
  outside:
  for (int i = 0; i < dim; ++i) {
    delete[]copy_bucket[i];
  }
  delete[] copy_bucket;
  cout << "Half Backward labeling: num_labels= " << num_backward_labels_in_enu << endl;
  if (status)return status;

  for (int i = 1; i < dim; ++i) {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      std::stable_sort(label_array_in_backward_sense[i][b].first.begin(),
                       label_array_in_backward_sense[i][b].first.begin()
                           + label_array_in_backward_sense[i][b].second,
                       CmpLabelRCLess);
    }
  }
  return 0;
}
#endif

int CVRP::generateColumnsByInspection(BbNode *node,
                                   bool if_only_need_value) {
  if (node->size_enumeration_col_pool == 0) return 0;
  int size_pool = node->size_enumeration_col_pool;
  auto &deleted_columns_in_enumeration_pool = node->deleted_columns_in_enumeration_pool;
  auto &mat = node->matrix_in_enumeration;

  safe_solver(node->solver.getDual(0, num_row, pi))
  safe_solver(node->solver.getObjVal(&lp_val))

  RowVectorXd rc = node->cost_for_columns_in_enumeration_column_pool;

  int num = 0;
  for (auto &it : mat) {
    RowVectorXd dual(it.rows());
    for (int i = 0; i < it.rows(); ++i)dual(i) = pi[num++];
    rc -= dual * it;
  }

  vector<int> col_to_be_added;
  col_to_be_added.reserve(size_pool);

  int cnt = 0;
  if (if_only_need_value) {
    for (int i = 0; i < size_pool; ++i) {
      if (deleted_columns_in_enumeration_pool[i]) continue;
      if (rc(i) < RC_TOLERANCE) {
        col_to_be_added.emplace_back(i);
        ++cnt;
        if (cnt == Config::MaxNumRoutesInExact) break;
      }
    }
  } else {
    for (int i = 0; i < size_pool; ++i) {
      if (deleted_columns_in_enumeration_pool[i]) continue;
      if (rc(i) < RC_TOLERANCE) {
        col_to_be_added.emplace_back(i);
        deleted_columns_in_enumeration_pool[i] = true;
        ++cnt;
        if (cnt == Config::MaxNumRoutesInExact) break;
      }
    }
  }

  if (col_to_be_added.empty()) {
    if (!if_only_need_value) {
      node->value = lp_val;
      opt_gap = calculateOptimalGap(node);
      for (int i = 0; i < size_pool; ++i) {
        if (rc(i) > opt_gap) {
          deleted_columns_in_enumeration_pool[i] = true;
        }
      }
      cleanColumnsRCGreaterThanOptimalGapInLP(node);
      regenerateEnumMat(node, nullptr);
    }
  } else {
    addColumnsByInspection(node, col_to_be_added);
  }
  return cnt;
}

void CVRP::generateVertex2IndexColsAndEdge2IndexCols(BbNode *node) {
  int size_pool = node->size_enumeration_col_pool, curr_node;
  if (size_pool == 0) return;
  auto &ptr = node->index_columns_in_enumeration_column_pool;
  sparseRowMatrixXd mat(num_row, size_pool);

  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(num_row) * double(size_pool) * 0.1));

  auto &deleted_columns_in_enumeration_pool = node->deleted_columns_in_enumeration_pool;
  deleted_columns_in_enumeration_pool.resize(size_pool, false);

  auto &map = node->column_pool_mapping;
  map.clear();
  for (int i = 0; i < size_pool; ++i) {
    int past_node = 0;
    for (auto j = ptr[i] + 1;; ++j) {
      curr_node = col_pool4_pricing[j];
      if (past_node < curr_node)
        map[{past_node, curr_node}].emplace_back(i);
      else map[{curr_node, past_node}].emplace_back(i);
      if (!curr_node) break;
      triplets.emplace_back(curr_node - 1, i, 1);
      past_node = curr_node;
    }
  }

  sparseRowMatrixXd tmpMat(real_dim, size_pool);
  tmpMat.setFromTriplets(triplets.begin(), triplets.end());

  for (int i = 0; i < size_pool; ++i) triplets.emplace_back(real_dim, i, 1);

  sparseRowMatrixXd sum(1, size_pool);
  unordered_map<int, double> tmp;
  tmp.reserve(size_pool);

  for (auto &rcc : node->rccs) {
    tmp.clear();
    if (rcc.form_rcc) {
      auto &info = rcc.info_rcc_customer;
      for (auto iter = info.begin(); iter != info.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != info.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
    } else {
      auto &infoRccCustomer = rcc.info_rcc_customer;
      auto &infoRccOutsideCustomer = rcc.info_rcc_outside_customer;
      for (auto iter = infoRccOutsideCustomer.begin(); iter != infoRccOutsideCustomer.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != infoRccOutsideCustomer.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
      for (auto customer_it : infoRccOutsideCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] += 0.5;
      }
      for (auto customer_it : infoRccCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] -= 0.5;
      }
    }
    int row = rcc.idx_rcc;
    for (auto &it : tmp) {
      if (abs(it.second) > TOLERANCE) {
        triplets.emplace_back(row, it.first, it.second);
      }
    }
  }

  for (auto &r1c : node->r1cs) {
    sum.setZero();
    auto &info = r1c.info_r1c;
    for (auto j : info) {
      sum += tmpMat.row(j - 1);
    }
    sum /= 2;
    int row = r1c.idx_r1c;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }

  for (auto &r1c : node->r1cs_multi) {
    sum.setZero();
    auto &info = r1c.info_r1c;
    const auto &plan = map_rank1_multiplier[(int) info.first.size()][info.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto &j : info.first) {
      sum += tmpMat.row(j - 1) * multi[count++];
    }
    sum /= denominator;
    int row = r1c.idx_r1c;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }

  mat.setFromTriplets(triplets.begin(), triplets.end());
  node->matrix_in_enumeration.push_back((mat));

  vector<bool> bad_cols(size_pool);
  for (auto &brc : node->brcs) {
    if (brc.br_dir) {//must use: 1. use and only use one edge .2 use two but not next to each other
      int ai = brc.edge.first, aj = brc.edge.second;
      if (ai) {
        sum = mat.row(ai - 1) + mat.row(aj - 1);
        std::fill(bad_cols.begin(), bad_cols.end(), true);
        for (int it : map[{ai, aj}]) bad_cols[it] = false;
        for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
          if (it.value() > 0.5) {
            if (it.value() < 1.5) {//==1
              deleted_columns_in_enumeration_pool[it.col()] = true;
            } else {//==2 further test
              if (bad_cols[it.col()]) deleted_columns_in_enumeration_pool[it.col()] = true;
            }
          }
        }
      } else {
        std::fill(bad_cols.begin(), bad_cols.end(), true);
        for (int it : map[{0, aj}]) bad_cols[it] = false;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, aj - 1); it; ++it) {
          if (it.value() > 0.5 && bad_cols[it.col()]) {
            deleted_columns_in_enumeration_pool[it.col()] = true;
          }
        }
      }
    }
  }
  regenerateEnumMat(node, nullptr, false);
  if (size_pool != node->size_enumeration_col_pool)
    cout << "now size= " << node->size_enumeration_col_pool << endl;
}

void CVRP::addColumnsByInspection(BbNode *const node, const vector<int> &Col_added) {
  if (Col_added.empty()) return;
  auto &ptr = node->index_columns_in_enumeration_column_pool;
  auto &ptr_cost = node->cost_for_columns_in_enumeration_column_pool;
  auto &mat = node->matrix_in_enumeration;
  auto &col_idx = node->index_columns;
  int index = num_col;

  vector<int> if_col_added(node->size_enumeration_col_pool, -1);
  int cnt = 0;
  for (auto col : Col_added) {
    if_col_added[col] = cnt++;
  }

  double nonZeros = 0;
  SparseMatrix<double, ColMajor> tmp_mat(num_row, (int) Col_added.size());
  vector<Triplet<double>> triplets;
  triplets.reserve(int((double)Col_added.size() * num_row * 0.1));
  size_t num = 0;
  for (auto &it : mat) {
    nonZeros += (double) it.nonZeros();
    for (int i = 0; i < it.rows(); ++i) {
      for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator inner_it(it, i); inner_it; ++inner_it) {
        if (if_col_added[inner_it.col()] != -1) {
          triplets.emplace_back(num, if_col_added[inner_it.col()], inner_it.value());
        }
      }
      ++num;
    }
  }

  tmp_mat.setFromTriplets(triplets.begin(), triplets.end());

  int _add= index+(int)Col_added.size();
  if (_add >= col_idx.size()) col_idx.resize(_add );
  vector<double> solver_obj(Col_added.size()), solver_val;
  vector<size_t> solver_beg(Col_added.size() + 1);
  vector<int> solver_ind;
  size_t numReserve= Col_added.size() * aver_route_length;
  solver_ind.reserve(numReserve);
  solver_val.reserve(numReserve);
  int ccnt = 0;
  for (auto col : Col_added) {
    col_idx[index++] = ptr[col];
    solver_obj[ccnt] = ptr_cost[col];
    solver_beg[ccnt] = solver_ind.size();
    for (Eigen::SparseMatrix<double, ColMajor>::InnerIterator it(tmp_mat, ccnt); it; ++it) {
	  solver_ind.emplace_back((int) it.row());
	  solver_val.emplace_back(it.value());
    }
    ++ccnt;
  }
  solver_beg[ccnt] = solver_ind.size();


  safe_solver(node->solver.XaddVars(ccnt,
                                          solver_ind.size(),
                                          solver_beg.data(),
                                          solver_ind.data(),
                                          solver_val.data(),
                                          solver_obj.data(),
                                          nullptr,
                                          nullptr,
                                          nullptr,
                                          nullptr))
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getNumCol(&num_col))
}

void CVRP::regenerateEnumMat(BbNode *node, BbNode *node2, bool if_force) {
	  auto & deleted_columns_in_enumeration_pool= node2? node2->deleted_columns_in_enumeration_pool:node->deleted_columns_in_enumeration_pool;
  int size_pool = (int) node->size_enumeration_col_pool;
  if (!size_pool) {
    node->deleted_columns_in_enumeration_pool.clear();
    return;
  }
  BbNode *out_node;
  int del_size = 0;
  if (!node2) {
    out_node = node;
    for (int i = 0; i < size_pool; ++i) {
      if (deleted_columns_in_enumeration_pool[i]) ++del_size;
    }
    node->valid_size = size_pool - del_size;
    if (!if_force) {
      auto left = double(node->size_enumeration_col_pool - del_size) / node->size_enumeration_col_pool;
      if (left > Config::LeftThresholdRCFixing4EnumerationPool) {
#if VERBOSE_MODE==1
        cout << "stash_size= " << del_size << endl;
#endif
        return;
      } else {
#if VERBOSE_MODE==1
        cout << "del_size= " << del_size << endl;
#endif
      }
    }
  } else {
    out_node = node2;
    for (int i = 0; i < size_pool; ++i) {
      if (deleted_columns_in_enumeration_pool[i]) ++del_size;
    }
    node2->valid_size = size_pool - del_size;
  }

  int colIndex = 0;
  vector<int> new_col_map(size_pool, -1);
  for (int i = 0; i < size_pool; ++i) {
    if (!deleted_columns_in_enumeration_pool[i]) {
      new_col_map[i] = colIndex++;
    }
  }

  if (cstr_index.empty()) {
    safe_solver(node->solver.optimize())
    findNonActiveCuts(node);
    if (cstr_index.empty()) {
      cstr_index.resize(num_row);
      iota(cstr_index.begin(), cstr_index.end(), 0);
    }
  }

  int new_size_pool = size_pool - del_size;
  sparseRowMatrixXd tmpMatrix(num_row, new_size_pool);

  int rowIndex = 0;
  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(num_row) * double(new_size_pool) * 0.1));

  for (auto &it : node->matrix_in_enumeration) {
    for (int i = 0; i < it.rows(); ++i) {
      if (cstr_index[rowIndex] != -1) {
        int j = cstr_index[rowIndex];
        for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator inner_it(it, i); inner_it; ++inner_it) {
          if (new_col_map[inner_it.col()] != -1) {
            triplets.emplace_back(j, new_col_map[inner_it.col()], inner_it.value());
          }
        }
      }
      ++rowIndex;
    }
  }

  tmpMatrix.setFromTriplets(triplets.begin(), triplets.end());

  out_node->matrix_in_enumeration.clear();
  out_node->matrix_in_enumeration.push_back(std::move(tmpMatrix));

  colIndex = 0;
  for (int i = 0; i < out_node->size_enumeration_col_pool; ++i) {
    if (!out_node->deleted_columns_in_enumeration_pool[i]) {
      out_node->cost_for_columns_in_enumeration_column_pool[colIndex] = out_node->cost_for_columns_in_enumeration_column_pool[i];
      out_node->index_columns_in_enumeration_column_pool[colIndex++] = out_node->index_columns_in_enumeration_column_pool[i];
    }
  }
  out_node->size_enumeration_col_pool = new_size_pool;
  out_node->cost_for_columns_in_enumeration_column_pool.conservativeResize(out_node->size_enumeration_col_pool);
  out_node->index_columns_in_enumeration_column_pool.conservativeResize(out_node->size_enumeration_col_pool);
  out_node->deleted_columns_in_enumeration_pool.assign(out_node->size_enumeration_col_pool, false);

  auto &column_pool_mapping = out_node->column_pool_mapping;
  column_pool_mapping.clear();
  auto ptr_col = out_node->index_columns_in_enumeration_column_pool;
  for (int i = 0; i < out_node->size_enumeration_col_pool; ++i) {
    int past_node = 0;
    for (auto j = ptr_col[i] + 1;; ++j) {
      int curr_node = col_pool4_pricing[j];
      if (past_node < curr_node) column_pool_mapping[{past_node, curr_node}].emplace_back(i);
      else column_pool_mapping[{curr_node, past_node}].emplace_back(i);
      if (!curr_node) break;
      past_node = curr_node;
    }
  }
  cstr_index.clear();
}

void CVRP::solveLPByInspection(BbNode *const node, bool if_only_need_value,
                               bool if_heuristic, bool if_record_sol) {
  node->is_integer = false;
  int ccnt = 0, iter_exact = 0, old_ncol = num_col;

  time_point<high_resolution_clock, duration<long long, ratio<1L, 1000000000LL>>> mt_beg,
      mt_end, spt_beg, spt_end;
  duration<long long int, ratio<1LL, 1000000000LL>> mt_elap{0}, spt_elap{0};
  spt_beg = high_resolution_clock::now();

  optimizeLPForOneIterationInEnum(node);

  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.getEnvMethod(&env_method))
  if (env_method != SOLVER_PRIMAL_SIMPLEX) {
    safe_solver(node->solver.setEnvMethod(SOLVER_PRIMAL_SIMPLEX))
    if_changed = true;
  }

  if (if_heuristic)
    cout << "Heuristic phase begin...\n";
  else
    cout << "Exact phase begin...\n";




  while (true) {

    ++iter_exact;

    spt_beg = high_resolution_clock::now();

    ccnt = generateColumnsByInspection(node, if_only_need_value);

    if (!ccnt) {
      --iter_exact;
      break;
    }

    spt_end = high_resolution_clock::now();
    spt_elap += duration_cast<milliseconds>(spt_end - spt_beg);
    mt_beg = high_resolution_clock::now();

    optimizeLPForOneIterationInEnum(node);

    mt_end = high_resolution_clock::now();
    mt_elap += duration_cast<milliseconds>(mt_end - mt_beg);

    if (!(iter_exact % PRINT_LABELING_STEP_SIZE)) {
      glo_end = high_resolution_clock::now();
      glo_eps = duration<double>(glo_end - glo_beg).count();
      printInfoLabeling(iter_exact, num_col - old_ncol, num_col, num_row, double(mt_elap.count()) * 1e-9,
                        double(spt_elap.count()) * 1e-9, glo_eps,
                        lp_val, lb, ub);
      if (!if_only_need_value) {
        cout << "col_pool= " << node->size_enumeration_col_pool << "  remain "
             << (max_num_enu_col_pool ? (double(node->size_enumeration_col_pool) / max_num_enu_col_pool * 100) : 0) << "%\n";
      }
      spt_elap = duration_values<milliseconds>::zero();
      mt_elap = duration_values<milliseconds>::zero();
      old_ncol = num_col;
    }
  }




  if (if_record_sol) {
    recordOptimalColumn(node);

  }

  if (if_changed) safe_solver(node->solver.setEnvMethod(env_method))
  if (iter_exact % PRINT_LABELING_STEP_SIZE || !iter_exact) {
    glo_end = high_resolution_clock::now();
    glo_eps = duration<double>(glo_end - glo_beg).count();
    printInfoLabeling(iter_exact, num_col - old_ncol, num_col, num_row, double(mt_elap.count()) * 1e-9,
                      double(spt_elap.count()) * 1e-9, glo_eps,
                      lp_val, lb, ub);
    if (!if_only_need_value) {
      cout << "col_pool= " << node->size_enumeration_col_pool << "  remain "
           << (max_num_enu_col_pool ? (double(node->size_enumeration_col_pool) / max_num_enu_col_pool * 100) : 0) << "%\n";
    }
  }
}

void CVRP::cleanColumnsNonElement(BbNode *const node) {
  int len = 0, keep = 1, current_node;
  bool if_break;
  yzzLong local_pi;
  vector<int> solver_ind(num_col);
  for (int i = keep; i < num_col; ++i) {
    if_break = false;
	local_pi = 0;
    for (size_t j = node->index_columns[i] + 1;; ++j) {
      current_node = col_pool4_pricing[j];
      if (!current_node) break;
      if (local_pi[current_node]) {
        solver_ind[len++] = i;
        if_break = true;
        break;
      }
	  local_pi.set(current_node);
    }
    if (!if_break) {
      node->index_columns[keep++] = node->index_columns[i];
    }
  }

  if (len) {
    safe_solver(node->solver.delVars(len, solver_ind.data()))
    safe_solver(node->solver.reoptimize())
    safe_solver(node->solver.getNumCol(&num_col))
  }

  safe_solver(node->solver.getObjVal(&lp_val))
  cout << "after clean non-ele routes lpval= " << lp_val << endl;
}

void CVRP::deleteBranchCutsAndR1C1s(BbNode *const node) {
  if (node->brcs.empty() && node->r1cs.empty()) return;
  if (!node->brcs.empty()) {
	vector<int> solver_ind(num_col);
	vector<size_t> solver_beg(num_col + 1);
	vector<double> solver_val(num_col);
    set<int> delete_col;
    int *ai_col = new int[num_col];
    int *aj_col = new int[num_col];
    int tmp;
    size_t numnzP;
    int ai, aj;
    for (auto &brc : node->brcs) {
      ai = brc.edge.first;
      aj = brc.edge.second;
      tmp = ai * dim + aj;
      if (brc.br_dir) {
        for (int i = 0; i < num_col; ++i) {
          ai_col[i] = 0;
          aj_col[i] = 0;
        }
        if (ai) {
          safe_solver(node->solver.XgetConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), ai - 1, 1))
          for (size_t i = 0; i < numnzP; ++i)ai_col[solver_ind[i]] = 1;
        }
        safe_solver(node->solver.XgetConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), aj - 1, 1))
        for (size_t i = 0; i < numnzP; ++i)aj_col[solver_ind[i]] = 1;
        safe_solver(node->solver.XgetConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), brc.idx_br_c, 1))
        for (int j = 0; j < solver_ind[0]; ++j)if (aj_col[j] || ai_col[j]) delete_col.insert(j);
        for (size_t i = 1; i < numnzP; ++i)
          for (int j = solver_ind[i - 1] + 1; j < solver_ind[i]; ++j)
            if (aj_col[j] || ai_col[j])delete_col.insert(j);
        for (int j = solver_ind[numnzP - 1] + 1; j < num_col; ++j)if (aj_col[j] || ai_col[j])delete_col.insert(j);
      } else {//delete all cstrs that coefficient == 1
        safe_solver(node->solver.XgetConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), brc.idx_br_c, 1))
        for (size_t i = 0; i < numnzP; ++i) delete_col.insert(solver_ind[i]);
      }
    }

    for (int i = 0; i < num_col; ++i) ai_col[i] = 0;
    for (auto col : delete_col) ai_col[col] = 1;
    ai_col[0] = 0;

    int len = 0, keep = 0;
    for (int i = keep; i < num_col; ++i) {
      if (ai_col[i]) solver_ind[len++] = i;
      else node->index_columns[keep++] = node->index_columns[i];
    }
    if (len) {
      safe_solver(node->solver.delVars(len, solver_ind.data()))
    }
    delete[]ai_col;
    delete[]aj_col;
  }

  vector<int> solver_ind(num_col);
  int len = 0, keep;
  auto local_cstr_index = new int[num_row];
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(num_row);
  iota(local_cstr_index, local_cstr_index + num_row, 0);
  for (auto &brc : node->brcs) {
    if (brc.idx_br_c == -1) continue;
    keep = brc.idx_br_c;
    solver_ind[len++] = keep;
	local_cstr_index[keep] = -1;
    deleted_cstrs.emplace_back(keep);
  }
  for (auto &r1c : node->r1cs) {
    if (r1c.info_r1c.size() == 1) {
      keep = r1c.idx_r1c;
      solver_ind[len++] = keep;
	  local_cstr_index[keep] = -1;
      deleted_cstrs.emplace_back(keep);
    }
  }
  if (deleted_cstrs.empty()) {
    delete[]local_cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;

  for (auto &rcc : node->rccs) rcc.idx_rcc = local_cstr_index[rcc.idx_rcc];

  for (auto i = node->r1cs.begin(); i < node->r1cs.end();) {
    if (local_cstr_index[i->idx_r1c] == -1) {
      i = node->r1cs.erase(i);
    } else {
      i->idx_r1c = local_cstr_index[i->idx_r1c];
      ++i;
    }
  }

  for (auto &r1c : node->r1cs_multi) r1c.idx_r1c = local_cstr_index[r1c.idx_r1c];

  for (auto &brc : node->brcs) brc.idx_br_c = -1;

  safe_solver(node->solver.delConstraints(len, solver_ind.data()))
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumCol(&num_col))
  safe_solver(node->solver.getNumRow(&num_row))
  safe_solver(node->solver.getObjVal(&lp_val))
  delete[]local_cstr_index;
}

void CVRP::cleanColumnsRCGreaterThanOptimalGapInLP(BbNode *const node) {
  opt_gap = calculateOptimalGap(node);
  vector<double> rc(num_col);
  vector<int> solver_ind(num_col);
  safe_solver(node->solver.getRC(0, num_col, rc.data()))
  int len = 0, keep = 1;
  for (int i = keep; i < num_col; ++i) {
    if (rc[i] > opt_gap) {
      solver_ind[len++] = i;
    } else {
      node->index_columns[keep++] = node->index_columns[i];
    }
  }
  safe_solver(node->solver.delVars(len, solver_ind.data()))
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumCol(&num_col))
}

void CVRP::organizeColumnsInMemoryToPricing(BbNode *const node) {
  int num_parent_cols = node->num_parent_cols;
  size_t start;
  pool_beg4_pricing = 0;
  for (int i = node->num_parent_cols; i < num_col; ++i) {
    start = pool_beg4_pricing;
    col_pool4_pricing[pool_beg4_pricing] = 0;
    ++pool_beg4_pricing;
    for (size_t j = node->index_columns[i] + 1;; ++j) {
      if (!(col_pool4_pricing[j]))break;
      col_pool4_pricing[pool_beg4_pricing] = col_pool4_pricing[j];
      ++pool_beg4_pricing;
    }
    col_pool4_pricing[pool_beg4_pricing] = 0;
    ++pool_beg4_pricing;
    node->index_columns[i] = start;
  }
  for (int i = 0; i < num_parent_cols; ++i) {
    start = pool_beg4_pricing;
    col_pool4_pricing[pool_beg4_pricing] = 0;
    ++pool_beg4_pricing;
    for (size_t j = node->index_columns[i] + 1;; ++j) {
      if (!(col_pool4_mem[j]))break;
      col_pool4_pricing[pool_beg4_pricing] = col_pool4_mem[j];
      ++pool_beg4_pricing;
    }
    col_pool4_pricing[pool_beg4_pricing] = 0;
    ++pool_beg4_pricing;
    node->index_columns[i] = start;
  }
}

void CVRP::recoverR1CsInEnum(BbNode *const node) {
  if (node->r1cs.empty() && node->r1cs_multi.empty()) return;
  int *ai_col = new int[num_col];
  auto *sum_col = new double[num_col];
  vector<double> val_(num_col);
  int index;
  size_t numnzP;
  vector<size_t> solver_beg(2);
  vector<int> solver_ind, solver_ind2;
  vector<double> solver_val;
  size_t numReserve= num_col*aver_route_length;
  solver_ind.reserve(numReserve);
  solver_ind2.reserve(numReserve);
  solver_val.reserve(numReserve);
  for (auto &r1c : node->r1cs) {
    index = r1c.idx_r1c;
    memset(sum_col, 0, sizeof(double) * num_col);
    for (auto i : r1c.info_r1c) {
      safe_solver(node->solver.XgetConstraints(
          &numnzP, solver_beg.data(), ai_col, val_.data(), i - 1, 1))
      for (size_t j = 0; j < numnzP; ++j) ++sum_col[ai_col[j]];
    }
    for (int i = 0; i < num_col; ++i) {
      int val = int(sum_col[i] / 2 + TOLERANCE);
      if (val) {
		solver_ind.emplace_back(index);
		solver_ind2.emplace_back(i);
		solver_val.emplace_back(val);
      }
    }
  }

  for (auto &r1c : node->r1cs_multi) {
    index = r1c.idx_r1c;
    memset(sum_col, 0, sizeof(double) * num_col);
    const auto &plan = map_rank1_multiplier[(int) r1c.info_r1c.first.size()][r1c.info_r1c.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto &i : r1c.info_r1c.first) {
      safe_solver(node->solver.XgetConstraints(
          &numnzP, solver_beg.data(), ai_col, val_.data(), i - 1, 1))
      for (size_t j = 0; j < numnzP; ++j) sum_col[ai_col[j]] += multi[count];
      ++count;
    }
    for (int i = 0; i < num_col; ++i) {
      int val = int(sum_col[i] / denominator + TOLERANCE);
      if (val) {
		solver_ind.emplace_back(index);
		solver_ind2.emplace_back(i);
		solver_val.emplace_back(val);
      }
    }
  }

  safe_solver(node->solver.XchangeCoeffs(solver_ind.size(), solver_ind.data(), solver_ind2.data(), solver_val.data()))
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getObjVal(&lp_val))
  cout << "after recover rank1c lpval= " << lp_val << endl;
  delete[]ai_col;
  delete[]sum_col;
}

void CVRP::optimizeLPForOneIterationInEnum(BbNode *const node) {

  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumCol(&num_col))
  safe_solver(node->solver.getObjVal(&lp_val))
  vector<double> X(num_col);
  safe_solver(node->solver.getX(0, num_col, X.data()))

  node->is_integer = true;
  for (int i = 0; i < num_col; ++i)
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      node->is_integer = false;
      break;
    }

  if (node->is_integer) {
    if (ceilTransformedNumberRelated(lp_val - TOLERANCE) + TOLERANCE < ub) {
      ub = ceilTransformedNumberRelated(lp_val - TOLERANCE);

      ip_opt_sol.clear();

      for (int i = 0; i < num_col; ++i)
        if (X[i] > TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(dim);
          tmp.emplace_back(0);
          for (auto j = node->index_columns[i] + 1;; ++j) {
            if (!col_pool4_pricing[j])
              break;
            tmp.emplace_back(col_pool4_pricing[j]);
          }
          tmp.emplace_back(0);
          ip_opt_sol.emplace_back(std::move(tmp));
        }

#ifdef MASTER_VALVE_ML
      updateUpperBoundEdgeSolution();
#endif

      if (ceilTransformedNumberRelated(lb_transformed - TOLERANCE) + TOLERANCE >= ub) {
        node->value = lp_val;
        node->is_terminated = true;
        cout << TERMINATED_MESSAGE_PROMISING_UPDATE_UB;
        return;
      }
    }
  }
}

void CVRP::terminateByMIP(BbNode *node) {
  cout << "b4_col= " << num_col << " ";
  regenerateEnumMat(node, nullptr, true);
  auto &deleted_columns_in_enumeration_pool = node->deleted_columns_in_enumeration_pool;
  vector<int> added_cols;
  added_cols.reserve(node->size_enumeration_col_pool);
  for (int i = 0; i < node->size_enumeration_col_pool; ++i) {
    if (deleted_columns_in_enumeration_pool[i])continue;
    added_cols.emplace_back(i);
  }
  addColumnsByInspection(node, added_cols);
  cout << "after_cols= " << num_col << endl;

  solveMIP(node, true);
  if (if_mip_enumeration_suc) {
    node->is_terminated = true;
  } else {
    cout << "MIP enumeration failed, seek for branch!\n";
    node->is_terminated = false;
    node->column_pool_mapping.clear();
    node->size_enumeration_col_pool = 0;
    node->matrix_in_enumeration.clear();
    node->cost_for_columns_in_enumeration_column_pool.resize(0);
    node->index_columns_in_enumeration_column_pool.resize(0);
  }
}



void CVRP::terminateNode(BbNode *&root_node) {
  if_in_enu_state = true;
  int root_index = root_node->index;

  cout << "To terminate the current node, " <<
       "chg node selection strategy by exploring the children nodes of the node first!\n";

  sub_bbt.push(root_node);
  --num_explored_nodes;
  for (auto &r1c : root_node->r1cs) {
    r1c.mem.clear();
  }
  for (auto &r1c : root_node->r1cs_multi) {
    r1c.mem.clear();
  }

#ifdef USE_M_DYNAMICS
  bool if_reset{false};
#endif

  while (!sub_bbt.empty()
      ) {

    auto node = sub_bbt.top();
    sub_bbt.pop();

    ++num_explored_nodes;

    safe_solver(node->solver.updateModel())
    safe_solver(node->solver.getNumRow(&num_row))
    safe_solver(node->solver.getNumCol(&num_col))

    if (ceilTransformedNumberRelated(node->value - TOLERANCE) + TOLERANCE >= ub) {
      cout << "The subtree rooted by " << root_index <<
           " has been terminated for tree_lb= " << node->value << " but ub= " << ub << endl;
      delete node;
      goto DELETE;
    }

    cout << BIG_PHASE_SEPARATION;
    glo_end = std::chrono::high_resolution_clock::now();
    glo_eps = duration<double>(glo_end - glo_beg).count();
    cout << "nd_ind= " << node->index << "  nd_col= " << num_col << "  nd_val= " << node->value << "  nd_dep= "
         << node->tree_level << "  et= " << glo_eps << "  lb= " << lb << "  ub= "
         << ub << "  sub_nd_rmn= " << sub_bbt.size() << "  ColPool= " << node->size_enumeration_col_pool << endl;
    for (auto &brc : node->brcs) {
      cout << (brc.br_dir ? "true" : "false") << "(" << brc.edge.first << "," << brc.edge.second << ")" << " ";
    }
	cout << endl;


    auto old_val = node->value;


#ifdef USE_M_DYNAMICS
    auto beg_c = std::chrono::high_resolution_clock::now();
#endif

    max_num_enu_col_pool = node->size_enumeration_col_pool;

    auto beg = high_resolution_clock::now();

    solveLPByInspection(node, false, false, true);

    auto end = high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();

    cout << "time for inspection is " << eps << endl;

    if (node->index) {
      auto &brc = node->brcs.back();
      if (brc.br_dir) {
        real_improvement_up[brc.edge].first += node->value - old_val;
        ++real_improvement_up[brc.edge].second;
      } else {
        real_improvement_down[brc.edge].first += node->value - old_val;
        ++real_improvement_down[brc.edge].second;
      }
    }

    if (node->is_terminated) {
      delete node;
      continue;
    } else if (num_col + node->valid_size <= max_num_route4_mip) {
      terminateByMIP(node);
      if (node->is_terminated) {
        delete node;
        continue;
      }
    }

	separateHybridCuts(node);
	if (node->is_terminated) {
	  delete node;
	  continue;
	} else if (num_col + node->valid_size <= max_num_route4_mip) {
	  terminateByMIP(node);
	  if (node->is_terminated) {
		delete node;
		continue;
	  }
	}

#ifdef DELUXING_APPLIED
	if(num_col + node->valid_size < 100000 && num_col + node->valid_size> 50000){
	  cout << "apply DELUXING..." << endl;
	  applyRCF(node, DELUXING_APPLIED, true);
	  solveLPByInspection(node, false, false, true);
	  if (num_col + node->valid_size <= max_num_route4_mip) {
		terminateByMIP(node);
		if (node->is_terminated) {
		  delete node;
		  continue;
		}
	  }
	}
#endif

#ifdef WRITE_ENUMERATION_TREES
    writeEnuTree(node);
    delete node;
    continue;
#endif

#ifdef USE_M_DYNAMICS
	if (node->index) {
	  auto end_c = std::chrono::high_resolution_clock::now();
	  eps = duration<double>(end_c - beg_c).count();
	  if (if_reset)
		BbNode::updateState(eps, node->c, (int)node->brcs.size() - 1);
	  else {
		node->c = eps;
		if_reset = true;
	  }
	  double new_r;
	  node->calculateRStar(node->value - old_val, new_r, this);
	  BbNode::updateState(new_r, node->geo_r_star, (int)node->brcs.size() - 1);
#if VERBOSE_MODE==1
	  cout << "node->c= " << node->c << " node->geo_r_star= " << node->geo_r_star << endl;
#endif
	}
#endif

    cout << "Begin branching...\n";

    regenerateEnumMat(node, nullptr, true);

    recordOptimalColumn(node);

    getEdgeInfo(node, true);

    writeMapEdgeColIndexInEnum(node);

    pair<int, int> info;

    doSB(node, info);

    addBranchCutToUnsolvedInEnum(node, info);

    ++num_br;
  }

  DELETE:
  while (!sub_bbt.empty()) {
    BbNode *extra_node = sub_bbt.top();
    sub_bbt.pop();
    delete extra_node;
  }

  root_node = nullptr;
  if_in_enu_state = false;
}

void
CVRP::addBranchCutToUnsolvedInEnum(BbNode *const node, const std::pair<int, int> &info) {
  if (node->is_terminated) return;
  Brc bf;
  bf.edge = info;
  ++node->tree_level;

  bf.br_dir = true;
  bf.idx_br_c = -1;

  cstr_index.resize(num_row);
  iota(cstr_index.begin(), cstr_index.end(), 0);

  auto node2 = new BbNode(node, num_col, idx_node + 2, bf);
  reviseEnumColInfoByBrC(node, node2, bf);
  safe_solver(node->solver.getNumCol(&num_col))
  bf.br_dir = false;
  node->brcs.emplace_back(bf);

  cstr_index.resize(num_row);
  iota(cstr_index.begin(), cstr_index.end(), 0);
  reviseEnumColInfoByBrC(node, node, bf);
  node->index = ++idx_node;
  ++idx_node;

  sub_bbt.push(node2);
  sub_bbt.push(node);
}

void CVRP::reviseEnumColInfoByBrC(BbNode *node, BbNode *out_node, const Brc &bf) {
  int ai = bf.edge.first, aj = bf.edge.second;
  if (ai >= aj) throw runtime_error("Wrong in reviseEnumColInfoByBrC!");
  int col_idx;
  int size_pool = int(node->size_enumeration_col_pool);
  auto &deleted_columns_in_enumeration_pool = out_node->deleted_columns_in_enumeration_pool;
  cleanColumnsInLPByNewBranchCut(out_node, bf);
  if (!node->size_enumeration_col_pool) return;
  auto &mat0 = *node->matrix_in_enumeration.begin();

  if (bf.br_dir) {//must use: 1. use and only use one edge .2 use two but not next to each other
    vector<bool> bad_cols(size_pool, true);
    for (int it : node->column_pool_mapping[{ai, aj}]) bad_cols[it] = false;
    if (ai) {
      sparseRowMatrixXd tmp = mat0.row(ai - 1) + mat0.row(aj - 1);
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
        if (it.value() > 0.5) {
          if (it.value() < 1.5) {//==1
            deleted_columns_in_enumeration_pool[it.col()] = true;
          } else {//==2 further test
            if (bad_cols[it.col()]) deleted_columns_in_enumeration_pool[it.col()] = true;
          }
        }
      }
    } else {
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat0, aj - 1); it;
           ++it) {
        if (it.value() > 0.5 && bad_cols[it.col()]) deleted_columns_in_enumeration_pool[it.col()] = true;
      }
    }
  } else {//use two and two adjacent
    for (int it : node->column_pool_mapping[{ai, aj}]) deleted_columns_in_enumeration_pool[it] = true;
  }
  regenerateEnumMat(node, out_node);
}

void CVRP::cleanColumnsInLPByNewBranchCut(BbNode *const node, const Brc &bf) {
  vector<int> delete_col;
  int *sum_col = new int[num_col]();
  int col, curr_node;
  int numnzP;
  bool if_keep;
  int ai = bf.edge.first, aj = bf.edge.second;
  if (ai >= aj) throw runtime_error("Wrong in cleanColumnsInLPByNewBranchCut!");
  vector<int> solver_beg(2);
  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
  if (!ai) {//=0
    safe_solver(node->solver.getConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), aj - 1, 1))
    if (bf.br_dir) {//一定要使用 // 找到有aj的col，然后判断其sequence，如果不是在开头或者结尾，就去除
      for (size_t i = 0; i < numnzP; ++i) {
        col = solver_ind[i];
        auto j = node->index_columns[col] + 1;
        curr_node = col_pool4_pricing[j];
        if (curr_node == aj)continue;
        ++j;
        if_keep = true;
        for (;; ++j) {
          curr_node = col_pool4_pricing[j];
          if (!curr_node) {
            if (col_pool4_pricing[j - 1] == aj) {
              if_keep = false;
            }
            break;
          }
        }
        if (!if_keep) continue;
        delete_col.emplace_back(col);
      }
    } else {//find col that traverse aj, check its sequence, if it's in the beginning or end, remove it
      for (size_t i = 0; i < numnzP; ++i) {
        col = solver_ind[i];
        auto j = node->index_columns[col] + 1;
        curr_node = col_pool4_pricing[j];
        if (curr_node == aj) {
          delete_col.emplace_back(col);
          continue;
        }
        ++j;
        if_keep = false;
        for (;; ++j) {
          curr_node = col_pool4_pricing[j];
          if (!curr_node) {
            if (col_pool4_pricing[j - 1] == aj) {
              if_keep = true;
            }
            break;
          }
        }
        if (!if_keep) continue;
        delete_col.emplace_back(col);
      }
    }
  } else {
    safe_solver(node->solver.getConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), ai - 1, 1))
    for (size_t i = 0; i < numnzP; ++i)++sum_col[solver_ind[i]];
    safe_solver(node->solver.getConstraints(&numnzP, solver_beg.data(), solver_ind.data(), solver_val.data(), aj - 1, 1))
    for (size_t i = 0; i < numnzP; ++i)++sum_col[solver_ind[i]];
    if (bf.br_dir) {//if 1 exists, remove it directly, if 2 exists, check if connect together
      for (int i = 1; i < num_col; ++i) {
        if (sum_col[i]) {
          if (sum_col[i] == 1) {
            delete_col.emplace_back(i);
          } else {//==2
            if_keep = true;
            for (auto j = node->index_columns[i] + 1;; ++j) {
              curr_node = col_pool4_pricing[j];
              if (!curr_node) break;
              if (curr_node == ai) {
                if (col_pool4_pricing[j + 1] == aj) {
                  if_keep = false;
                  break;
                }
              } else if (curr_node == aj) {
                if (col_pool4_pricing[j + 1] == ai) {
                  if_keep = false;
                  break;
                }
              }
            }
            if (!if_keep) continue;
            delete_col.emplace_back(i);
          }
        }
      }
    } else {//cstr only check if 2 exists, check if connect together
      for (int i = 1; i < num_col; ++i) {
        if (sum_col[i] == 2) {
          if_keep = false;
          for (auto j = node->index_columns[i] + 1;; ++j) {
            curr_node = col_pool4_pricing[j];
            if (!curr_node) break;
            if (curr_node == ai) {
              if (col_pool4_pricing[j + 1] == aj) {
                if_keep = true;
                break;
              }
            } else if (curr_node == aj) {
              if (col_pool4_pricing[j + 1] == ai) {
                if_keep = true;
                break;
              }
            }
          }
          if (if_keep)delete_col.emplace_back(i);
        }
      }
    }
  }

  for (int i = 0; i < num_col; ++i) sum_col[i] = 0;
  for (auto c : delete_col) sum_col[c] = 1;
  sum_col[0] = 0;//keep the first col

  int len = 0, keep = 1;
  for (int i = 1; i < num_col; ++i) {
    if (sum_col[i]) solver_ind[len++] = i;
    else node->index_columns[keep++] = node->index_columns[i];
  }

  if (len) {
    safe_solver(node->solver.delVars(len, solver_ind.data()))
    safe_solver(node->solver.reoptimize())
    safe_solver(node->solver.getNumCol(&num_col))
  }

  delete[]sum_col;
}

void CVRP::generateRCCsInEnum(BbNode *const node) {
  int numnz, cnt = 0;
  auto &ptr_col = node->index_columns_in_enumeration_column_pool;
  int curr_node;
  set<int> tmp_cutinfo;
  vector<int> solver_ind(num_col);//one time thing
  vector<double> solver_val(num_col);

  CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

  CMGR_CreateCMgr(&MyCutsCMP, 100);
  CMGR_CreateCMgr(&MyOldCutsCMP, 100);

  getEdgeInfo(node, false);
  for (int i = 1; i <= node->num_edges; ++i) if (!node->edge_tail[i]) node->edge_tail[i] = dim; else break;

  CAPSEP_SeparateCapCuts(real_dim, demand, cap, node->num_edges, node->edge_tail.data(),
                         node->edge_head.data(), node->edge_value.data(), MyOldCutsCMP,
                         MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                         &if_int_n_feasible, &max_vio, MyCutsCMP);
  for (int i = 1; i <= node->num_edges; ++i) if (node->edge_tail[i] == dim) node->edge_tail[i] = 0; else break;
  unordered_map<pair<int, int>, vector<int>, PairHasher> map_node_lp;
  int old_num_rcc = (int) node->rccs.size();
  if (!MyCutsCMP->Size) goto QUIT;

  for (int col = 0; col < num_col; ++col) {
    int past_node = 0;
    for (auto j = node->index_columns[col] + 1;; ++j) {
      curr_node = col_pool4_pricing[j];
      if (past_node < curr_node) {
        map_node_lp[{past_node, curr_node}].emplace_back(col);
      } else {
        map_node_lp[{curr_node, past_node}].emplace_back(col);
      }
      if (!curr_node) break;
      past_node = curr_node;
    }
  }

  for (int i = 0; i < MyCutsCMP->Size; ++i) {
    Rcc rcc;
    auto &tmp_customerInfo = rcc.info_rcc_customer;
    for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
      tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
    }
    if (tmp_customerInfo.size() <= dim / 2) {
      rcc.form_rcc = true;
      rcc.rhs = MyCutsCMP->CPL[i]->RHS;

      if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) {
        continue;
      }

      ++cnt;

      vector<double> col_val(num_col, 0);
      for (auto iter = tmp_customerInfo.begin(); iter != tmp_customerInfo.end(); ++iter) {
        int ai = *iter;
        auto inner_iter = iter;
        ++inner_iter;
        for (; inner_iter != tmp_customerInfo.end(); ++inner_iter) {
          int aj = *inner_iter;
          if (ai < aj) {
            for (auto it : map_node_lp[{ai, aj}]) {
              ++col_val[it];
            }
          } else {
            for (auto it : map_node_lp[{aj, ai}]) {
              ++col_val[it];
            }
          }
        }
      }
      col_val[0] = rcc.rhs;
	  numnz=0;
      for (int j = 0; j < num_col; ++j) {
        if (col_val[j] != 0) {
          solver_ind[numnz] = j;
          solver_val[numnz++] = col_val[j];
        }
      }
    } else {
      rcc.form_rcc = false;
      rcc.rhs = real_dim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;

      auto &tmp_NoCustomerInfo = rcc.info_rcc_outside_customer;
      vector<bool> tmp_supp(dim, false);
      for (int j : tmp_customerInfo) {
        tmp_supp[j] = true;
      }
      for (int j = 1; j < dim; ++j) {
        if (!tmp_supp[j]) {
          tmp_NoCustomerInfo.emplace_back(j);
        }
      }

      if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) {
        continue;
      }

      ++cnt;

      vector<double> col_val(num_col, 0);
      for (auto iter = tmp_NoCustomerInfo.begin(); iter != tmp_NoCustomerInfo.end(); ++iter) {
        int ai = *iter;
        auto inner_iter = iter;
        ++inner_iter;
        for (; inner_iter != tmp_NoCustomerInfo.end(); ++inner_iter) {
          int aj = *inner_iter;
          if (ai < aj) {
            for (auto it : map_node_lp[{ai, aj}]) {
              ++col_val[it];
            }
          } else {
            for (auto it : map_node_lp[{aj, ai}]) {
              ++col_val[it];
            }
          }
        }
      }
      for (auto customer_it : tmp_NoCustomerInfo) {
        for (auto it : map_node_lp[{0, customer_it}]) col_val[it] += 0.5;
      }
      for (auto customer_it : tmp_customerInfo) {
        for (auto it : map_node_lp[{0, customer_it}]) col_val[it] -= 0.5;
      }
      col_val[0] = rcc.rhs;
      numnz = 0;
      for (int j = 0; j < num_col; ++j) {
        if (col_val[j] != 0) {
          solver_ind[numnz] = j;
          solver_val[numnz++] = col_val[j];
        }
      }
    }
    rcc.idx_rcc = num_row;
    node->rccs.emplace_back(rcc);
    safe_solver(node->solver.addConstraint(numnz, solver_ind.data(), solver_val.data(), SOLVER_LESS_EQUAL, rcc.rhs, nullptr))
    safe_solver(node->solver.updateModel())
    safe_solver(node->solver.getNumRow(&num_row))
  }

  QUIT:
  for (int i = 0; i < MyCutsCMP->Size; ++i) CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
  MyCutsCMP->Size = 0;

  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  CMGR_FreeMemCMgr(&MyCutsCMP);
}

void CVRP::deleteNonActiveCutsSafely(BbNode *const node, int old_num) {
  /**
   * since this mode only delete cuts by slack and only new cuts will be tested,
   * therefore, the rcc can be safely deleted
   */
  int cnt = 0;
  safe_solver(node->solver.getSlack(0, num_row, slack))
  vector<int> local_cstr_index(num_row);
  iota(local_cstr_index.begin(), local_cstr_index.end(), 0);
  int delta = 0;
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(num_row);
  decltype(deleted_cstrs.begin()) stop_sign;
  vector<int> solver_ind(num_row);
  for (int i = old_num; i < num_row; ++i) {
    if (abs(slack[i]) > TOLERANCE) {
      solver_ind[cnt++] = i;
	  local_cstr_index[i] = -1;
      deleted_cstrs.emplace_back(i);
    }
  }
  if (deleted_cstrs.empty())  goto QUIT;
  std::sort(deleted_cstrs.begin(), deleted_cstrs.end());
  stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;

  safe_solver(node->solver.delConstraints(cnt, solver_ind.data()))
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumRow(&num_row))//cannot be deleted

  for (auto i = node->rccs.begin(); i < node->rccs.end();) {
    if (local_cstr_index[i->idx_rcc] == -1) {
      i = node->rccs.erase(i);
    } else {
      i->idx_rcc = local_cstr_index[i->idx_rcc];
      ++i;
    }
  }

  for (auto i = node->r1cs.begin(); i < node->r1cs.end();) {
    if (local_cstr_index[i->idx_r1c] == -1) {
      i = node->r1cs.erase(i);
    } else {
      i->idx_r1c = local_cstr_index[i->idx_r1c];
      ++i;
    }
  }
  for (auto i = node->r1cs_multi.begin(); i < node->r1cs_multi.end();) {
    if (local_cstr_index[i->idx_r1c] == -1) {
      i = node->r1cs_multi.erase(i);
    } else {
      i->idx_r1c = local_cstr_index[i->idx_r1c];
      ++i;
    }
  }

  QUIT:
  if (if_in_enu_state) {
    if (node->size_enumeration_col_pool) {
      changeEnumMatByCuts(node);
    }
  }
}

void CVRP::changeEnumMatByCuts(BbNode *node) {
  if (!node->size_enumeration_col_pool) return;

  auto &mat0 = node->matrix_in_enumeration.front();
  auto &map = node->column_pool_mapping;
  int oldNum = 0;
  for (auto &it : node->matrix_in_enumeration) oldNum += (int)it.rows();

  sparseRowMatrixXd mat(num_row - oldNum, node->size_enumeration_col_pool);
  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(num_row) * double(node->size_enumeration_col_pool) * 0.1));
  auto &deleted_ColsInEnuPool = node->deleted_columns_in_enumeration_pool;

  sparseRowMatrixXd sum(1, mat.cols());
  unordered_map<int, double> tmp;
  tmp.reserve(mat.cols());

  for (auto &rcc : node->rccs) {
    if (rcc.idx_rcc < oldNum) continue;
    tmp.clear();
    if (rcc.form_rcc) {
      auto &info = rcc.info_rcc_customer;
      for (auto iter = info.begin(); iter != info.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != info.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
    } else {
      auto &infoRccCustomer = rcc.info_rcc_customer;
      auto &infoRccOutsideCustomer = rcc.info_rcc_outside_customer;
      for (auto iter = infoRccOutsideCustomer.begin(); iter != infoRccOutsideCustomer.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != infoRccOutsideCustomer.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
      for (auto customer_it : infoRccOutsideCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] += 0.5;
      }
      for (auto customer_it : infoRccCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] -= 0.5;
      }
    }
    int row = rcc.idx_rcc - oldNum;
    for (auto &it : tmp) {
      if (abs(it.second) > TOLERANCE && !deleted_ColsInEnuPool[it.first]) {
        triplets.emplace_back(row, it.first, it.second);
      }
    }
  }

  for (auto &r1c : node->r1cs) {
    if (r1c.idx_r1c < oldNum)continue;
    sum.setZero();
    auto &info = r1c.info_r1c;
    for (auto j : info) {
      sum += mat0.row(j - 1);
    }
    sum /= 2;
    int row = r1c.idx_r1c - oldNum;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_ && !deleted_ColsInEnuPool[it.col()]) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }

  for (auto &r1c : node->r1cs_multi) {
    if (r1c.idx_r1c < oldNum)continue;
    sum.setZero();
    auto &info = r1c.info_r1c;
    const auto &plan = map_rank1_multiplier[(int) info.first.size()][info.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto &j : info.first) {
      sum += mat0.row(j - 1) * multi[count++];
    }
    sum /= denominator;
    int row = r1c.idx_r1c - oldNum;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_ && !deleted_ColsInEnuPool[it.col()]) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }
  mat.setFromTriplets(triplets.begin(), triplets.end());
  node->matrix_in_enumeration.push_back(std::move(mat));
}


void CVRP::findNonActiveCuts(BbNode *node) {
  safe_solver(node->solver.getDual(0, num_row, pi))
  vector<int> nonactive_cuts;
  nonactive_cuts.reserve(num_row);
  for (auto &rcc : node->rccs) {
#ifdef SOLVER_VRPTW
	if(rcc.if_keep) continue;
#endif
    int idx = rcc.idx_rcc;
    if (abs(pi[idx]) < TOLERANCE) {
      nonactive_cuts.emplace_back(idx);
    }
  }
  for (auto &r1c : node->r1cs) {
    int idx = r1c.idx_r1c;
    if (abs(pi[idx]) < TOLERANCE) {
      nonactive_cuts.emplace_back(idx);
    }
  }
  for (auto &r1c_mul : node->r1cs_multi) {
    int idx = r1c_mul.idx_r1c;
    if (abs(pi[idx]) < TOLERANCE) {
      nonactive_cuts.emplace_back(idx);
    }
  }
  deleteNonactiveCuts(node, nonactive_cuts);
}

void CVRP::deleteNonactiveCuts(BbNode *node, std::vector<int> &nonactive_cuts) {
  if (nonactive_cuts.empty()) return;
  int old_num_row = num_row;
  sort(nonactive_cuts.begin(), nonactive_cuts.end());
  vector<int> local_cstr_index(num_row);
  iota(local_cstr_index.begin(), local_cstr_index.end(), 0);
  int cnt = 0;
  vector<int> solver_ind(num_row);
  for (auto &cut : nonactive_cuts) {
    solver_ind[cnt++] = cut;
	local_cstr_index[cut] = -1;
  }
  int delta = 0;
  auto stop_sign = nonactive_cuts.end() - 1;
  for (auto i = nonactive_cuts.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;
  safe_solver(node->solver.delConstraints(cnt, solver_ind.data()))
  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.getEnvMethod(&env_method))
  if (env_method != SOLVER_DUAL_SIMPLEX) {
    safe_solver(node->solver.setEnvMethod(SOLVER_DUAL_SIMPLEX))
    if_changed = true;
  }
  safe_solver(node->solver.reoptimize())
  if (if_changed) {
    safe_solver(node->solver.setEnvMethod(env_method))
  }
  safe_solver(node->solver.getNumRow(&num_row))


  for (auto i = node->rccs.begin(); i < node->rccs.end();) {
    if (local_cstr_index[i->idx_rcc] == -1) {
      i = node->rccs.erase(i);
    } else {
      i->idx_rcc = local_cstr_index[i->idx_rcc];
      ++i;
    }
  }

  for (auto i = node->r1cs.begin(); i < node->r1cs.end();) {
    if (local_cstr_index[i->idx_r1c] == -1) {
      i = node->r1cs.erase(i);
    } else {
      i->idx_r1c = local_cstr_index[i->idx_r1c];
      ++i;
    }
  }

  for (auto i = node->r1cs_multi.begin(); i < node->r1cs_multi.end();) {
    if (local_cstr_index[i->idx_r1c] == -1) {
      i = node->r1cs_multi.erase(i);
    } else {
      i->idx_r1c = local_cstr_index[i->idx_r1c];
      ++i;
    }
  }

  if (if_in_enu_state) {
    cstr_index = std::move(local_cstr_index);
  } else {
    for (auto &i : node->brcs) i.idx_br_c = local_cstr_index[i.idx_br_c];
  }
}
