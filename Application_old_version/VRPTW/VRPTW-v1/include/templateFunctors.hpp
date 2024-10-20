
#ifndef CVRP_TEMPLATEFUNCTORS_HPP
#define CVRP_TEMPLATEFUNCTORS_HPP
#include "CVRP.hpp"

template<typename T, bool dir, bool if_last_half, bool if_symmetry, int heuristic_level>
int CVRP::extendKernel4Exact(BbNode*node,Label *ki,
                             int i,
							 ResTuple res,
                             const std::vector<T> &arc) {
  int state = 0;
  int bj;
  bool if_suc;

  for (auto &pair : arc) {
    int j;
    if constexpr (std::is_same<T, int>::value) {
      j = pair;
    } else {
      j = pair.second;
      res.first_res = pair.first;
    }

    updateLabel<dir, if_last_half, if_symmetry, false, false>(node,res, ki, i, j, bj,
                                                              if_suc);

    if (!if_suc) continue;

    doDominance<dir, heuristic_level>(ki, j, bj,  if_suc);

    if (!if_suc) continue;

    if constexpr (!if_last_half) {
	  if (if_symmetry){
		if(dir? node->canLeaveDepot_forward.test(j) : node->canLeaveDepot_backward.test(j)) {
		  addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
					  all_label + idx_glo,
					  nullptr,
					  Config::MaxNumRoutesInExact);
		}
	  }else{
		if(dir ? node->canLeaveDepot_backward.test(j) : node->canLeaveDepot_forward.test(j)){//should be reversed
		  addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
					  all_label + idx_glo,
					  nullptr,
					  Config::MaxNumRoutesInExact);
		}
	  }
    }

    ++idx_glo;
    if (idx_glo == label_assign) {
      rollback = ROLL_BACK_LACK_MEMORY;
      state = 2;
      goto QUIT;
    }
  }
  QUIT:

  return state;
}

template<bool dir, bool if_last_half, bool if_symmetry, bool if_std_optgap, bool if_res_updated>
void CVRP::updateLabel(BbNode* node,const ResTuple& res, Label *ki, int i, int j, int &bj,
                       bool &if_suc) {
  if_suc = false;
  if (ki->pi[j]) return;
  double &which_rc = (if_std_optgap ? opt_gap : rc_std);
  auto new_label= all_label+idx_glo;
  auto &tmp_res= new_label->res;
  auto &tmp_rc = new_label->rc;
  tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];//real rc
  if constexpr (!if_res_updated) {
    if constexpr (dir) {
	  if (!increaseMainResourceConsumption(res, tmp_res, i, j)) return;
      if constexpr (!if_last_half) {
        if (tmp_res.first_res > meet_point_resource_in_bi_dir) {
          concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_res);
          return;
        }
      }
    } else {
      if (!decreaseMainResourceConsumption(res, tmp_res, i, j)) return;
      if constexpr (!if_last_half) {
        if (tmp_res.first_res < meet_point_resource_in_bi_dir) {
          concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_res);
          return;
        }
      }
    }
  }
  if constexpr (if_last_half) {
    /**
     * here we only test the RC2TillThisBinInXXXwardSense, we do not get into details of the labels in the bins
     */
    int concate_bj;
    int if_state;
    if constexpr (if_symmetry) {
      concate_bj = int((resource.first_res - tmp_res.first_res) / step_size);
    } else {
      concate_bj = int(tmp_res.first_res / step_size);
    }
    constexpr bool if_dif = dir ^ if_symmetry;
    for (; if_dif ? concate_bj < num_buckets_per_vertex : concate_bj >= 0;
         if_dif ? ++concate_bj : --concate_bj) {
      concatenateOneLabelWithOtherLabels<dir, if_symmetry, true>(ki,
																  j,
																  concate_bj,
																  tmp_rc,
																  tmp_res,
																  if_state);
      if (if_state == -2) break;
      else if (if_state >= 0) goto outside;
    }
    outside:
    if (if_state < 0) {
      return;
    }
  }

  if_suc = true;
  bj = int(tmp_res.first_res / step_size);
  auto &tmp_PI = new_label->pi;
  tmp_PI = (ki->pi) & (ng_mem4_vertex[j]);
  tmp_PI.set(j);

  updateR1CStates(tmp_rc, new_label->r1c, ki->r1c,i ,j);
}


template<bool dir, int heuristic_level>
void CVRP::doDominance(Label *ki, int j, int bj, bool &if_suc) {
  auto &labelList_j = dir ? label_array_in_forward_sense[j][bj] : label_array_in_backward_sense[j][bj];
  auto new_label = all_label + idx_glo;
  auto &tmp_Resource = new_label->res;
  auto &tmp_rc = new_label->rc;
  auto &tmp_PI = new_label->pi;
  auto &r1c= new_label->r1c;

  double tmp_rc_add = tmp_rc + RC_TOLERANCE
  , tmp_rc_sub = tmp_rc - RC_TOLERANCE;
  if_suc = true;
#ifdef CHECK_PRICING_LABELS
  double len=0;
#endif
  auto it=labelList_j.begin();
  for(;it!=labelList_j.end();++it){
	++num_dominance_checks;
#ifdef CHECK_PRICING_LABELS
	++len;
#endif
	auto kj = *it;
	if(kj->rc < tmp_rc_add) {
	  if(dominanceCore<dir, heuristic_level>(kj,new_label)){
		HERE1:
		labelList_j.splice(labelList_j.begin(), labelList_j, it);
#ifdef CHECK_PRICING_LABELS
		inner_bin_len.first+=len;
		inner_bin_len.second++;
#endif
		if_suc = false;
		return;
	  }
	}else if(kj->rc > tmp_rc_sub){
	  if (dominanceCore<dir, heuristic_level>(new_label, kj)) {
		HERE2:
		kj->is_extended = true;
		it=labelList_j.erase(it);
		break;
	  }
	}else{
	  if(dominanceCore<dir, heuristic_level>(kj, new_label)){
		goto HERE1;
	  }else if(dominanceCore<dir, heuristic_level>(new_label, kj)){
		goto HERE2;
	  }
	}
  }

  for (; it != labelList_j.end();) {
	auto kj = *it;
	if (kj->rc < tmp_rc_add) {
	  ++it;
	  continue;
	}
	if (dominanceCore<dir, heuristic_level>(new_label, kj)) {
	  kj->is_extended = true;
	  it = labelList_j.erase(it);
	} else ++it;
  }

  labelList_j.push_front(new_label);

  new_label->p_label = ki;
  new_label->end_vertex = j;
  new_label->is_extended = false;
  auto &bucket = dir ? if_exist_extra_labels_in_forward_sense[j][bj] : if_exist_extra_labels_in_backward_sense[j][bj];
  bucket.first[bucket.second++] = new_label;
  if (bucket.second == bucket.first.size()) {
	bucket.first.resize(bucket.first.size() * 2);
  }
}



template<bool dir, int heuristic_level>
bool CVRP::dominanceCore(Label *ki, Label *kj){
  if constexpr (heuristic_level==0) {//exact
	if (dir ? tellResTupleRelations<'p'>(ki->res, kj->res) : tellResTupleRelations<'p'>(kj->res, ki->res)) return false;
	if (((ki->pi & kj->pi) ^ (ki->pi)).any())  return false;
	if (!doR1CDominance(ki->rc, kj->rc, ki->r1c, kj->r1c))  return false;
  }else if (heuristic_level==1){
	if (dir ? tellResTupleRelations<'p'>(ki->res, kj->res) : tellResTupleRelations<'p'>(kj->res, ki->res)) return false;
  }else{

  }
  return true;
}

template<bool dir, int heuristic_level>
void CVRP::checkIfDominated(Label *&ki, int i, int b,
                            bool &if_suc) {
  if_suc = true;
  double dif, rc;
  double tmp_ki_rc_sub = ki->rc - RC_TOLERANCE;
#ifdef CHECK_PRICING_LABELS
  int len = 0;
#endif
  for (int b4_b = (dir ? b - 1 : b + 1); dir ? b4_b >= 0 : b4_b < num_buckets_per_vertex; dir ? --b4_b : ++b4_b) {
	auto &b4_label_list = dir ? label_array_in_forward_sense[i][b4_b] : label_array_in_backward_sense[i][b4_b];
	if ((dir ? rc2_till_this_bin_in_forward_sense[i][b4_b] : rc2_till_this_bin_in_backward_sense[i][b4_b])
		> tmp_ki_rc_sub)
	  break;
	++num_dominance_checks;
	for (auto &p : b4_label_list) {
#ifdef CHECK_PRICING_LABELS
	  ++len;
#endif
	  if (p->rc > tmp_ki_rc_sub) break;
	  if (dominanceCore<dir, heuristic_level>(p, ki)) {
		if_suc = false;
#ifdef CHECK_PRICING_LABELS
		outer_bin_len.first += len;
		outer_bin_len.second++;
#endif
		return;
	  }
	}
  }
#ifdef CHECK_PRICING_LABELS
  outer_bin_but_keep_len.first += len;
  outer_bin_but_keep_len.second++;
#endif
}

template<bool dir>
void CVRP::sortLabelsInBinByRC(int i, int b){
  auto &label_bin = (dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b]);
  label_bin.sort(CmpLabelRCLess);
  auto &rc2_till_this_bin = (dir ? rc2_till_this_bin_in_forward_sense[i][b] : rc2_till_this_bin_in_backward_sense[i][b]);
  auto &rc2_bin = (dir ? rc2_bin_in_forward_sense[i][b] : rc2_bin_in_backward_sense[i][b]);
  rc2_bin= label_bin.empty()? LARGE_FLOAT : label_bin.front()->rc;
  if constexpr (dir) {
	rc2_till_this_bin = (b ? std::min(rc2_till_this_bin_in_forward_sense[i][b - 1], rc2_bin)
						   : rc2_bin);
  } else {
	rc2_till_this_bin = (b < num_buckets_per_vertex - 1 ? std::min(rc2_till_this_bin_in_backward_sense[i][b + 1],
																   rc2_bin)
														: rc2_bin);
  }
}

template<char type>
[[nodiscard]] bool CVRP::tellResTupleRelations(const ResTuple &res1, const ResTuple &res2) const {
  /**
   * second is priority since the first resource is sorted to some extent
   */
  if constexpr (type=='c') {
	/**
	 * concatenate: at least one is larger than max_res, return true
	 */
#ifdef USE_TWO_RESOURCE
	if(res1.second_res + res2.second_res> resource.second_res) return true;
#endif
	 if(res1.first_res + res2.first_res> resource.first_res) return true;
  }else if(type=='p') {
	/**
	 * compare: at least one is larger, return true
	 */
#ifdef USE_TWO_RESOURCE
	if(res1.second_res> res2.second_res) return true;
#endif
	 if(res1.first_res> res2.first_res) return true;
  }else if(type=='e') {
	/**
	 * equal: both are equal, return true
	 */
#ifdef USE_TWO_RESOURCE
	 if(res1.second_res==res2.second_res && res1.first_res==res2.first_res) return true;
#else
	 if(res1.first_res==res2.first_res) return true;
#endif
  }
  return false;
}

template<bool dir, bool if_symmetry, bool if_std_optgap>
void CVRP::concatenateOneLabelWithOtherLabels(Label *ki, int j, int arr_bj, double tmp_rc,const ResTuple &tmp_res,
                                              int &if_state) {
  double path_rc;
  double &which_rc = if_std_optgap ? opt_gap : rc_std;
  auto ptr_rc_till_this_bin = &rc2_till_this_bin_in_forward_sense[j][arr_bj];
  auto ptr_rc_bin = &rc2_bin_in_forward_sense[j][arr_bj];
  if constexpr ((dir && !if_symmetry) || (!dir && if_symmetry)) {
    ptr_rc_till_this_bin = &rc2_till_this_bin_in_backward_sense[j][arr_bj];
    ptr_rc_bin = &rc2_bin_in_backward_sense[j][arr_bj];
  }

  if (*ptr_rc_till_this_bin + tmp_rc > which_rc) {
    if_state = -2;//no need to continue
    return;
  }
  if_state = -1;//do some operations
  if (*ptr_rc_bin + tmp_rc < which_rc) {//most_negative_rc_in_this_bin
    auto &label_arr =
        ((!dir && !if_symmetry) || (dir && if_symmetry)) ? label_array_in_forward_sense[j][arr_bj]
                                                         : label_array_in_backward_sense[j][arr_bj];
    for (auto &kj : label_arr) {
      path_rc = kj->rc + tmp_rc;
      if (path_rc > which_rc) break;
	  if constexpr (if_symmetry) {
		if(tellResTupleRelations<'c'>(tmp_res, kj->res)) continue;
	  } else {
		if(dir? tellResTupleRelations<'p'>(tmp_res, kj->res) : tellResTupleRelations<'p'>(kj->res, tmp_res)) continue;
	  }

	  if ((ki->pi & kj->pi).any()) continue;

	  if(!concatenateR1CStates(path_rc, which_rc, ki->r1c, kj->r1c, ki->end_vertex, j)) continue;

      if constexpr (!if_std_optgap) {
        addPathByRC(path_rc, ki, kj, Config::MaxNumRoutesInExact);
      } else {
        if (path_rc < opt_gap) {
          if_state = arr_bj;//what bin to concatenate
          break;
        }
      }
    }
  }
}

template<bool dir, bool if_last_half, bool if_symmetry, int heuristic_level>
void CVRP::runLabeling(BbNode *node) {
  rollback=ROLL_BACK_INIT;
  bool if_suc;
  if constexpr (!if_last_half) {
    if constexpr (dir) {
      initializeLabels(node, 1, true, {true, 1, true});
    } else {
      initializeLabels(node, 2, false, {true, 2, false});
    }
  }
  auto beg = std::chrono::high_resolution_clock::now();
  auto end = beg;

#ifdef CHECK_PRICING_LABELS
  inner_bin_len={0,0};
  outer_bin_len={0,0};
  outer_bin_but_keep_len={0,0};
#endif

  safe_solver(node->solver.getObjVal(&lp_val))
  double gap_improvement= calculateGapImprovement(lp_val, node->value);

  double eps;
  int min_sorted_b = dir ? -1 : num_buckets_per_vertex;
  for (int b = (dir ? 0 : num_buckets_per_vertex - 1); (dir ? b < num_buckets_per_vertex : b >= 0); (dir ? ++b : --b)) {
	for(auto &comp: dir? node->topological_order_forward[b]: node->topological_order_backward[b]){
	  int index=0;
	  STILL_EXIST:
	  for(;index<comp.size();++index){
		int i=comp[index];
		auto &valid_num =
			dir ? if_exist_extra_labels_in_forward_sense[i][b].second
				: if_exist_extra_labels_in_backward_sense[i][b].second;
		if (!valid_num) continue;
		auto &label_array =
			dir ? if_exist_extra_labels_in_forward_sense[i][b].first
				: if_exist_extra_labels_in_backward_sense[i][b].first;
		for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
		  auto &ki = label_array[vec_index];
		  if (ki->is_extended) continue;
		  checkIfDominated<dir, heuristic_level>(ki, i, b, if_suc);
		  ki->is_extended = true;
		  if (!if_suc)continue;

		  auto sig = extendKernel4Exact<int,
										dir, if_last_half, if_symmetry, heuristic_level>(node, ki, i, ki->res,
																		dir
																		? node->all_forward_buckets[i][b].bucket_arcs
																		: node->all_backward_buckets[i][b].bucket_arcs);

		  if (sig == 2) goto populateBin;

		  sig = extendKernel4Exact<std::pair<double, int>,
								   dir, if_last_half, if_symmetry, heuristic_level>(node, ki, i, ki->res,
																   dir ? node->all_forward_buckets[i][b].jump_arcs
																	   : node->all_backward_buckets[i][b].jump_arcs);

		  if (sig == 2) goto populateBin;
		}
		valid_num = 0;
	  }
	  for(index=0;index<comp.size();++index){
		int i=comp[index];
		if(dir? if_exist_extra_labels_in_forward_sense[i][b].second : if_exist_extra_labels_in_backward_sense[i][b].second) {
		  goto STILL_EXIST;
		}
	  }
	}
	for(int i=1;i<dim;++i) sortLabelsInBinByRC<dir>(i,b);
	min_sorted_b= b;
	if(node->r1cs.empty()) continue;
	end = std::chrono::high_resolution_clock::now();
	eps = std::chrono::duration<double>(end - beg).count();
	if(heuristic_level==0){
	  if constexpr (if_last_half) {
		if (eps > arc_elimination_time) {
		  rollback = ROLL_BACK_LONG_TIME;
		  goto QUIT;
		}
	  } else {
		if (eps > soft_time) {
		  rollback = ROLL_BACK_TAIL_OFF;
		  if(!force_not_rollback){
			if (eps > Config::HardTimeThresholdInPricing ||
				((last_max_time_labeling> eps) && (gap_improvement/ (last_max_time_labeling-eps) < Config::RatioGapImprovedVSTimeIncreased_LB))) {
			  rollback = ROLL_BACK_LONG_TIME;
			  goto QUIT;
			}
		  }
		}
	  }
	}
  }

  populateBin:
  if constexpr (!if_last_half) {
    end = std::chrono::high_resolution_clock::now();
    eps = std::chrono::duration<double>(end - beg).count();//overall time
    last_max_time_labeling = std::max(last_max_time_labeling, eps);
  }

  for(int i=1;i<dim;++i){
	for (int b = dir ? min_sorted_b + 1 : min_sorted_b - 1; dir ? b < num_buckets_per_vertex : b >= 0;
		 dir ? ++b : --b) {
	  sortLabelsInBinByRC<dir>(i, b);
	}
  }

  QUIT:
  if (rollback == ROLL_BACK_LONG_TIME) {
    std::cout << "rollback to original states!" << std::endl;
  } else if (rollback == ROLL_BACK_LACK_MEMORY) {
    std::cout << "rollback with larger mem!" << std::endl;
  }
#ifdef CHECK_PRICING_LABELS
  using  namespace  std;
  if(if_exact_labeling_cg) {
	cout << "inner= " << inner_bin_len.first / inner_bin_len.second << " abs= " << inner_bin_len.second << " " <<
		 "outer= " << outer_bin_len.first / outer_bin_len.second << " abs= " << outer_bin_len.second << " " <<
		 "outer2= " << outer_bin_but_keep_len.first / outer_bin_but_keep_len.second << " abs= "
		 << outer_bin_but_keep_len.second << endl;
	cout<<"num_dominance_checks= "<<num_dominance_checks<<" eps= "<<eps<<endl;
  }
#endif
}

template<bool dir, bool if_symmetry>
void CVRP::concatenatePhaseInArcElimination(BbNode*node) {
  bool if_find;
  int concate_bj;
  int if_state;
  int bj;
  bool if_suc;
  auto beg = std::chrono::high_resolution_clock::now();
  auto end = beg;
  double eps;
  constexpr bool if_dif = dir ^ if_symmetry;
  for (auto &label_list : dir ? concatenate_labels_in_forward_cg : concatenate_labels_in_backward_cg) {
    int i = label_list.first.first;
    int j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
	  auto &tmp_Resource = all_label[idx_glo].res;
	  tmp_Resource = pr.second;//do not change here! since res will not be updated later
      updateLabel<dir, true, if_symmetry, true, true>(node, tmp_Resource, ki, i, j, bj,if_suc);
      if (!if_suc) continue;

      doDominance<dir, 0>(ki, j, bj,  if_suc);
      if (!if_suc) continue;

      ++idx_glo;
      if (idx_glo == label_assign) {
        rollback = ROLL_BACK_LACK_MEMORY;
        goto QUIT;
      }
    }
  }
  if(!node->r1cs.empty()) {
	end = std::chrono::high_resolution_clock::now();
	eps = std::chrono::duration<double>(end - beg).count();
	if (eps > Config::HardTimeThresholdInArcEliminationMidConcatenate) {
	  rollback = ROLL_BACK_LONG_TIME;
#if VERBOSE_MODE == 1
	  std::cout << "concatenatePhaseInArcElimination time= " << eps << " failed!" << std::endl;
#endif
	}
  }
  QUIT:
  return;
}

template<bool dir, bool if_symmetry>
void CVRP::eliminateBucketArcs(BbNode *node,
                              int dim_sq,
                              bool *stateBetween2Buckets,
                              int *latest_bucket) {
  memset(stateBetween2Buckets, 0, dim_sq * num_buckets_per_vertex * sizeof(bool));
  constexpr bool if_dif = dir ^ if_symmetry;
  if constexpr (if_dif) {
    memset(latest_bucket, -1, dim_sq * sizeof(int));
  } else {
    std::fill_n(latest_bucket, dim_sq, num_buckets_per_vertex);
  }
  int num_bucket_arcs = 0;
  for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
    for (int i = 1; i < dim; ++i) {
      concatenateTestKernelInArcElimination<std::pair<double, int>, dir, if_symmetry>(i,
                                                                                      b,
                                                                                      (dir
                                                                                       ? node->all_forward_buckets[i][b].jump_arcs
                                                                                       : node->all_backward_buckets[i][b].jump_arcs),
                                                                                      dim_sq,
                                                                                      stateBetween2Buckets,
                                                                                      latest_bucket);
      concatenateTestKernelInArcElimination<int, dir, if_symmetry>(i,
                                                                   b,
                                                                   (dir
                                                                    ? node->all_forward_buckets[i][b].bucket_arcs
                                                                    : node->all_backward_buckets[i][b].bucket_arcs),
                                                                   dim_sq,
                                                                   stateBetween2Buckets,
                                                                   latest_bucket);
    }
  }

  std::vector<int> tmp_vec;
  tmp_vec.reserve(dim);
  int map1 = 0;
  for (int b = 0; b < num_buckets_per_vertex; ++b, map1 += dim_sq) {
    int map2 = map1 + dim;
    for (int i = 1; i < dim; ++i, map2 += dim) {
      tmp_vec.clear();
      for (int j : dir ? node->all_forward_buckets[i][b].bucket_arcs : node->all_backward_buckets[i][b].bucket_arcs) {
        if (stateBetween2Buckets[map2 + j]) {
          tmp_vec.emplace_back(j);
          ++num_bucket_arcs;
        }
      }
      (dir ? node->all_forward_buckets[i][b].bucket_arcs : node->all_backward_buckets[i][b].bucket_arcs) = tmp_vec;
    }
  }

  eliminatebuketArc4Depot<dir, if_symmetry>(node);
#if VERBOSE_MODE==1
  std::cout << "Num of " << (dir ? "Forward" : "Backward") << " bucket_arcs= " << num_bucket_arcs << " prev.= "
            << double(num_bucket_arcs) / (dir ? node->num_forward_bucket_arcs : node->num_backward_bucket_arcs) * 100 << "%"
            << " max.= " << double(num_bucket_arcs) / max_num_forward_graph_arc * 100 << "%"
            << std::endl;
#endif
  (dir ? node->num_forward_bucket_arcs : node->num_backward_bucket_arcs) = num_bucket_arcs;
  if constexpr (dir){
	node->canLeaveDepot_forward.reset();
	for (auto j : node->all_forward_buckets[0][0].bucket_arcs)  node->canLeaveDepot_forward.set(j);
  }else{
	node->canLeaveDepot_backward.reset();
	for (auto j : node->all_backward_buckets[0][0].bucket_arcs)  node->canLeaveDepot_backward.set(j);
  }
  getTopologicalOrder(node);
}


template<bool dir, bool if_symmetry>
void CVRP::eliminatebuketArc4Depot(BbNode *node) {
  bool constexpr if_dif = dir ^ if_symmetry;
  auto &allBuckets = dir ? node->all_forward_buckets : node->all_backward_buckets;
  auto &RC2TillThisBin = if_dif ? rc2_till_this_bin_in_backward_sense : rc2_till_this_bin_in_forward_sense;
  auto &RC2Bin = if_dif ? rc2_bin_in_backward_sense : rc2_bin_in_forward_sense;
  auto &labelArray = if_dif ? label_array_in_backward_sense : label_array_in_forward_sense;

  std::unordered_set<int> depot_set;
  for (auto &arc : allBuckets[0][0].bucket_arcs) depot_set.emplace(arc);

  for (int i = 1; i < dim; ++i) {
    if (depot_set.find(i) == depot_set.end()) continue;
    bool if_delete = true;

    for (int b = (if_dif ? 0 : (num_buckets_per_vertex - 1)); if_dif ? (b < num_buckets_per_vertex) : (b >= 0);
         if_dif ? ++b : --b) {
      if (RC2TillThisBin[i][b] + chg_cost_mat4_vertex[i][0] > opt_gap) break;
      if (RC2Bin[i][b] + chg_cost_mat4_vertex[i][0] > opt_gap) continue;

      if (!labelArray[i][b].empty()) {
        auto &labels = labelArray[i][b];
		for(auto &label: labels){
          if (label->rc + chg_cost_mat4_vertex[i][0] > opt_gap) break;
          auto tmp_res = label->res;//cannot use &
          if (if_dif ? (decreaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) : (
              increaseMainResourceConsumption(tmp_res, tmp_res, i, 0))) {
            if_delete = false;
            break;
          }
        }
      }
    }

    if (if_delete) {
      depot_set.erase(i);
    }
  }

  if (depot_set.size() != allBuckets[0][0].bucket_arcs.size()) {
    allBuckets[0][0].bucket_arcs.clear();
    allBuckets[0][0].bucket_arcs.assign(depot_set.begin(), depot_set.end());
    std::sort(allBuckets[0][0].bucket_arcs.begin(), allBuckets[0][0].bucket_arcs.end());
  }
}

template<typename T, bool dir, bool if_symmetry>
void CVRP::concatenateTestKernelInArcElimination(int i,
                                                 int b,
                                                 const std::vector<T> &arc,
                                                 int dim_sq,
                                                 bool *stateBetween2Buckets,
                                                 int *latest_bucket) {
  ResTuple tmp_Resource;
  double tmp_rc;
  int if_state, arr_bj;
  for (auto &pr : arc) {
    int j;
    if constexpr (std::is_same<T, int>::value) {
      j = pr;
    } else {
      j = pr.second;
	  ResTuple res;
	  res.first_res=pr.first;
#ifdef USE_TWO_RESOURCE
	  dir? res.second_res=0 : res.second_res=resource.second_res;
#endif
      if (dir ? !increaseMainResourceConsumption(res, tmp_Resource, i, j) :
          !decreaseMainResourceConsumption(res, tmp_Resource, i, j))
        continue;
    }
    int map = i * dim + j;
    int &latest = latest_bucket[map];
    int old_latest = latest;
    auto &label_array = dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b];
    constexpr bool if_dif = dir ^ if_symmetry;
    for (auto ki: label_array) {
      if constexpr (std::is_same<T, int>::value) {
        if (dir ? !increaseMainResourceConsumption(ki->res, tmp_Resource, i, j) :
            !decreaseMainResourceConsumption(ki->res, tmp_Resource, i, j))
          continue;
      }
      tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
      if constexpr (dir) {
        arr_bj =
            if_symmetry ? std::min(int((resource.first_res - tmp_Resource.first_res) / step_size), latest - 1) : std::max(int(
                (tmp_Resource.first_res) / step_size), latest + 1);
      } else {
        arr_bj =
            if_symmetry ? std::max(int((resource.first_res - tmp_Resource.first_res) / step_size), latest + 1) : std::min(int(
                (tmp_Resource.first_res) / step_size), latest - 1);
      }

      for (int bj = if_dif ? num_buckets_per_vertex - 1 : 0; if_dif ? bj >= arr_bj : bj <= arr_bj;
           if_dif ? --bj : ++bj) {
        concatenateOneLabelWithOtherLabels<dir, if_symmetry, true>(ki,
																   j,
																   bj,
																   tmp_rc,
																   tmp_Resource,
																   if_state);
        if (if_state >= 0) {
          latest = if_state;
          goto outside;
        }
      }
      outside:
      continue;
    }
    if (latest != old_latest) {
      int bi, chg_latest;
      if constexpr (if_symmetry) {
        chg_latest = int((resource.first_res - latest * step_size) / step_size);
        if (chg_latest < 0) chg_latest = 0;
        else if (chg_latest >= num_buckets_per_vertex) chg_latest = num_buckets_per_vertex - 1;
      } else chg_latest = latest;
      bi = (dir ? tell_which_bin4_arc_elimination_in_forward_sense[map + chg_latest * dim_sq]
                : tell_which_bin4_arc_elimination_in_backward_sense[map + chg_latest * dim_sq]);
      int map2 = map + b * dim_sq;
      for (int k = b; dir ? k <= bi : k >= bi; dir ? ++k : --k, dir ? map2 += dim_sq : map2 -= dim_sq)
        stateBetween2Buckets[map2] = true;
    }
  }
}

template<bool dir>
void CVRP::populateRC2TillThisBinNRC2Bin() {
  for (int i = 1; i < dim; ++i) {
	for(int b=0;b<num_buckets_per_vertex;++b){
	  sortLabelsInBinByRC<dir>(i,b);
	}
  }
}

template<bool if_symmetry>
int CVRP::generateColsByBidir(BbNode *node) {
  if_exact_labeling_cg=true;
  priceLabeling(node, pi4_labeling);
  int ccnt=0;
  double NumExistedLabels = 0;
  double NumExistedLabel_back = 0;
  RE_TRY:
  runLabeling<true, false, if_symmetry,0>(node);

  if (rollback == ROLL_BACK_LACK_MEMORY){
    reallocateLabel();
    goto RE_TRY;
  } else if (rollback == ROLL_BACK_LONG_TIME){
    goto QUIT;
  }

  if (!if_symmetry) {
    RE_TRY2:
    runLabeling<false, false, if_symmetry,0>(node);//ccnt can only be applied for forward case

    if (rollback == ROLL_BACK_LACK_MEMORY){
      reallocateLabel();
      goto RE_TRY2;
    } else if (rollback == ROLL_BACK_LONG_TIME){
	  goto QUIT;
    }
  }

  ccnt = concatenateCols_prior_forward<if_symmetry>(node);

  addColumns(node, ccnt);

  for (int i = 1; i < dim; ++i) {
	for (int b = 0; b < num_buckets_per_vertex; ++b) {
	  NumExistedLabels += (double ) label_array_in_forward_sense[i][b].size();
	  if constexpr (!if_symmetry) {
		NumExistedLabel_back += (double ) label_array_in_backward_sense[i][b].size();
	  }
	}
  }
  ratio_dominance_checks_non_dominant.first +=
	  if_symmetry ? num_dominance_checks / NumExistedLabels : num_dominance_checks
		  / (NumExistedLabels + NumExistedLabel_back);
  ++ratio_dominance_checks_non_dominant.second;
  if constexpr (!if_symmetry) {
	double dif = abs(NumExistedLabels - NumExistedLabel_back);
	double over = dif / std::min(NumExistedLabels, NumExistedLabel_back);

#ifdef FIXED_MEET_POINT_FOR_RESOURCE
	if(node->r1cs.empty()) {
#endif
	  if (over > Config::NumberOfOverLabelsInMeetPoint) {
		if (NumExistedLabels > NumExistedLabel_back) {
		  meet_point_resource_in_bi_dir *= (1 - Config::MeetPointFactor);
		} else {
		  meet_point_resource_in_bi_dir *= (1 + Config::MeetPointFactor);
		}
	  }
#ifdef FIXED_MEET_POINT_FOR_RESOURCE
	}
#endif
  }

  if (ccnt) {
    smallest_rc = std::get<2>(negative_rc_label_tuple[0]);
	if(if_exact_labeling_cg){
	  if (ceilTransformedNumberRelated(smallest_rc * max_num_vehicle + lp_val + RC_TOLERANCE) + TOLERANCE >= ub) {
		node->value = lp_val;
		node->is_terminated = true;
		std::cout << TERMINATED_MESSAGE_PROMISING_VEHICLES;
		goto QUIT;
	  }
	}
  } else {
    smallest_rc = 0;
  }
  if (!ccnt) {
	optimal_dual_vector = pi4_labeling;
	goto QUIT;
  }

  QUIT:
  return ccnt;
}
template<bool if_symmetry>
int CVRP::concatenateCols_prior_forward(BbNode *node) {
  int index = num_col - 1;
  double tmp_rc;
  int i, j, arr_bj;
  Label *p;
  int if_state;//==0 means work

  for (auto &label_list : concatenate_labels_in_forward_cg) {
    i = label_list.first.first;
    j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      auto &tmp_Resource = pr.second;
      arr_bj =
          if_symmetry ? int((resource.first_res - tmp_Resource.first_res) / step_size) : (int) (tmp_Resource.first_res / step_size);
      tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];

      for (; if_symmetry ? arr_bj >= 0 : arr_bj < num_buckets_per_vertex;
           if_symmetry ? --arr_bj : ++arr_bj) {
        concatenateOneLabelWithOtherLabels<true, if_symmetry, false>(ki,
																	  j,
																	  arr_bj,
																	  tmp_rc,
																	  tmp_Resource,
																	  if_state);
        if (if_state == -2)break;
      }
    }
  }

  writeColumnsInPricingPool(node, index);

  return index - num_col + 1;
}

template<bool dir>
void CVRP::populateTellWhichBin4ArcElimination() {
  /**
   * one time calculation except for regenerate the bucket graph
   */
  if constexpr (dir) {
    if (!tell_which_bin4_arc_elimination_in_forward_sense.empty()) return;
  } else {
    if (!tell_which_bin4_arc_elimination_in_backward_sense.empty()) return;
  }
#if VERBOSE_MODE==1
  std::cout << "populateTellWhichBin4ArcElimination" << std::endl;
#endif
  size_t size = dim * dim * num_buckets_per_vertex;
  int dim_sq = dim * dim;
  if constexpr (dir) tell_which_bin4_arc_elimination_in_forward_sense.reserve(size);
  else tell_which_bin4_arc_elimination_in_backward_sense.reserve(size);
  std::unordered_map<int, int> map_bj_B;
  std::vector<std::pair<int, int>> vec_bj_B(num_buckets_per_vertex);
  map_bj_B.reserve(num_buckets_per_vertex + 1);
  ResTuple tmp_Resource, res_tuple;
  for (int b = 0; b <= num_buckets_per_vertex; ++b) map_bj_B[b] = dir ? -1 : num_buckets_per_vertex;
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      int base = i * dim + j;
      for (auto &key_value : map_bj_B) key_value.second = dir ? -1 : num_buckets_per_vertex;
      for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
        auto res = dir ? b * step_size : std::min((b + 1) * step_size, resource.first_res);
		res_tuple.first_res=res;
#ifdef USE_TWO_RESOURCE
		dir? res_tuple.second_res=0 : res_tuple.second_res=resource.second_res;
#endif
        if (dir ? !increaseMainResourceConsumption(res_tuple, tmp_Resource, i, j) :
            !decreaseMainResourceConsumption(res_tuple, tmp_Resource, i, j))
          break;
        int bj = int(tmp_Resource.first_res / step_size);
        if (dir ? map_bj_B[bj] < b : map_bj_B[bj] > b) map_bj_B[bj] = b;
      }
      for (int bj = dir ? num_buckets_per_vertex - 1 : 0; dir ? bj >= 0 : bj < num_buckets_per_vertex; dir ? --bj : ++bj) {
        if (map_bj_B[bj] != (dir ? -1 : num_buckets_per_vertex)) {
          for (int bj2 = dir ? bj + 1 : bj - 1; dir ? bj2 < num_buckets_per_vertex : bj2 >= 0; dir ? ++bj2 : --bj2)
            map_bj_B[bj2] = map_bj_B[bj];
          break;
        }
      }
      for (int bj = 0; bj < num_buckets_per_vertex; ++bj) {
        (dir ? tell_which_bin4_arc_elimination_in_forward_sense[base + bj * dim_sq] :
		 tell_which_bin4_arc_elimination_in_backward_sense[base + bj * dim_sq]) = map_bj_B[bj];
      }
    }
  }
}

template<bool dir>
void CVRP::obtainjumpArcs(BbNode *node, std::bitset<2> **bitMap) const {
  int num_jump_arcs = 0;
  bool if_used;

  for (int i = 1; i < dim; ++i) {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      if constexpr (dir) {
        node->all_forward_buckets[i][b].jump_arcs.clear();
      } else {
        node->all_backward_buckets[i][b].jump_arcs.clear();
      }
      for (int j = 1; j < dim; ++j) bitMap[j][b] = 2;
      bitMap[i][b] = 1;
      for (int j : (dir ? node->all_forward_buckets[i][b].bucket_arcs :
                    node->all_backward_buckets[i][b].bucket_arcs))
        bitMap[j][b] = 0;
    }
    for (int b = (dir ? 0 : num_buckets_per_vertex - 1); dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
      for (int j = 1; j < dim; ++j) {
        if (bitMap[j][b] == 2) {//need a jump arc
          if_used = false;
          for (int b4_i = (dir ? b + 1 : b - 1); dir ? b4_i < num_buckets_per_vertex : b4_i >= 0;
               dir ? ++b4_i : --b4_i) {
            if (bitMap[j][b4_i] == 0) {//find a going arc
              std::pair<int, int> map = (dir ? std::make_pair(b4_i * step_size, j) :
                                         std::make_pair((b4_i + 1) * step_size, j));
              for (int tmp_b = b; dir ? tmp_b < b4_i : tmp_b > b4_i; dir ? ++tmp_b : --tmp_b) {//不取等号< qi
                bitMap[j][tmp_b] = 1;
                if constexpr (dir) {
                  node->all_forward_buckets[i][tmp_b].jump_arcs.emplace_back(map);
                } else {
                  node->all_backward_buckets[i][tmp_b].jump_arcs.emplace_back(map);
                }
                ++num_jump_arcs;
              }
              if_used = true;
              break;
            }
          }
          if (!if_used) {
            for (int tmp_b = b; dir ? tmp_b < num_buckets_per_vertex : tmp_b >= 0; dir ? ++tmp_b : --tmp_b)
              bitMap[j][tmp_b] = 1;
          }
        }
      }
    }
  }

  (dir ? node->num_forward_jump_arcs : node->num_backward_jump_arcs) = num_jump_arcs;
#if VERBOSE_MODE==1
  std::cout << "Obtain" << (dir ? "Forward" : "Backward") << " Jump Arcs= " << num_jump_arcs << std::endl;
#endif
}

template<bool dir, bool if_symmetry>
void CVRP::extendKernel4Enumeration(BbNode*node, int i, int b, Label*ki, std::vector<Label *> **copy_bucket,
									std::unordered_map<yzzLong, std::tuple<Label *, Label *, double>> &Tags,
									int &num_routes_now, int &status) {
  status=0;
  constexpr bool  if_dif = dir ^ if_symmetry;
  for (int j : (dir ? node->all_forward_buckets[i][b].bucket_arcs :
                      node->all_backward_buckets[i][b].bucket_arcs)) {
	if (ki->pi[j]) continue;
	auto new_label = all_label+idx_glo;
	auto &tmp_Resource = new_label->res;
	if (dir ? !increaseMainResourceConsumption(ki->res, tmp_Resource, i, j) :
		!decreaseMainResourceConsumption(ki->res, tmp_Resource, i, j))
	  continue;
	auto &tmp_rc = new_label->rc;
	tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];//real rc
	auto if_keep = false;
	int arr_bj = (if_symmetry ? int((resource.first_res- tmp_Resource.first_res) / step_size) :
				  int(tmp_Resource.first_res / step_size));
	if constexpr (dir) {
	  if (tmp_Resource.first_res > meet_point_resource_in_bi_dir_enu) {
		if constexpr (if_symmetry) {
		  if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap)
			concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_Resource);
		} else {
		  if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap)
			concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_Resource);
		}
		continue;
	  }
	} else {
	  if (tmp_Resource.first_res < meet_point_resource_in_bi_dir_enu) {
		if constexpr (if_symmetry) {
		  if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap)
			concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_Resource);
		} else {
		  if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap)
			concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_Resource);
		}
		continue;
	  }
	}

	for (; if_dif ? arr_bj < num_buckets_per_vertex : arr_bj >= 0;
		 if_dif ? ++arr_bj : --arr_bj) {
	  if (tmp_rc + (if_dif ? rc2_till_this_bin_in_backward_sense[j][arr_bj] : rc2_till_this_bin_in_forward_sense[j][arr_bj])
		  > opt_gap)
		break;
	  if (tmp_rc + (if_dif ? rc2_bin_in_backward_sense[j][arr_bj] : rc2_bin_in_forward_sense[j][arr_bj])
		  > opt_gap)
		continue;
	  for (auto &kkj : copy_bucket[j][arr_bj]) {
		if constexpr (if_symmetry) {
		  if(tellResTupleRelations<'c'>(tmp_Resource, kkj->res)) continue;
		} else {
		  if constexpr (dir) {
			if(tellResTupleRelations<'p'>(tmp_Resource, kkj->res)) continue;
		  } else {
			if(tellResTupleRelations<'p'>(kkj->res, tmp_Resource)) continue;
		  }
		}
		double path_rc = tmp_rc + kkj->rc;

		if (path_rc > opt_gap) break;

		if ((ki->pi & kkj->pi).any()) continue;
		if(!concatenateR1CStates(path_rc, opt_gap, ki->r1c, kkj->r1c, i, j)) continue;

		if_keep = true;
		goto OUT;
	  }
	}

	OUT:
	if(!if_keep) continue;

	int bj;
	updateEnumerationLabel<dir>(ki, i, j, bj);

	doDominanceEnumerationLabel<dir>(ki, i, j, bj, if_keep);

	if(!if_keep) continue;

	updateR1CStates(new_label->rc, new_label->r1c, ki->r1c, i, j);
	new_label->p_label= ki;
	new_label->end_vertex = j;
	new_label->is_extended = false;
	auto &bucket = (dir ? if_exist_extra_labels_in_forward_sense[j][bj] : if_exist_extra_labels_in_backward_sense[j][bj]);
	bucket.first[bucket.second++] = new_label;
	if (bucket.second == bucket.first.size()) {
	  bucket.first.resize(bucket.first.size() * 2);
	}

	if constexpr (dir) {
	  if (tmp_rc + chg_cost_mat4_vertex[j][0] < opt_gap) {
		auto path_cost = new_label->cost + cost_mat4_vertex[j][0];
		auto &tmp_PI = new_label->pi;
		if (Tags.find(tmp_PI) == Tags.end()) {
		  Tags[tmp_PI] = {all_label + idx_glo, nullptr, path_cost};
		  ++num_routes_now;
		  if (num_routes_now > Config::MaxNumRouteInEnumeration_half) {
			status = 2;
			goto HERE;
		  }
		} else if (std::get<2>(Tags[tmp_PI]) > path_cost) {
		  Tags[tmp_PI] = {all_label + idx_glo, nullptr, path_cost};
		}
	  }
	}

	if ((dir ? ++num_forward_labels_in_enu : ++num_backward_labels_in_enu) > Config::MaxNumLabelInEnumeration) {
	  status = 3;
	  goto HERE;
	}
	++idx_glo;//can be put here, because once go outside, the function will end
	if (idx_glo == label_assign) {
	  status = 1;
	  rollback = ROLL_BACK_LACK_MEMORY;
	  goto HERE;
	}
  }
  HERE:;
}

template<bool dir>
bool CVRP::dominanceCoreInEnumeration(Label *ki, Label *kj) {
  if(dir? tellResTupleRelations<'p'>(ki->res, kj->res) : tellResTupleRelations<'p'>(kj->res, ki->res)) return false;
  if((ki->pi ^ kj->pi).any()) return false;
  return true;
}

template<bool dir>
void CVRP::doDominanceEnumerationLabel(Label*ki, int i, int j, int bj, bool &if_suc) {
  auto new_label= all_label+idx_glo;
  auto &tmp_Resource = new_label->res;
  auto &tmp_PI = new_label->pi;
  auto &tmp_Cost = new_label->cost;
  double tmp_Cost_sub = tmp_Cost - TOLERANCE, tmp_Cost_add = tmp_Cost + TOLERANCE;
  auto &labelList_j = dir ? label_array_in_forward_sense[j][bj] : label_array_in_backward_sense[j][bj];

  if_suc=true;
  auto it= labelList_j.begin();
  for (; it!=labelList_j.end(); ++it) {
	auto kj= *it;
	if(kj->cost < tmp_Cost_sub){
	  if(dominanceCoreInEnumeration<dir>(kj, new_label)){
		HERE1:
		if_suc = false;
		labelList_j.splice(labelList_j.begin(), labelList_j, it);
		return;
	  }
	} else if(kj->cost > tmp_Cost_add) {
	  if (dominanceCoreInEnumeration<dir>(new_label, kj)) {
		HERE2:
		kj->is_extended = true;
		it = labelList_j.erase(it);
		dir ? --num_forward_labels_in_enu : --num_backward_labels_in_enu;
		break;
	  }
	}else{
	  if(dominanceCoreInEnumeration<dir>(kj, new_label)){
		goto HERE1;
	  }else if(dominanceCoreInEnumeration<dir>(new_label, kj)){
		goto HERE2;
	  }
	}
  }

  for(;it!=labelList_j.end();++it){
	auto kj= *it;
	if(kj->cost < tmp_Cost_sub) continue;
	if(dominanceCoreInEnumeration<dir>(new_label, kj)){
	  kj->is_extended = true;
	  it = labelList_j.erase(it);
	  dir ? --num_forward_labels_in_enu : --num_backward_labels_in_enu;
	}
  }

  labelList_j.push_front(new_label);
}

template<bool dir>
void CVRP::updateEnumerationLabel(Label *ki, int i, int j, int &bj){
  auto new_label= all_label+idx_glo;
  auto &tmp_Resource= new_label->res;
  bj = int(tmp_Resource.first_res / step_size);
  auto &labelList_j = dir ? label_array_in_forward_sense[j][bj] : label_array_in_backward_sense[j][bj];
  auto &tmp_PI = new_label->pi;
  auto &tmp_Cost = new_label->cost;
  tmp_PI = ki->pi;
  tmp_PI.set(j);
  tmp_Cost = ki->cost + cost_mat4_vertex[i][j];
};

template<bool dir, bool if_symmetry>
int CVRP::enumerateHalfwardRoutes(BbNode *node,
                                  std::unordered_map<yzzLong, std::tuple<Label *, Label *, double>> &Tags,
                                  std::vector<Label *> **copy_bucket,
                                  int &num_routes_now) {
  using namespace std::chrono;
  int status;
  (dir ? num_forward_labels_in_enu : num_backward_labels_in_enu) = 0;
  constexpr bool if_dif = dir ^ if_symmetry;
  if constexpr (dir) {
    initializeLabels(node, 1, false, {true, 1, true});
  } else {
    initializeLabels(node, 2, false, {true, 2, false});
  }

  double left_time = Config::HardTimeThresholdInAllEnumeration;
  auto initial_time= high_resolution_clock::now();
  time_point<high_resolution_clock> beg, end;

  int stop_b= (int)meet_point_resource_in_bi_dir_enu/step_size;

  for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
	double time_4_b;
	if(dir? b< stop_b : b>stop_b) {
	  beg=high_resolution_clock::now();
	  left_time= Config::HardTimeThresholdInAllEnumeration-(duration<double>(beg-initial_time).count());
	  time_4_b= left_time /abs(stop_b-b);
	}else time_4_b= std::numeric_limits<double>::max();

   for(auto &comp: dir? node->topological_order_forward[b]: node->topological_order_backward[b]){
	  int index=0;
	  STILL_EXIST:
	 end=high_resolution_clock::now();
	 if(duration<double>(end-beg).count()>time_4_b) {
	   status = 2;
	   goto outside;
	 }
	  for(;index<comp.size();++index){
		int i=comp[index];
      auto &valid_num = (dir ? if_exist_extra_labels_in_forward_sense[i][b].second :
                         if_exist_extra_labels_in_backward_sense[i][b].second);
      if (!valid_num) continue;
      auto &label_array = (dir ? if_exist_extra_labels_in_forward_sense[i][b].first :
                           if_exist_extra_labels_in_backward_sense[i][b].first);
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->is_extended) continue;
        ki->is_extended = true;
		extendKernel4Enumeration<dir, if_symmetry>(node, i, b, ki, copy_bucket, Tags, num_routes_now, status);
		if(status) goto outside;
      }
      valid_num = 0;
	  }
	 for(index=0;index<comp.size();++index){
	   int i=comp[index];
	   if(dir? if_exist_extra_labels_in_forward_sense[i][b].second : if_exist_extra_labels_in_backward_sense[i][b].second) {
		 goto STILL_EXIST;
	   }
	 }
	}
  }
  outside:
  for (int i = 0; i < dim; ++i) {
    delete[]copy_bucket[i];
  }
  delete[] copy_bucket;
#if VERBOSE_MODE==1
  std::cout << "Half" << (dir ? "Forward" : "Backward") << " labeling: num_labels= "
            << (dir ? num_forward_labels_in_enu : num_backward_labels_in_enu)
            << " num_routes= " << num_routes_now <<
            std::endl;
#endif
  if (status)return status;
  return 0;
}

template<bool dir, bool if_symmetry>
int CVRP::inspectLabeling(BbNode *node) {
  using namespace std;
  auto num_r1c = node->r1cs.size();
  vector<double> pi4r1c(num_r1c);
  for (int i = 0; i < num_r1c; ++i)pi4r1c[i] = pi4_labeling[node->r1cs[i].idx_r1c];

  for(auto &label_bin: label_concatenate_bin){
	auto ki= label_bin.label;
	auto &aux_label = ki->aux_label;
	int num = aux_label->sparse_num;
	auto &sparse_rep = aux_label->sparse_lp_states;
	auto &states = aux_label->states;
	double tmp_rc = ki->rc + chg_cost_mat4_vertex[ki->end_vertex][label_bin.bin.front()->end_vertex];
	for(auto kj: label_bin.bin){
	  auto &state_back = kj->aux_label->states;
	  double path_rc = kj->rc + tmp_rc;
	  if(path_rc > rc_std) continue;
	  for (int k = 0; k < num; ++k) {
		int idx = sparse_rep[k];
		if (states[idx] + state_back[idx] >= lp_r1c_denominator[idx]) {
		  path_rc -= pi4r1c[idx];
		  if (path_rc > rc_std) goto CONTINUE;
		}
	  }
	  addPathByRC(path_rc, ki, kj, Config::MaxNumRoutesInExact);
	  CONTINUE:;
	}
  }

  int index = num_col - 1;
  writeColumnsInPricingPool(node, index);
  return index - num_col + 1;
}

template<bool if_symmetry>
int CVRP::retrieveLabeling(BbNode *node,
						   const sparseColMatrixXd &rc_mat,
						   const RowVectorXd &cost) {
  using namespace std;
  if_exact_labeling_cg = false;
  priceLabeling(node, pi4_labeling);

  Eigen::Map<const RowVectorXd> pi(pi4_labeling.data(), num_row);
  RowVectorXd reduced_cost = cost - pi * rc_mat;

  for (auto &label : useful_label_vec) label->rc = reduced_cost(label->aux_label->index4_rc_matrix);

  rc_std = RC_TOLERANCE;
  negative_rc_label_tuple.clear();
  map4_each_negative_rc_route.clear();

  int ccnt = inspectLabeling<true, if_symmetry>(node);
  addColumns(node, ccnt);
  return ccnt;
}

template<bool dir, bool if_symmetry>
void CVRP::takeOutUsefulLabels(){
  using namespace std;
  unordered_map<int, unordered_set<yzzLong>> map_end_vertex_ng;
  for(auto &label: negative_rc_label_tuple){
	auto ki= get<0>(label);
	map_end_vertex_ng[ki->end_vertex].emplace(ki->pi);
  }
  label_concatenate_bin.clear();
  label_concatenate_bin.resize(Config::MaxNumLabelReCalculateRC);

  unordered_set<Label*> useful_label_set;
  useful_label_set.reserve(Config::MaxNumLabelReCalculateRC);

  int cnt=0;
  for (auto &label_list : concatenate_labels_in_forward_cg) {
	int i = label_list.first.first;
	if(map_end_vertex_ng.find(i)==map_end_vertex_ng.end()) continue;
	int j=label_list.first.second;
	for (auto &pr : label_list.second) {
	  auto ki = pr.first;
	  bool if_suc = false;
	  for(auto &pi: map_end_vertex_ng[i]){
		if((((ki->pi & pi)^ ki->pi).none())) {
		  if_suc=true;
		  break;
		}
	  }
	  if(!if_suc) continue;
	  auto &label_bin=label_concatenate_bin[cnt];
	  label_bin.label=ki;
	  useful_label_set.emplace(ki);
	  auto &tmp_Resource = pr.second;
	  int arr_bj = int((resource.first_res - tmp_Resource.first_res) / step_size);
	  auto &label_arr =
		  ((!dir && !if_symmetry) || (dir && if_symmetry)) ? label_array_in_forward_sense[j][arr_bj]
														   : label_array_in_backward_sense[j][arr_bj];
	  double tmp_rc= ki->rc+ chg_cost_mat4_vertex[i][j];
	   for (auto &kj : label_arr) {
		 if(tmp_rc+ kj->rc>RC_TOLERANCE) break;
		if constexpr (if_symmetry) {
		  if (tellResTupleRelations<'c'>(tmp_Resource, kj->res)) continue;
		} else {
		  if (dir ? tellResTupleRelations<'p'>(tmp_Resource, kj->res) : tellResTupleRelations<'p'>(kj->res,
																								   tmp_Resource))
			continue;
		}

		if ((ki->pi & kj->pi).any()) continue;
		label_bin.bin.emplace_back(kj);
		useful_label_set.emplace(kj);
	  }
	   if(!label_bin.bin.empty()){
		 ++cnt;
		 if(cnt==label_concatenate_bin.size()){
		   label_concatenate_bin.resize(label_concatenate_bin.size()*2);
		 }
	   }
	}
  }

  label_concatenate_bin.resize(cnt);

 if(useful_label_set.size()>Config::MaxNumLabelReCalculateRC) {
   vector<Label*> vec(useful_label_set.begin(), useful_label_set.end());
   sort(vec.begin(), vec.end(), [](const Label* a, const Label* b){
	 		 return a->rc < b->rc;
	   });
	useful_label_set.clear();
	for(int i=0;i<Config::MaxNumLabelReCalculateRC;++i) useful_label_set.emplace(vec[i]);
   	int count = 0;
   	unordered_set<Label*> tmp_set;
	tmp_set.reserve(vec.size());
	for(auto &label_bin: label_concatenate_bin){
	  auto &bin= label_bin.bin;
	  tmp_set.emplace(label_bin.label);
	  int i=1;
	  for(;i<bin.size();++i){
		if(useful_label_set.find(bin[i])==useful_label_set.end()) {
		  break;
		}
	  }
	  if(i!=bin.size())bin.resize(i);
	  for(int j=0;j<i;++j){
		tmp_set.emplace(bin[j]);
	  }
	}
	useful_label_vec.assign(tmp_set.begin(), tmp_set.end());
 }else{
   useful_label_vec.assign(useful_label_set.begin(), useful_label_set.end());
 }
}




#endif //CVRP_TEMPLATEFUNCTORS_HPP