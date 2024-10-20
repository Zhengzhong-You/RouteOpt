//
// Created by You, Zhengzhong on 4/3/24.
//

#include "heuristic.hpp"

using namespace std;
using namespace chrono;
using namespace Eigen;

double Heuristic::fake_ub = 0.0;
bool Heuristic::if_fail = false;
double Heuristic::heuristic_mip_time = 0.0;
std::vector<SequenceInfo> Heuristic::heuristic_col_info;

void Heuristic::cleanHeuristicMIPCols() {
  heuristic_col_info.clear();
}

void Heuristic::addHeuristicMIPCols(CVRP *cvrp, BbNode *node) {
  auto num_col = cvrp->getNumCol();
  auto old_size = heuristic_col_info.size();
  heuristic_col_info.resize(heuristic_col_info.size() + node->getCols().size() - num_col);
  transform(node->getCols().begin() + num_col, node->getCols().end(), heuristic_col_info.begin() + (int)old_size,
			[](const SequenceInfo &info) {
			  return info;
			});
}

template<bool if_symmetry>
void Heuristic::completeHeuristic(CVRP *cvrp, BbNode *node) {
  cvrp->opt_gap = fake_ub;
  cvrp->priceLabeling(node, cvrp->optimal_dual_vector);

  if (if_symmetry) {
	cvrp->concatenatePhaseInArcElimination<true, true>(node);
	if (cvrp->if_roll_back || cvrp->if_short_memory) goto QUIT;
	cvrp->runLabeling<true, true, false, true, 0>(node);
	if (cvrp->if_roll_back || cvrp->if_short_memory) goto QUIT;
  } else {
	cvrp->concatenatePhaseInArcElimination<true, false>(node);
	if (cvrp->if_roll_back || cvrp->if_short_memory) goto QUIT;
	cvrp->runLabeling<true, true, false, false, 0>(node); // heuristic
	if (cvrp->if_roll_back || cvrp->if_short_memory) goto QUIT;
	cvrp->concatenatePhaseInArcElimination<false, false>(node);
	if (cvrp->if_roll_back || cvrp->if_short_memory) goto QUIT;
	cvrp->runLabeling<false, true, false, false, 0>(node); // heuristic
	if (cvrp->if_roll_back || cvrp->if_short_memory) goto QUIT;
  }

QUIT:
  if (cvrp->if_roll_back) {
	if_fail = true;
	cout << "the labeling process is too long!" << endl;
  }
}

void Heuristic::enumerateHeuristicMIP(CVRP *cvrp, BbNode *node) {
#ifdef SYMMETRY_PROHIBIT
  completeHeuristic<false>(cvrp, node);
#else
  completeHeuristic<true>(cvrp, node);
#endif

  if (abs(cvrp->meet_point_resource_in_bi_dir_enu) < TOLERANCE)
	cvrp->meet_point_resource_in_bi_dir_enu = cvrp->meet_point_resource_in_bi_dir;

  int Max_routes_all = Config::MaxNumRouteInEnumeration;
  int num_routes_now = 0;

  int index;
  cvrp->opt_gap = fake_ub;
  unordered_map<yzzLong, tuple<Label *, Label *, double>> Tags;
  Tags.reserve(Max_routes_all);
  int dim = cvrp->dim;
  int num_buckets_per_vertex = cvrp->num_buckets_per_vertex;
  auto &label_array_in_forward_sense = cvrp->label_array_in_forward_sense;
  auto &label_array_in_backward_sense = cvrp->label_array_in_backward_sense;

  auto copy_Forward_bucket = new vector < Label * > *[dim];
  symmetry_prohibit_call(auto
  copy_Backward_bucket =
	  new vector < Label * > *[dim];//these memory will be freed inside the forward function
  )
  for (int i = 0; i < dim; ++i) {
	copy_Forward_bucket[i] = new vector<Label *>[num_buckets_per_vertex];
	symmetry_prohibit_call(
		copy_Backward_bucket[i] = new vector<Label *>[num_buckets_per_vertex];)
	for (int b = 0; b < num_buckets_per_vertex; ++b) {
	  copy_Forward_bucket[i][b].assign(label_array_in_forward_sense[i][b].begin(),
									   label_array_in_forward_sense[i][b].end());
	  symmetry_prohibit_call(
		  copy_Backward_bucket[i][b].assign(label_array_in_backward_sense[i][b].begin(),
											label_array_in_backward_sense[i][b].end());)
	}
  }

  auto beg = high_resolution_clock::now();

  int status =
#ifdef SYMMETRY_PROHIBIT
	  cvrp->enumerateHalfwardRoutes<true, false>(node, Tags,
												 copy_Backward_bucket, num_routes_now);
#else
	  cvrp->enumerateHalfwardRoutes<true, true>(node, Tags, copy_Forward_bucket, num_routes_now);
#endif

  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  verbose_call(cout << "Half Forward time= " << eps << "s" << endl;)

  if (status || cvrp->if_roll_back || cvrp->if_short_memory) {
	if (status == 1) {
	  if_fail = true;
	  cout << "the number of labels in Forward reached its limit!" << endl;
	}
  }

  symmetry_prohibit_call(beg = high_resolution_clock::now();

  status = cvrp->enumerateHalfwardRoutes<false, false>(node,
													   Tags,
													   copy_Forward_bucket,
													   num_routes_now);

  end = high_resolution_clock::now();
  eps = duration<double>(end - beg).count();
  verbose_call(cout << "Half Backward time= " << eps << "s" << endl;
  )

  if (status || cvrp->if_roll_back || cvrp->if_short_memory) {
	if (status == 1) {
	  if_fail = true;
	  cout << "the number of labels in Backward reached its limit!" << endl;
	}
  })

  beg = high_resolution_clock::now();

  status = cvrp->concatenateRoutesPriorForwardInEnumeration(node, Tags, num_routes_now);

  end = high_resolution_clock::now();
  eps = duration<double>(end - beg).count();

  verbose_call(cout << "Concatenate time= " << eps << "s" << endl;
  )

  if (status) {
	if_fail = true;
	cout << "the number of routes reached its limit!" << endl;
  }

  index = (int)heuristic_col_info.size();
  heuristic_col_info.resize(heuristic_col_info.size() + Tags.size());
  Label *ki, *li, *p;
  vector<int> seq(dim + 1);
  vector<res_int> res_vec(dim + 1);
  for (auto &tag : Tags) {
	ki = get<0>(tag.second);
	li = get<1>(tag.second);
	int cnt = 0;
	p = ki;
	while (p->p_label) {
	  res_vec[cnt] = p->res.first_res;
	  seq[cnt++] = p->end_vertex;
	  p = p->p_label;
	}
	auto &col = heuristic_col_info[index].col_seq;
	auto &res = heuristic_col_info[index].main_res;
	col.resize(cnt);
	res.resize(cnt);
	for (int k = 0; k < cnt; ++k) {
	  col[k] = seq[cnt - k - 1];
	  res[k] = res_vec[cnt - k - 1];
	}

	heuristic_col_info[index].forward_concatenate_pos = cnt - 1;

	if (li) {
	  p = li;
	  while (p->p_label) {
		col.emplace_back(p->end_vertex);
		res.emplace_back(p->res.first_res);
		p = p->p_label;
	  }
	}
	++index;
  }
}

void Heuristic::heuristicApplyRCF(Solver &local_solver,
								  int round,
								  bool if_verbose,
								  double ub,
								  int real_dim,
								  vector<int> &col_map) {
  cout << "WARNING: RCF can only be applied by GRB model so far!" << endl;
#if SOLVER_TYPE == 0
  col_map.clear();
  vector<int> idx;
  Solver tmp_solver{};
  tmp_solver.model = local_solver.copyModel();
  deLuxing(tmp_solver.model,
		   ub,
		   // real_dim,
		   round,
		   HEURISTIC_DELUXING_BETA1,
		   HEURISTIC_DELUXING_BETA2,
		   idx,
		   HEURISTIC_DELUXING_TIME_LIMIT,
		   MIP_GAP_TOLERANCE,
		   HEURISTIC_DELUXING_VERBOSE);
  tmp_solver.freeModel();
  int num_col;
  safe_solver(local_solver.getNumCol(&num_col))
  if (idx.back() > num_col) {
	cout << idx.back() << " " << num_col << endl;
	throw runtime_error("idx.back()>num_col");
  }
  sort(idx.begin(), idx.end());
  if (idx.front() == 0) idx.erase(idx.begin());
  safe_solver(local_solver.delVars(idx.size(), idx.data()))
  safe_solver(local_solver.updateModel())
  col_map.resize(num_col);
  int delta = 0;
  auto stop_sign = idx.end() - 1;
  for (auto i = idx.begin(); i < stop_sign; ++i) {
	++delta;
	for (int j = *i + 1; j < *(i + 1); ++j) col_map[j - delta] = j;
  }
  ++delta;
  safe_solver(local_solver.getNumCol(&num_col))
  for (int j = *stop_sign + 1; j < num_col; ++j) col_map[j - delta] = j;
#endif
}

void Heuristic::rmNonCapCols(const double *demand, double cap) {
  if (heuristic_col_info.empty()) return;
  for (auto iter = heuristic_col_info.begin() + 1; iter != heuristic_col_info.end();) {
	double local_cap = 0;
	bool if_erase = false;
	for (auto &i : iter->col_seq) {
	  local_cap += demand[i];
	  if (local_cap > cap + TOLERANCE) {
		if_erase = true;
		iter = heuristic_col_info.erase(iter);
		break;
	  }
	}
	if (!if_erase)++iter;
  }
}

void Heuristic::rmNonEleCols() {
  if (heuristic_col_info.empty()) return;
  for (auto iter = heuristic_col_info.begin() + 1; iter != heuristic_col_info.end();) {
	yzzLong tmp = 0;
	bool if_erase = false;
	for (auto &i : iter->col_seq) {
	  if (tmp.test(i)) {
		if_erase = true;
		iter = heuristic_col_info.erase(iter);
		break;
	  }
	  tmp.set(i);
	}
	if (!if_erase)++iter;
  }
}

void Heuristic::heuristicMIP(CVRP *cvrp, BbNode *&node) {
#ifndef HEURISTIC
  cout << "heuristic is shut down!" << endl;
  return;
#endif
  if (node->index) return;

  auto beg = high_resolution_clock::now();
  Solver tmp_node_solver{};
  tmp_node_solver.model = node->getSolver().copyModel();
  std::swap(node->getSolver().model, tmp_node_solver.model);
  auto cols = node->getCols();

  cvrp->force_not_rollback = true;
  cvrp->solveLPInLabeling(node);
  cvrp->force_not_rollback = false;
  fake_ub = min(node->getCurrentNodeVal() * GUESSED_OPT_GAP_PER, BaseBranching::ub);
  cout << "fake_ub is " << fake_ub << endl;
  rmNonEleCols();
  heuristic_enumeration_call(enumerateHeuristicMIP(cvrp, node);)
  solver_vrptw_call(rmNonCapCols(cvrp->demand, cvrp->cap);)
  if (heuristic_col_info.empty()) return;
  if (!node->getBrCs().empty()) {
	throw runtime_error("brcs is not empty!");
  }
  int ccnt = (int)heuristic_col_info.size();
  int dim = cvrp->dim;
  int num_row = cvrp->getNumRow();
  auto &cost_mat4_vertex = cvrp->cost_mat4_vertex;
  int ccnt_cnt;
  int past_node;
  double cost_sum;
  unordered_map<pair<int, int>, vector<int>, PairHasher> edge_map;
  edge_map.reserve(dim * dim);
  sparseColMatrixXd mat(num_row, ccnt);
  mat.setZero();
  Eigen::RowVectorXd cost(ccnt);

  unordered_map<pair<int, int>, double, PairHasher> coeff_map;
  coeff_map.reserve(dim * dim);

  ccnt_cnt = 0;
  for (int i = 0; i < ccnt; ++i) {
	auto &col = heuristic_col_info[i].col_seq;
	past_node = 0;
	cost_sum = 0;
	for (int curr_node : col) {
	  cost_sum += cost_mat4_vertex[past_node][curr_node];
	  ++coeff_map[{curr_node - 1, ccnt_cnt}];
	  auto pr = past_node < curr_node ? make_pair(past_node, curr_node) : make_pair(curr_node, past_node);
	  edge_map[pr].emplace_back(ccnt_cnt);
	  past_node = curr_node;
	}
	edge_map[make_pair(0, past_node)].emplace_back(ccnt_cnt);
	cost_sum += cost_mat4_vertex[past_node][0];
	cost(ccnt_cnt) = cost_sum;
	++ccnt_cnt;
  }

  safe_Hyperparameter(ccnt != ccnt_cnt)

  for (auto &rcc : node->rccs) {
	if (rcc.form_rcc) {
	  auto &info = rcc.info_rcc_customer;
	  int idx = rcc.idx_rcc;
	  for (auto it = info.begin(); it != info.end(); ++it) {
		int ai = *it;
		auto it_inner = it;
		++it_inner;
		for (; it_inner != info.end(); ++it_inner) {
		  int aj = *it_inner;
		  auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
		  for (auto col : edge_map[pr]) ++coeff_map[{idx, col}];
		}
	  }
	} else {
	  auto &customer_info = rcc.info_rcc_customer;
	  auto &outside_customer_info = rcc.info_rcc_outside_customer;
	  int idx = rcc.idx_rcc;
	  for (auto it = outside_customer_info.begin(); it != outside_customer_info.end(); ++it) {
		int ai = *it;
		auto it_inner = it;
		++it_inner;
		for (; it_inner != outside_customer_info.end(); ++it_inner) {
		  int aj = *it_inner;
		  auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
		  for (auto it_map : edge_map[pr]) ++coeff_map[{idx, it_map}];
		}
	  }
	  for (int aj : outside_customer_info) {
		for (auto it_map : edge_map[make_pair(0, aj)]) coeff_map[{idx, it_map}] += 0.5;
	  }
	  for (int aj : customer_info) {
		for (auto it_map : edge_map[make_pair(0, aj)]) coeff_map[{idx, it_map}] -= 0.5;
	  }
	}
  }

  for (auto &br : node->getBrCs()) {
	int ai = br.edge.first;
	int aj = br.edge.second;
	int idx = br.idx_br_c;
	auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
	for (auto it_map : edge_map[pr]) ++coeff_map[{idx, it_map}];
  }

  sparseColMatrixXd mat_r1c;
  cvrp->getLimitedR1CCoeffs(heuristic_col_info, mat_r1c);

  size_t n = coeff_map.size() + mat_r1c.nonZeros() + ccnt;
  vector<Eigen::Triplet<double>> triplet(n);
  n = 0;

  for (auto &it : coeff_map) {
	triplet[n++] = {it.first.first, it.first.second, it.second};
  }

  for (int i = 0; i < mat_r1c.outerSize(); ++i) {
	for (sparseColMatrixXd::InnerIterator it(mat_r1c, i); it; ++it) {
	  triplet[n++] = {cvrp->lp_r1c_map[it.row()], (int)it.col(), it.value()};
	}
  }
  for (int i = 0; i < ccnt; ++i) {
	triplet[n++] = {cvrp->real_dim, i, 1};
  }

  mat.setFromTriplets(triplet.begin(), triplet.end());

  {
	int num_col = cvrp->getNumCol();
	vector<size_t> solver_beg;
	vector<int> solver_ind;
	vector<double> solver_val;
	vector<double> solver_obj;
	Solver local_solver{};
	local_solver.model = node->getSolver().copyModel();
	vector<int> idx(num_col);
	iota(idx.begin(), idx.end(), 0);
	safe_solver(local_solver.delVars(num_col, idx.data()))
	safe_solver(local_solver.updateModel())
	for (int i = 0; i < ccnt; i++) {
	  solver_obj.emplace_back(cost[i]);
	  solver_beg.emplace_back(solver_ind.size());
	  for (sparseColMatrixXd::InnerIterator it(mat, i); it; ++it) {
		solver_ind.emplace_back(it.row());
		solver_val.emplace_back(it.value());
	  }
	}
	solver_beg.emplace_back(solver_ind.size());
	safe_solver(local_solver.XaddVars(solver_beg.size() - 1,
									  solver_ind.size(),
									  solver_beg.data(),
									  solver_ind.data(),
									  solver_val.data(),
									  solver_obj.data(),
									  nullptr,
									  nullptr,
									  nullptr,
									  nullptr))
	safe_solver(local_solver.updateModel())
	vector<int> col_map;
#ifdef HEURISTIC_DELUXING_ON
	heuristicApplyRCF(local_solver, HEURISTIC_DELUXING_ROUND, true, node->getCurrentNodeVal() + fake_ub, cvrp->real_dim, col_map);
#else
	col_map.resize(ccnt);
	iota(col_map.begin(), col_map.end(), 0);
#endif

	Map<RowVectorXd> local_pi(cvrp->optimal_dual_vector.data(), num_row);
	RowVectorXd rc = cost - local_pi * mat;

	int local_num_col;
	safe_solver(local_solver.getNumCol(&local_num_col))
	vector<double> rc_cp(local_num_col);
	for (int i = 0; i < local_num_col; ++i) {
	  rc_cp[i] = rc[col_map[i]];
	}

	int num_columns = min(NUM_COL_HEURISTIC_MIP, local_num_col);
	if (num_columns != local_num_col) {
	  auto rc_cp2 = rc_cp;
	  nth_element(rc_cp2.begin(), rc_cp2.begin() + num_columns, rc_cp2.end(), less<>());
	  double threshold = rc_cp2[num_columns];
	  vector<int> delete_cols;
	  for (int i = 1; i < local_num_col; ++i) {
		if (rc_cp[i] > threshold) {
		  delete_cols.emplace_back(i);
		}
	  }

	  safe_solver(local_solver.delVars(delete_cols.size(), delete_cols.data()))
	  safe_solver(local_solver.updateModel())

	  int delta = 0;
	  auto stop_sign = delete_cols.end() - 1;
	  for (auto i = delete_cols.begin(); i < stop_sign; ++i) {
		++delta;
		for (int j = *i + 1; j < *(i + 1); ++j) col_map[j - delta] = col_map[j];
	  }
	  ++delta;
	  for (int j = *stop_sign + 1; j < local_num_col; ++j) col_map[j - delta] = col_map[j];
	  local_num_col -= delta;
	  col_map.resize(local_num_col);
	}
	vector<char> vtype(local_num_col, GRB_BINARY);
	safe_solver(local_solver.setVTypeArray(0, local_num_col, vtype.data()))
	safe_solver(local_solver.setEnvCutoff(BaseBranching::ub))
	safe_solver(local_solver.setEnvThreads(HEURISTIC_MIP_THREADS_NUM, false))
	safe_solver(local_solver.optimize())
	safe_solver(local_solver.setEnvThreads(1, false))
	int status;
	safe_solver(local_solver.getStatus(&status))
	if (false && status == SOLVER_OPTIMAL) {
	  double value;
	  safe_solver(local_solver.getObjVal(&value))
	  BaseBranching::ub = min(BaseBranching::ub, value);
	  cout << "ub is updated to " << BaseBranching::ub << endl;
	  vector<double> x(local_num_col);
	  safe_solver(local_solver.getX(0, local_num_col, x.data()))
	  cvrp->ip_opt_sol.clear();
	  for (int i = 0; i < local_num_col; ++i) {
		if (x[i] > 0.5) {
		  cvrp->ip_opt_sol.emplace_back(heuristic_col_info[col_map[i]].col_seq);
		}
	  }
	} else {
	  cout << "fail to find any better solutions!" << endl;
	  goto FREE;
	  cout << "we terminate this solution process!" << endl;
	  node->getIfTerminated() = true;
	  if_fail = true;
	}
FREE:
	safe_solver(local_solver.setEnvCutoff(GRB_INFINITY))
	local_solver.freeModel();
  }

  heuristic_mip_time = duration<double>(high_resolution_clock::now() - beg).count();
  cleanHeuristicMIPCols();
  node->getSolver().freeModel();
  std::swap(node->getSolver().model, tmp_node_solver.model);
  tmp_node_solver.model = nullptr;
  node->getCols() = cols;
  safe_solver(node->getSolver().getNumCol(&cvrp->getNumCol()))
  if (node->getIfTerminated()) {
	delete node;
	node = nullptr;
  }
}