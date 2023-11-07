

#include <utility>
#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::testCG(BbNode *node, bool if_exact_CG, bool if_record_product_val,
				  bool if_record_improvement, bool if_force_complete) {
  if (branch_pair.size() == 1) {
	cout << "testCG: branch_pair.size() == 1, return!" << endl;
#ifdef USE_M_DYNAMICS
	node->objective_change[branch_pair[0]] = {0, 0, opt_k};
#endif
	return;
  }
  auto beg = high_resolution_clock::now();
  int BeforeNumRow = num_row;
  auto org_val = node->value;
  int cnt = 0;
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(branch_pair.size());
  Brc bf;
  bf.idx_br_c = num_row;
  node->brcs.emplace_back(bf);
  int col_start = num_col;
  auto pseudo_down = &heuristic_improvement_down;
  auto pseudo_up = &heuristic_improvement_up;
  if (if_exact_CG) {
	if (!if_in_enu_state) {
	  force_not_rollback = true;
	  if_force_not_regenerate_bucket_graph = true;
	}
	pseudo_down = &real_improvement_down;
	pseudo_up = &real_improvement_up;
  }
  if (if_in_enu_state) {
	sparseRowMatrixXd mat(1, node->size_enumeration_col_pool);
	node->matrix_in_enumeration.push_back(std::move(mat));
  }
  vector<int> solver_ind(num_col), solver_ind1(num_col), solver_ind2(num_col, BeforeNumRow);
  iota(solver_ind1.begin(), solver_ind1.end(), 0);
  vector<double> solver_val(num_col), solver_val3(num_col);
  vector<int> const_for_branching(num_col);
  iota(const_for_branching.begin(), const_for_branching.end(), num_col);
  bool if_changed = false;

  ++num_row;
  safe_Hyperparameter(checkCSTLimit())

#ifdef USE_M_DYNAMICS
  unordered_map<pair<int, int>, pair<double, double>, PairHasher> objective_change;
#endif
  double max_product = 0;
  for (auto &edge : branch_pair) {
	int ai = edge.first;
	int aj = edge.second;
	if (!if_force_complete) if (max_product > branch_lp[edge]) continue;
	cout << MID_PHASE_SEPARATION;
	cout << "Evaluate on ( " << ai << " , " << aj << " )...\n";
	int numnz;
	double temp_val;
	getNewConstraintCoefficientByEdge(node, edge, solver_ind.data(), solver_val.data(), numnz);
	if (!if_changed) {
	  safe_solver(addBranchConstraint(numnz,
									  solver_ind.data(),
									  solver_val.data(),
									  SOLVER_LESS_EQUAL,
									  0,
									  nullptr,
									  node->solver))
	  if_changed = true;
	} else {
	  changeBranchConstraint(solver_val3.data(),
							 solver_ind2.data(),
							 solver_ind1.data(),
							 numnz,
							 solver_ind.data(),
							 solver_val.data(),
							 SOLVER_LESS_EQUAL,
							 0,
							 node->solver);
	}
	safe_solver(node->solver.updateModel())
	node->brcs.back().edge = {ai, aj};
	node->brcs.back().br_dir = false;
	if (if_in_enu_state) {
	  addBranchConstraint2ColPoolInEnumByColMap(node, edge);
	  solveLPByInspection(node, true, !if_exact_CG, false);
	} else {
	  pool_beg4_pricing = 0;
	  solveLPInLabeling(node, true, if_exact_CG, false);
	}
	node->is_integer = false;
	node->is_terminated = false;
	safe_solver(node->solver.getObjVal(&temp_val))
	auto dif1 = calculateDifference(temp_val, org_val);
	cout << SMALL_PHASE_SEPARATION;
	int len = num_col - col_start;
	if (const_for_branching.size() < len) {
	  int arr_beg = (int)const_for_branching.size();
	  int arr_val = const_for_branching.back() + 1;
	  const_for_branching.resize(len);
	  iota(const_for_branching.begin() + arr_beg, const_for_branching.end(), arr_val);
	}
	safe_solver(node->solver.delVars(num_col - col_start, const_for_branching.data()))
	safe_solver(inverseLastBranchConstraint(SOLVER_GREATER_EQUAL, 1, node->solver))
	node->brcs.back().br_dir = true;
	safe_solver(node->solver.updateModel())
	safe_solver(node->solver.getNumCol(&num_col))
	if (if_in_enu_state) {
	  solveLPByInspection(node, true, !if_exact_CG, false);
	} else {
	  pool_beg4_pricing = 0;
	  solveLPInLabeling(node, true, if_exact_CG, false);
	}
	node->is_integer = false;
	node->is_terminated = false;
	safe_solver(node->solver.getObjVal(&temp_val))
	auto dif2 = calculateDifference(temp_val, org_val);
	auto product = dif1 * dif2;
	cout << "ldf= " << setw(6) << left << dif1 << "  rdf= " << setw(6) << left << dif2 << "  pd= " << setw(6)
		 << left << product << endl;
	len = num_col - col_start;
	if (const_for_branching.size() < len) {
	  int arr_beg = (int)const_for_branching.size();
	  int arr_val = const_for_branching.back() + 1;
	  const_for_branching.resize(len);
	  iota(const_for_branching.begin() + arr_beg, const_for_branching.end(), arr_val);
	}
	safe_solver(node->solver.delVars(num_col - col_start, const_for_branching.data()))
	safe_solver(node->solver.updateModel())
	safe_solver(node->solver.getNumCol(&num_col))
#ifdef USE_M_DYNAMICS
	objective_change[edge] = {dif1, dif2};
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
  safe_solver(node->solver.delConstraints(1, &BeforeNumRow))
  safe_solver(node->solver.updateModel())
  --num_row;

  if (if_exact_CG && !if_in_enu_state) {
	force_not_rollback = false;
	if_force_not_regenerate_bucket_graph = false;
  }
  if (if_in_enu_state) {
	node->matrix_in_enumeration.pop_back();
  }
  safe_solver(node->solver.reoptimize())

  if (if_record_product_val) {
	sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
	  return a.second > b.second;
	});
	branch_pair_val = Branch_Val;
	branch_pair = {Branch_Val[0].first};
  } else {
	auto it = max_element(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
	  return a.second < b.second;
	});
	branch_pair = {it->first};
  }
  node->brcs.pop_back();
#ifdef USE_M_DYNAMICS
  node->objective_change[branch_pair[0]] =
	  {objective_change[branch_pair[0]].first, objective_change[branch_pair[0]].second, opt_k};
#endif
  /**
   * is_terminated, value. these two cannot be deleted! since when exact cg could determine if this node is no promising to be solved
   */
  node->is_terminated = false;
  node->value = org_val;
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "testCG spent= " << eps << "s" << endl;
}

#ifdef MASTER_VALVE_ML

void CVRP::useModelInPhase2(BbNode *node, int num) {
  if (branch_pair.size() == 1) {
	cout << "useModelInPhase2: only one edge, no need to use model" << endl;
	return;
  }
  auto beg = high_resolution_clock::now();
  getTrainingDataInPhase2(node);
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(branch_pair.size());
  transform(branch_pair.begin(),
			branch_pair.end(),
			Branch_Val.begin(),
			[&](const auto &edge) {
			  return std::make_pair(edge, 0.0);
			});
  ml.predict(Branch_Val, 2);

  branch_pair.resize(num);
  transform(Branch_Val.begin(),
			Branch_Val.begin() + num,
			branch_pair.begin(),
			[&](const auto &val) {
			  return val.first;
			});
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();

#ifdef USE_M_DYNAMICS
  double average_t = eps / 2 / (double)Branch_Val.size();
  BbNode::updateState(average_t, node->t_for_one_lp, (int)node->brcs.size());
  cout << "eps= " << eps << " average_t= " << average_t << " t_for_one_lp= " << node->t_for_one_lp << endl;
#endif
  cout << "useModelInPhase2 spent= " << eps << "s" << endl;
}

void CVRP::getTrainingDataInPhase2(BbNode *node) {
  if (branch_pair.size() == 1) {
	cout << "getTrainingDataInPhase2: only one edge, no need to get this data!" << endl;
	return;
  }
  int BeforeNumRow = num_row;
  auto org_val = node->value;

  Brc bf;
  bf.idx_br_c = num_row;
  node->brcs.emplace_back(bf);

  vector<int> solver_ind(num_col), solver_ind2(num_col, BeforeNumRow);
  vector<double> solver_val(num_col), solver_val3(num_col);
  vector<int> solver_ind1(num_col);
  iota(solver_ind1.begin(), solver_ind1.end(), 0);
  bool if_changed = false;

  branch_lp.clear();
  ++num_row;
  safe_Hyperparameter(checkCSTLimit())
  for (auto &edge : branch_pair) {
	int numnz;
	double tmp_val;
	getNewConstraintCoefficientByEdge(node, edge, solver_ind.data(), solver_val.data(), numnz);
	if (!if_changed) {
	  safe_solver(addBranchConstraint(numnz,
									  solver_ind.data(),
									  solver_val.data(),
									  SOLVER_LESS_EQUAL,
									  0,
									  nullptr,
									  node->solver))
	  if_changed = true;
	} else {
	  changeBranchConstraint(solver_val3.data(),
							 solver_ind2.data(),
							 solver_ind1.data(),
							 numnz,
							 solver_ind.data(),
							 solver_val.data(),
							 SOLVER_LESS_EQUAL,
							 0,
							 node->solver);
	}
	safe_solver(node->solver.reoptimize())
	safe_solver(node->solver.getObjVal(&tmp_val))
	auto dif1 = calculateDifference(tmp_val, org_val);
	ml.collect_resolving_features(node, edge, BeforeNumRow, tmp_val, org_val, numnz, false);
	safe_solver(inverseLastBranchConstraint(SOLVER_GREATER_EQUAL, 1, node->solver))
	safe_solver(node->solver.reoptimize())
	safe_solver(node->solver.getObjVal(&tmp_val))
	auto dif2 = calculateDifference(tmp_val, org_val);
	ml.collect_resolving_features(node, edge, BeforeNumRow, tmp_val, org_val, numnz, true);
	branch_lp[edge] = dif1 * dif2;
	lp_testing_improvement_down[edge].first += dif1;
	++lp_testing_improvement_down[edge].second;
	lp_testing_improvement_up[edge].first += dif2;
	++lp_testing_improvement_up[edge].second;
  }
  --num_row;
  safe_solver(node->solver.delConstraints(1, &BeforeNumRow))
  safe_solver(node->solver.reoptimize())
  node->brcs.pop_back();
  cout << "getTrainingDataInPhase2!" << endl;
}

#endif

void CVRP::addBranchConstraint2ColPoolInEnumByColMap(BbNode *node,
													 const pair<int, int> &edge) {
  auto &mat = node->matrix_in_enumeration;
  int size = node->size_enumeration_col_pool;
  if (!size) return;
  auto &colmap = node->column_pool_mapping;
  auto &mat_last = mat.back();
  mat_last.setZero();
  for (auto i : colmap[edge]) {
	mat_last.insert(0, i) = 1;
  }
}

