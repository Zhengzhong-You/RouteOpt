

#include <utility>
#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;
using namespace Eigen;

void CVRP::testCG(BbNode *node, bool if_exact_CG, bool if_record_product_val,
				  bool if_record_improvement, bool if_force_complete) {
  if (branch_pair.size() == 1) {
#if VERBOSE_MODE == 1
	cout << "testCG: branch_pair.size() == 1, return!" << endl;
#endif
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
  vector<int> solver_ind;
  vector<double> solver_val;
  vector<int> const_for_branching(num_col);
  iota(const_for_branching.begin(), const_for_branching.end(), num_col);
  bool if_changed = false;

  ++num_row;

#ifdef USE_M_DYNAMICS
  unordered_map<pair<int, int>, pair<double, double>, PairHasher> objective_change;
#endif
  double max_product = 0;
  for (auto &edge : branch_pair) {
	if (!if_force_complete) if (max_product > branch_lp[edge]) continue;
	cout << MID_PHASE_SEPARATION;
	cout << "Evaluate on ( " << edge.first << " , " << edge.second << " )...\n";
	double temp_val;
	getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
	if (!if_changed) {
	  safe_solver(addBranchConstraint(solver_ind,
									  solver_val,
									  SOLVER_LESS_EQUAL,
									  0,
									  nullptr,
									  node->solver))
	  if_changed = true;
	} else {
	  changeBranchConstraint(solver_ind,
							 solver_val,
							 SOLVER_LESS_EQUAL,
							 0,
							 BeforeNumRow,
							 node->solver);
	}
	safe_solver(node->solver.updateModel())
	node->brcs.back().edge = edge;
	node->brcs.back().br_dir = false;
#ifdef TEST_CG_TIME
	time_point<high_resolution_clock> start, end;
#endif
	if (if_in_enu_state) {
#ifdef TEST_CG_TIME
	  start = high_resolution_clock::now();
#endif
	  addBranchConstraint2ColPoolInEnumByColMap(node, edge);
#ifdef TEST_CG_TIME
	  end = high_resolution_clock::now();
	  cout << "addBranchConstraint2ColPoolInEnumByColMap spent= " << duration<double>(end - start).count() << "s"
		   << endl;
	  beg = high_resolution_clock::now();
#endif
	  solveLPByInspection(node, true, !if_exact_CG, false);
#ifdef TEST_CG_TIME
	  end = high_resolution_clock::now();
	  cout << "solveLPByInspection spent= " << duration<double>(end - start).count() << "s" << endl;
#endif
	} else {
#ifdef TEST_CG_TIME
	  start = high_resolution_clock::now();
#endif
	  solveLPInLabeling(node, true, if_exact_CG, false);
#ifdef TEST_CG_TIME
	  end = high_resolution_clock::now();
	  cout << "solveLPInLabeling spent= " << duration<double>(end - start).count() << "s" << endl;
#endif
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
	rmLPCols(node, vector<int>(const_for_branching.begin(), const_for_branching.begin() + len));
	safe_solver(inverseLastBranchConstraint(SOLVER_GREATER_EQUAL, 1, node->solver))
	node->brcs.back().br_dir = true;
	safe_solver(node->solver.updateModel())
	safe_solver(node->solver.getNumCol(&num_col))
	if (if_in_enu_state) {
#ifdef TEST_CG_TIME
	  start = high_resolution_clock::now();
#endif
	  solveLPByInspection(node, true, !if_exact_CG, false);
#ifdef    TEST_CG_TIME
	  end = high_resolution_clock::now();
	  cout << "solve LPByInspection spent= " << duration<double>(end - start).count() << "s" << endl;
#endif
	} else {
#ifdef TEST_CG_TIME
	  start = high_resolution_clock::now();
#endif
	  solveLPInLabeling(node, true, if_exact_CG, false);
#ifdef TEST_CG_TIME
	  end = high_resolution_clock::now();
	  cout << "solveLPInLabeling spent= " << duration<double>(end - start).count() << "s" << endl;
#endif
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
	rmLPCols(node, vector<int>(const_for_branching.begin(), const_for_branching.begin() + len));
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
#ifdef TEST_CG_TIME
	time_point<high_resolution_clock> beg2, end2;
	beg2 = high_resolution_clock::now();
#endif
  safe_solver(node->solver.reoptimize())
#ifdef TEST_CG_TIME
  end2 = high_resolution_clock::now();
  cout << "reoptimize spent= " << duration<double>(end2 - beg2).count() << "s" << endl;
#endif

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
#if VERBOSE_MODE == 1
  cout << "testCG spent= " << eps << "s" << endl;
#endif
}

#ifdef MASTER_VALVE_ML

void CVRP::useModelInPhase2(BbNode *node, int num) {
  if (branch_pair.size() == 1) {
#if VERBOSE_MODE == 1
	cout << "useModelInPhase2: only one edge, no need to use model" << endl;
#endif
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
  getAverageT4LP(node, eps, (int)Branch_Val.size());
#endif
}

void CVRP::getTrainingDataInPhase2(BbNode *node) {
  if (branch_pair.size() == 1) {
#if VERBOSE_MODE == 1
	cout << "getTrainingDataInPhase2: only one edge, no need to get this data!" << endl;
#endif
	return;
  }
  int BeforeNumRow = num_row;
  auto org_val = node->value;

  Brc bf;
  bf.idx_br_c = num_row;
  node->brcs.emplace_back(bf);

  vector<int> solver_ind;
  vector<double> solver_val;
  bool if_changed = false;

  branch_lp.clear();
  ++num_row;
  int method;
  safe_solver(node->solver.getEnvMethod(&method))
  safe_solver(node->solver.setEnvMethod(SOLVER_BARRIER))
  safe_solver(node->solver.setEnvCrossOver(SOLVER_CROSSOVER_CLOSE))
  for (auto &edge : branch_pair) {
	double tmp_val;
	getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
	if (!if_changed) {
	  safe_solver(addBranchConstraint(solver_ind,
									  solver_val,
									  SOLVER_LESS_EQUAL,
									  0,
									  nullptr,
									  node->solver))
	  if_changed = true;
	} else {
	  changeBranchConstraint(solver_ind,
							 solver_val,
							 SOLVER_LESS_EQUAL,
							 0,
							 BeforeNumRow,
							 node->solver);
	}
	safe_solver(node->solver.reoptimize())
	safe_solver(node->solver.getObjVal(&tmp_val))
	auto dif1 = calculateDifference(tmp_val, org_val);
	ml.collect_resolving_features(node, edge, BeforeNumRow, tmp_val, org_val);
	safe_solver(inverseLastBranchConstraint(SOLVER_GREATER_EQUAL, 1, node->solver))
	safe_solver(node->solver.reoptimize())
	safe_solver(node->solver.getObjVal(&tmp_val))
	auto dif2 = calculateDifference(tmp_val, org_val);
	ml.collect_resolving_features(node, edge, BeforeNumRow, tmp_val, org_val);
	branch_lp[edge] = dif1 * dif2;
	lp_testing_improvement_down[edge].first += dif1;
	++lp_testing_improvement_down[edge].second;
	lp_testing_improvement_up[edge].first += dif2;
	++lp_testing_improvement_up[edge].second;
  }
  safe_solver(node->solver.setEnvMethod(method))
  --num_row;
  safe_solver(node->solver.delConstraints(1, &BeforeNumRow))
  safe_solver(node->solver.reoptimize())
  node->brcs.pop_back();
#if VERBOSE_MODE == 1
  cout << "getTrainingDataInPhase2!" << endl;
#endif
}

#endif

void CVRP::addBranchConstraint2ColPoolInEnumByColMap(BbNode *node,
													 const pair<int, int> &edge) const {
  auto &mat = node->matrix_in_enumeration;
  int size = node->size_enumeration_col_pool;
  if (!size) return;
  auto &mat_last = mat.back();
  mat_last.setZero();
  vector<Eigen::Triplet<double>> triplets(size);
  sparseRowMatrixXd tmp(1, size);
  int cnt = 0;
  int ai = edge.first, aj = edge.second;
  auto &mat0 = mat.front();
  if (ai) {
	tmp = mat0.row(ai - 1) + mat0.row(aj - 1);
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
	  if (it.value() > 1.5) {
		for (auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;; ++j) {
		  int current_node = col_pool4_pricing[j];
		  if (!current_node) break;
		  if (current_node == ai) {
			if (col_pool4_pricing[j + 1] == aj || col_pool4_pricing[j - 1] == aj)
			  triplets[cnt++] = {0, (int)it.col(), 1};
		  }
		}
	  }
	}
  } else {//now ai ==0
	for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat0, aj - 1); it; ++it) {
	  if (it.value() > 0.5) {//further test
		auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;
		if (col_pool4_pricing[j] == aj) {
		  triplets[cnt++] = {0, (int)it.col(), 1};
		  continue;
		}
		for (;; ++j) {
		  int current_node = col_pool4_pricing[j];
		  if (!current_node) break;
		}
		if (col_pool4_pricing[j - 1] == aj) {
		  triplets[cnt++] = {0, (int)it.col(), 1};
		}
	  }
	}
  }
  mat_last.setFromTriplets(triplets.begin(), triplets.end());
}

