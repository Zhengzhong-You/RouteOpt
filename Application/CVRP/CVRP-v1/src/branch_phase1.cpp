
#include <utility>

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;


void getM(double &m, double &alpha, const std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val);

void CVRP::testLP(BbNode *node, int num, bool if_writeBranch_pair, bool if_record_sb_scores,
				  bool if_record_LP_improvement) {
  if (branch_pair.size() == 1) {
	cout << "No LP Testing: a single candidate" << endl;
	return;
  }
  time_point<high_resolution_clock> beg, end;
  double eps;
  beg = high_resolution_clock::now();
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(branch_pair.size());
  double org_val = node->value;
  int BeforeNumRow = num_row;
  int cnt = 0;
  vector<int> solver_ind;
  vector<double> solver_val;
  bool if_changed = false;
  branch_lp.clear();
  int method;
  safe_solver(node->solver.getEnvMethod(&method))
  safe_solver(node->solver.setEnvMethod(SOLVER_BARRIER))
  safe_solver(node->solver.setEnvCrossOver(SOLVER_CROSSOVER_CLOSE))
  for (auto &edge : branch_pair) {
	int numnz;
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
#ifdef TEST_LP_TIME
	time_point<high_resolution_clock> beg1, end1;
	beg1 = high_resolution_clock::now();
#endif
	safe_solver(node->solver.reoptimize())
#ifdef TEST_LP_TIME
	end1 = high_resolution_clock::now();
	cout << "time for reoptimize: " << duration<double>(end1 - beg1).count() << endl;
#endif
	safe_solver(node->solver.getObjVal(&tmp_val))
	auto dif1 = calculateDifference(tmp_val, org_val);
	safe_solver(inverseLastBranchConstraint(SOLVER_GREATER_EQUAL, 1, node->solver))
#ifdef TEST_LP_TIME
	beg1 = high_resolution_clock::now();
#endif
	safe_solver(node->solver.reoptimize())
#ifdef TEST_LP_TIME
	end1 = high_resolution_clock::now();
	cout << "time for reoptimize2: " << duration<double>(end1 - beg1).count() << endl;
#endif
	safe_solver(node->solver.getObjVal(&tmp_val))
	auto dif2 = calculateDifference(tmp_val, org_val);
	auto product = dif1 * dif2;
	Branch_Val[cnt++] = {edge, product};
	branch_lp[edge] = product;
	if (if_record_LP_improvement) {
	  lp_testing_improvement_down[edge].first += dif1;
	  ++lp_testing_improvement_down[edge].second;
	  lp_testing_improvement_up[edge].first += dif2;
	  ++lp_testing_improvement_up[edge].second;
	}
  }
  safe_solver(node->solver.delConstraints(1, &BeforeNumRow))
  safe_solver(node->solver.updateModel())

#ifdef USE_M_DYNAMICS
  getM(est_m, alpha, Branch_Val);//cannot be sorted!
  end = high_resolution_clock::now();
  eps = duration<double>(end - beg).count();
  getAverageT4LP(node, eps, (int)Branch_Val.size());
#endif

  if (if_writeBranch_pair) {
	int size = min(num, int(Branch_Val.size()));
	sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
	  return a.second > b.second;
	});
	branch_pair.resize(size);
	transform(Branch_Val.begin(),
			  Branch_Val.begin() + size,
			  branch_pair.begin(),
			  [](const auto &a) {
				return a.first;
			  });
  }

  if (if_record_sb_scores) {
	auto is_sorted = std::is_sorted(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
	  return a.second > b.second;
	});
	if (!is_sorted) {
	  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
		return a.second > b.second;
	  });
	}
	branch_pair_val = Branch_Val;
  }

  std::unordered_map<std::pair<int, int>, double, PairHasher> branch_LP_tmp;
  branch_LP_tmp.reserve(branch_pair.size());
  for (auto &edge : branch_pair) {
	branch_LP_tmp[edge] = branch_lp[edge];
  }
  branch_lp = branch_LP_tmp;
  safe_solver(node->solver.setEnvMethod(method))
  safe_solver(node->solver.reoptimize())
  end = high_resolution_clock::now();
  eps = duration<double>(end - beg).count();
  cout << "testLP: ";
  for (auto &edge : branch_pair) {
	cout << edge.first << "-" << edge.second << " " << branch_lp[edge] << ",";
  }
  cout << ", time= " << eps << endl;
}

#ifdef USE_M_DYNAMICS
void getM(double &m, double &alpha, const std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val) {
  if (abs(m) > TOLERANCE) return;
  int num = min((int)Branch_Val.size(), Config::BranPhase1);
  auto upper_m = int(num * MAX_M_COEFF);
  auto lower_m = int(max(num * MIN_M_COEFF, 1.));
  double average = accumulate(Branch_Val.begin(), Branch_Val.end(), 0.0, [](double sum, const auto &a) {
	return sum + a.second;
  }) / (double)Branch_Val.size();
  double step_size = 2 * average / (double)Branch_Val.size();
  vector<double> diff;
  diff.reserve(Branch_Val.size());
  for (int i = lower_m; i < upper_m; ++i) {
	double
		real_average = accumulate(Branch_Val.begin(), Branch_Val.begin() + i, 0.0, [](double sum, const auto &a) {
	  return sum + a.second;
	}) / i;
	double est_aver = step_size / 2 * (2 * (double)Branch_Val.size() - i) * (i + 1) / i;
	diff.emplace_back(abs(real_average - est_aver));
  }
  auto it = min_element(diff.begin(), diff.end()) - diff.begin() + lower_m;
  m = (double)it;
  alpha = min(m / (double)Branch_Val.size(), MAX_M_COEFF);
  cout << "m= " << m << " alpha= " << alpha << endl;
}

void CVRP::getAverageT4LP(BbNode *node, double t, int num) {
  double average_t = t / 2 / num;
  BbNode::updateState(average_t, node->t_for_one_lp, (int)node->brcs.size());
#if VERBOSE_MODE == 1
  cout << "eps= " << t << " average_t= " << average_t << " t_for_one_lp= " << node->t_for_one_lp << endl;
  cout << "useModelInPhase2 spent= " << t << "s" << endl;
#endif
}
#endif

#ifdef MASTER_VALVE_ML

void CVRP::useModelPhase1Separate(BbNode *node, int num) {
  if (branch_pair.size() == 1) {
	cout << "No Model Phase1: a single candidate" << endl;
	return;
  }
  getTrainingDataPhase1(node);
  auto local_ratio_pseudo =
	  ((double)branch_pair_from_pseudo.size()
		  / (double)(branch_pair_from_pseudo.size() + branch_pair_from_fractional.size()));
  int pseudo_size, frac_size;
  if (num == 1) {
	pseudo_size = branch_pair_from_pseudo.empty() ? 0 : 1;
	frac_size = 1 - pseudo_size;
  } else {
	pseudo_size = min((int)(num * local_ratio_pseudo), (int)branch_pair_from_pseudo.size());
	frac_size = min(num - pseudo_size, (int)branch_pair_from_fractional.size());
  }
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_pseudo(branch_pair_from_pseudo.size());
  transform(branch_pair_from_pseudo.begin(),
			branch_pair_from_pseudo.end(),
			Branch_Val_model_pseudo.begin(),
			[](const auto &edge) {
			  return make_pair(edge, 0.0);
			});
  ml.predict(Branch_Val_model_pseudo, 1);
  transform(Branch_Val_model_pseudo.begin(),
			Branch_Val_model_pseudo.begin() + pseudo_size,
			branch_pair.begin(),
			[](const auto &a) {
			  return a.first;
			});
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_fractional(
	  branch_pair_from_fractional.size());
  transform(branch_pair_from_fractional.begin(),
			branch_pair_from_fractional.end(),
			Branch_Val_model_fractional.begin(),
			[](const auto &edge) {
			  return make_pair(edge, 0.0);
			});
  ml.predict(Branch_Val_model_fractional, 1);
  transform(Branch_Val_model_fractional.begin(),
			Branch_Val_model_fractional.begin() + frac_size,
			branch_pair.begin() + pseudo_size,
			[](const auto &a) {
			  return a.first;
			});
  branch_pair.resize(num);
#if VERBOSE_MODE==1
  cout << "candidates in phase1: " << branch_pair.size() << endl;
#endif
}

void CVRP::getTrainingDataPhase1(BbNode *node) {
  if (branch_pair.size() == 1) {
	cout << "No Training Data Phase1: a single candidate" << endl;
	return;
  }
  double org_val = node->value;
  int BeforeNumRow = num_row;
  ml.collect_edge_related_features(node, org_val);
  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
  for (auto &edge : branch_pair) {
	getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
	ml.collect_variable_related_features(node, edge, solver_ind.data(), BeforeNumRow, (int)solver_ind.size(), org_val);
  }
}
#endif
