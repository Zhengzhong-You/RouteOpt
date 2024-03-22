
#include "CVRP.hpp"

using namespace std;
using namespace chrono;


#if  defined(FIND_DIF_LP) or defined(FIND_INDICATOR)
void CVRP::copyModel2GetDual(BbNode *node, bool dir) {
  if (node->r1cs.empty()) return;
  Solver local_solver{};
  local_solver.model = node->solver.copyModel();
  vector<double> new_var_coeff(num_row);
  safe_solver(node->solver.getRhs(0, num_row, new_var_coeff.data()))
  for (auto &i : new_var_coeff) i *= -1;
  double new_var_obj;
  safe_solver(node->solver.getObjVal(&new_var_obj))
  new_var_obj -= TOLERANCE;
  new_var_obj *= -1;
  vector<int> new_var_ind;
  vector<double> new_var_val;
  for (int i = 0; i < num_row; ++i) {
	if (abs(new_var_coeff[i]) < TOLERANCE) continue;
	new_var_ind.emplace_back(i);
	new_var_val.emplace_back(new_var_coeff[i]);
  }

  safe_solver(local_solver.addVar(new_var_ind.size(),
								  new_var_ind.data(),
								  new_var_val.data(),
								  new_var_obj,
								  0.,
								  SOLVER_INFINITY,
								  SOLVER_CONTINUOUS,
								  nullptr))
  vector<double> new_model_rhs(num_row, 0);
  auto max_rhs = int(Config::MaxRowRank1 / 2);

#if LIMITED_MEMORY_TYPE == 1
  for (auto &r1c : node->r1cs) {
	new_model_rhs[r1c.idx_r1c] = 1000 * (max_rhs - r1c.rhs) + 1 + (double)r1c.mem.size();
  }
#elif LIMITED_MEMORY_TYPE == 2
  for (auto &r1c : node->r1cs) {
	new_model_rhs[r1c.idx_r1c] = 1000 * (max_rhs - r1c.rhs) + 1 + pow((double)r1c.arc_mem.size(), 2);
  }
#endif

  if (!dir) {
	for (int i = 0; i < num_row; ++i) {
	  new_model_rhs[i] = -new_model_rhs[i];
	}
  }

  safe_solver(local_solver.setRhs(0, num_row, new_model_rhs.data()))
  int env_method;
  bool if_changed = false;
  safe_solver(local_solver.getEnvMethod(&env_method))
  if (env_method != SOLVER_DUAL_SIMPLEX) {
	safe_solver(local_solver.setEnvMethod(SOLVER_DUAL_SIMPLEX))
	if_changed = true;
  }
  time_point<high_resolution_clock> start, end;
  start = high_resolution_clock::now();
  safe_solver(local_solver.reoptimize())
  end = high_resolution_clock::now();
  good.emplace_back(duration<double>(end - start).count());
  int status;
  safe_solver(node->solver.getStatus(&status))
  if (status != SOLVER_UNBOUNDED && status != SOLVER_INF_OR_UNBD) {
	safe_solver(local_solver.getDual(0, num_row, pi4_labeling.data()))
  } else {
	throw runtime_error("dual is unbounded or infeasible");
  }

  if (if_changed) {
	safe_solver(local_solver.setEnvMethod(env_method))
  }
  safe_solver(local_solver.freeModel())
}
#endif



void CVRP::reviseSubPricingModel(BbNode *node) {
  if (sub_pricing_solver.model) sub_pricing_solver.freeModel();
  if (node->r1cs.empty()) return;
  vector<double> old_rhs(num_row);
  safe_solver(node->solver.getRhs(0, num_row, old_rhs.data()))
#ifdef CHECK_FIRST_COLUMN
  vector<double> first_column_coeff(num_row);
  for (int i = 0; i < num_row; ++i) {
	double coeff;
	safe_solver(node->solver.getCoeff(i, 0, coeff))
	if (abs(coeff - old_rhs[i]) > TOLERANCE) {
	  cout << "coeff= " << coeff << ", rhs= " << old_rhs[i] << endl;
	  cout << "i= " << i << endl;
	  throw runtime_error("first column coefficient is not equal to rhs");
	}
  }
#endif
  sub_pricing_solver.model = node->solver.copyModel();
  vector<int> cind;
  vector<double> val;
  for (int i = 0; i < num_row; ++i) {
	if (abs(old_rhs[i]) < TOLERANCE) continue;
	cind.emplace_back(i);
	val.emplace_back(-old_rhs[i]);
  }
  vector<int> vind(cind.size(), 0);
  safe_solver(sub_pricing_solver.changeCoeffs(cind.size(), cind.data(), vind.data(), val.data()))
  vector<double> new_model_rhs(num_row, 0);
  auto max_rhs = int(Config::MaxRowRank1 / 2);

#if LIMITED_MEMORY_TYPE == 1
  for (auto &r1c : node->r1cs) {
	new_model_rhs[r1c.idx_r1c] = 1000 * (max_rhs - r1c.rhs) + 1 + (double)r1c.mem.size();
  }
#elif LIMITED_MEMORY_TYPE == 2
  for (auto &r1c : node->r1cs) {
	new_model_rhs[r1c.idx_r1c] = 1000 * (max_rhs - r1c.rhs) + 1 + pow((double)r1c.arc_mem.size(), 2);
  }
#endif

  safe_solver(sub_pricing_solver.setRhs(0, num_row, new_model_rhs.data()))
  safe_solver(sub_pricing_solver.updateModel())
}

void CVRP::changeModelForBetterDual(BbNode *node) {
  if (node->r1cs.empty() || !if_exact_labeling_cg || !sub_pricing_solver.model) {
	safe_solver(node->solver.getDual(0, num_row, pi4_labeling.data()))
	return;
  }
#ifdef CHANGE_DUAL
  double new_obj;
  double old_obj;
  safe_solver(node->solver.getObjVal(&new_obj))
  safe_solver(sub_pricing_solver.getObj(0, 1, &old_obj))
  old_obj *= -1;
  old_obj += TOLERANCE;
  if (abs(new_obj - old_obj) > TOLERANCE) {
	new_obj -= TOLERANCE;
	new_obj *= -1;
	safe_solver(sub_pricing_solver.changeObj(0, 1, &new_obj))
  }
  int env_method;
  bool if_changed = false;
  safe_solver(sub_pricing_solver.getEnvMethod(&env_method))
  if (env_method != SOLVER_PRIMAL_SIMPLEX) {
	safe_solver(sub_pricing_solver.setEnvMethod(SOLVER_PRIMAL_SIMPLEX))
	if_changed = true;
  }

#ifdef FIND_DIF_LP
	time_point<high_resolution_clock> start, end;
	start = high_resolution_clock::now();
#endif
  safe_solver(sub_pricing_solver.reoptimize())
#ifdef FIND_DIF_LP
  end = high_resolution_clock::now();
  better.emplace_back(duration<double>(end - start).count());
#endif

  int status;
  safe_solver(sub_pricing_solver.getStatus(&status))
  if (status != SOLVER_UNBOUNDED && status != SOLVER_INF_OR_UNBD) {
	safe_solver(sub_pricing_solver.getDual(0, num_row, pi4_labeling.data()))
  } else {
	throw runtime_error("dual is unbounded or infeasible");
  }

  if (if_changed) {
	safe_solver(sub_pricing_solver.setEnvMethod(env_method))
  }
#ifdef FIND_DIF_LP
  copyModel2GetDual(node, true);

  cout << "good= " << good.back() << ", better= " << better.back() << endl;
#endif

#ifdef CHECK_DUAL
  cout << "check if the dual is optimal dual in the original model" << endl;
  sparseColMatrixXd A(num_row, num_col);
  size_t numnz;
  safe_solver(node->solver.XgetVars(&numnz, nullptr, nullptr, nullptr, 0, num_col))
  vector<size_t> vbeg(num_col + 1);
  vector<int> vind(numnz);
  vector<double> vval(numnz);
  vector<Eigen::Triplet<double>> triplet(numnz);
  safe_solver(node->solver.XgetVars(&numnz, vbeg.data(), vind.data(), vval.data(), 0, num_col))
  vbeg[num_col] = numnz;
  numnz = 0;
  for (size_t i = 0; i < num_col; ++i) {
	for (size_t j = vbeg[i]; j < vbeg[i + 1]; ++j) {
	  triplet[numnz++] = Eigen::Triplet<double>(vind[j], i, vval[j]);
	}
  }
  A.setFromTriplets(triplet.begin(), triplet.end());
  vector<double> cost(num_col);
  safe_solver(node->solver.getObj(0, num_col, cost.data()))
  Eigen::Map<Eigen::RowVectorXd> eigen_cost(cost.data(), num_col);
  Eigen::Map<Eigen::RowVectorXd> pi(pi4_labeling.data(), num_row);
  RowVectorXd rc = eigen_cost - pi * A;
  for (int i = 0; i < num_col; ++i) {
	if (rc[i] < -TOLERANCE) {
	  cout << "rc[" << i << "]= " << rc[i] << endl;
	  cout << "dual is not optimal" << endl;
	  throw runtime_error("dual is not optimal");
	}
  }
  cout << "sub_dual= ";
  for (int i = 0; i < num_row; ++i) {
	cout << pi4_labeling[i] << ", ";
  }
  cout << endl;
  vector<double> old_pi(num_row);
  safe_solver(node->solver.getDual(0, num_row, old_pi.data()))
  cout << "old_dual= ";
  for (int i = 0; i < num_row; ++i) {
	cout << old_pi[i] << ", ";
  }
  cout << endl;
#endif
#else
  safe_solver(node->solver.getDual(0, num_row, pi4_labeling.data()))
#endif
}

