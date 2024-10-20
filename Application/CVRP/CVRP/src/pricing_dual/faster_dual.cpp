//
// Created by You, Zhengzhong on 3/8/24.
//

#include "faster_dual.hpp"
#ifdef COMBINE_FASTER_SMOOTHING
#include "combine_faster_smoothing.hpp"
#endif

using namespace std;
using namespace chrono;


CVRP *FasterDual::cvrp{};
BbNode *FasterDual::node{};
Solver FasterDual::sub_pricing_solver{};

void checkFirstColumn(BbNode *node, int num_row, const vector<double> &old_rhs);

void checkDual(CVRP *cvrp, BbNode *node);

void FasterDual::changeModel4BetterDual() {
  auto num_row = cvrp->getNumRow();
  auto &pi4_labeling = cvrp->pi4_labeling;
  if (!sub_pricing_solver.model) {
	safe_solver(node->getSolver().getDual(0, num_row, pi4_labeling.data()))
	return;
  }
  double new_obj;
  double old_obj;
  safe_solver(node->getSolver().getObjVal(&new_obj))
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
  safe_solver(sub_pricing_solver.reoptimize(SOLVER_PRIMAL_SIMPLEX))
  int status;
  safe_solver(sub_pricing_solver.getStatus(&status))
  if (status == SOLVER_OPTIMAL) {
	safe_solver(sub_pricing_solver.getDual(0, num_row, pi4_labeling.data()))
	combine_faster_call({
						  double val_;
						  safe_solver(sub_pricing_solver.getObjVal(&val_))
						  CombineFasterSmoothing::getLP2Value(val_);
						})
  } else {
	cout << "sub_pricing_solver status= " << status << endl;
  }
  if (if_changed) {
	safe_solver(sub_pricing_solver.setEnvMethod(env_method))
  }
#ifdef CHECK_DUAL
  checkDual(cvrp, node);
#endif
}

void FasterDual::reviseSubPricingModel() {
  if (sub_pricing_solver.model) return;
  if (cvrp->stopChangingDualCondition(node)) {
	return;
  }
  int num_row = cvrp->getNumRow();
  vector<double> old_rhs(num_row);
  safe_solver(node->getSolver().getRhs(0, num_row, old_rhs.data()))
#ifdef CHECK_FIRST_COLUMN
  checkFirstColumn(node, num_row, old_rhs);
#endif
  sub_pricing_solver.model = node->getSolver().copyModel();
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

  cvrp->getNewRhs4ChangingDual(node, new_model_rhs);

  combine_faster_call(CombineFasterSmoothing::getLP2Obj(new_model_rhs))
  safe_solver(sub_pricing_solver.setRhs(0, num_row, new_model_rhs.data()))
  safe_solver(sub_pricing_solver.updateModel())
}

void FasterDual::addCols2SubPricingModel(int ccnt_cnt, size_t nzcnt, vector<size_t> &solver_beg,
										 vector<int> &solver_ind, vector<double> &solver_val,
										 vector<double> &solver_obj) {
  if (sub_pricing_solver.model) {
	safe_solver(sub_pricing_solver.XaddVars(ccnt_cnt,
											nzcnt,
											solver_beg.data(),
											solver_ind.data(),
											solver_val.data(),
											solver_obj.data(),
											nullptr,
											nullptr,
											nullptr,
											nullptr))
	safe_solver(sub_pricing_solver.updateModel())
  }
}

void FasterDual::init(CVRP *pr_cvrp) {
  cvrp = pr_cvrp;
}

void FasterDual::updateNode(BbNode *pr_node) {
  node = pr_node;
}

void FasterDual::freeSubPricingModel() {
  if (sub_pricing_solver.model) sub_pricing_solver.freeModel();
}

void FasterDual::deleteCols4SubPricingModel(std::vector<int> &col_idx) {
  if (sub_pricing_solver.model) {
	cout << "revise sub_pricing_solver as well!" << endl;
	safe_solver(sub_pricing_solver.delVars(col_idx.size(), col_idx.data()))
	safe_solver(sub_pricing_solver.updateModel())
	int local_col;
	safe_solver(sub_pricing_solver.getNumCol(&local_col))
	if (local_col != cvrp->getNumCol()) {
	  cout << "local_col= " << local_col << " num_col= " << cvrp->getNumCol() << endl;
	  throw runtime_error("Error in sub_pricing_solver");
	}
  }
}

void checkFirstColumn(BbNode *node, int num_row, const vector<double> &old_rhs) {
  vector<double> first_column_coeff(num_row);
  for (int i = 0; i < num_row; ++i) {
	double coeff;
	safe_solver(node->getSolver().getCoeff(i, 0, coeff))
	if (abs(coeff - old_rhs[i]) > TOLERANCE) {
	  cout << "coeff= " << coeff << ", rhs= " << old_rhs[i] << endl;
	  cout << "i= " << i << endl;
	  throw runtime_error("first column coefficient is not equal to rhs");
	}
  }
}

void checkDual(CVRP *cvrp, BbNode *node) {
  cout << "check if the dual is optimal dual in the original model" << endl;
  int num_row = cvrp->getNumRow(), num_col = cvrp->getNumCol();
  sparseColMatrixXd A(num_row, num_col);
  size_t numnz;
  safe_solver(node->getSolver().XgetVars(&numnz, nullptr, nullptr, nullptr, 0, num_col))
  vector<size_t> vbeg(num_col + 1);
  vector<int> vind(numnz);
  vector<double> vval(numnz);
  vector<Eigen::Triplet<double>> triplet(numnz);
  safe_solver(node->getSolver().XgetVars(&numnz, vbeg.data(), vind.data(), vval.data(), 0, num_col))
  vbeg[num_col] = numnz;
  numnz = 0;
  for (int i = 0; i < num_col; ++i) {
	for (size_t j = vbeg[i]; j < vbeg[i + 1]; ++j) {
	  triplet[numnz++] = Eigen::Triplet<double>(vind[j], i, vval[j]);
	}
  }
  safe_eigen(A.setFromTriplets(triplet.begin(), triplet.end());)
  vector<double> cost(num_col);
  safe_solver(node->getSolver().getObj(0, num_col, cost.data()))
  Eigen::Map<Eigen::RowVectorXd> eigen_cost(cost.data(), num_col);
  auto &pi4_labeling = cvrp->pi4_labeling;
  Eigen::Map<Eigen::RowVectorXd> pi(pi4_labeling.data(), num_row);
  RowVectorXd rc;
  safe_eigen(rc = eigen_cost - pi * A;)
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
  safe_solver(node->getSolver().getDual(0, num_row, old_pi.data()))
  cout << "old_dual= ";
  for (int i = 0; i < num_row; ++i) {
	cout << old_pi[i] << ", ";
  }
  cout << endl;
}