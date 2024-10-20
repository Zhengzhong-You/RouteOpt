//
// Created by You, Zhengzhong on 5/2/24.
//


#include "dual_smoothing.hpp"

using namespace std;

CVRP *DualSmoothing::cvrp{};
double DualSmoothing::alpha{};
int DualSmoothing::num_k{};
BbNode *DualSmoothing::node{};
double DualSmoothing::lp_value{};
double DualSmoothing::best_lagrange_bound{};
double DualSmoothing::most_negative_reduced_cost_in_lp{};
denseColMatrixXd DualSmoothing::lp_matrix{};//easy for adding and deleting columns
int DualSmoothing::most_negative_reduced_cost_col_idx{};
std::vector<double> DualSmoothing::most_negative_reduced_cost_col{};
std::vector<double> DualSmoothing::rhs_lp_matrix{};
RowVectorXd DualSmoothing::lp_obj{};
std::vector<double> DualSmoothing::best_lagrange_multipliers{};
std::vector<double> DualSmoothing::pi_in{};
std::vector<double> DualSmoothing::pi_out{};
std::vector<double> DualSmoothing::pi_price{};
bool DualSmoothing::if_mis_pricing{};
int DualSmoothing::counter_mis_pricing{};

void decreaseAlpha(double &alpha);

void increaseAlpha(double &alpha);

void DualSmoothing::reset() {
  alpha = 0.5;
  best_lagrange_bound = -numeric_limits<double>::max();
  best_lagrange_multipliers.clear();
  if_mis_pricing = false;
  counter_mis_pricing = 0;
}

void DualSmoothing::buildLPMatrix() {
  int num_row = cvrp->getNumRow();
  int num_col = cvrp->getNumCol();
  lp_obj.resize(num_col);
  safe_solver(node->getSolver().getObj(0, num_col, lp_obj.data()))

  lp_matrix = denseColMatrixXd::Zero(num_row, num_col);
  int numnzP;
  vector<int> vbeg(num_col + 1);
  vector<int> vind;
  vector<double> vval;
  safe_solver(node->getSolver().getVars(&numnzP, nullptr, nullptr, nullptr, 0, num_col))
  vind.resize(numnzP);
  vval.resize(numnzP);
  safe_solver(node->getSolver().getVars(&numnzP, vbeg.data(), vind.data(), vval.data(), 0, num_col))
  vbeg[num_col] = numnzP;
  for (int i = 0; i < num_col; ++i) {
	for (int j = vbeg[i]; j < vbeg[i + 1]; ++j) {
	  lp_matrix(vind[j], i) = vval[j];
	}
  }
  rhs_lp_matrix.resize(num_row);
  safe_solver(node->getSolver().getRhs(0, num_row, rhs_lp_matrix.data()))
}

void DualSmoothing::updateBestLagrangeBound() {
  double new_bound = lp_value + most_negative_reduced_cost_in_lp * num_k;
  if (new_bound > best_lagrange_bound) {
	best_lagrange_bound = new_bound;
	best_lagrange_multipliers = pi_price;
  }
}

void DualSmoothing::recordLPValue() {
  safe_solver(node->getSolver().getObjVal(&lp_value))
  initNumK();
}

void DualSmoothing::updatePiIn(int mode) {
  if (mode == MODE_IN_BEST_LAGRANGE_BOUND) {
	pi_in = best_lagrange_multipliers;
  } else if (mode == MODE_IN_LAST_ITERATION) {
	pi_in = pi_price;
  }
}

void DualSmoothing::updatePiOut() {
  pi_out = cvrp->pi4_labeling;
}

void DualSmoothing::updatePiPrice() {
  if (pi_in.empty()) {
	pi_in = pi_out;
	pi_price = pi_out;
	return;
  }
  pi_price.resize(pi_in.size());
  transform(pi_in.begin(),
			pi_in.end(),
			pi_out.begin(),
			pi_price.begin(),
			[](double a, double b) { return alpha * a + (1 - alpha) * b; });
  cvrp->pi4_labeling = pi_price;
}

void DualSmoothing::init(CVRP *pr_cvrp) {
  cvrp = pr_cvrp;
}

void DualSmoothing::addCol2LPMatrix(int ccnt_cnt,
									std::vector<size_t> &solver_beg,
									std::vector<int> &solver_ind,
									std::vector<double> &solver_val,
									std::vector<double> &solver_obj) {
  int old_num_col = (int)lp_matrix.cols();
  int num_col = old_num_col + ccnt_cnt;
  lp_matrix.conservativeResize(Eigen::NoChange, num_col);
  lp_matrix.block(0, old_num_col, lp_matrix.rows(), ccnt_cnt).setZero();
  lp_obj.conservativeResize(num_col);
  for (int i = 0; i < ccnt_cnt; ++i) {
	int col_idx = old_num_col + i;
	for (size_t j = solver_beg[i]; j < solver_beg[i + 1]; ++j) {
	  lp_matrix(solver_ind[j], col_idx) = solver_val[j];
	}
	lp_obj[col_idx] = solver_obj[i];
  }
  calculateMostNegativeReducedCostInLP();
  auto ptr = lp_matrix.col(most_negative_reduced_cost_col_idx).data();
  most_negative_reduced_cost_col.assign(ptr, ptr + lp_matrix.rows());
  updateBestLagrangeBound();
  cvrp->if_exact_labeling_cg ? updatePiIn(MODE_IN_BEST_LAGRANGE_BOUND) : updatePiIn(MODE_IN_LAST_ITERATION);
  adjustAlpha();
}

void DualSmoothing::deleteCol4LPMatrix(std::vector<int> &col_idx) {
  if (!lp_matrix.size()) return;
  if (!is_sorted(col_idx.begin(), col_idx.end())) {
	sort(col_idx.begin(), col_idx.end());
  }

  int old_num_col = (int)lp_matrix.cols();
  int delta = 0;
  auto stop_sign = col_idx.end() - 1;
  for (auto i = col_idx.begin(); i < stop_sign; ++i) {
	++delta;
	for (int j = *i + 1; j < *(i + 1); ++j) {
	  lp_matrix.col(j - delta) = lp_matrix.col(j);
	  lp_obj[j - delta] = lp_obj[j];
	}
  }
  ++delta;
  for (int j = *stop_sign + 1; j < old_num_col; ++j) {
	lp_matrix.col(j - delta) = lp_matrix.col(j);
	lp_obj[j - delta] = lp_obj[j];
  }
  int num_col = old_num_col - delta;
  lp_matrix.conservativeResize(Eigen::NoChange, num_col);
  lp_obj.conservativeResize(num_col);
}

void DualSmoothing::calculateMostNegativeReducedCostInLP() {
  Eigen::Map<Eigen::RowVectorXd> pi(pi_price.data(), (int)pi_price.size());
  auto rc = lp_obj - pi * lp_matrix;
  Eigen::Index min_idx;
  safe_eigen(most_negative_reduced_cost_in_lp = rc.minCoeff(&min_idx);
				 most_negative_reduced_cost_col_idx = (int)min_idx;)
}

void DualSmoothing::adjustAlpha() {
  if (if_mis_pricing) {
	alpha = max(0., 1 - counter_mis_pricing * (1 - alpha));
	return;
  }
  if (std::sqrt(std::inner_product(pi_out.begin(), pi_out.end(), pi_in.begin(), 0.0,
								   std::plus<>(),
								   [](double a, double b) {
									 return pow(a - b, 2);
								   })) / lp_value < NORM_TOLERANCE) {
	alpha = 0;
  }

  vector<double> pi_out_in(pi_out.size());
  transform(pi_out.begin(),
			pi_out.end(),
			pi_in.begin(),
			pi_out_in.begin(),
			[](double a, double b) { return a - b; });
  vector<double> subgradient(pi_out.size());
  transform(rhs_lp_matrix.begin(),
			rhs_lp_matrix.end(),
			most_negative_reduced_cost_col.begin(),
			subgradient.begin(),
			[](double a, double b) { return a - num_k * b; });//num_vehicle
  std::inner_product(pi_out_in.begin(),
					 pi_out_in.end(),
					 subgradient.begin(),
					 0.0,
					 plus<>(),
					 [](double a, double b) { return a * b; }) < 0 ? decreaseAlpha(alpha) : increaseAlpha(alpha);
}

void DualSmoothing::tellIfMisPrice(int &ccnt) {
  bool mark;
  ((ccnt == 0) && (most_negative_reduced_cost_in_lp < RC_TOLERANCE)) ? mark = true : mark = false;
  if_mis_pricing = mark;
  if (mark) {
	++counter_mis_pricing;
	++ccnt;
	adjustAlpha();
  } else counter_mis_pricing = 0;
}

void DualSmoothing::updateNode(BbNode *pr_node) {
  node = pr_node;
}

void DualSmoothing::initNumK() {
  vector<double> x(cvrp->getNumCol());
  safe_solver(node->getSolver().getX(0, cvrp->getNumCol(), x.data()))
  num_k = ceil(accumulate(x.begin(), x.end(), 0.0) - TOLERANCE);
}

void DualSmoothing::freeLPMatrix() {
  lp_matrix.resize(0, 0);
  lp_obj.resize(0);
  most_negative_reduced_cost_col.clear();
  rhs_lp_matrix.clear();
  best_lagrange_multipliers.clear();
  pi_in.clear();
  pi_out.clear();
  pi_price.clear();
}

void DualSmoothing::changeRCStd() {
  cvrp->rc_std = most_negative_reduced_cost_in_lp;
}

void DualSmoothing::printInfo() {
  cout << "alpha: " << alpha << endl;
  cout << "most negative reduced cost in LP: " << most_negative_reduced_cost_in_lp << endl;
  cout << "best lagrange bound: " << best_lagrange_bound << endl;
  cout << "most negative reduced cost column index: " << most_negative_reduced_cost_col_idx << endl;
  cout << "pi out: " << endl;
  for (auto &i : pi_out) {
	cout << i << " ";
  }
  cout << endl;
  cout << "pi price: " << endl;
  for (auto &i : pi_price) {
	cout << i << " ";
  }
  cout << endl;
  cout << "if mis pricing: " << if_mis_pricing << endl;
  cout << "counter mis pricing: " << counter_mis_pricing << endl;
  cout << "lp value: " << lp_value << endl;
  cout << BIG_PHASE_SEPARATION;
}

void decreaseAlpha(double &alpha) {
  alpha = max(0.0, alpha - 0.1);
}

void increaseAlpha(double &alpha) {
  alpha = 0.1 + 0.9 * alpha;
}
