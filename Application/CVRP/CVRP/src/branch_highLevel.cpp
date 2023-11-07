
#include "CVRP.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <functional>

using namespace std;
using namespace std::chrono;
using namespace boost::math::tools;
using namespace std::placeholders;
using namespace boost::math;

void CVRP::doSB(BbNode *node, pair<int, int> &info) {
  if (node->is_terminated) return;
  branch_pair.clear();
#ifdef MASTER_VALVE_ML
  ml.read_enu_state(if_in_enu_state);
  switch (ML_STATE) {
	case 1:generateModelInPhase1(node);
	  break;
	case 2:generateModelInPhase2(node);
	  break;
	case 3: useMLInGeneralFramework(node);
//	  findAccML(node);
	  break;
#ifdef USE_M_DYNAMICS
	case 4:useMLVaryLPTesting(node);
	  break;
#endif
	default:throw runtime_error("doSB: ML state not found");
  }
  ml.edge_tmp_info.clear();
#else
  useDefaultSB(node);
#endif
  info = branch_pair[0];
  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++branch_choice[info];
  ++branch_times;
}

void CVRP::useDefaultSB(BbNode *node) {
  int num1, num2;
  double frac;

  if (if_in_enu_state) {
	num1 = Config::BranPhase1InEnu;
	num2 = Config::BranPhase2InEnu;
  } else {
	num1 = Config::BranPhase1;
	num2 = Config::BranPhase2;
  }
  frac = Config::Frac4sudoCostBranPhase0;

  initialScreen(node, true, true, num1, frac);
  testLP(node, num2, true);
  testCG(node, false, false, true);
}

#ifdef MASTER_VALVE_ML

void CVRP::findAccML(BbNode *node) {
  if (if_in_enu_state) throw runtime_error("findAccInitialML: not in enu state");
  initialScreen(node, true, true, ml.give_initial_screening_num_candidates(), Config::Frac4sudoCostBranPhase0);
  int num = min(ml.give_testing_num_candidates(node->tree_level), int(branch_pair.size()));
  auto commonSet = branch_pair;
  useModelPhase1Separate(node, num);
  useModelInPhase2(node, Config::BranPhase2);
  auto initialS = branch_pair;
  branch_pair.clear();
  branch_pair = commonSet;
  testCG(node, true, true, false, true);
  auto overallS = branch_pair_val;
  self_mkdir("testAcc");
  ofstream out("testAcc/" + file_name + ".txt", ios::app);
  out << "---------------------" << endl;
  out << "MLS: " << endl;
  for (auto &edge : initialS) {
	out << edge.first << " " << edge.second << ", ";
  }
  out << endl;
  out << "overallS: " << endl;
  for (auto &edge : overallS) {
	out << edge.first.first << " " << edge.first.second << " " << edge.second << ", ";
  }
  out << endl;
  out.close();
}

void CVRP::generateModelInPhase1(BbNode *node) {
  initialScreen(node, true, true, ml.give_initial_screening_num_candidates(), Config::Frac4sudoCostBranPhase0);
  getTrainingDataPhase1(node);

  simulateWriteLPPseudoCost(node);
  testCG(node, true, true, false, true);
  ml.write_training_lp_file();
}

void CVRP::generateModelInPhase2(BbNode *node) {
  initialScreen(node, true, true, ml.give_initial_screening_num_candidates(), Config::Frac4sudoCostBranPhase0);
  int num = min(ml.give_testing_num_candidates(node->tree_level), int(branch_pair.size()));
  useModelPhase1Separate(node, num);
  getTrainingDataInPhase2(node);
  testCG(node, true, true, false, true);
  ml.write_training_exact_file();
}

void CVRP::useMLInGeneralFramework(BbNode *node) {
  initialScreen(node, true, true, ml.give_initial_screening_num_candidates(), Config::Frac4sudoCostBranPhase0);
  int num = min(ml.give_testing_num_candidates(node->tree_level), int(branch_pair.size()));
  useModelPhase1Separate(node, num);
  num = if_in_enu_state ? Config::BranPhase2InEnu : Config::BranPhase2;
  useModelInPhase2(node, num);
  testCG(node, false, false, true);
}

void CVRP::simulateWriteLPPseudoCost(BbNode *node) {
  auto tmp_pair = branch_pair;
  branch_pair.clear();
  auto ratio_pseudo =
	  ((double)branch_pair_from_pseudo.size()
		  / (double)(branch_pair_from_pseudo.size() + branch_pair_from_fractional.size()));
  int num = min(ml.give_testing_num_candidates(node->tree_level), int(tmp_pair.size()));
  int pseudo_size = (int)(num * ratio_pseudo);
  int frac_size = num - pseudo_size;
  transform(branch_pair_from_pseudo.begin(), branch_pair_from_pseudo.begin() + pseudo_size,
			back_inserter(branch_pair),
			[](pair<int, int> &p) { return make_pair(p.first, p.second); });
  transform(branch_pair_from_fractional.begin(),
			branch_pair_from_fractional.begin() + frac_size,
			back_inserter(branch_pair),
			[](pair<int, int> &p) { return make_pair(p.first, p.second); });
  testLP(node, 0, false);
  branch_pair = tmp_pair;
}

#ifdef USE_M_DYNAMICS
void CVRP::useMLVaryLPTesting(BbNode *node) {
  initialScreen(node, true, true, ml.give_initial_screening_num_candidates(), Config::Frac4sudoCostBranPhase0);
  int num;
  giveDynamicM(node, num);
  useModelPhase1Separate(node, num);
  if (node->index == 0) cp_branch_pair = branch_pair;
  num = if_in_enu_state ? Config::BranPhase2InEnu : Config::BranPhase2;
  useModelInPhase2(node, num);
  testCG(node, false, false, true);
  if (node->index == 0) evaluateM1();
}

double getLBTk(double omega, double alpha, double B, double k) {
  return (omega + k) * pow(B, 1 / (1 - alpha / (k + 1)));
}

double getUBTk(double omega, double alpha, double B, double k) {
  double exponent = alpha / (1.0 - alpha) * log(B);
  double omega_term = boost::math::gamma_p(k, exponent);
  return (omega + k) * k * pow(B, 1 / (1 - alpha)) * pow(exponent, -k) * omega_term;
}

double getRealTk(double omega, double alpha, double B, double k) {
  using namespace boost::math::quadrature;

  // Define the integrand directly within the getRealTk function
  auto integrand = [&](double x) {
	return k * std::pow(x, k - 1) * std::pow(B, 1.0 / (alpha * x + 1 - alpha));
  };

  double error;
  double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
	  integrand,
	  0.0,
	  1.0,
	  15,
	  std::numeric_limits<double>::epsilon(),
	  &error
  );

  return result * (omega + k);
}

void getRange(double omega,
			  double alpha,
			  double B,
			  double target,
			  double k1,
			  double est_m,
			  int &lb_k,
			  int &ub_k) {
  std::function<double(double)> bound_equation =
	  [omega, alpha, B, target](double k) {
		return getLBTk(omega, alpha, B, k) - target;
	  };

  // Define the search ranges
  std::array<std::pair<double, double>, 2> search_ranges = {
	  std::make_pair(1., k1),
	  std::make_pair(k1, est_m)
  };

  boost::uintmax_t max_iterations = 1000; // Maximum number of iterations

  int cnt = 0;
  for (const auto &range : search_ranges) {
	try {
	  auto result = boost::math::tools::toms748_solve(
		  bound_equation,
		  range.first,
		  range.second,
		  boost::math::tools::eps_tolerance<double>(6),
		  max_iterations
	  );
	  if (cnt == 0) {
		lb_k = int(ceil(result.second));
	  } else {
		ub_k = int(floor(result.first));
	  }
	} catch (const std::exception &e) {
	  if (cnt == 0) {
		lb_k = 1;
	  } else {
		ub_k = int(est_m);
	  }
	}
	++cnt;
  }
}

void CVRP::giveDynamicM(BbNode *node, int &num) {
  if (node->index == 0) {
	num = min(ml.give_initial_screening_num_candidates(), int(branch_pair.size()));
	opt_k = num;
	cout << "we directly use the max n= " << num << endl;
	return;
  }
  double beta = (ub - node->value - f) / node->geo_r_star;
  if (is_use_full_k || beta > 16) {
	num = min(int(est_m), int(branch_pair.size()));
	opt_k = num;
	cout << "we directly use the full m= " << num << endl;
	return;
  }

  double omega = node->c / node->t_for_one_lp;
  double B = pow(2, beta);

  double alpha_beta = alpha * beta * log(sqrt(2));
  double k1 = sqrt(alpha_beta * (alpha_beta + 2 * alpha + 2 * omega - 2)) + alpha_beta + alpha - 1;
  int k1_up = (int)max(min(ceil(k1), double(est_m)), 1.);
  int k1_down = (int)max(min(floor(k1), double(est_m)), 1.);
  cout << "alpha= " << alpha << " beta= " << beta << " omega= " << omega;

  if (k1_up == k1_down) {
	k1 = k1_up;
  } else {
	double k1_up_val = getLBTk(omega, alpha, B, k1_up);
	double k1_down_val = getLBTk(omega, alpha, B, k1_down);
	if (k1_up_val < k1_down_val) {
	  k1 = k1_up;
	} else {
	  k1 = k1_down;
	}
  }
  cout << " k1= " << k1 << endl;

  //calculate the UB
  double ub_tk = getUBTk(omega, alpha, B, k1);
  int lb_k, ub_k;
  getRange(omega, alpha, B, ub_tk, k1, est_m, lb_k, ub_k);
  cout << "lb_k= " << lb_k << " ub_k= " << ub_k << endl;
  if (lb_k == ub_k) {
	num = lb_k;
  } else {
	double real_tk = getRealTk(omega, alpha, B, k1);
	getRange(omega, alpha, B, real_tk, k1, est_m, lb_k, ub_k);
	cout << "revise: lb_k= " << lb_k << " ub_k= " << ub_k << endl;
	if (lb_k == ub_k) {
	  num = lb_k;
	} else {
	  vector<double> tmp;
	  for (int i = lb_k; i <= ub_k; ++i) {
		tmp.emplace_back(getRealTk(omega, alpha, B, i));
	  }
	  auto it = min_element(tmp.begin(), tmp.end()) - tmp.begin() + lb_k;
	  num = (int)it;
	}
  }

  num = min(num, int(est_m));
  num = min(num, int(branch_pair.size()));
  opt_k = num;
  cout << "opt_k= " << opt_k << endl;
}

void CVRP::evaluateM1() {
  est_m = double(find(cp_branch_pair.begin(), cp_branch_pair.end(), branch_pair[0]) - cp_branch_pair.begin()) + 1;
#ifdef SOLVER_VRPTW
  est_m = max(21., est_m);
#else
  est_m = max(20., est_m);
#endif
  est_m = min(est_m, double(ml.give_initial_screening_num_candidates()));
  cp_branch_pair.clear();
  alpha = min(est_m / double(ml.give_initial_screening_num_candidates()), 0.9);
  cout << "m= " << est_m << " alpha= " << alpha << endl;
}
#endif
#endif