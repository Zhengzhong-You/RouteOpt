
#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

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

void CVRP::giveDynamicM(BbNode *node, int &num) {
  if (node->index == 0) {
	num = min(ml.give_initial_screening_num_candidates(), int(branch_pair.size()));
	esti_k = num;
	cout << "we directly use the max n= " << num << endl;
	return;
  }
  if (is_use_full_k) {
	num = min(int(esti_m), int(branch_pair.size()));
	esti_k = num;
	cout << "we directly use the full n= " << num << endl;
	return;
  }
  double t_for_one_lp = node->t_for_one_lp;
  double r_star = node->geo_r_star;
  double c = node->c;
  double n = ml.give_initial_screening_num_candidates();
  double G = ub - node->value;

  double A = r_star * t_for_one_lp * pow(n + 1, 2);
  double B = t_for_one_lp * (log(2) * (f - G) * n * (esti_m + 1) + 2 * r_star * (n - esti_m) * (n + 1));
  double C = c * log(2) * (f - G) * n * (esti_m + 1) + r_star * t_for_one_lp * pow(esti_m - n, 2);

  double k1 = (-B + sqrt(pow(B, 2) - 4 * A * C)) / (2 * A);
  int k1_up = (int)max(min(ceil(k1), double(esti_m)), 1.);
  int k1_down = (int)max(min(floor(k1), double(esti_m)), 1.);

  cout << "___________________" << endl;
  cout << "G= " << G << " f= " << f << " c= " << c << " t_for_one_lp= " << t_for_one_lp << " r_star= " << r_star
	   << " n= " << n
	   << " m= " << esti_m << " k1= " << k1 << endl;

  unordered_set<int> tmp{k1_up, k1_down};

  double least_T_N_ratio = numeric_limits<double>::max();
  for (auto i : tmp) {
	double
		T_N_ratio =
		(c + i * t_for_one_lp) * pow(2, (G - f) / (r_star * ((esti_m + 1) / n * i / (i + 1) + 1 - esti_m / n)));
	if (T_N_ratio == numeric_limits<double>::infinity()) {
	  num = k1_up;
	  break;
	}
	if (T_N_ratio < least_T_N_ratio) {
	  least_T_N_ratio = T_N_ratio;
	  num = i;
	}
	cout << "i= " << i << " T_N_ratio= " << T_N_ratio << ", ";
  }
  num = min(num, int(esti_m));
  num = min(num, int(branch_pair.size()));
  esti_k = num;
  cout << " num= " << num << endl;
}

void CVRP::evaluateM1() {
  esti_m = double(find(cp_branch_pair.begin(), cp_branch_pair.end(), branch_pair[0]) - cp_branch_pair.begin()) + 1;
#ifdef SOLVER_VRPTW
  esti_m = max(21., esti_m);
#else
  esti_m = max(20., esti_m);
#endif
  esti_m = min(esti_m, double(ml.give_initial_screening_num_candidates()));
  cp_branch_pair.clear();
  cout << "m= " << esti_m << endl;
}

#endif
#endif