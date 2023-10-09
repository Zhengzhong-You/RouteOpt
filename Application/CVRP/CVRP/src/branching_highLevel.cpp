//
// Created by Zhengzhong You on 3/28/23.
//

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

#if do_SB_mode == 1
void CVRP::do_SB(BBNODE *node, pair<int, int> &info) {
  if (node->if_terminated) return;
  Branch_pair.clear();
#ifdef MASTER_VALVE_ML
  ml.readEnuState(If_in_Enu_State);
//  bool if_sep = true;
  //although generate data can be sep, but when use the model, the sep must be on!
  switch (ML_state) {
	case 1:generateModelInPhase1(node, true);
	  break;
	case 2:generateModelInPhase2(node, true);
	  break;
	case 3: useMLInGeneralFramework(node, true);
//	  findAccML(node);
	  break;
	case 4:onlyUseStageOneModel(node, true);
	  break;
#ifdef useM_dynamically
	  case 5:useMLVaryLPTesting(node);//sep is determined
		break;
#endif
	default:throw runtime_error("do_SB: ML state not found");
  }
  ml.EdgeTmpInfo.clear();
#elif defined(random_selection)
  randomSelection(node);
#else
  useDefaultSB(node);
//  findAccInitial(node);
//  randomPickLP(node, 20);
//  initialPickLP(node, 20);
#endif
  info = Branch_pair[0];
  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++BranchChoice[info];
  ++BranchTimes;
}
#endif

void CVRP::useDefaultSB(BBNODE *node) {
  int num1, num2;
  double frac;

  if (If_in_Enu_State) {
	num1 = CONFIG::BranPhase1InEnu;
	num2 = CONFIG::BranPhase2InEnu;
  } else {
	num1 = CONFIG::BranPhase1;
	num2 = CONFIG::BranPhase2;
  }
  frac = CONFIG::Frac4sudoCostBranPhase0;

  InitialScreening(node, true, true, num1, frac);
  LPTesting(node, num2, true);
  CGTesting(node, false, false, true);
}

void CVRP::findAccInitial(BBNODE *node) {
  int num1, num2;
  double frac;

  num1 = 100;
  frac = CONFIG::Frac4sudoCostBranPhase0;

  InitialScreening(node, true, true, num1, frac);
  auto initialS = Branch_pair;
  Branch_pair.clear();
  constructBranchingSet(node);
//  LPTesting(node, num2, true);
  CGTesting(node, false, true, true, true);
  auto overallS = Branch_Pair_Val;
  self_mkdir("testAcc");
  ofstream out("testAcc/" + FileName + ".txt", ios::app);
  out << "---------------------" << endl;
  out << "initialS: " << endl;
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

#ifdef MASTER_VALVE_ML

void CVRP::findAccML(BBNODE *node) {
  if (If_in_Enu_State) throw runtime_error("findAccInitialML: not in enu state");
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
  useModelInPhase1_sep(node, num);
  useModelInPhase2(node, CONFIG::BranPhase2);
  auto initialS = Branch_pair;
  Branch_pair.clear();
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
  CGTesting(node, false, true, true, true);
  auto overallS = Branch_Pair_Val;
  self_mkdir("testAcc");
  ofstream out("testAcc/" + FileName + ".txt", ios::app);
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

void CVRP::generateModelInPhase1(BBNODE *node, bool if_sep) {
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
  getTrainingDataInPhase1(node);

  simulateWriteLPPseudocost(node);
  CGTesting(node, true, true, false, true);

  for (auto &edge : Branch_Pair_Val) {
	ml.EdgeTmpInfo[edge.first].SB_scores = edge.second;
  }

  if_sep ? ml.writeTrainingLPFile(false) : ml.writeTrainingLPFile_combinedFashion(false);
}

void CVRP::generateModelInPhase2(BBNODE *node, bool if_sep) {
#ifdef NoInitialScreening
  constructBranchingSet(node);
#else
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
#endif
#if defined(getOneStageModel) && getOneStageModel == 1
  int num = min(ml.giveInitialScreeningNumCandidates(), int(Branch_pair.size()));
#else
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
#endif
  if_sep ? useModelInPhase1_sep(node, num) : useModelInPhase1(node, num);
  getTrainingDataInPhase2(node);
  CGTesting(node, true, true, false, true);
  for (auto &edge : Branch_Pair_Val) {
	ml.EdgeTmpInfo[edge.first].SB_scores = edge.second;
  }
  ml.writeTrainingExactFile(false);
}

void CVRP::useMLInGeneralFramework(BBNODE *node, bool if_sep) {
#ifdef NoInitialScreening
  constructBranchingSet(node);
#else
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
#endif
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
  if_sep ? useModelInPhase1_sep(node, num) : useModelInPhase1(node, num);
  If_in_Enu_State ? useModelInPhase2(node, CONFIG::BranPhase2InEnu) : useModelInPhase2(node, CONFIG::BranPhase2);
  CGTesting(node, false, false, true);
}

void CVRP::useMLInStageOne(BBNODE *node) {
#ifdef NoInitialScreening
  constructBranchingSet(node);
#else
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
#endif
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
  useModelInPhase1_sep(node, num);
  LPTesting(node, CONFIG::BranPhase2, true);
  CGTesting(node, false, false, true);
}

void CVRP::simulateWriteLPPseudocost(BBNODE *node) {
  auto tmp_pair = Branch_pair;
  Branch_pair.clear();
  auto ratio_pseudo =
	  ((double)Branch_pair_from_pseudo.size()
		  / (double)(Branch_pair_from_pseudo.size() + Branch_pair_from_fractional.size()));
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(tmp_pair.size()));
  int pseudo_size = (int)(num * ratio_pseudo);
  int frac_size = num - pseudo_size;
  transform(Branch_pair_from_pseudo.begin(), Branch_pair_from_pseudo.begin() + pseudo_size,
			back_inserter(Branch_pair),
			[](pair<int, int> &p) { return make_pair(p.first, p.second); });
  transform(Branch_pair_from_fractional.begin(),
			Branch_pair_from_fractional.begin() + frac_size,
			back_inserter(Branch_pair),
			[](pair<int, int> &p) { return make_pair(p.first, p.second); });
  LPTesting(node, 0, false);
  Branch_pair = tmp_pair;
}

void CVRP::onlyUseStageOneModel(BBNODE *node, bool if_sep) {
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
  if_sep ? useModelInPhase1_sep(node, num) : useModelInPhase1(node, num);
  LPTesting(node, 0, false);//write record!
  CGTesting(node, true, true, false);
}

#ifdef useM_dynamically
void CVRP::useMLVaryLPTesting(BBNODE *node) {
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
  int num;
  giveDynamicM(node, num);
  useModelInPhase1_sep(node, num);
  if (node->Idx == 0) cp_Branch_pair = Branch_pair;
  If_in_Enu_State ? useModelInPhase2(node, CONFIG::BranPhase2InEnu) : useModelInPhase2(node, CONFIG::BranPhase2);
  CGTesting(node, false, false, true);
  if (node->Idx == 0) evaluateM1();
}

void CVRP::giveDynamicM(BBNODE *node, int &num) {
  if (node->Idx == 0) {
	num = min(ml.giveInitialScreeningNumCandidates(), int(Branch_pair.size()));
	k = num;
	cout << "we directly use the max n= " << num << endl;
	return;
  }
  if (if_use_full_k) {
	num = min(int(m), int(Branch_pair.size()));
	k = num;
	cout << "we directly use the full n= " << num << endl;
	return;
  }
  double t4oneLP = node->t4oneLP;
  double r_star = node->geo_r_star;
  double c = node->c;
  double n = ml.giveInitialScreeningNumCandidates();
  double G = UB - node->Val;

  double A = r_star * t4oneLP * pow(n + 1, 2);
  double B = t4oneLP * (log(2) * (F - G) * n * (m + 1) + 2 * r_star * (n - m) * (n + 1));
  double C = c * log(2) * (F - G) * n * (m + 1) + r_star * t4oneLP * pow(m - n, 2);

  double k1 = (-B + sqrt(pow(B, 2) - 4 * A * C)) / (2 * A);
  int k1_up = min(ceil(k1), double(m));
  int k1_down = min(floor(k1), double(m));

  cout << "___________________" << endl;
  cout << "G= " << G << " F= " << F << " c= " << c << " t4oneLP= " << t4oneLP << " r_star= " << r_star << " n= " << n
	   << " m= " << m << endl;
  cout << "k1= " << k1 << endl;

  unordered_set<double> tmp{1};
  if (k1_up > 1) tmp.emplace(k1_up);
  if (k1_down > 1) tmp.emplace(k1_down);

  cout << "---test---" << endl;
  cout << "k1_up= " << k1_up << endl;
  cout << "k1_down= " << k1_down << endl;

  if (tmp.size() == 1) {
	num = 1;
  } else {
	double least_T_N_ratio = numeric_limits<double>::max();
	for (auto i : tmp) {
	  double T_N_ratio = (c + i * t4oneLP) * pow(2, (G - F) / (r_star * ((m + 1) / n * i / (i + 1) + 1 - m / n)));
	  if (T_N_ratio == numeric_limits<double>::infinity()) {
		num = k1_up;
		break;
	  }
	  if (T_N_ratio < least_T_N_ratio) {
		least_T_N_ratio = T_N_ratio;
		num = i;
	  }
	  cout << "i= " << i << " T_N_ratio= " << T_N_ratio << endl;
	}
  }
  num = min(num, int(m));
  num = min(num, int(Branch_pair.size()));
  k = num;
  cout << "num= " << num << endl;
}

//void CVRP::evaluateM1() {
//  /**
//   * there is no pseudo candidates
//   */
//
//  int tmp_k;
//  for (int i = 0; i < cp_Branch_pair.size(); ++i) {
//    if (cp_Branch_pair[i] == Branch_pair[0]) {
//      m = i;
//      tmp_k = i;
//      break;
//    }
//  }
//  cp_Branch_pair.clear();
//
//  ++m;
//  m *= 2;
//  m = min(int(m), ml.giveInitialScreeningNumCandidates());
//  cout << "m= " << m << endl;
//
//
//  /**
//   * now test another calculation method
//   */
//
//  int n = tested_Branch_pair.size();
//  vector<int> rank(n);
//
//  for (int i = 0; i < n; ++i) {
//    rank[i] =
//        std::find(tested_Branch_pair.begin(), tested_Branch_pair.end(), cp_Branch_pair[i]) - tested_Branch_pair.begin()
//            + 1;
//  }
//
//  cout << "rank= ";
//  for (auto i : rank) {
//    cout << i << " ";
//  }
//  cout << endl;
//
//  //use the whole testing data!
//  for (int i = tmp_k + 1; i < n; ++i) {
//    //3rd min of i to n of rank vector
//    int min_range;
//    vector<int> tmp(rank.begin() + i, rank.end());
//    sort(tmp.begin(), tmp.end());
//    min_range = tmp[2];
//    cout << "min_range= " << min_range << " i= " << i << endl;
//    if (min_range > i) {
//      m = i;
//      cout << "m= " << m << endl;
//      break;
//    }
//  }
//
//  exit(0);
//
//}

void CVRP::evaluateM1() {
  int tmp_k = find(cp_Branch_pair.begin(), cp_Branch_pair.end(), Branch_pair[0]) - cp_Branch_pair.begin() + 1;
  vector<int> rank(tested_Branch_pair.size());
  for (int i = 0; i < tested_Branch_pair.size(); ++i) {
	rank[i] =
		std::find(tested_Branch_pair.begin(), tested_Branch_pair.end(), cp_Branch_pair[i]) - tested_Branch_pair.begin()
			+ 1;
  }
  int n = std::find(rank.begin(), rank.end(), 1) - rank.begin() + 1;
//  m = *max_element(rank.begin(), rank.begin() + n);
  m = n;
#ifdef SOLVER_VRPTW
  m = max(21., m);
#else
  m = max(20., m);
#endif
  m = min(m, double(ml.giveInitialScreeningNumCandidates()));

  cout << "m= " << m << endl;

  tested_Branch_pair.clear();
  cp_Branch_pair.clear();
}

#endif

#endif

void CVRP::exactCGTesting(BBNODE *node) {
#if defined(if_tell_model_works) || defined(if_tell_LPTesting_works)
  auto tmp_LPTestingImprovement_up = LPTestingImprovement_up;
  auto tmp_LPTestingImprovement_down = LPTestingImprovement_down;
  LPTestingImprovement_up.clear();
  LPTestingImprovement_down.clear();
#endif
  InitialScreening(node, true, true, CONFIG::ML_BranchPhase0, CONFIG::Frac4sudoCostBranPhase0);
#ifdef if_tell_LPTesting_works
  LPTestingImprovement_up = tmp_LPTestingImprovement_up;
  LPTestingImprovement_down = tmp_LPTestingImprovement_down;
  CGTesting(node, false, true, true);//get map
  int max_rank = 3;
  vector<pair<int, int>> top_n(max_rank);
  for (int i = 0; i < max_rank; ++i) {
	top_n[i] = Branch_Pair_Val[i].first;
  }
  Branch_pair.clear();
  //use the old pseudo cost
  auto tmp_RealImprovement_up = RealImprovement_up;
  auto tmp_RealImprovement_down = RealImprovement_down;
  RealImprovement_up.clear();
  RealImprovement_down.clear();
  InitialScreening(node, true, true, CONFIG::ML_NumBrCandiInPhase1, CONFIG::ML_Fracsudo);
  RealImprovement_up = tmp_RealImprovement_up;
  RealImprovement_down = tmp_RealImprovement_down;
  int num = min(ML3::giveTestingNumCandidates(node->TreeLevel, false), int(Branch_pair.size()));
  LPTesting(node, num, false);
  calculate_acc(top_n, Branch_pair, top_n_percentage_stage1);
  ml.EdgeTmpInfo.clear();
  Branch_pair = {Branch_Pair_Val[0].first};
#endif
#ifdef if_tell_model_works
  LPTestingImprovement_up = tmp_LPTestingImprovement_up;
  LPTestingImprovement_down = tmp_LPTestingImprovement_down;
  CGTesting(node, false, true, true);//get map
  int max_rank = 3;
  vector<pair<int, int>> top_n(max_rank);
  for (int i = 0; i < max_rank; ++i) {
	top_n[i] = Branch_Pair_Val[i].first;
  }
  Branch_pair.clear();
  //use the old pseudo cost
  auto tmp_RealImprovement_up = RealImprovement_up;
  auto tmp_RealImprovement_down = RealImprovement_down;
  RealImprovement_up.clear();
  RealImprovement_down.clear();
  InitialScreening(node, true, true, CONFIG::ML_NumBrCandiInPhase1, CONFIG::ML_Fracsudo);
  RealImprovement_up = tmp_RealImprovement_up;
  RealImprovement_down = tmp_RealImprovement_down;
  int num = min(ML3::giveTestingNumCandidates(node->TreeLevel, false), int(Branch_pair.size()));
  useModelInPhase1(node, num);//tell the top-n of all
  calculate_acc(top_n, Branch_pair, top_n_percentage_stage1);
  int take_out = 3;
  useModelInPhase2(node, take_out);//tell the top-n of all
  if (top_n_percentage_stage2_takeout_n.size() != take_out) {
	top_n_percentage_stage2_takeout_n.resize(take_out);
  }
  for (int i = 0; i < take_out; ++i) {
	calculate_acc(top_n, Branch_pair, top_n_percentage_stage2_takeout_n[take_out - i - 1]);
	Branch_pair.pop_back();
  }
  ml.EdgeTmpInfo.clear();
  Branch_pair = {Branch_Pair_Val[0].first};
#else
  CGTesting(node, false, false, true);
#endif
}

void CVRP::calculate_acc(const vector<pair<int, int>> &top_n, const vector<pair<int, int>> &branch_pair,
						 std::unordered_map<int, std::pair<int, double>> &top_n_per) {
  int max_rank = (int)top_n.size();
  int cnt = 0;
  for (int i = 1; i <= max_rank; ++i) {
	++top_n_per[i].second;
  }
  for (auto &edge : top_n) {
	auto if_find = std::find(branch_pair.begin(), branch_pair.end(), edge);
	if (if_find == branch_pair.end()) {
	  ++cnt;
	  continue;
	} else {
	  for (int i = cnt + 1; i <= max_rank; ++i) {
		++top_n_per[i].first;
	  }
	  break;
	}
  }
}

void CVRP::takeTheOneWithMaxColumns(BBNODE *node) {
  if (!If_in_Enu_State) throw runtime_error("takeTheOneWithMaxColumns: not in enu state");
  InitialScreening(node, true, true, numeric_limits<int>::max(), CONFIG::Frac4sudoCostBranPhase0);
  if (node->MatInEnu.size() == 0) {
	Branch_pair.resize(1);
	return;
  }
  auto &colmap = node->map_col_pool;
  size_t max_col = 0;
  auto dec = make_pair(0, 0);
  for (auto &edge : Branch_pair) {
	auto col_left = colmap[edge].size();
	size_t col_right = 0;
	int ai = edge.first, aj = edge.second;
	for (int i = 0; i < Dim; ++i) {
	  if (i == ai || i == aj) {
		continue;
	  }
	  auto pr = i > aj ? make_pair(aj, i) : make_pair(i, aj);
	  if (colmap.find(pr) != colmap.end()) {
		col_right += colmap[pr].size();
	  }
	  if (ai != 0) {
		pr = i > ai ? make_pair(ai, i) : make_pair(i, ai);
		if (colmap.find(pr) != colmap.end()) {
		  col_right += colmap[pr].size();
		}
	  }
	}
	if (max_col < col_left + col_right) {
	  max_col = col_left + col_right;
	  dec = edge;
	}
  }
  Branch_pair = {dec};
}


