//
// Created by Zhengzhong You on 3/27/23.
//

#include <utility>

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

#define AccuracyTest_outputDir std::string ("AccuracyTest")

void CVRP::LPTesting(BBNODE *node, int num, bool if_writeBranch_pair, bool if_record_sb_scores,
					 bool if_record_LP_improvement) {
  if (Branch_pair.size() == 1) {
	cout << "LPTesting: Branch_pair.size() == 1, return!" << endl;
	return;
  }
  auto begin = high_resolution_clock::now();
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  double org_val = node->Val;
  int BeforeNumRow = NumRow;
  int cnt = 0;

  fill(solver_ind2, solver_ind2 + NumCol, BeforeNumRow);
  vector<double> solver_val3(NumCol);
  vector<int> solver_ind1(NumCol);
  iota(solver_ind1.begin(), solver_ind1.end(), 0);
  bool if_changed = false;
  branch_LP.clear();
  for (auto &edge : Branch_pair) {
	int numnz;
	double tmp_val;
	getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
	if (!if_changed) {
	  safe_solver(addBranchconstr(numnz,
								  solver_ind,
								  solver_val,
								  SOLVER_LESS_EQUAL,
								  0,
								  nullptr,
								  node->solver))
	  if_changed = true;
	} else {
	  chgBranchconstr(solver_val3.data(),
					  solver_ind2,
					  solver_ind1.data(),
					  numnz,
					  solver_ind,
					  solver_val,
					  SOLVER_LESS_EQUAL,
					  0,
					  node->solver);
	}
	safe_solver(node->solver.SOLVERreoptimize())
	safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
	auto dif1 = calculateDif(tmp_val, org_val);
	safe_solver(inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
	safe_solver(node->solver.SOLVERreoptimize())
	safe_solver(node->solver.SOLVERgetObjVal(&tmp_val))
	auto dif2 = calculateDif(tmp_val, org_val);
	auto product = dif1 * dif2;
	Branch_Val[cnt++] = {edge, product};
	branch_LP[edge] = product;
	cout << "( " << edge.first << " , " << edge.second << " ): " << product << endl;
//    cout << product << " ";
	if (if_record_LP_improvement) {
	  LPTestingImprovement_down[edge].first += dif1;
	  ++LPTestingImprovement_down[edge].second;
	  LPTestingImprovement_up[edge].first += dif2;
	  ++LPTestingImprovement_up[edge].second;
	}
  }
  safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
  safe_solver(node->solver.SOLVERupdatemodel())

  if (if_writeBranch_pair) {
	int size = min(num, int(Branch_Val.size()));
	sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
	  return a.second > b.second;
	});
	Branch_pair.resize(size);
	transform(Branch_Val.begin(),
			  Branch_Val.begin() + size,
			  Branch_pair.begin(),
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
	Branch_Pair_Val = Branch_Val;
  }

  std::unordered_map<std::pair<int, int>, double, PairHasher> branch_LP_tmp;
  branch_LP_tmp.reserve(Branch_pair.size());
  for (auto &edge : Branch_pair) {
	branch_LP_tmp[edge] = branch_LP[edge];
  }
  branch_LP = branch_LP_tmp;

  safe_solver(node->solver.SOLVERreoptimize())
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - begin).count();
  cout << "LPTesting Branch_pair size= " << Branch_pair.size() << " took " << eps << " s!" << endl;
}

#ifdef MASTER_VALVE_ML

void CVRP::useModelInPhase1(BBNODE *node, int num, int if_force_frac) {
  if (Branch_pair.size() == 1) {
	cout << "useModelInPhase1: Branch_pair.size() == 1, return!" << endl;
	return;
  }
  getTrainingDataInPhase1(node);
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val(Branch_pair.size());
  transform(Branch_pair.begin(),
			Branch_pair.end(),
			Branch_Val.begin(),
			[](const auto &edge) {
			  return make_pair(edge, 0.0);
			});
  ml.predict(Branch_Val, 1);

  Branch_pair.resize(num);

  if (if_force_frac) {
	cout << "force_frac" << endl;
	auto ratio_pseudo =
		((double) Branch_pair_from_pseudo.size()
			/ (double) (Branch_pair_from_pseudo.size() + Branch_pair_from_fractional.size()));
	unordered_set < pair<int, int>, PairHasher > p_set;
	p_set.reserve(Branch_pair_from_pseudo.size());
	for (auto &edge : Branch_pair_from_pseudo) p_set.emplace(edge);
	int pseudo_size = (int) (num * ratio_pseudo);
	int frac_size = num - pseudo_size;
	int cnt_p = 0, cnt_f = 0, cnt = 0;
	for (auto &edge : Branch_Val) {
	  if (p_set.find(edge.first) != p_set.end()) {
		if (cnt_p < pseudo_size) {
		  Branch_pair[cnt++] = edge.first;
		  if (cnt == num) break;
		  ++cnt_p;
		}
	  } else {
		if (cnt_f < frac_size) {
		  Branch_pair[cnt++] = edge.first;
		  if (cnt == num) break;
		  ++cnt_f;
		}
	  }
	}
  } else {
	transform(Branch_Val.begin(),
			  Branch_Val.begin() + num,
			  Branch_pair.begin(),
			  [](const auto &a) {
				return a.first;
			  });
  }
  cout << "useModelInPhase1 Branch_pair size= " << Branch_pair.size() << endl;
}

void CVRP::getTrainingDataInPhase1(BBNODE *node) {
  if (Branch_pair.size() == 1) {
	cout << "getTrainingDataInPhase1: Branch_pair.size() == 1, return!" << endl;
	return;
  }
  auto beg = high_resolution_clock::now();
  double org_val = node->Val;
  int BeforeNumRow = NumRow;
  ml.collectEdgeRelatedFeatures(this, node, org_val);
  for (auto &edge : Branch_pair) {
	int numnz;
	getNewCstrCoeffByEdge(node, edge, solver_ind, solver_val, numnz);
	ml.collectVariableRelatedFeatures(this, node, edge, solver_ind, BeforeNumRow, numnz, org_val);
  }
  if (If_in_Enu_State) {
	ml.collectExtraEdgeFeatures(this, node);
  }
  auto end = high_resolution_clock::now();
  auto eps = duration<double>(end - beg).count();
  cout << "getTrainingDataInPhase1 took " << eps << " s!" << endl;
}

#endif

void CVRP::writeSB_scores2Files(int node_dep) {
  // if AccuracyTest_outputDir does not exist, create it
  if (!std::experimental::filesystem::exists(AccuracyTest_outputDir)) {
	std::experimental::filesystem::create_directory(AccuracyTest_outputDir);
  }
  //change Branch_Pair_Val into map
  unordered_map<pair<int, int>, double, PairHasher> Map_Edge_SB_scores;
  for (auto &edge : Branch_Pair_Val) {
	Map_Edge_SB_scores[edge.first] = edge.second;
  }
  ofstream fout;
  string out_path = AccuracyTest_outputDir + "/" + FileName + ".txt";
  // 2 models, check candidates selected by model 1 and model 2, print their rank and ratio, more detail, better.
  // later information can be checked by python
  fout.open(out_path, ios::app);
  fout << "node_dep= " << node_dep << endl;
  fout << "model_1= " << endl;
  for (auto &edge : Branch_selected_from_model1) {
	fout << Map_Edge_SB_scores[edge] << " ";
  }
  fout << endl;
  fout << "model_2= " << endl;
  for (auto &edge : Branch_selected_from_model2) {
	fout << Map_Edge_SB_scores[edge] << " ";
  }
  fout << endl;
  fout.close();
}