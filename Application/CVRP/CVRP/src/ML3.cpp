//
// Created by Zhengzhong You on 3/27/23.
//

#include "ML3.hpp"
#include "CVRP.hpp"
//#include "fastforest.h"
using namespace std;

#define safe_xgboost(call) {  int Xerr = (call);\
if (Xerr != 0) { \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call + ":" + XGBGetLastError()); \
}\
}

#define PseudoMark std::string("pseudo-")

namespace fs = std::experimental::filesystem;

void self_mkdir(const std::string &path) {
  // Check if directory exists
  if (fs::exists(path) && fs::is_directory(path)) {
	cout << path + " already exists" << endl;
  } else {
	// Directory does not exist, create it
	fs::create_directory(path);
	cout << path + " created" << endl;
  }
}

#ifdef MASTER_VALVE_ML

void ML3::debugInputData(const std::pair<std::string, double> &fs) {
  if (isnan(fs.second)) {
	cout << fs.first << " is nan" << endl;
  } else if (isinf(fs.second)) {
	cout << fs.first << " is inf" << endl;
  }
}

#ifdef testAcc
void ML3::writeTrainingLPFile(bool if_debug_features) {
  auto &path1 = if_in_enu ? enum_lp_path : lp_Output_path;
  auto &path2 = if_in_enu ? enum_lp_path2 : lp_Output_path2;
  ofstream trainingData1, trainingData2;
  trainingData1.open(path1, ios::app);
  trainingData2.open(path2, ios::app);
  vector<pair<int, int>> record;
  record.reserve(EdgeTmpInfo.size());
  for (auto &tmp_info : EdgeTmpInfo) {
	//check when tmp_info.second.BasicFeatures.first is PseudoMark + "ever_lp_find" and second is 1
	auto find = find_if(tmp_info.second.BasicFeatures.begin(), tmp_info.second.BasicFeatures.end(),
						[](const pair<string, double> &p) {
						  return p.first == PseudoMark + "ever_lp_find";
						});
	if (abs(find->second - 1) < TOLERANCE) {
	  trainingData1 << tmp_info.second.SB_scores;
	  trainingData1 << " qid:" << QID;
	  int cnt = 0;
	  for (auto &feature : tmp_info.second.BasicFeatures) {
		trainingData1 << " " << cnt << ":" << (float) feature.second;
		++cnt;
	  }
	  trainingData1 << endl;
	} else {
	  trainingData2 << tmp_info.second.SB_scores;
	  trainingData2 << " qid:" << QID;
	  int cnt = 0;
	  for (auto &feature : tmp_info.second.BasicFeatures) {
		trainingData2 << " " << cnt << ":" << (float) feature.second;
		++cnt;
	  }
	  trainingData2 << endl;
	}
  }

  ++QID;
  trainingData1.close();
  trainingData2.close();
  if (if_debug_features) {
	printFeatures();
  }
}
#else
void ML3::writeTrainingLPFile(bool if_debug_features) {
  auto &path = if_in_enu ? enum_lp_path : lp_Output_path;
  ofstream trainingData;
  trainingData.open(path, ios::app);
  vector<pair<int, int>> record;
  record.reserve(EdgeTmpInfo.size());
  for (auto &tmp_info : EdgeTmpInfo) {
	//check when tmp_info.second.BasicFeatures.first is PseudoMark + "ever_lp_find" and second is 1
	auto find = find_if(tmp_info.second.BasicFeatures.begin(), tmp_info.second.BasicFeatures.end(),
						[](const pair<string, double> &p) {
						  return p.first == PseudoMark + "ever_lp_find";
						});
	if (abs(find->second - 1) < TOLERANCE) {
	  trainingData << tmp_info.second.SB_scores;
	  trainingData << " qid:" << QID;
	  int cnt = 0;
	  for (auto &feature : tmp_info.second.BasicFeatures) {
		trainingData << " " << cnt << ":" << (float)feature.second;
		++cnt;
	  }
	  trainingData << endl;
	} else {
	  record.emplace_back(tmp_info.first);
	}
  }

  ++QID;

  for (auto &pr : record) {
	auto &tmp_info = EdgeTmpInfo[pr];
	trainingData << tmp_info.SB_scores;
	trainingData << " qid:" << QID;
	int cnt = 0;
	for (auto &feature : tmp_info.BasicFeatures) {
	  trainingData << " " << cnt << ":" << (float)feature.second;
	  ++cnt;
	}
	trainingData << endl;
  }

  ++QID;

  trainingData.close();
  if (if_debug_features) {
	printFeatures();
  }
}
#endif

void ML3::writeTrainingExactFile(bool if_debug_features) {
  auto &path = if_in_enu ? enum_exact_path : exact_Output_path;
  ofstream trainingData;
  trainingData.open(path, ios::app);
  for (auto &tmp_info : EdgeTmpInfo) {
	if (tmp_info.second.ResolvingLPFeatures.empty())
	  continue;
	trainingData << tmp_info.second.SB_scores;
	cout << tmp_info.second.SB_scores << " | ";
	trainingData << " qid:" << QID;
	int cnt = 0;
	for (auto &feature : tmp_info.second.BasicFeatures) {
	  trainingData << " " << cnt << ":" << (float)feature.second;
	  cout << feature.second << " ";
	  ++cnt;
	}
	for (auto &feature : tmp_info.second.ResolvingLPFeatures) {
	  trainingData << " " << cnt << ":" << (float)feature.second;
	  cout << feature.second << " ";
	  ++cnt;
	}
	trainingData << endl;
	cout << endl;
  }

  ++QID;

  trainingData.close();
  if (if_debug_features) {
	printFeatures();
  }
}

void ML3::loadModel(const std::string &model_path, int phase) {
  //load model

#ifdef SOLVER_VRPTW
  string m1 = "vrptw_model_1.bin";
  string m2 = "vrptw_model_2.bin";
#ifdef NoDifEnu_b4_af
  string m3 = m1;
  string m4 = m2;
#else
  string m3 = "vrptw_model_enu_1.bin";
  string m4 = "vrptw_model_enu_2.bin";
#endif
#else
  string m1 = "cvrp_model_1.bin";
  string m2 = "cvrp_model_2.bin";
#ifdef NoDifEnu_b4_af
  string m3 = m1;
  string m4 = m2;
#else
  string m3 = "cvrp_model_enu_1.bin";
  string m4 = "cvrp_model_enu_2.bin";
#endif
#endif

  if (phase == 1) {
	auto path1 = model_path + "/" + m1;
	safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_1))
	safe_xgboost(XGBoosterLoadModel(booster_1, path1.c_str()))
	cout << m1 + " loaded successfully" << endl;
	auto path2 = model_path + "/" + m3;
	safe_xgboost(XGBoosterCreate(nullptr, 0, &enum_booster1))
	safe_xgboost(XGBoosterLoadModel(enum_booster1, path2.c_str()))
	cout << m3 + " loaded successfully" << endl;
  } else if (phase == 2) {
	auto path1 = model_path + "/" + m2;
	safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_2))
	safe_xgboost(XGBoosterLoadModel(booster_2, path1.c_str()))
	cout << m2 + " loaded successfully" << endl;
	auto path2 = model_path + "/" + m4;
	safe_xgboost(XGBoosterCreate(nullptr, 0, &enum_booster2))
	safe_xgboost(XGBoosterLoadModel(enum_booster2, path2.c_str()))
	cout << m4 + " loaded successfully" << endl;
  } else {
	throw std::runtime_error("phase is not valid");
  }
}

//void ML3::predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx) {
//  if (Branch_Val.empty())
//    return;
//  DMatrixHandle test;
//
//  int numFeatures;
//  BoosterHandle booster;
//  bst_ulong output_length;
//  const float *output_result;
//
//  auto &edge = EdgeTmpInfo[Branch_Val[0].first];
//  if (model_idx == 1) {
//    numFeatures = (int) (edge.BasicFeatures.size());
//    booster = if_in_enu ? enum_booster1 : booster_1;
//  } else if (model_idx == 2) {
//    numFeatures = (int) (edge.BasicFeatures.size() + edge.ResolvingLPFeatures.size());
//    booster = if_in_enu ? enum_booster2 : booster_2;
//  } else {
//    throw std::runtime_error("model_idx is not valid");
//  }
//  auto data = new float[Branch_Val.size() * numFeatures];
//  if (model_idx == 1) {
//    for (int i = 0; i < Branch_Val.size(); i++) {
//      auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
//      int j = 0;
//      for (auto &fs : tmp_edge.BasicFeatures) {
//#ifdef debugMLInputData
//        debugInputData(fs);
//#endif
//        data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//    }
//  } else {
//    for (int i = 0; i < Branch_Val.size(); i++) {
//      auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
//      int j = 0;
//      for (auto &fs : tmp_edge.BasicFeatures) {
//#ifdef debugMLInputData
//        debugInputData(fs);
//#endif
//        data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//      for (auto &fs : tmp_edge.ResolvingLPFeatures) {
//#ifdef debugMLInputData
//        debugInputData(fs);
//#endif
//        data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//    }
//  }
//#ifdef debugMLInputData
//  cout << "check in MLInputData" << endl;
//  size_t tmp_len = (int) Branch_Val.size() * numFeatures;
//  for (int i = 0; i < tmp_len; i++) {
//    if (isnan(data[i]))
//      throw std::runtime_error("nan in data");
//    else if (isinf(data[i]))
//      throw std::runtime_error("inf in data");
//  }
//#endif
//
//  safe_xgboost(XGDMatrixCreateFromMat(data,
//                                      (int) Branch_Val.size(),
//                                      numFeatures,
//                                      numeric_limits<float>::quiet_NaN(),
//                                      &test))
//
//  safe_xgboost(XGBoosterPredict(booster, test, 1, 0, 0, &output_length, &output_result))
//
//  for (unsigned int i = 0; i < output_length; i++) {
//    Branch_Val[i].second = output_result[i];
//  }
//  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
//    return a.second > b.second;
//  });
//
//  safe_xgboost(XGDMatrixFree(test))
//  delete[] data;
//}

void ML3::predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx) {
  if (Branch_Val.empty())
	return;
  DMatrixHandle test;

  int numFeatures;
  BoosterHandle booster;
  bst_ulong output_length;
  const float *output_result;

  auto &edge = EdgeTmpInfo[Branch_Val[0].first];
  if (model_idx == 1) {
	numFeatures = (int)(edge.BasicFeatures.size());
	booster = if_in_enu ? enum_booster1 : booster_1;
  } else if (model_idx == 2) {
	numFeatures = (int)(edge.BasicFeatures.size() + edge.ResolvingLPFeatures.size());
	booster = if_in_enu ? enum_booster2 : booster_2;
  } else {
	throw std::runtime_error("model_idx is not valid");
  }
  auto data = new float[Branch_Val.size() * numFeatures];
  if (model_idx == 1) {
	for (int i = 0; i < Branch_Val.size(); i++) {
	  auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
	  int j = 0;
	  for (auto &fs : tmp_edge.BasicFeatures) {
#ifdef debugMLInputData
		debugInputData(fs);
#endif
		data[i * numFeatures + j] = (float)fs.second;
//		cout << "name= :" << fs.first << " value= " << fs.second << " j= " << j << endl;
		++j;
	  }
	}
  } else {
	for (int i = 0; i < Branch_Val.size(); i++) {
	  auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
	  int j = 0;
	  for (auto &fs : tmp_edge.BasicFeatures) {
#ifdef debugMLInputData
		debugInputData(fs);
#endif
		data[i * numFeatures + j] = (float)fs.second;
		++j;
	  }
	  for (auto &fs : tmp_edge.ResolvingLPFeatures) {
#ifdef debugMLInputData
		debugInputData(fs);
#endif
		data[i * numFeatures + j] = (float)fs.second;
		++j;
	  }
	}
  }
#ifdef debugMLInputData
  cout << "check in MLInputData" << endl;
  size_t tmp_len = (int) Branch_Val.size() * numFeatures;
  for (int i = 0; i < tmp_len; i++) {
	if (isnan(data[i]))
	  throw std::runtime_error("nan in data");
	else if (isinf(data[i]))
	  throw std::runtime_error("inf in data");
  }
#endif

  safe_xgboost(XGDMatrixCreateFromMat(data,
									  (int)Branch_Val.size(),
									  numFeatures,
									  numeric_limits<float>::quiet_NaN(),
									  &test))

  safe_xgboost(XGBoosterPredict(booster, test, 1, 0, 0, &output_length, &output_result))

  for (unsigned int i = 0; i < output_length; i++) {
	Branch_Val[i].second = output_result[i];
  }
  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
	return a.second > b.second;
  });

//  for (int i = 0; i < Branch_Val.size(); i++) {
//    cout << "(" << Branch_Val[i].first.first << "," << Branch_Val[i].first.second << "):" << Branch_Val[i].second
//         << endl;
//  }
//
//  cout << MID_PHASE_SEPARATION;
//
//  vector<string> features(numFeatures);
//  for (int i = 0; i < numFeatures; i++) {
//    features[i] = "f" + to_string(i);
//  }

//  string m1 = "cvrp_model_1.txt";
//  string m2 = "cvrp_model_2.txt";
////  string m3 = "cvrp_model_enu_1.bin";
////  string m4 = "cvrp_model_enu_2.bin";
//
//  string m;
//
//  if (model_idx == 1) {
//    m = m1;
//  } else if (model_idx == 2) {
//    m = m2;
//  }

//  const auto fastforest = fastforest::load_txt("model/" + m, features);
//
//  for (int i = 0; i < Branch_Val.size(); i++) {
//    auto beg = chrono::high_resolution_clock::now();
//    const auto predictions = fastforest(data + i * numFeatures);
//    auto end = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::microseconds>(end - beg);
//    cout << "fastforest time: " << duration.count() << endl;
//    Branch_Val[i].second = predictions;
//  }
//
//  sort(Branch_Val.begin(), Branch_Val.end(), [](const auto &a, const auto &b) {
//    return a.second > b.second;
//  });
//
//  for (int i = 0; i < Branch_Val.size(); i++) {
//    cout << "(" << Branch_Val[i].first.first << "," << Branch_Val[i].first.second << "):" << Branch_Val[i].second
//         << endl;
//  }

  safe_xgboost(XGDMatrixFree(test))
  delete[] data;
}

void ML3::freeModel() {
  if (booster_1) {
	safe_xgboost(XGBoosterFree(booster_1))
	booster_1 = nullptr;
  }
  if (booster_2) {
	safe_xgboost(XGBoosterFree(booster_2))
	booster_2 = nullptr;
  }
  if (enum_booster1) {
	safe_xgboost(XGBoosterFree(enum_booster1))
	enum_booster1 = nullptr;
  }
  if (enum_booster2) {
	safe_xgboost(XGBoosterFree(enum_booster2))
	enum_booster2 = nullptr;
  }
}

int ML3::giveTestingNumCandidates(int node_dep) const {
//  std::vector<std::pair<int, int>> *p;
//  if (if_in_enu) p = &CONFIG::TreeLevel_Num_Vec_enu;
//  else p = &CONFIG::TreeLevel_Num_Vec;
//  for (auto &it : *p) {
//    if (node_dep <= it.first)
//      return it.second;
//  }
  if (if_in_enu) return CONFIG::BranPhase1InEnu;
  else return CONFIG::BranPhase1;
}

void ML3::collectEdgeRelatedFeatures(CVRP *cvrp, BBNODE *node, double org_val) {
  auto NumCol = cvrp->NumCol;
  auto Dim = cvrp->Dim;
  auto &Branch_pair = cvrp->Branch_pair;
  auto &CostMat4Vertex = cvrp->CostMat4Vertex;
  auto MaxMainResource = cvrp->MaxMainResource;
  auto &MainResourceAcrossArcsInForwardSense = cvrp->MainResourceAcrossArcsInForwardSense;
#ifdef SYMMETRY_PROHIBIT
  auto &MainResourceAcrossArcsInBackwardSense = cvrp->MainResourceAcrossArcsInBackwardSense;
#endif
  if (!EdgeTmpInfo.empty()) throw std::runtime_error("EdgeTmpInfo is not empty");

  edge_val.clear();
  for (int i = 1; i <= node->NumEdges; ++i) {
	edge_val[make_pair(node->EdgeTail[i], node->EdgeHead[i])] = node->EdgeVal[i];
  }

  //depth of the node
  for (auto &edge : Branch_pair) {
	auto &e = EdgeTmpInfo[edge];
	e.BasicFeatures.emplace_back("TreeLevel", node->TreeLevel);
	e.BasicFeatures.emplace_back("EdgeCost_ratio", CostMat4Vertex[edge.first][edge.second] / MaxEdgeCost);
#ifdef SYMMETRY_PROHIBIT
	auto aver_res = (MainResourceAcrossArcsInForwardSense[edge.first][edge.second]
		+ MainResourceAcrossArcsInBackwardSense[edge.first][edge.second]) / 2;
	e.BasicFeatures.emplace_back("EdgeRes_ratio", aver_res / MaxMainResource);
#else
	e.BasicFeatures.emplace_back("EdgeRes_ratio",
								 MainResourceAcrossArcsInForwardSense[edge.first][edge.second] / MaxMainResource);
#endif
	e.BasicFeatures.emplace_back("EdgeDis2Depot_ratio",
								 MidPointEdgeCord_2_depot[edge.first][edge.second] / MaxMidPointEdgeCord_2_depot);

	e.BasicFeatures.emplace_back("NodeDensity_in_std_dis_vec_form_ratio",
								 (double)pow(NodeDensity_in_std_dis_vec_form[edge.first][edge.second].size(), 2) / Dim
									 / Dim);
	e.BasicFeatures.emplace_back("Edge_2_other_convert_dis", Edge_2_other_convert_dis[edge.first][edge.second]);
  }
  //find the index of fractional variable
  if_in_solution.resize(NumCol);
  fill(if_in_solution.begin(), if_in_solution.end(), 0);
  unordered_map<size_t, int> idx_in_solution;
  idx_in_solution.reserve(NumCol);
  for (int i = 0; i < NumCol; ++i) {
	idx_in_solution[node->IdxCols[i]] = i;
  }
  for (auto &i : node->Idx4LPSolsInColPool) {
	if_in_solution[idx_in_solution[i.first]] = i.second;
  }
}

void ML3::collectVariableRelatedFeatures(CVRP *cvrp,
										 BBNODE *node,
										 pair<int, int> edge,
										 const int *solver_ind,
										 int BeforeNumRow,
										 int numnz,
										 double org_val) {

  auto &e = EdgeTmpInfo[edge];
  //if 0, no related to this edge, if 1, only related to i, if 2, only related to j, if 3, related to both i and j
  vector<double> frac_up;
  frac_up.reserve(numnz);
  for (int i = 0; i < numnz; ++i) {
	int col_idx = solver_ind[i];
	if (abs(if_in_solution[col_idx]) > TOLERANCE) {
	  frac_up.emplace_back(1 - if_in_solution[col_idx]);
	}
  }
  if (frac_up.empty()) throw runtime_error("frac_up is empty");
  //mean row_density
  e.BasicFeatures.emplace_back("Mean_frac_up", accumulate(frac_up.begin(), frac_up.end(), 0.0) / (int)frac_up.size());
  e.BasicFeatures.emplace_back("Min_frac_up", *min_element(frac_up.begin(), frac_up.end()));
  e.BasicFeatures.emplace_back("Max_frac_up", *max_element(frac_up.begin(), frac_up.end()));
  e.BasicFeatures.emplace_back("Frac_edge", edge_val[edge]);
  auto &RealImprovement_up = cvrp->RealImprovement_up;
  auto &RealImprovement_down = cvrp->RealImprovement_down;
  auto frac_edge_down = edge_val[edge];
  auto frac_edge_up = 1 - frac_edge_down;
  auto if_up = RealImprovement_up.find(edge) == RealImprovement_up.end() ? false : true;
  auto if_down = RealImprovement_down.find(edge) == RealImprovement_down.end() ? false : true;
  double improvement_up = if_up ? (RealImprovement_up[edge].first / RealImprovement_up[edge].second) : 0;
  double improvement_down = if_down ? (RealImprovement_down[edge].first / RealImprovement_down[edge].second) : 0;
  double pseudo_cost_up = max(improvement_up * frac_edge_up, 0.0);
  double pseudo_cost_down = max(improvement_down * frac_edge_down, 0.0);
  double pseudo_cost_mean = sqrt(pseudo_cost_up * pseudo_cost_down);
  e.BasicFeatures.emplace_back("pseudo_cost_geomean_ratio", pseudo_cost_mean / org_val);
  e.BasicFeatures.emplace_back("ever_geomean", if_up && if_down);

  auto &LPTestingImprovement_up = cvrp->LPTestingImprovement_up;
  auto &LPTestingImprovement_down = cvrp->LPTestingImprovement_down;
  bool if_find = LPTestingImprovement_down.find(edge) == LPTestingImprovement_down.end() ? false : true;
  double improvement_lp_up, improvement_lp_down;
  if (if_find) {
	improvement_lp_up = LPTestingImprovement_up[edge].first / LPTestingImprovement_up[edge].second;
	improvement_lp_down = LPTestingImprovement_down[edge].first / LPTestingImprovement_down[edge].second;
  } else {
	improvement_lp_up = improvement_lp_down = 0;
  }
  double pseudo_cost_lp_up = max(improvement_lp_up * frac_edge_up, 0.0);
  double pseudo_cost_lp_down = max(improvement_lp_down * frac_edge_down, 0.0);
  double pseudo_cost_lp_mean = sqrt(pseudo_cost_lp_up * pseudo_cost_lp_down);
  e.BasicFeatures.emplace_back(PseudoMark + "pseudo_cost_lp_up_ratio", pseudo_cost_lp_up / org_val);
  e.BasicFeatures.emplace_back(PseudoMark + "pseudo_cost_lp_down_ratio", pseudo_cost_lp_down / org_val);
  e.BasicFeatures.emplace_back(PseudoMark + "pseudo_cost_lp_geomean_ratio", pseudo_cost_lp_mean / org_val);
  e.BasicFeatures.emplace_back(PseudoMark + "improvement_lp_up_ratio", improvement_lp_up / org_val);
  e.BasicFeatures.emplace_back(PseudoMark + "improvement_lp_down_ratio", improvement_lp_down / org_val);
  e.BasicFeatures.emplace_back(PseudoMark + "ever_lp_find", if_find);
  e.BasicFeatures.emplace_back("branch_times", cvrp->BranchChoice[edge]);
  auto &lp = EdgeLongInfo[edge].AverEdgeLP;
  e.BasicFeatures.emplace_back("aver_edge_lp", lp.first / lp.second);
}

void ML3::collectResolvingFeatures(CVRP *cvrp,
								   BBNODE *node,
								   std::pair<int, int> edge,
								   int BeforeNumRow,
								   double tmp_val,
								   double org_val,
								   int numnz,
								   bool dir) {
  auto &e = EdgeTmpInfo[edge];
  double dual;
  auto dif = cvrp->calculateDif(tmp_val, org_val);
  safe_solver(node->solver.SOLVERgetDual(BeforeNumRow, 1, &dual))
  auto
	  if_find = find_if(e.ResolvingLPFeatures.begin(), e.ResolvingLPFeatures.end(), [&](const pair<string, double> &p) {
	return p.first == "dif / (tmp_val + org_val)";
  });
  bool if_use = (if_find != e.ResolvingLPFeatures.end());
  auto where = if_find - e.ResolvingLPFeatures.begin();
  e.ResolvingLPFeatures.emplace_back("dual/tmp_val", dual / tmp_val);
  e.ResolvingLPFeatures.emplace_back("dif / (tmp_val + org_val)", dif / (tmp_val + org_val));
  // if in the ResolvingLPFeatures, we find two dif / (tmp_val + org_val), then we use to create a new feature
  if (if_use) {
	auto if_find2 =
		find_if(e.ResolvingLPFeatures.begin(), e.ResolvingLPFeatures.end(), [&](const pair<string, double> &p) {
		  return p.first == "product";
		});
	if (if_find2 == e.ResolvingLPFeatures.end()) {
	  e.ResolvingLPFeatures.emplace_back("product",
										 e.ResolvingLPFeatures[where].second * dif / (tmp_val + org_val) * 1e6);
	}
  }
}

void ML3::printFeatures() {
  for (auto &tmp_info : EdgeTmpInfo) {
	int cnt = 0;
	for (auto &feature : tmp_info.second.BasicFeatures) {
	  cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
	  ++cnt;
	}
	for (auto &feature : tmp_info.second.ResolvingLPFeatures) {
	  cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
	  ++cnt;
	}
	cout << BIG_PHASE_SEPARATION;
  }
  auto &edge = *EdgeTmpInfo.begin();
  int cnt = 0;
  vector<int> f_set;
  for (auto &feature : edge.second.BasicFeatures) {
	if (feature.first.find(PseudoMark) != string::npos) {
	  cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
	  f_set.emplace_back(cnt);
	}
	++cnt;
  }
  cout << "f_set: " << endl;
  cout << "[";
  for (auto &tmp : f_set) {
	cout << "," << tmp;
  }
  cout << "]";
}

void ML3::calculatePrerequisites(CVRP *cvrp) {
  auto Dim = cvrp->Dim;
  auto &InfoVertex = cvrp->InfoVertex;
  auto &CostMat4Vertex = cvrp->CostMat4Vertex;
  MaxEdgeCost = 0;
  int RealDim = Dim - 1;
  for (int i = 0; i < RealDim; ++i) {
	double pre_cost = *max_element(CostMat4Vertex[i].begin() + i + 1, CostMat4Vertex[i].end());
	MaxEdgeCost = max(MaxEdgeCost, pre_cost);
  }
  MidPointEdgeCord.resize(Dim, vector<pair<double, double >>(Dim));
  for (int i = 0; i < Dim; ++i) {
	for (int j = i; j < Dim; ++j) {
	  MidPointEdgeCord[i][j].first = (InfoVertex[i][1] + InfoVertex[j][1]) / 2;
	  MidPointEdgeCord[i][j].second = (InfoVertex[i][2] + InfoVertex[j][2]) / 2;
	  MidPointEdgeCord[j][i] = MidPointEdgeCord[i][j];
	}
  }
  MidPointEdgeCord_2_depot.resize(Dim, vector<double>(Dim));
  for (int i = 0; i < Dim; ++i) {
	for (int j = i; j < Dim; ++j) {
	  MidPointEdgeCord_2_depot[i][j] =
		  sqrt_self(float(
			  (MidPointEdgeCord[i][j].first - InfoVertex[0][1]) * (MidPointEdgeCord[i][j].first - InfoVertex[0][1]) +
				  (MidPointEdgeCord[i][j].second - InfoVertex[0][2])
					  * (MidPointEdgeCord[i][j].second - InfoVertex[0][2])));
	  MidPointEdgeCord_2_depot[j][i] = MidPointEdgeCord_2_depot[i][j];
	}
  }
  MaxMidPointEdgeCord_2_depot = 0;
  for (int i = 0; i < RealDim; ++i) {
	double max_dis = *max_element(MidPointEdgeCord_2_depot[i].begin() + i + 1, MidPointEdgeCord_2_depot[i].end());
	MaxMidPointEdgeCord_2_depot = max(MaxMidPointEdgeCord_2_depot, max_dis);
  }
  //calculate the geoMean distance of the graph
  double geo_dis = exp(accumulate(CostMat4Vertex[0].begin() + 1,
								  CostMat4Vertex[0].end(),
								  0.0,
								  [](double a, double b) { return a + log(b); }) / (Dim - 1));
  vector<yzzLong> v_neighbor(Dim);
  for (int i = 0; i < Dim; ++i) {//need 0 included
	for (int j = 0; j < Dim; ++j) {
	  if (CostMat4Vertex[i][j] <= geo_dis) {
		v_neighbor[i].set(j);
	  }
	}
  }
  vector<vector<yzzLong>> density_std_dis(Dim, vector<yzzLong>(Dim));
  for (int i = 0; i < Dim; ++i) {
	for (int j = i + 1; j < Dim; ++j) {
	  density_std_dis[i][j] = v_neighbor[i] & v_neighbor[j];
	  density_std_dis[j][i] = density_std_dis[i][j];
	}
  }
  NodeDensity_in_std_dis_vec_form.resize(Dim, vector<vector<int >>(Dim));
  for (int i = 0; i < Dim; ++i) {
	for (int j = i + 1; j < Dim; ++j) {
	  for (int k = 0; k < Dim; ++k) {
		if (density_std_dis[i][j].test(k)) {
		  NodeDensity_in_std_dis_vec_form[i][j].emplace_back(k);
		}
	  }
	  NodeDensity_in_std_dis_vec_form[j][i] = NodeDensity_in_std_dis_vec_form[i][j];
	}
  }
  Edge_2_other_convert_dis.resize(Dim, vector<double>(Dim));
  for (int i = 0; i < Dim; ++i) {
	for (int j = i + 1; j < Dim; ++j) {
	  double aver_dis = 0;
	  for (int k = 0; k < NodeDensity_in_std_dis_vec_form[i][j].size(); ++k) {
		for (int l = k + 1; l < NodeDensity_in_std_dis_vec_form[i][j].size(); ++l) {
		  double dif_x =
			  MidPointEdgeCord[NodeDensity_in_std_dis_vec_form[i][j][k]][NodeDensity_in_std_dis_vec_form[i][j][l]].first
				  - MidPointEdgeCord[i][j].first;
		  double dif_y =
			  MidPointEdgeCord[NodeDensity_in_std_dis_vec_form[i][j][k]][NodeDensity_in_std_dis_vec_form[i][j][l]].second
				  - MidPointEdgeCord[i][j].second;
		  aver_dis += sqrt_self(float(dif_x * dif_x + dif_y * dif_y));
		}
	  }
	  if (NodeDensity_in_std_dis_vec_form[i][j].empty())
		aver_dis = 0;
	  else
		aver_dis /=
			(double)pow(NodeDensity_in_std_dis_vec_form[i][j].size(), 2) / 2;
	  Edge_2_other_convert_dis[i][j] = aver_dis / geo_dis;
	  Edge_2_other_convert_dis[j][i] = Edge_2_other_convert_dis[i][j];
	}
  }
}

void ML3::getInfo(CVRP *cvrp) {
#ifdef ML_state
  switch (ML_state) {
	case 1:self_mkdir("train_enu_lp");
	  self_mkdir("train_lp");
#ifdef testAcc
	  self_mkdir("train_enu_lp2");
	  self_mkdir("train_lp2");
	  lp_Output_path2 = "train_lp2/" + cvrp->FileName + ".txt";
	  enum_lp_path2 = "train_enu_lp2/" + cvrp->FileName + ".txt";
#endif
	  enum_lp_path = "train_enu_lp/" + cvrp->FileName + ".txt";
	  lp_Output_path = "train_lp/" + cvrp->FileName + ".txt";
	  break;
	case 2:loadModel("model", 1);
	  self_mkdir("train_enu_exact");
	  self_mkdir("train_exact");
#ifdef testAcc
	  enum_exact_path2 = "train_enu_exact2/" + cvrp->FileName + ".txt";
	  exact_Output_path2 = "train_exact2/" + cvrp->FileName + ".txt";
#endif
	  enum_exact_path = "train_enu_exact/" + cvrp->FileName + ".txt";
	  exact_Output_path = "train_exact/" + cvrp->FileName + ".txt";
	  break;
	case 3:
	case 5:loadModel("model", 1);
	  loadModel("model", 2);
	  break;
	case 4:loadModel("model", 1);
	  break;
	default:throw runtime_error("ML_state error");
  }
#endif
}

void ML3::collectResolvingFeatures_runCG(CVRP *cvrp,
										 BBNODE *node,
										 pair<int, int> edge,
										 double min_rc,
										 double mean_rc,
										 int add,
										 int BeforeNumRow,
										 double tmp_val,
										 double org_val) {
  auto &e = EdgeTmpInfo[edge];
  double dual;
  auto dif = cvrp->calculateDif(tmp_val, org_val);
  safe_solver(node->solver.SOLVERgetDual(BeforeNumRow, 1, &dual))
  auto
	  if_find = find_if(e.ResolvingLPFeatures.begin(), e.ResolvingLPFeatures.end(), [&](const pair<string, double> &p) {
	return p.first == "dif_CG";
  });
  bool if_use = if_find != e.ResolvingLPFeatures.end();
  auto where = if_find - e.ResolvingLPFeatures.begin();
  e.ResolvingLPFeatures.emplace_back("dual_CG", dual);
  e.ResolvingLPFeatures.emplace_back("dif_CG", dif);
  e.ResolvingLPFeatures.emplace_back("min_rc_CG", min_rc);
  e.ResolvingLPFeatures.emplace_back("mean_rc_CG", mean_rc);
  e.ResolvingLPFeatures.emplace_back("add_CG", add);
  if (if_use) {
	auto if_find2 =
		find_if(e.ResolvingLPFeatures.begin(), e.ResolvingLPFeatures.end(), [&](const pair<string, double> &p) {
		  return p.first == "product_CG";
		});
	if (if_find2 == e.ResolvingLPFeatures.end()) {
	  e.ResolvingLPFeatures.emplace_back("product_CG", e.ResolvingLPFeatures[where].second * dif);
	}
  }
}

void ML3::collectExtraEdgeFeatures(CVRP *cvrp, BBNODE *node) {
#ifdef NoDifEnu_b4_af
  return;
#endif
  auto &map_col_pool = node->map_col_pool;
  auto &Branch_pair = cvrp->Branch_pair;
  auto &Map_Edge_ColIdx_in_Enu = cvrp->Map_Edge_ColIdx_in_Enu;//without record the first one!
  auto Dim = cvrp->Dim;
  auto &Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;

  //for each edge, collect the rc of the cols in the edge
  constexpr int inner_size = 3;
  vector<vector<double>> val_(Branch_pair.size(), vector<double>(inner_size));
  int Br_cnt = 0;
  for (auto &edge : Branch_pair) {
	int RCs_cnt = 0;
	int ai = edge.first, aj = edge.second;
	RCs_cnt += (int)Map_Edge_ColIdx_in_Enu[edge].size();
	for (auto &col_idx : map_col_pool[edge]) {
	  if (Deleted_ColsInEnuPool[col_idx]) continue;
	  ++RCs_cnt;
	}
	val_[Br_cnt][0] = RCs_cnt;

	RCs_cnt = 0;
	if (ai) {
	  for (int i = 0; i < Dim; ++i) {
		if (i == ai || i == aj)continue;
		auto pr = i < ai ? make_pair(i, ai) : make_pair(ai, i);
		RCs_cnt += (int)Map_Edge_ColIdx_in_Enu[pr].size();
		for (auto &col_idx : map_col_pool[pr]) {
		  if (Deleted_ColsInEnuPool[col_idx]) continue;
		  ++RCs_cnt;
		}
	  }
	}

	val_[Br_cnt][1] = RCs_cnt;

	RCs_cnt = 0;

	for (int i = 0; i < Dim; ++i) {
	  if (i == ai || i == aj)continue;
	  auto pr = i < aj ? make_pair(i, aj) : make_pair(aj, i);
	  RCs_cnt += (int)Map_Edge_ColIdx_in_Enu[pr].size();
	  for (auto &col_idx : map_col_pool[pr]) {
		if (Deleted_ColsInEnuPool[col_idx]) continue;
		++RCs_cnt;
	  }
	}
	val_[Br_cnt][2] = RCs_cnt;
	++Br_cnt;
  }

  //normalize val
  for (int i = 0; i < inner_size; ++i) {
	//print val_
	double max_ = (*max_element(val_.begin(), val_.end(),
								[i](vector<double> &a, vector<double> &b) {
								  return a[i] < b[i];
								}))[i];
	double min_ = (*min_element(val_.begin(), val_.end(),
								[i](vector<double> &a, vector<double> &b) {
								  return a[i] < b[i];
								}))[i];
	if (max_ == min_) {
	  if (min_ != 0) {
		transform(val_.begin(), val_.end(), val_.begin(),
				  [i](vector<double> &a) {
					a[i] = 1;
					return a;
				  });
	  }
	  continue;
	}
	transform(val_.begin(), val_.end(), val_.begin(),
			  [i, max_, min_](vector<double> &a) {
				a[i] = (a[i] - min_) / (max_ - min_);
				return a;
			  });
  }

  Br_cnt = 0;
  for (auto &edge : Branch_pair) {
	auto &e = EdgeTmpInfo[edge];
	for (int i = 0; i < inner_size; ++i) {
	  e.BasicFeatures.emplace_back("enu_" + to_string(i), val_[Br_cnt][i]);
	}
	++Br_cnt;
  }
}

#endif