
#include "ML.hpp"
#include "CVRP.hpp"
using namespace std;

#define safe_xgboost(call) {  int Xerr = (call);\
if (Xerr != 0) { \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call + ":" + XGBGetLastError()); \
}\
}

#define PseudoMark std::string("pseudo-")

namespace fs = std::filesystem;

void self_mkdir(const std::string &path1) {
  if (fs::exists(path1) && fs::is_directory(path1)) {
	cout << path1 + " already exists" << endl;
  } else {
	fs::create_directory(path1);
	cout << path1 + " created" << endl;
  }
}

#ifdef MASTER_VALVE_ML

void ML::debug_input_data(const std::pair<std::string, double> &fs) {
  if (isnan(fs.second)) {
	cout << fs.first << " is nan" << endl;
  } else if (isinf(fs.second)) {
	cout << fs.first << " is inf" << endl;
  }
}

void ML::write_training_lp_file() {
  for (auto &edge : cvrp->branch_pair_val) {
	edge_tmp_info[edge.first].sb_scores = edge.second;
  }
  auto &path = lp_output_path;
  ofstream trainingData;
  trainingData.open(path, ios::app);
  vector<pair<int, int>> record;
  record.reserve(edge_tmp_info.size());
  for (auto &tmp_info : edge_tmp_info) {
	auto find = find_if(tmp_info.second.basic_features.begin(), tmp_info.second.basic_features.end(),
						[](const pair<string, double> &p) {
						  return p.first == PseudoMark + "ever_lp_find";
						});
	if (abs(find->second - 1) < TOLERANCE) {
	  trainingData << tmp_info.second.sb_scores;
	  trainingData << " qid:" << qid;
	  int cnt = 0;
	  for (auto &feature : tmp_info.second.basic_features) {
		trainingData << " " << cnt << ":" << (float)feature.second;
		++cnt;
	  }
	  trainingData << endl;
	} else {
	  record.emplace_back(tmp_info.first);
	}
  }

  ++qid;

  for (auto &pr : record) {
	auto &tmp_info = edge_tmp_info[pr];
	trainingData << tmp_info.sb_scores;
	trainingData << " qid:" << qid;
	int cnt = 0;
	for (auto &feature : tmp_info.basic_features) {
	  trainingData << " " << cnt << ":" << (float)feature.second;
	  ++cnt;
	}
	trainingData << endl;
  }

  ++qid;

  trainingData.close();

}

void ML::write_training_exact_file() {
  for (auto &edge : cvrp->branch_pair_val) {
	edge_tmp_info[edge.first].sb_scores = edge.second;
  }
  auto &path = exact_output_path;
  ofstream trainingData;
  trainingData.open(path, ios::app);
  for (auto &tmp_info : edge_tmp_info) {
	if (tmp_info.second.resolving_lp_features.empty())
	  continue;
	trainingData << tmp_info.second.sb_scores;
	cout << tmp_info.second.sb_scores << " | ";
	trainingData << " qid:" << qid;
	int cnt = 0;
	for (auto &feature : tmp_info.second.basic_features) {
	  trainingData << " " << cnt << ":" << (float)feature.second;
	  cout << feature.second << " ";
	  ++cnt;
	}
	for (auto &feature : tmp_info.second.resolving_lp_features) {
	  trainingData << " " << cnt << ":" << (float)feature.second;
	  cout << feature.second << " ";
	  ++cnt;
	}
	trainingData << endl;
	cout << endl;
  }

  ++qid;

  trainingData.close();

  print_features();
}

void ML::loadModel(const std::string &model_path, int phase) {

#ifdef SOLVER_VRPTW
  string m1 = "vrptw_model_1.bin";
  string m2 = "vrptw_model_2.bin";
#else
  string m1 = "cvrp_model_1.bin";
  string m2 = "cvrp_model_2.bin";
#endif

  if (phase == 1) {
	auto path1 = model_path + "/" + m1;
	safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_1))
	safe_xgboost(XGBoosterLoadModel(booster_1, path1.c_str()))
	cout << m1 + " loaded successfully" << endl;
  } else if (phase == 2) {
	auto path1 = model_path + "/" + m2;
	safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_2))
	safe_xgboost(XGBoosterLoadModel(booster_2, path1.c_str()))
	cout << m2 + " loaded successfully" << endl;
  } else {
	throw std::runtime_error("phase is not valid");
  }
}

void ML::predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx) {
  if (Branch_Val.empty())
	return;
  DMatrixHandle test;

  int numFeatures;
  BoosterHandle booster;
  bst_ulong output_length;
  const float *output_result;

  auto &edge = edge_tmp_info[Branch_Val[0].first];
  if (model_idx == 1) {
	numFeatures = (int)(edge.basic_features.size());
	booster = booster_1;
  } else if (model_idx == 2) {
	numFeatures = (int)(edge.basic_features.size() + edge.resolving_lp_features.size());
	booster = booster_2;
  } else {
	throw std::runtime_error("model_idx is not valid");
  }
  auto data = new float[Branch_Val.size() * numFeatures];
  if (model_idx == 1) {
	for (int i = 0; i < Branch_Val.size(); i++) {
	  auto &tmp_edge = edge_tmp_info[Branch_Val[i].first];
	  int j = 0;
	  for (auto &fs : tmp_edge.basic_features) {
#ifdef debugMLInputData
		debug_input_data(fs);
#endif
		data[i * numFeatures + j] = (float)fs.second;
		++j;
	  }
	}
  } else {
	for (int i = 0; i < Branch_Val.size(); i++) {
	  auto &tmp_edge = edge_tmp_info[Branch_Val[i].first];
	  int j = 0;
	  for (auto &fs : tmp_edge.basic_features) {
#ifdef debugMLInputData
		debug_input_data(fs);
#endif
		data[i * numFeatures + j] = (float)fs.second;
		++j;
	  }
	  for (auto &fs : tmp_edge.resolving_lp_features) {
#ifdef debugMLInputData
		debug_input_data(fs);
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

  safe_xgboost(XGDMatrixFree(test))
  delete[] data;
}

void ML::freeModel() {
  if (booster_1) {
	safe_xgboost(XGBoosterFree(booster_1))
	booster_1 = nullptr;
  }
  if (booster_2) {
	safe_xgboost(XGBoosterFree(booster_2))
	booster_2 = nullptr;
  }
}

int ML::give_testing_num_candidates(int node_dep) const {
  if (is_in_enu) return Config::BranPhase1InEnu;
  else return Config::BranPhase1;
}

void ML::collect_edge_related_features(BbNode *node, double org_val) {
  auto num_col = cvrp->num_col;
  auto dim = cvrp->dim;
  auto &branch_pair = cvrp->branch_pair;
  auto &cost_mat4_vertex = cvrp->cost_mat4_vertex;
  auto max_main_resource = cvrp->max_main_resource;
  auto &main_resource_across_arcs_in_forward_sense = cvrp->main_resource_across_arcs_in_forward_sense;
#ifdef SYMMETRY_PROHIBIT
  auto &MainResourceAcrossArcsInBackwardSense = cvrp->main_resource_across_arcs_in_backward_sense;
#endif
  if (!edge_tmp_info.empty()) throw std::runtime_error("edge_tmp_info is not empty");

  edge_val.clear();
  for (int i = 1; i <= node->num_edges; ++i) {
	edge_val[make_pair(node->edge_tail[i], node->edge_head[i])] = node->edge_value[i];
  }

  for (auto &edge : branch_pair) {
	auto &e = edge_tmp_info[edge];
	e.basic_features.emplace_back("tree_level", node->tree_level);
	e.basic_features.emplace_back("EdgeCost_ratio", cost_mat4_vertex[edge.first][edge.second] / max_edge_cost);
#ifdef SYMMETRY_PROHIBIT
	auto aver_res = (main_resource_across_arcs_in_forward_sense[edge.first][edge.second]
		+ MainResourceAcrossArcsInBackwardSense[edge.first][edge.second]) / 2;
	e.basic_features.emplace_back("EdgeRes_ratio", aver_res / max_main_resource);
#else
	e.basic_features.emplace_back("EdgeRes_ratio",
								  main_resource_across_arcs_in_forward_sense[edge.first][edge.second]
									  / max_main_resource);
#endif
	e.basic_features.emplace_back("EdgeDis2Depot_ratio",
								  mid_point_edge_cord_2_depot[edge.first][edge.second]
									  / max_mid_point_edge_cord_2_depot);

	e.basic_features.emplace_back("NodeDensity_in_std_dis_vec_form_ratio",
								  (double)pow(node_density_in_std_dis_vec_form[edge.first][edge.second].size(), 2) / dim
									  / dim);
	e.basic_features.emplace_back("edge_2_other_convert_dis", edge_2_other_convert_dis[edge.first][edge.second]);
  }
  is_in_solution.resize(num_col);
  fill(is_in_solution.begin(), is_in_solution.end(), 0);
  unordered_map<size_t, int> idx_in_solution;
  idx_in_solution.reserve(num_col);
  for (int i = 0; i < num_col; ++i) {
	idx_in_solution[node->index_columns[i]] = i;
  }
  for (auto &i : node->index_for_lp_solutions_in_column_pool) {
	is_in_solution[idx_in_solution[i.first]] = i.second;
  }
}

void ML::collect_variable_related_features(
	BbNode *node,
	pair<int, int> edge,
	const int *solver_ind,
	int BeforeNumRow,
	int numnz,
	double org_val) {

  auto &e = edge_tmp_info[edge];
  vector<double> frac_up;
  frac_up.reserve(numnz);
  for (int i = 0; i < numnz; ++i) {
	int col_idx = solver_ind[i];
	if (abs(is_in_solution[col_idx]) > TOLERANCE) {
	  frac_up.emplace_back(1 - is_in_solution[col_idx]);
	}
  }
  if (frac_up.empty()) throw runtime_error("frac_up is empty");
  e.basic_features.emplace_back("Mean_frac_up", accumulate(frac_up.begin(), frac_up.end(), 0.0) / (int)frac_up.size());
  e.basic_features.emplace_back("Min_frac_up", *min_element(frac_up.begin(), frac_up.end()));
  e.basic_features.emplace_back("Max_frac_up", *max_element(frac_up.begin(), frac_up.end()));
  e.basic_features.emplace_back("Frac_edge", edge_val[edge]);
  auto &real_improvement_up = cvrp->real_improvement_up;
  auto &real_improvement_down = cvrp->real_improvement_down;
  auto frac_edge_down = edge_val[edge];
  auto frac_edge_up = 1 - frac_edge_down;
  auto if_up = real_improvement_up.find(edge) == real_improvement_up.end() ? false : true;
  auto if_down = real_improvement_down.find(edge) == real_improvement_down.end() ? false : true;
  double improvement_up = if_up ? (real_improvement_up[edge].first / real_improvement_up[edge].second) : 0;
  double improvement_down = if_down ? (real_improvement_down[edge].first / real_improvement_down[edge].second) : 0;
  double pseudo_cost_up = max(improvement_up * frac_edge_up, 0.0);
  double pseudo_cost_down = max(improvement_down * frac_edge_down, 0.0);
  double pseudo_cost_mean = sqrt(pseudo_cost_up * pseudo_cost_down);
  e.basic_features.emplace_back("pseudo_cost_geomean_ratio", pseudo_cost_mean / org_val);
  e.basic_features.emplace_back("ever_geomean", if_up && if_down);

  auto &lp_testing_improvement_up = cvrp->lp_testing_improvement_up;
  auto &lp_testing_improvement_down = cvrp->lp_testing_improvement_down;
  bool if_find = lp_testing_improvement_down.find(edge) == lp_testing_improvement_down.end() ? false : true;
  double improvement_lp_up, improvement_lp_down;
  if (if_find) {
	improvement_lp_up = lp_testing_improvement_up[edge].first / lp_testing_improvement_up[edge].second;
	improvement_lp_down = lp_testing_improvement_down[edge].first / lp_testing_improvement_down[edge].second;
  } else {
	improvement_lp_up = improvement_lp_down = 0;
  }
  double pseudo_cost_lp_up = max(improvement_lp_up * frac_edge_up, 0.0);
  double pseudo_cost_lp_down = max(improvement_lp_down * frac_edge_down, 0.0);
  double pseudo_cost_lp_mean = sqrt(pseudo_cost_lp_up * pseudo_cost_lp_down);
  e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_up_ratio", pseudo_cost_lp_up / org_val);
  e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_down_ratio", pseudo_cost_lp_down / org_val);
  e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_geomean_ratio", pseudo_cost_lp_mean / org_val);
  e.basic_features.emplace_back(PseudoMark + "improvement_lp_up_ratio", improvement_lp_up / org_val);
  e.basic_features.emplace_back(PseudoMark + "improvement_lp_down_ratio", improvement_lp_down / org_val);
  e.basic_features.emplace_back(PseudoMark + "ever_lp_find", if_find);
  e.basic_features.emplace_back("branch_times", cvrp->branch_choice[edge]);
  auto &lp = edge_long_info[edge].aver_edge_lp;
  e.basic_features.emplace_back("aver_edge_lp", lp.first / lp.second);
}

void ML::collect_resolving_features(
	BbNode *node,
	std::pair<int, int> edge,
	int BeforeNumRow,
	double tmp_val,
	double org_val,
	int numnz,
	bool dir) {
  auto &e = edge_tmp_info[edge];
  double dual;
  auto dif = cvrp->calculateDifference(tmp_val, org_val);
  safe_solver(node->solver.getDual(BeforeNumRow, 1, &dual))
  auto
	  if_find =
	  find_if(e.resolving_lp_features.begin(), e.resolving_lp_features.end(), [&](const pair<string, double> &p) {
		return p.first == "dif / (tmp_val + org_val)";
	  });
  bool if_use = (if_find != e.resolving_lp_features.end());
  auto where = if_find - e.resolving_lp_features.begin();
  e.resolving_lp_features.emplace_back("dual/tmp_val", dual / tmp_val);
  e.resolving_lp_features.emplace_back("dif / (tmp_val + org_val)", dif / (tmp_val + org_val));
  if (if_use) {
	auto if_find2 =
		find_if(e.resolving_lp_features.begin(), e.resolving_lp_features.end(), [&](const pair<string, double> &p) {
		  return p.first == "product";
		});
	if (if_find2 == e.resolving_lp_features.end()) {
	  e.resolving_lp_features.emplace_back("product",
										   e.resolving_lp_features[where].second * dif / (tmp_val + org_val) * 1e6);
	}
  }
}

void ML::print_features() {
  for (auto &tmp_info : edge_tmp_info) {
	int cnt = 0;
	for (auto &feature : tmp_info.second.basic_features) {
	  cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
	  ++cnt;
	}
	for (auto &feature : tmp_info.second.resolving_lp_features) {
	  cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
	  ++cnt;
	}
	cout << BIG_PHASE_SEPARATION;
  }
  auto &edge = *edge_tmp_info.begin();
  int cnt = 0;
  vector<int> f_set;
  for (auto &feature : edge.second.basic_features) {
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

void ML::calculate_prerequisites() {
  auto dim = cvrp->dim;
  auto &info_vertex = cvrp->info_vertex;
  auto &cost_mat4_vertex = cvrp->cost_mat4_vertex;
  max_edge_cost = 0;
  int real_dim = dim - 1;
  for (int i = 0; i < real_dim; ++i) {
	double pre_cost = *max_element(cost_mat4_vertex[i].begin() + i + 1, cost_mat4_vertex[i].end());
	max_edge_cost = max(max_edge_cost, pre_cost);
  }
  mid_point_edge_cord.resize(dim, vector<pair<double, double >>(dim));
  for (int i = 0; i < dim; ++i) {
	for (int j = i; j < dim; ++j) {
	  mid_point_edge_cord[i][j].first = (info_vertex[i][1] + info_vertex[j][1]) / 2;
	  mid_point_edge_cord[i][j].second = (info_vertex[i][2] + info_vertex[j][2]) / 2;
	  mid_point_edge_cord[j][i] = mid_point_edge_cord[i][j];
	}
  }
  mid_point_edge_cord_2_depot.resize(dim, vector<double>(dim));
  for (int i = 0; i < dim; ++i) {
	for (int j = i; j < dim; ++j) {
	  mid_point_edge_cord_2_depot[i][j] =
		  sqrt_self(float(
			  (mid_point_edge_cord[i][j].first - info_vertex[0][1])
				  * (mid_point_edge_cord[i][j].first - info_vertex[0][1]) +
				  (mid_point_edge_cord[i][j].second - info_vertex[0][2])
					  * (mid_point_edge_cord[i][j].second - info_vertex[0][2])));
	  mid_point_edge_cord_2_depot[j][i] = mid_point_edge_cord_2_depot[i][j];
	}
  }
  max_mid_point_edge_cord_2_depot = 0;
  for (int i = 0; i < real_dim; ++i) {
	double max_dis = *max_element(mid_point_edge_cord_2_depot[i].begin() + i + 1, mid_point_edge_cord_2_depot[i].end());
	max_mid_point_edge_cord_2_depot = max(max_mid_point_edge_cord_2_depot, max_dis);
  }
  double geo_dis = exp(accumulate(cost_mat4_vertex[0].begin() + 1,
								  cost_mat4_vertex[0].end(),
								  0.0,
								  [](double a, double b) { return a + log(b); }) / (dim - 1));
  vector<yzzLong> v_neighbor(dim);
  for (int i = 0; i < dim; ++i) {//need 0 included
	for (int j = 0; j < dim; ++j) {
	  if (cost_mat4_vertex[i][j] <= geo_dis) {
		v_neighbor[i].set(j);
	  }
	}
  }
  vector<vector<yzzLong>> density_std_dis(dim, vector<yzzLong>(dim));
  for (int i = 0; i < dim; ++i) {
	for (int j = i + 1; j < dim; ++j) {
	  density_std_dis[i][j] = v_neighbor[i] & v_neighbor[j];
	  density_std_dis[j][i] = density_std_dis[i][j];
	}
  }
  node_density_in_std_dis_vec_form.resize(dim, vector<vector<int >>(dim));
  for (int i = 0; i < dim; ++i) {
	for (int j = i + 1; j < dim; ++j) {
	  for (int k = 0; k < dim; ++k) {
		if (density_std_dis[i][j].test(k)) {
		  node_density_in_std_dis_vec_form[i][j].emplace_back(k);
		}
	  }
	  node_density_in_std_dis_vec_form[j][i] = node_density_in_std_dis_vec_form[i][j];
	}
  }
  edge_2_other_convert_dis.resize(dim, vector<double>(dim));
  for (int i = 0; i < dim; ++i) {
	for (int j = i + 1; j < dim; ++j) {
	  double aver_dis = 0;
	  for (int k = 0; k < node_density_in_std_dis_vec_form[i][j].size(); ++k) {
		for (int l = k + 1; l < node_density_in_std_dis_vec_form[i][j].size(); ++l) {
		  double dif_x =
			  mid_point_edge_cord[node_density_in_std_dis_vec_form[i][j][k]][node_density_in_std_dis_vec_form[i][j][l]].first
				  - mid_point_edge_cord[i][j].first;
		  double dif_y =
			  mid_point_edge_cord[node_density_in_std_dis_vec_form[i][j][k]][node_density_in_std_dis_vec_form[i][j][l]].second
				  - mid_point_edge_cord[i][j].second;
		  aver_dis += sqrt_self(float(dif_x * dif_x + dif_y * dif_y));
		}
	  }
	  if (node_density_in_std_dis_vec_form[i][j].empty())
		aver_dis = 0;
	  else
		aver_dis /=
			(double)pow(node_density_in_std_dis_vec_form[i][j].size(), 2) / 2;
	  edge_2_other_convert_dis[i][j] = aver_dis / geo_dis;
	  edge_2_other_convert_dis[j][i] = edge_2_other_convert_dis[i][j];
	}
  }
}

void ML::get_info(CVRP *cvrp_ptr) {
  cvrp = cvrp_ptr;
  switch (ML_STATE) {
	case 1: self_mkdir("train_lp");
	  lp_output_path = "train_lp/" + cvrp->file_name + ".txt";
	  break;
	case 2:loadModel("model", 1);
	  self_mkdir("train_exact");
	  exact_output_path = "train_exact/" + cvrp->file_name + ".txt";
	  break;
	case 3:
	case 4:loadModel("model", 1);
	  loadModel("model", 2);
	  break;
	default:throw runtime_error("ML_state error");
  }
}
#endif