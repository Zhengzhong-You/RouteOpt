

#ifdef tmp
#include "learnArcMemory.hpp"
#include "CVRP.hpp"
#include "nlohmann/json.hpp"
using namespace std;
#endif



/**

### 静态特征 （当 图静态结构有削减的时候，会重新计算）

静态特征是指在整个分析期间不会改变的特征。它们通常是固有的属性，不会因为图的状态或外部环境而变化。例如：

#### 节点的静态特征：

2. **节点属性**：特定于节点的不变属性： 接近中心性， 特征向量中心性

#### 边的静态特征：

1. **边类型**：标准化长度
2. **边属性**：边两端节点共同的邻居数目

### 动态特征（当前 动态特征：

动态特征是指随着时间或网络状态的改变而变化的特征。这些特征可以反映网络的动态性质和节点之间的交互。

#### 节点的动态特征（您已提供的）：

1. **lp边数量**：节点的度，即与节点相连的边的数量。
2. **normalized lp边数量**：标准化的度，可能与网络中的最大度数或平均度数相关。
3. **normalized vertex cord**：节点的标准化坐标，可能与图的几何或空间结构有关。
4. **聚类系数**：衡量节点的邻居之间互相连接的程度。
5. **PageRank**：反映节点在整个网络中的影响力。

#### 边的动态特征（您已提供的）：

1. **normalized edge length**：边的标准化长度，可能与边的重要性或权重有关。
2. **edge value**：边的值或权重。
3. **边两端节点共同的邻居数目**：这反映了边连接的两个节点之间的紧密程度。

### 平均化的 动态特征：（上面的数据取平均）
 */



/**
 * input data：一个图，建立好这个图，然后 设定这个 图上面 和 点上面的 每个feature；
 * 最后label 是 ：每个点 是否是 关键点 以及 arc 是不是 关键arc；
 * 那我的c++的输出 training data 的文件是 json 文件；那我具体应该怎么输出呢？能否提供示范代码
 */

#ifdef tmp

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

void LeanArcMemory::predictUsingModel(std::vector<std::pair<int, int>> &pre_arcs) {
  namespace py = pybind11;

  pre_arcs.clear();
  auto beg = std::chrono::high_resolution_clock::now();
  py::scoped_interpreter guard{};
  auto end = std::chrono::high_resolution_clock::now();
  cout << "time for init python: " << std::chrono::duration<double>(end - beg).count() << endl;

  py::module_ sys = py::module_::import("sys");
  py::list sys_path = sys.attr("path");
  sys_path.append("../../../ML_python/LTA");

  beg = std::chrono::high_resolution_clock::now();
  py::module_ predict_module = py::module_::import("predict");
  end = std::chrono::high_resolution_clock::now();
  cout << "time for import predict: " << std::chrono::duration<double>(end - beg).count() << endl;

  auto data_path = OutFolder + "/" + cvrp->file_name + ".json";
  auto model_path = ModelPath + "/model.pt";

  beg = std::chrono::high_resolution_clock::now();
  py::object result = predict_module.attr("predict")(data_path, model_path);
  end = std::chrono::high_resolution_clock::now();
  cout << "time for predict: " << std::chrono::duration<double>(end - beg).count() << endl;
  vector<int> res;
  if (py::isinstance<py::list>(result)) {
	auto list = result.cast<py::list>();
	for (auto item : list) {
	  try {
		int value = item.cast<int>();
		res.emplace_back(value);
	  } catch (const py::cast_error &e) {
		throw runtime_error("Result is not a list of int.");
	  }
	}
  } else {
	throw runtime_error("Result is not a list.");
  }

  for (int i = 0; i < res.size(); ++i) {
	if (res[i] == 1) {
	  pre_arcs.emplace_back(considered_arcs[i]);
	}
  }

  if (res.size() != considered_arcs.size()) {
	throw runtime_error("res.size()!= considered_arcs.size()");
  }
}

void LeanArcMemory::getFeatureNames() {
  int dim = cvrp->dim;
  node_data.resize(dim);
  arc_data.resize(dim, vector<DataLeanArcMemory>(dim));
  node_feature_name.static_feature_name.emplace_back("normalized_demand");
  node_feature_name.current_features_name.emplace_back("node_in_cut_number");
  node_feature_name.current_features_name.emplace_back("current_node_in_mem_cut_number");
  for (int i = 0; i < dim; ++i) {
	node_data[i].static_feature.resize(1);
	node_data[i].current_features_normalized.resize(2);
  }
  arc_feature_name.static_feature_name.emplace_back("normalized_distance");
  arc_feature_name.current_features_name.emplace_back("current_arc_in_mem_cut_number");
  arc_feature_name.long_term_features_name.emplace_back("if_used");
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  arc_data[i][j].static_feature.resize(1);
	  arc_data[i][j].current_features_normalized.resize(1);
	  arc_data[i][j].long_term_features_normalized.resize(1, 0);
	}
  }
}

void LeanArcMemory::writeSupVectors(BbNode *node) {
  arc_value.clear();
  node_in_cut.clear();
  current_node_in_mem_cut.clear();
  current_arc_in_mem_cut.clear();

  auto &sol = node->all_lp_sol.first;//should use all
  auto &sol_value = node->all_lp_sol.second;
  int dim = cvrp->dim;
  vector<vector<double>> val_matrix(dim, vector<double>(dim, 0));
  for (int i = 0; i < sol.size(); ++i) {
	auto &seq = sol[i].col_seq;
	int b4 = 0;
	int break_point = sol[i].forward_concatenate_pos;
	for (int j = 0; j <= break_point; ++j) {
	  val_matrix[b4][seq[j]] += sol_value[i];
	  b4 = seq[j];
	}

	b4 = 0;
	for (int j = (int)seq.size() - 1; j >= break_point; --j) {
	  val_matrix[b4][seq[j]] += sol_value[i];
	  b4 = seq[j];
	}
  }

  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  if (val_matrix[i][j] > 0) {
		arc_value.emplace_back(make_pair(i, j), val_matrix[i][j]);
	  }
	}
  }

  auto &r1cs = node->r1cs;
  unordered_map<int, int> node_in_cut_map;
  unordered_map<int, int> node_in_mem_cut_map;
  unordered_map<pair<int, int>, int, PairHasher> arc_in_mem_cut_map;
  for (auto &r1c : r1cs) {
	for (auto &v : r1c.info_r1c.first) {
	  ++node_in_cut_map[v];
	}
	for (auto &v : r1c.arc_mem) {
	  ++node_in_mem_cut_map[v.second];
	  for (auto &vv : v.first) {
		++arc_in_mem_cut_map[make_pair(vv, v.second)];
	  }
	}
  }

  for (auto &i : node_in_cut_map) {
	node_in_cut.emplace_back(i.first, i.second);
  }
  for (auto &i : node_in_mem_cut_map) {
	current_node_in_mem_cut.emplace_back(i.first, i.second);
  }
  for (auto &i : arc_in_mem_cut_map) {
	current_arc_in_mem_cut.emplace_back(i.first, i.second);
  }

}

void LeanArcMemory::updateData() {

  int dim = cvrp->dim;

  for (int i = 0; i < dim; ++i) {
	std::fill(node_data[i].current_features_normalized.begin(), node_data[i].current_features_normalized.end(), 0);
  }
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  std::fill(arc_data[i][j].current_features_normalized.begin(), arc_data[i][j].current_features_normalized.end(),
				0);
	  arc_data[i][j].is_key = false;
	}
  }

  int node_cnt = 0, arc_cnt = 0;
  for (auto &i : arc_value) {
	arc_data[i.first.first][i.first.second].current_features_normalized[arc_cnt] = i.second;
  }
  for (auto &i : node_in_cut) {
	node_data[i.first].current_features_normalized[node_cnt] = i.second;
  }
  ++node_cnt;
  for (auto &i : current_node_in_mem_cut) {
	node_data[i.first].current_features_normalized[node_cnt] = i.second;
  }
  ++arc_cnt;
  for (auto &i : current_arc_in_mem_cut) {
	arc_data[i.first.first][i.first.second].is_key = true;
  }
}

void LeanArcMemory::collectLongTermDataEdgeOnceUsed(std::vector<SequenceInfo> &sol) {
  for (auto &i : sol) {
	int b4 = 0;
	for (int j = 0; j <= i.forward_concatenate_pos; ++j) {
	  ++arc_data[b4][i.col_seq[j]].long_term_features_normalized[0];
	  if (considered_arcs_set.find(make_pair(b4, i.col_seq[j])) == considered_arcs_set.end()) {
		considered_arcs.emplace_back(b4, i.col_seq[j]);
		considered_arcs_set.emplace(b4, i.col_seq[j]);
	  }
	  b4 = i.col_seq[j];
	}
	b4 = 0;
	for (int j = (int)i.col_seq.size() - 1; j > i.forward_concatenate_pos; --j) {
	  ++arc_data[b4][i.col_seq[j]].long_term_features_normalized[0];
	  if (considered_arcs_set.find(make_pair(b4, i.col_seq[j])) == considered_arcs_set.end()) {
		considered_arcs.emplace_back(b4, i.col_seq[j]);
		considered_arcs_set.emplace(b4, i.col_seq[j]);
	  }
	  b4 = i.col_seq[j];
	}
  }
}

void LeanArcMemory::write2Json() {
  using namespace nlohmann;
  self_mkdir(OutFolder);
  string file_name = OutFolder + "/" + cvrp->file_name + ".json";

  ofstream out(file_name);
  if (!out.is_open()) {
	throw runtime_error("cannot open " + file_name);
  }

  json info;

  int dim = cvrp->dim;
  for (int i = 0; i < dim; ++i) {
	info["node"].emplace_back(json{{"id", i},
								   {"is_key", node_data[i].is_key},
								   {"static_feature", node_data[i].static_feature},
								   {"current_features", node_data[i].current_features_normalized},
								   {"long_term_features", node_data[i].long_term_features_normalized}});
  }

  for (auto &i : considered_arcs) {
	info["arc"].emplace_back(json{{"source", i.first},
								  {"target", i.second},
								  {"is_key", arc_data[i.first][i.second].is_key},
								  {"static_feature", arc_data[i.first][i.second].static_feature},
								  {"current_features", arc_data[i.first][i.second].current_features_normalized},
								  {"long_term_features", arc_data[i.first][i.second].long_term_features_normalized}});
  }

  std::ofstream o(file_name);
  o.clear();
  o << info.dump(4) << std::endl;
  o.close();
}

void LeanArcMemory::readGraph(CVRP *ptr) {
  cvrp = ptr;
  int dim = cvrp->dim;
  auto &info_vertex = cvrp->info_vertex;
  auto &cost_mat4_vertex = cvrp->cost_mat4_vertex;
  normalized_distance.resize(dim, vector<double>(dim));
  normalized_demand.resize(dim);

  double max_dis = 0;
  for (int i = 0; i < dim; ++i) {
	max_dis = max(max_dis, *max_element(cost_mat4_vertex[i].begin(), cost_mat4_vertex[i].end()));
  }

  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  normalized_distance[i][j] = cost_mat4_vertex[i][j] / max_dis;
	}
  }

  double max_demand = *max_element(cvrp->demand, cvrp->demand + dim);
  for (int i = 0; i < dim; ++i) {
	normalized_demand[i] = cvrp->demand[i] / max_demand;
  }

  getFeatureNames();

  for (int i = 0; i < dim; ++i) {
	node_data[i].static_feature[0] = (normalized_demand[i]);
  }
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  arc_data[i][j].static_feature[0] = (normalized_distance[i][j]);
	}
  }
}

#endif


