//
// Created by You, Zhengzhong on 6/8/24.
//


#include "machine_learning.hpp"
#include "branching.hpp"

using namespace std;

std::unordered_map<std::pair<int, int>, std::pair<double, double>, PairHasher> MachineLearning::edge_lp_change{};

void MachineLearning::updateOneSideLPChange(const std::pair<int, int> &edge, double obj_change, bool if_left) {
  if_left ? (edge_lp_change[edge].first = obj_change) : (edge_lp_change[edge].second = obj_change);
}

void addEdgeInfo(
	vector<pair<string, double>> &edge_info,
	const std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> &history,
	const std::pair<int, int> &edge,
	const std::string &indicator_key,
	const std::string &change_key) {
  auto it = history.find(edge);
  if (it != history.end()) {
	const auto &[first, second] = it->second;  // structured binding for clarity
	edge_info.emplace_back(indicator_key, 1);
	edge_info.emplace_back(change_key, first / second);
  } else {
	edge_info.emplace_back(indicator_key, false);
	edge_info.emplace_back(change_key, 0);
  }
}

void MachineLearning::collectOneSideEdgeFeatures() {
  double current_gap = 1 - node->getCurrentNodeVal() / BaseBranching::ub;
  auto &lp_his_edge0 = BaseBranching::branching_history.lp_testing_improvement_down;
  auto &lp_his_edge1 = BaseBranching::branching_history.lp_testing_improvement_up;
  auto &heuristic_his_edge0 = BaseBranching::branching_history.heuristic_improvement_down;
  auto &heuristic_his_edge1 = BaseBranching::branching_history.heuristic_improvement_up;
  auto &exact_his_edge0 = BaseBranching::branching_history.exact_improvement_down;
  auto &exact_his_edge1 = BaseBranching::branching_history.exact_improvement_up;

  for (auto &edge : BaseBranching::current_branching_info.branch_pair) {
	auto &tmp_edge0_info = edge_tmp_info[edge].extra_features_edge0;
	auto &tmp_edge1_info = edge_tmp_info[edge].extra_features_edge1;
	tmp_edge0_info.emplace_back("current_gap", current_gap);
	tmp_edge1_info.emplace_back("current_gap", current_gap);

	addEdgeInfo(tmp_edge0_info, lp_his_edge0, edge, "indicator_lp_change", "historical_lp_obj_change");
	addEdgeInfo(tmp_edge1_info, lp_his_edge1, edge, "indicator_lp_change", "historical_lp_obj_change");

	addEdgeInfo(tmp_edge0_info,
				heuristic_his_edge0,
				edge,
				"indicator_heuristic_change",
				"historical_heuristic_obj_change");
	addEdgeInfo(tmp_edge1_info,
				heuristic_his_edge1,
				edge,
				"indicator_heuristic_change",
				"historical_heuristic_obj_change");

	addEdgeInfo(tmp_edge0_info, exact_his_edge0, edge, "indicator_exact_change", "historical_exact_obj_change");
	addEdgeInfo(tmp_edge1_info, exact_his_edge1, edge, "indicator_exact_change", "historical_exact_obj_change");
  }
}

