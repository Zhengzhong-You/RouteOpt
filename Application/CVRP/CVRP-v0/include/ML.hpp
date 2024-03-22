
#ifndef VRPTW_INCLUDE_ML_HPP_
#define VRPTW_INCLUDE_ML_HPP_

#include "MACRO.hpp"
#include "BbNode.hpp"
#include <string>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <bitset>

#ifdef MASTER_VALVE_ML
#include <xgboost/c_api.h>

class CVRP;

class BbNode;

struct TmpEdgeRelatedData {
  std::vector<std::pair<std::string, double>> basic_features;
  std::vector<std::pair<std::string, double>>
	  resolving_lp_features;// 1 is the left branch, 3 is the right branch, and will be used!
  double sb_scores{};
  TmpEdgeRelatedData() = default;
};

struct LongEdgeRelatedData {
  std::pair<double, int> aver_edge_lp;
  LongEdgeRelatedData() = default;
};

class ML {
 public:
  CVRP *cvrp{};
  double max_edge_cost{};
  double max_mid_point_edge_cord_2_depot{};
  std::vector<std::vector<std::pair<double, double>>> mid_point_edge_cord;
  std::vector<std::vector<double>> mid_point_edge_cord_2_depot;
  std::vector<std::vector<std::vector<int>>> node_density_in_std_dis_vec_form;//
  std::vector<std::vector<double>> edge_2_other_convert_dis;
  int qid{};
  std::string lp_output_path{}, exact_output_path{};
  std::unordered_map<std::pair<int, int>, TmpEdgeRelatedData, PairHasher> edge_tmp_info;
  std::unordered_map<std::pair<int, int>, LongEdgeRelatedData, PairHasher> edge_long_info;

  BoosterHandle booster_1{};
  BoosterHandle booster_2{};
  void freeModel();
  ~ML() { freeModel(); }
  void loadModel(const std::string &model_path, int phase);
  void predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx);
  bool is_in_enu{};
  [[nodiscard]] int give_testing_num_candidates(int node_dep) const;
  [[nodiscard]] int give_initial_screening_num_candidates() const {
	return Config::ML_BranchPhase0;
  }
  void write_training_lp_file();
  void write_training_exact_file();
  void print_features();

  void read_enu_state(bool If_in_enu) {
	is_in_enu = If_in_enu;
  }
  void collect_edge_related_features(BbNode *node, double org_val);

  void collect_variable_related_features(
	  BbNode *node,
	  std::pair<int, int> edge,
	  const int *solver_ind,
	  int BeforeNumRow,
	  int numnz,
	  double org_val);

  void collect_resolving_features(
	  BbNode *node,
	  std::pair<int, int> edge,
	  int BeforeNumRow,
	  double tmp_val,
	  double org_val,
	  int numnz,
	  bool dir);

  void calculate_prerequisites();

  void get_info(CVRP *cvrp_ptr);

  std::vector<double> is_in_solution;
  std::unordered_map<std::pair<int, int>, double, PairHasher> edge_val;
  static void debug_input_data(const std::pair<std::string, double> &fs);
};

#endif
#endif //VRPTW_INCLUDE_ML_HPP_
