
#ifndef CVRP_CVRP_HPP
#define CVRP_CVRP_HPP

#include "BbNode.hpp"
#include "MACRO.hpp"
#include "Config.hpp"
#include "Solver.hpp"
#include "capsep.h"
#include "ML.hpp"
#include "InstanceData.hpp"
#include <cstdio>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <tuple>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <random>
#include <algorithm>
#include <climits>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <list>
#include <Eigen/Sparse>
#include <limits>
#include <cstdio>
#include <filesystem>

#ifdef HGS_APPLIED
#include "my_hgs.hpp"
#endif

#ifdef DELUXING_APPLIED
#include "deLuxing.hpp"
#endif

using yzzLong = std::bitset<MAX_NUM_CUSTOMERS>;
using R1CMem = std::bitset<MAX_NUM_R1CS>;

struct Label {
  bool is_extended{};
  int end_vertex{};
  double sum_main_resource{};
  double rc{};
  double cost{};
  yzzLong pi{};
  Label *p_label{};
  Label *c_label{};
  R1CMem rank1_cut_mem{};
  int num_valid_rank1_cut{};
  int *valid_rank1_cut{};
  int *rank1_cut_mem_multi{};
  int *valid_rank1_cut_multi{};
  int num_valid_rank1_cut_multi{};
};

struct PtrAllR1CS {
  double *r1c_to_pi{};
  int num_r1c_nzero{};
  double *r1c_multi_to_pi{};
  int num_r1c_multi_nzero{};

  PtrAllR1CS(BbNode *node, CVRP *cvrp);
  ~PtrAllR1CS();
};

struct Rank1MultiLabel {
  std::vector<int> c;
  std::vector<int> w_no_c;
  int plan_idx{};
  double vio{};
  char search_dir{};
  Rank1MultiLabel(std::vector<int> c, std::vector<int> w_no_c, int plan_idx, double vio, char search_dir) :
	  c(std::move(c)),
	  w_no_c(std::move(w_no_c)),
	  plan_idx(plan_idx),
	  vio(vio),
	  search_dir(search_dir) {}
  Rank1MultiLabel() = default;
};

class CmpNodeValue {
 public:
  auto operator()(BbNode *a, BbNode *b) -> bool {
	return a->value > b->value;
  }
};

struct VectorHash {
  size_t operator()(const std::vector<int> &V) const {
	size_t hash = 0;
	for (auto &i : V) {
	  hash *= MAX_NUM_CUSTOMERS;
	  hash += i;
	}
	return hash;
  }
};

struct Rank1MultiPairHasher {
  size_t operator()(const std::pair<std::vector<int>, int> &V) const {
	size_t hash = 0;
	for (auto &i : V.first) {
	  hash *= MAX_NUM_CUSTOMERS;
	  hash += i;
	}
	hash *= MAX_NUM_CUSTOMERS;
	hash += V.second;
	return hash;
  }
};

using VecLabel = std::pair<std::vector<Label *>, int>;
using QueueTree = std::priority_queue<BbNode *, std::vector<BbNode *>, CmpNodeValue>;
using StackTree = std::stack<BbNode *>;

struct Bucket {
  std::vector<int> bucket_arcs{};
  std::vector<std::pair<double, int>> jump_arcs{};
  int i{};
};

struct MyDoubleHash {
  size_t operator()(double d) const {
	return std::hash<int>()(std::round(d * 100));
  }
};

struct MyDoubleEqual {
  bool operator()(double a, double b) const {
	return std::abs(a - b) < TOLERANCE;
  }
};

class CVRP {

 public:

  std::unordered_map<std::vector<int>,
					 std::vector<std::vector<std::vector<int>>>, VectorHash> pureMap;

  std::unordered_map<std::vector<int>, std::vector<std::vector<int>>, VectorHash> rank1_multi_mem_plan_map;
  int num_valid_r1c_multi_in_cg{};

  Solver solver{};

  VecLabel **label_array_in_forward_sense{};

  VecLabel **if_exist_extra_labels_in_forward_sense{};

  std::chrono::time_point<std::chrono::high_resolution_clock> glo_beg{},
	  glo_end{};
  double glo_eps{};

  int num_br{};
  int num_explored_nodes{};
  char if_int_n_feasible{};
  double max_vio{};
  double *demand{};

  int rank1_mem_size_limit{};

  int branch_times{};
  int real_dim{};
  int num_col{}, num_row{};

  void updateLowerBound(double val);
  void resetEnvironment(BbNode *node);
  double lp_val{}, lb{};
  int max_num_route4_mip{};
  bool if_mip_enumeration_suc{};
  double lb_transformed{}, ub{};
  double obj4_first_col{};
  int transformed_number{};
  int idx_node{};
  double *pi{}, *slack{};
  bool if_exact_labeling_cg{};
  std::vector<double> pi4_labeling{};//this pi is used for labeling and cutting, need special care!
  std::vector<double> slack4_labeling{};//this slack is used for labeling and cutting, need special care!
  int *col_pool4_pricing{};
  int *copy_col_pool4_pricing{};//used in heuristic find ub!
  size_t pool_beg4_pricing{};
  size_t pool_beg4_copy_pricing{};
  size_t prior_pool_beg4_pricing{};
  int *col_pool4_mem{};
  size_t pool_beg4_mem{};
  std::vector<std::vector<int>> ip_opt_sol{};

  std::vector<std::pair<int, int>> branch_pair{};
  std::unordered_map<std::pair<int, int>, double, PairHasher> branch_lp{};
  std::vector<std::pair<int, int>> branch_pair_from_pseudo{};
  std::vector<std::pair<int, int>> branch_pair_from_fractional{};
  double smallest_rc{};
  double gap_between_last_smallest_rc_and_rc_threshold{};
  double rc_std{};
  double cap{};
  double max_main_resource{}, meet_point_resource_in_bi_dir{}, meet_point_resource_in_bi_dir_enu{};
  int dim{}, num_vehicle{};
#ifdef SOLVER_VRPTW
  bool if_force_keep_rcc{};
#endif

  std::vector<std::vector<double>> info_vertex{};
  std::vector<std::vector<double>> cost_mat4_vertex{};
  std::vector<yzzLong> ng_mem4_vertex{};
  std::vector<yzzLong> rank1_sep_heur_mem4_vertex{};
  std::vector<int> size_ng_mem4_vertex{};
  std::vector<std::vector<double>> chg_cost_mat4_vertex{};
  std::vector<int> prior_mem_sets{};
  std::vector<double> lb4_vertex{};
  std::vector<double> ub4_vertex{};

  std::vector<std::tuple<Label *, Label *, double>> negative_rc_label_tuple{};
  std::unordered_map<yzzLong, std::tuple<Label *, Label *, double>> map4_each_negative_rc_route{};

  int idx_glo{};
  int size_ng_mem{};
  std::string file_name;

  Label *all_label{};
  int *all_valid_r1cs{};
  int *all_valid_r1c_multi{};
  int *all_r1c_multi_mem{};

  double **arc_graph{};
  double **arc_graph_revised{};
  int num_edge{};
  double num_forward_labels_in_enu{}, num_backward_labels_in_enu{};
  double **rc2_till_this_bin_in_forward_sense{};
  double **rc2_bin_in_forward_sense{};

  std::vector<std::tuple<std::vector<int>, R1CMem, std::vector<int>>> Vertex2ActiveInOnePricingR1Cs;
  std::vector<std::pair<std::vector<int>, R1CMem>> vertex2_all_in_one_lp_r1cs{};

  std::vector<std::tuple<std::vector<std::tuple<int, int, int>>, std::vector<int>, std::vector<int>>>
	  Vertex2ActiveInOnePricingR1C_multi;

  std::vector<std::pair<std::vector<std::pair<int, int>>, std::vector<int> >> Vertex2AllInOneLPR1C_multi;

  std::vector<int> r1c_multi_denominator_in_cg{};//valid in CG
  std::vector<int> r1c_multi_denominator_in_lp{};//valid in CG

  std::vector<int> cstr_index{};//used in enumeration, delete rows

  size_t mem4_pricing{};
  size_t mem4_mem{};
  int num_arti_vars{};

  int num_strong_arti_vars{};

  int seq_size_arti_vars{};

  int max_num_forward_graph_arc{};

  bool if_stop_arc_elimination{};

  std::unordered_map<size_t, int> tell_which_bin4_arc_elimination_in_forward_sense{};

  std::unordered_map<std::pair<int, int>, std::vector<std::pair<Label *, double>>, PairHasher>
	  concatenate_labels_in_backward_cg{};
  int max_num_backward_graph_arc{};
  std::vector<std::vector<double>> main_resource_across_arcs_in_backward_sense{};
  std::unordered_map<size_t, int> tell_which_bin4_arc_elimination_in_backward_sense{};
  double **rc2_till_this_bin_in_backward_sense{};
  double **rc2_bin_in_backward_sense{};
  VecLabel **label_array_in_backward_sense{};
  VecLabel **if_exist_extra_labels_in_backward_sense{};

  int max_num_enu_col_pool{};

  double max_gap_improved_in_r1c{TOLERANCE};

  double gap_tolerance4_arc_elimination_n_enumeration{};

  double cut_gen_time_threshold_in_pricing{};

  bool if_in_enu_state{};

  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> map_edge_col_idx_in_enu{};
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> real_improvement_up{};// sum and cnt
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> real_improvement_down{};// sum and cnt
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> heuristic_improvement_up{};
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> heuristic_improvement_down{};
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> lp_testing_improvement_up{};
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> lp_testing_improvement_down{};
  std::vector<std::pair<std::pair<int, int>, double >> branch_pair_val{};
  std::vector<std::pair<int, int>> branch_selected_from_model1{};
  std::vector<std::pair<int, int>> branch_selected_from_model2{};
  std::unordered_map<std::pair<int, int>, int, PairHasher> branch_choice{};

  std::vector<std::vector<double>> main_resource_across_arcs_in_forward_sense{};

  int num_buckets_per_vertex{};

  bool if_force_not_regenerate_bucket_graph{};

  double guessed_ub{};

  bool if_enumeration_suc{};

  double step_size{};

  double hard_rollback_factor{};

  double last_max_time_labeling{};

  bool if_populate_dump_idx4_br_col{};
  std::vector<size_t> dump_idx4_left_br_col{};//record the size_t idx in colPool
  std::vector<size_t> dump_idx4_right_br_col{};

  std::vector<std::tuple<bool, int, yzzLong>> reset_cut_mem{};//cut index in r1c & type & memory

  std::unordered_map<std::pair<int, int>, std::vector<std::pair<Label *, double>>, PairHasher>
	  concatenate_labels_in_forward_cg{};

  size_t label_assign{};
  size_t route_in_mem_assign{};
  size_t route_in_pricing_assign{};
  int aver_route_length{};
  bool if_arc_elimination_succeed{};
  bool if_arc_elimination_tried_but_failed{};

  int rollback{};

  bool force_not_rollback{true};

  double gap4_last_fail_due2_difficulty_in_arc_elimination{1};

  bool enumeration_mode{};//0: use info provided by arc elimination; 1: standard enumeration

  double last_enumeration_fail_gap{1};

  double time_resolve_lp4_iter{};

  double time_pricing4_iter{};

  bool if_ban_convert_dual{};

  double old_ub{};

  bool if_fill_mem{};

  double round_up_tolerance{};

  bool if_can_arc_elimination_by_exact_cg{};

  bool final_decision_4_arc_elimination{};

  bool final_decision_4_enumeration{};

  double global_gap{1};

  std::vector<std::vector<std::vector<std::vector<int>>>> record_map_rank1_combinations{};

  int count4_tolerance4_try_enumeration_when_arc_elimination_fails{};

  std::unordered_map<int, std::vector<std::vector<int>>> recorded_combinations4_r1c{};

  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, int>>> map_rank1_multiplier{};

  double num_dominance_checks{};

  std::unordered_map<yzzLong, std::unordered_set<int>> cut_record{};

  void chooseCuts(const std::vector<std::pair<std::vector<int>, int>> &tmp_cuts,
				  std::vector<std::pair<std::vector<int>, int>> &cuts,
				  int numCuts);

  double opt_gap{};

  std::unordered_map<int, std::unordered_set<int>> map_rank1_multiplier_dominance{};

#ifdef MASTER_VALVE_ML
  ML ml;
#endif

  std::pair<double, int> ratio_dominance_checks_non_dominant{};//sum, int


  void applyRCF(BbNode *node, int round, bool if_verbose);

  void writeRCFIP(double ub, const Solver &solver) const;

  void checkBucketIJ(BbNode *node, int i, int j) const;

  void deleteArcByFalseBranchConstraint(BbNode *node, std::pair<int, int> edge) const;

  void regenerateGraphBucket(BbNode *node);

  void writeMapEdgeColIndexInEnum(BbNode *node);

  void doDefaultSB(BbNode *node, std::pair<int, int> &info);

  [[nodiscard]] double calculateGapImprovement(double nowVal, double b4Val) const;

  void createBranchingSet(BbNode *node);

  void writeISAccount(const std::vector<std::pair<int, int>> &tmp_br);

  void writeUseModelAccount(const std::vector<std::pair<int, int>> &tmp_br, int top_n);

  void generateModelPhase1CombineFashion(BbNode *node);

  void testInitialScreen(BbNode *node);

  void initialScreen(BbNode *node,
					 bool if_record_source,
					 bool if_fill_candidates_gap,
					 int num,
					 double pseudo_frac);

  void readSolutionFile(bool if_force);

  void reviseBranchPair();

  void updateUpperBoundEdgeSolution();

  void simulateWriteLPPseudoCost(BbNode *node);

  std::unordered_set<std::pair<int, int>, PairHasher> edge_sol{};

  void testLP(BbNode *node, int num, bool if_writeBranch_pair, bool if_record_sb_scores = false,
			  bool if_record_LP_improvement = true);

  void randomPickLP(BbNode *node, int num);

  void initialPickLP(BbNode *node, int num);

  void addSearchCrazy(int plan_idx,
					  const std::vector<int> &c,
					  const std::vector<int> &w_no_c,
					  double &new_vio,
					  int &add_j,
					  const std::vector<std::unordered_map<int, int>> &v_r_map,
					  const std::vector<double> &frac_routes,
					  std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void constructCutsCrazy(
	  const std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio,
	  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
	  std::vector<std::pair<std::vector<int>, int>> &cuts
  );

  void removeSearchCrazy(int plan_idx,
						 const std::vector<int> &c,
						 double &new_vio,
						 int &remove_j,
						 const std::vector<std::unordered_map<int, int>> &v_r_map,
						 const std::vector<double> &frac_routes,
						 std::unordered_map<yzzLong,
											std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void swapSearchCrazy(int plan_idx,
					   const std::vector<int> &c,
					   const std::vector<int> &w_no_c,
					   double &new_vio,
					   std::pair<int, int> &swap_i_j,
					   const std::vector<std::unordered_map<int, int>> &v_r_map,
					   const std::vector<double> &frac_routes,
					   std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void operationsCrazy(
	  Rank1MultiLabel &label,
	  const std::vector<std::unordered_map<int, int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  int &i,
	  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
	  std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void startSeedCrazy(
	  const std::vector<std::unordered_map<int, int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  const std::vector<std::pair<std::vector<int>, std::vector<int>>> &c_N_noC,
	  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
	  std::vector<Rank1MultiLabel> &rank1_multi_label_pool,
	  std::unordered_map<
		  yzzLong, std::vector<
			  std::pair<std::vector<int>, double>>> &map_cut_plan_vio,
	  int &num_label
  );

  void addR1CAtOnceInEnum(BbNode *node,
						  const std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateN4ColpoolMapInEnum(BbNode *node,
								  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> &colmap);

  static void addBranchConstraint2ColPoolInEnumByColMap(BbNode *node,
														const std::pair<int, int> &edge);

  void changeModelForBetterDual(BbNode *node);

  void testCG(BbNode *node,
			  bool if_exact_CG,
			  bool if_record_product_val,
			  bool if_record_improvement,
			  bool if_force_complete = false);

  void testExactCG(BbNode *node);

  void useModelPhase1Separate(BbNode *node, int num);

  std::unordered_map<int, std::pair<int, double>> top_n_percentage_stage1{};
  std::vector<std::unordered_map<int, std::pair<int, double>>> top_n_percentage_stage2_takeout_n{};

  std::unordered_map<int, int> map_level_num_candidates{};

  void doSB(BbNode *node, std::pair<int, int> &info);

  void searchForParameterPhase1(BbNode *node);

  void fillMemoryFirst(BbNode *node,
					   const std::vector<std::vector<int>> &routes,
					   const std::vector<double> &frac_routes,
					   std::vector<std::pair<std::vector<int>, int>> &cuts);

  void useDefaultSB(BbNode *node);

  void findNonGuillotineMemorySets(BbNode *node, bool &if_empty);

  void deleteColumnByNGMemory(BbNode *node);

  void augmentNonGuillotineRound(BbNode *node);

  void changeEnumMatByCuts(BbNode *node);

  void addR1CMultiInEnum(BbNode *node, const std::pair<std::vector<int>, int> &cut);

  void addR1CMulti(BbNode *node, const std::pair<std::vector<int>, int> &cut, const std::set<int> &mem, int cut_index);

  void findMemoryForRank1Multi(const std::vector<std::vector<int>> &routes,
							   const std::vector<std::unordered_map<int, int>> &v_r_map,
							   const std::pair<std::vector<int>, int> &cut_pair,
							   std::set<int> &mem,
							   bool &if_suc);

  void printCutsInformation(BbNode *node) const;

  void findPlanForRank1Multi(const std::vector<int> &vis, int denominator, yzzLong &mem,
							 std::vector<std::set<int>> &segment, std::vector<std::vector<int>> &plan);

  void combinationUtilAddOne(const std::vector<int> &arr,
							 std::vector<int> &tmp,
							 std::vector<std::vector<int>> &data,
							 int start,
							 int end,
							 int index,
							 int r);
  void findCombinationsForRankCutSetToGetBestMultiplier(int plan_idx,
														const std::vector<int> &cset,
														std::vector<int> &new_cset,//the c will be sorted!
														int c,
														double &new_vio,
														const std::vector<std::vector<int>> &v_r_map,
														const std::vector<double> &frac_routes,
														bool &if_succeed);
  void generatePureMap(std::vector<std::vector<std::vector<int>>> &map, const std::vector<int> &sorted_info);

  void combineAll(int n, int r, std::vector<std::vector<int>> &data);

  void combinationUtilOuter(std::vector<std::vector<int>> &data, int n);

  void addInRank1HeuristicSeparate(int plan_idx,
								   const std::vector<int> &c,
								   const std::vector<int> &w_no_c,
								   double &new_vio,
								   int &choice_j,
								   const std::vector<std::vector<int>> &v_r_map,
								   const std::vector<double> &frac_routes);

  void removeInRank1HeuristicSeparate(int plan_idx,
									  const std::vector<int> &c,
									  double &new_vio,
									  int &choice_j,
									  const std::vector<std::vector<int>> &v_r_map,
									  const std::vector<double> &frac_routes);

  void generatePermutations(std::unordered_map<int, int> &count_map,
							std::vector<int> &result,
							std::vector<std::vector<int>> &results,
							int remaining);

  void exactFindBestPermutationForOnePlan(const std::vector<std::unordered_map<int, int>> &v_r_map,
										  const std::vector<double> &frac_routes,
										  const std::vector<int> &cut,
										  int plan_idx,
										  double &vio,
										  std::unordered_map<
											  yzzLong, std::vector<
												  std::pair<std::vector<int>, double>>> &map_cut_plan_vio);

  void swapInRank1HeuristicSeparate(int plan_idx,
									const std::vector<int> &c,
									const std::vector<int> &w_no_c,
									double &new_vio,
									std::pair<int, int> &choice_i_j,
									const std::vector<std::vector<int>> &v_r_map,
									const std::vector<double> &frac_routes);

  void constructOperationsAggressive(
	  Rank1MultiLabel &label,
	  const std::vector<std::vector<int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  int &i,
	  std::unordered_map<int,
						 std::vector<std::tuple<std::vector<int>,
												int,
												double>>> &generated_rank1_multi_pool);

  void checkR1CTotally(BbNode *node);

  void constructCutsAggressive(std::unordered_map<int,
												  std::vector<std::tuple<std::vector<int>,
																		 int,
																		 double>>> &generated_rank1_multi_pool,
							   std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateR1C3sAggressive(const std::vector<std::vector<int>> &ele_routes,
							   const std::vector<std::vector<int>> &non_ele_routes,
							   const std::vector<double> &frac_ele_routes,
							   const std::vector<double> &frac_non_ele_routes,
							   std::vector<std::pair<std::vector<int>, int>> &cuts) const;

  void constructVRVecAggressive(const std::vector<std::vector<int>> &routes,
								std::vector<std::vector<int>> &v_r_map,
								std::vector<std::unordered_map<int, int>> &v_r_map2,
								std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> &c_n_w_no_c,
								std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map);

  void findR1CMultiAggressive(const std::vector<std::vector<int>> &routes,
							  const std::vector<double> &frac_routes,
							  std::vector<std::unordered_map<int, int>> &v_r_map,
							  std::vector<std::pair<std::vector<int>, int>> &cuts);

  void constructSeedAggressive(
	  const std::vector<std::vector<int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool,
	  std::vector<Rank1MultiLabel> &rank1_multi_label_pool,
	  std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map,
	  std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> &c_n_w_no_c,
	  int &num_label);

  void generateOptimalMultiplierForR1CMulti();

  virtual void getLowerBoundofMinimumNumberCars();

  void findR1CMulti(const std::vector<std::vector<int>> &routes,
					const std::vector<double> &frac_routes,
					const std::vector<std::vector<int>> &vbasis,
					std::vector<std::pair<std::vector<int>, double>> &cut_info,
					std::vector<std::tuple<std::vector<int>, int, double>> &multi_cut_info);

  void newFindR1CMulti(const std::vector<std::vector<int>> &routes,
					   const std::vector<double> &frac_routes,
					   std::vector<std::unordered_map<int, int>> &v_r_map,
					   std::vector<std::pair<std::vector<int>, int>> &cuts);

  void searchCrazy(const std::vector<std::vector<int>> &routes,
				   const std::vector<double> &frac_routes,
				   std::vector<std::unordered_map<int, int>> &v_r_map,
				   std::vector<std::pair<std::vector<int>, int>> &cuts);

  void constructVRMapAndSeedCrazy(const std::vector<std::vector<int>> &routes,
								  std::vector<std::unordered_map<int, int>> &v_r_map,
								  std::vector<std::pair<std::vector<int>,
														std::vector<int>>> &c_N_noC);

  void constructVRVec(const std::vector<std::vector<int>> &routes,
					  std::vector<std::vector<int>> &v_r_vec,
					  std::vector<std::unordered_map<int, int>> &v_r_map) const;

  void constructCuts(const std::vector<std::unordered_map<int, int>> &v_r_map,
					 const std::vector<std::vector<int>> &seed,
					 const std::vector<double> &frac_routes,
					 std::vector<std::pair<std::vector<int>, int>> &cuts// cut self, plan idx
  );

  void constructCutsByResourceGapOnly(const std::vector<std::unordered_map<int, int>> &v_r_map,
									  const std::vector<std::vector<int>> &seed,
									  const std::vector<double> &frac_routes,
									  std::vector<std::pair<std::vector<int>, int>> &cuts// cut self, plan idx
  );

  void constructMemory(BbNode *node, const std::vector<std::vector<int>> &routes,
					   const std::vector<std::unordered_map<int, int>> &v_r_map,
					   const std::vector<std::pair<std::vector<int>, int>> &cuts,
					   std::vector<std::tuple<std::vector<int>, int, int, std::set<int>>> &full_cuts);

  template<bool if_symmetry>
  void construct_v_resourceGap(const std::vector<std::vector<int>> &routes, std::vector<double> &v_resourceGap);

  void constructSeed(const std::vector<std::vector<int>> &routes,
					 const std::vector<std::vector<int>> &v_r_vec,
					 const std::vector<double> &v_resourceGap,
					 std::vector<std::vector<int>> &seed);

  static void checkR1CMultiLook(BbNode *node, int start);

  template<bool dir, bool if_symmetry>
  int enumerateHalfwardRoutes(BbNode *node,
							  const double *r1c_to_pi,
							  const double *r1c_multi_to_pi,
							  std::unordered_map<yzzLong, std::tuple<Label *, Label *, double>> &Tags,
							  std::vector<Label *> **copy_bucket,
							  int &num_routes_now);

  void generateAllR1Cs(BbNode *node, bool if_use_MIP);

  void newGenerateAllR1Cs(BbNode *node, int aggressive_level);

  void addR1C(BbNode *node, const std::vector<int> &cut, const std::set<int> &mem, int cut_index);

  void addR1CInEnum(BbNode *node, const std::vector<int> &cut);

  void heuristicFindBestPermutationForOnePlan(
	  const std::vector<std::unordered_map<int, int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  const std::vector<int> &cut, int plan_idx,
	  std::vector<int> &new_cut,
	  double &vio);

  void addR1CAtOnce(BbNode *node,
					const std::vector<std::tuple<std::vector<int>, int, int, std::set<int>>> &full_cuts);

  void runCertainNumberOfHeuristicCGs(BbNode *node, int iteration);

  void generateAllR1CsInEnum(BbNode *node, bool if_use_MIP);

  void refineArcsHeuristic(BbNode *node) const;

  explicit CVRP(const InstanceData &instanceData);

  virtual ~CVRP();

  void buildModel();

  bool addBranchCutToUnsolved(BbNode *node, const std::pair<int, int> &info);

  void separateHybridCuts(BbNode *&node);

  void findMemoryForR1CsMode3ByEnumerationAndMIP(yzzLong v_comb,
												 std::set<int> &mem,
												 const std::vector<std::vector<int>> &routes,
												 bool &if_suc);

  void getMemoryByMIP(const std::vector<std::vector<std::vector<int>>> &array,
					  const std::vector<std::vector<std::set<int>>> &vec_segment,
					  std::set<int> &mem, bool &if_suc);

  static void combinations(const std::vector<std::vector<std::vector<int>>> &array,
						   const std::vector<std::vector<std::set<int>>> &vec_segment,
						   int i,
						   const std::vector<int> &accum,
						   const std::set<int> &mem,
						   int &record_min,
						   std::set<int>
						   &new_mem);

  static void combinationUtil(const std::vector<int> &arr, std::vector<int> &tmp,
							  std::vector<std::vector<int>> &data,
							  int start, int end,
							  int index, int r);
  void separateRCCs(BbNode *&node);

  void generateRCCs(BbNode *node);

  void randomSelection(BbNode *node);

  void solveLPInLabeling(BbNode *node, bool if_open_heur = true, bool if_open_exact = true, bool if_record_sol = true);

  void priceLabeling(BbNode *node);

  void rollbackEasyWay(BbNode *node, int old_num);

  void assignMemory();

  bool solveSBModel();

  void recordOptimalColumn(BbNode *node, bool if_force_rewrite = false);

  int optimizeLPForOneIteration(BbNode *node, double prior_value);

  void reformulateIntegerProgramSolution(BbNode *node);

  void constructSeedMode3(const std::vector<std::vector<int>> &routes,
						  const std::vector<std::vector<int>> &v_r_vec,
						  const std::vector<double> &v_resourceGap,
						  std::vector<std::vector<int>> &seed);

  void newFindR1CMultiMode3(const std::vector<std::vector<int>> &routes,
							const std::vector<double> &frac_routes,
							std::vector<std::unordered_map<int, int>> &v_r_map,
							std::vector<std::pair<std::vector<int>, int>> &cuts);

  void cleanIndexColForNode(BbNode *node, int beg, bool if_only_rcfixing = false);

  [[nodiscard]] double calculateDifference(double tmp_val, double prior_val, bool if_chg = false) const;

  void writeColumnToMemory(BbNode *node);

  void constructMap(BbNode *node, int beg) const;

  void getNewConstraintCoefficientByEdge(BbNode *node,
										 std::pair<int, int> edge,
										 int *c_ind,
										 double *c_val,
										 int &num_nz);

  void getCoefficientRCC(BbNode *node, Rcc &rcc, int *c_ind, double *c_val, int &num_nz) const;

  void convertVertexToR1CsInPricing(BbNode *node);

  void convertVertexToR1CsInOneLP(BbNode *node);

  template<bool if_symmetry>
  int generateColsByBidir(BbNode *node);

  void addColumns(BbNode *node, int &ccnt);//直接从pri一直遍历到curr吧

  int generateColumnsByLighterHeuristic(BbNode *node);

  int generateColumnsByHeavierHeuristic(BbNode *node);

  double calculateOptimalGap(BbNode *node) const;

  void eliminateArcs(BbNode *node);

  template<bool dir>
  void populateRC2TillThisBinNRC2Bin(BbNode *node) const;

  void initializeLabels(BbNode *node, int mode,
						bool if_resetLabelPoint, std::tuple<bool, int, bool> control_cleanAllPtr);
  void deleteNonActiveCutsBySlack(BbNode *node, bool if_rcc, bool if_r1c);

  void terminateByMIP(BbNode *node);

  void enumerateMIP(BbNode *&node);

  void cleanAllPointers(BbNode *node, int mode, bool if_clear_concatenate);

  void cleanColumnsNonElement(BbNode *node);

  void cleanColumnsRCGreaterThanOptimalGapInLP(BbNode *node);

  void organizeColumnsInMemoryToPricing(BbNode *node);

  void solveLPByInspection(BbNode *node, bool if_only_need_value,
						   bool if_heuristic, bool if_record_sol);

  void deleteBranchCutsAndR1C1s(BbNode *node);

  int generateColumnsByInspection(BbNode *node,
								  bool if_only_need_value);

  void addColumnsByInspection(BbNode *node, const std::vector<int> &Col_added);

  void regenerateEnumMat(BbNode *node, BbNode *node2, bool if_force = false);

  void generateVertex2IndexColsAndEdge2IndexCols(BbNode *node);

  void recoverR1CsInEnum(BbNode *node);

  void optimizeLPForOneIterationInEnum(BbNode *node);

  void terminateNode(BbNode *&root_node);

  void addBranchCutToUnsolvedInEnum(BbNode *node, const std::pair<int, int> &info);

  void reviseEnumColInfoByBrC(BbNode *node, BbNode *out_node, const Brc &bf);

  void cleanColumnsInLPByNewBranchCut(BbNode *node, const Brc &bf);

  void generateRCCsInEnum(BbNode *node);

  void deleteNonActiveCutsSafely(BbNode *node, int old_num);

  void writeColumnsInPricingPool(BbNode *node, int &index);

  void determineIfArcElimination(BbNode *node);

  template<bool dir, bool if_last_half, bool if_symmetry, bool if_std_optgap, bool if_res_updated>
  void updateLabel(double res, Label *&ki, int i, int j, int &bj,
				   const double *r1c_to_pi,
				   const double *r1c_multi_to_pi,
				   int min_sort_b,
				   bool &if_suc);

  template<bool dir>
  void doDominance(Label *&ki, int j, int bj, const double *r1c_to_pi,
				   const double *r1c_multi_to_pi, bool &if_suc);

  template<bool dir>
  bool define_dir_resource(double a, double b) {
	if constexpr (dir) return a < b;
	else return a > b;
  }
  template<typename T, bool dir, bool if_symmetry>
  void concatenateTestKernelInArcElimination(int i,
											 int b,
											 const std::vector<T> &arc,
											 int dim_sq,
											 bool *stateBetween2Buckets,
											 int *latest_bucket,
											 const double *r1c_to_pi,
											 const double *r1c_multi_to_pi);
  template<bool dir, bool if_symmetry, bool if_check_res, bool if_std_optgap>
  void concatenateOneLabelWithOtherLabels(Label *ki, int j, int arr_bj, double tmp_rc, double tmp_mainResource,
										  const double *r1c_to_pi,
										  const double *r1c_multi_to_pi,
										  int &if_state);

  template<bool dir, bool if_symmetry>
  void eliminateBucketArcs(BbNode *node, const double *r1c_to_pi,
						   const double *r1c_multi_to_pi,
						   int dim_sq,
						   bool *stateBetween2Buckets,
						   int *latest_bucket);

  template<bool dir, bool if_symmetry>
  void eliminatebuketArc4Depot(BbNode *node);

  template<bool dir, bool if_symmetry>
  void concatenatePhaseInArcElimination(const double *r1c_to_pi,
										const double *r1c_multi_to_pi);

  void addPathByRC(double path_rc, Label *ki, Label *kj, int num);

  void determineIfEnumeration(BbNode *node);

  template<bool dir>
  void populateTellWhichBin4ArcElimination();

  void getEdgeInfo(BbNode *node, bool if_br) const;

  void deleteNonActiveCutsByDual(BbNode *node, bool if_rcc_by_slack);

  void deleteNewAddedNonActiveCutsBySlack(BbNode *node, int olNum, bool if_keep_rcc);

  void deleteNewAddedNonActiveCutsBySlackInEnum(BbNode *node, int olNum, bool if_delete_rcc);

  void initializeBucketGraphForNode(BbNode *node);

  [[nodiscard]] int checkMaximumNumberCustomers() const;

  [[nodiscard]] static int checkMaximumNumberR1Cs(int num_r1cs);

  [[nodiscard]] int checkCSTLimit() const;

  [[nodiscard]] static int checkIfAssignGreaterMAXINT(const size_t &value, const std::string &str);

  void reallocateLabel();

  void reallocateMemoryPool();

  [[nodiscard]]  int checkMemoryPool() const;

  void reallocatePricingPool();

  [[nodiscard]]  int checkPricingPool() const;

  void rollbackToPreviousState(BbNode *node, const std::vector<Rcc> &old_rcc,
							   const std::vector<R1c> &old_r1c,
							   const std::vector<R1cMulti> &old_r1c_multi,
							   const std::vector<Brc> &old_brc, const std::vector<size_t> &old_col_idx);

  static int inverseLastBranchConstraint(char sense, double rhs, Solver &solver);

  static int addBranchConstraint(int numnz,
								 int *cind,
								 double *cval,
								 char sense,
								 double rhs,
								 const char *constrname, Solver &solver);

  void changeBranchConstraint(double *val,
							  int *cind,
							  int *vind,
							  int numnz,
							  const int *old_ind,
							  const double *old_val,
							  char sense,
							  double rhs,
							  Solver &solver);

  void initializeBucketGraph();

  void checkIfCutsLegal(BbNode *node) const;

  void lateProcessing();

  virtual double transformCost(double x);

  void solveMIP(BbNode *node, bool if_inEnu);

  [[nodiscard]] double ceilTransformedNumberRelated(double x) const;

  void initializeLabels();

  template<typename T, bool dir, bool if_last_half, bool if_symmetry>
  int extendKernel4Exact(Label *&ki,
						 int i,
						 double res,
						 const std::vector<T> &arc,
						 const double *r1c_to_pi,
						 const double *r1c_multi_to_pi,
						 int min_sort_b);

  template<bool dir>
  void obtainjumpArcs(BbNode *node, std::bitset<2> **bitMap) const;

  template<typename T>
  int extendKernel4Exact_Backward(Label *&ki,
								  int i,
								  double res,
								  const std::vector<T> &arc,
								  const double *r1c_to_pi,
								  const double *r1c_multi_to_pi);

  void runHalfBackwardLabeling(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  int concatenateColumnsNoneSymmetry(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  int generateColumnsByBidirectionalNoneSymmetry(BbNode *node);

  bool decreaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end);
  template<typename T>
  int extendKernel4ArcElimination_last_half_Backward(Label *&ki,
													 int i,
													 double res,
													 const std::vector<T> &arc,
													 const double *r1c_to_pi,
													 const double *r1c_multi_to_pi);

  template<bool dir>
  void checkIfDominated(Label *&ki, int i, int b, const double *r1c_to_pi,
						const double *r1c_multi_to_pi,
						bool &if_suc);

  template<bool dir>
  void checkIfNoLabelsLeft(int &i, int b, int &min_sorted_b, bool &if_suc);

  template<bool dir, bool if_last_half, bool if_symmetry>
  void runLabeling(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  int extendKernelForArcEliminationInnerBackward(Label *&ki,
												 int i, int j, double &tmp_mainResource,
												 const double *r1c_to_pi,
												 const double *r1c_multi_to_pi);

  int lastHalfBackwardInArcElimination(BbNode *node, const double *r1c_to_pi,
									   const double *r1c_multi_to_pi);

  int backwardConcatenateInArcElimination(const double *r1c_to_pi,
										  const double *r1c_multi_to_pi);

  void eliminateBackwardArcs(BbNode *node, const double *r1c_to_pi,
							 const double *r1c_multi_to_pi,
							 int dim_sq,
							 bool *stateBetween2Buckets,
							 int *latest_bucket);

  void obtainBackwardJumpArcs(BbNode *node, std::bitset<2> **bitMap) const;

  int enumerateHalfBackwardRoutes(BbNode *node,
								  const double *r1c_to_pi,
								  const double *r1c_multi_to_pi,
								  std::vector<Label *> **copy_bucket);

  void setTailOffStandardAndRollBackStandard() const;

  int concatenateRoutesPriorForwardInEnumeration(BbNode *node,
												 const double *r1c_to_pi,
												 const double *r1c_multi_to_pi,
												 std::unordered_map<yzzLong,
																	std::tuple<Label *, Label *, double>> &Tags,
												 int &num_routes_now
  );
  template<typename T>
  void extendKernel4LightHeur(Label *&ki,
							  int i,
							  double res,
							  const std::vector<T> &arc,
							  const double *r1c_to_pi,
							  const double *r1c_multi_to_pi,
							  bool if_return_to_depot);
  template<typename T>
  void extendKernel4HeavierHeur(Label *&ki,
								int i,
								double res,
								const std::vector<T> &arc,
								const double *r1c_to_pi,
								const double *r1c_multi_to_pi,
								bool if_return_to_depot);

  int lastHalfForwardInArcElimination(BbNode *node,
									  const double *r1c_to_pi,
									  const double *r1c_multi_to_pi);

  int forwardConcatenateInArcElimination(const double *r1c_to_pi, const double *r1c_multi_to_pi);

  void findAccML(BbNode *node);

  void eliminateBucketArcs(BbNode *node,
						   const double *r1c_to_pi,
						   const double *r1c_multi_to_pi,
						   int dim_sq,
						   bool *stateBetween2Buckets,
						   int *latest_bucket);

  template<typename T>
  int extendKernel4ArcElimination_last_half_Forward(Label *&ki,
													int i,
													double res,
													const std::vector<T> &arc,
													const double *r1c_to_pi,
													const double *r1c_multi_to_pi);
  int extendKernelForArcEliminationInnerForward(Label *&ki,
												int i, int j, double &tmp_mainResource,
												const double *r1c_to_pi,
												const double *r1c_multi_to_pi);

  bool runColumnGenerationType(BbNode *node, int mode);

  void findNonActiveCuts(BbNode *node);

  void deleteNonactiveCuts(BbNode *node, std::vector<int> &nonactive_cuts);

  int runLighterHeuristicLabeling(BbNode *node);

  int runHeavierHeuristicLabeling(BbNode *node);

  void runHalfForwardLabeling(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  bool increaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end);

  template<bool if_symmetry>
  int concatenateCols_prior_forward(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  void runLabelingForArcElimination(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  void eliminateBucketArcs(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  void obtainJumpArcs(BbNode *node) const;

  bool enumerateRoutes(BbNode *node, const PtrAllR1CS &ptrAllR1Cs);

  virtual void setResourceInBucketGraph();

  virtual void checkSolutionFeasibleByCapacity(bool &feasible) {};

  virtual void popArcGraph(BbNode *node) {};

  virtual void cleanColumnsNonFeasible(BbNode *node) {};

  void generateR1C3s(const std::vector<std::vector<int>> &ele_routes,
					 const std::vector<std::vector<int>> &non_ele_routes,
					 const std::vector<double> &frac_ele_routes,
					 const std::vector<double> &frac_non_ele_routes,
					 std::vector<std::pair<std::vector<int>, double>> &cut_info) const;

  void newGenerateR1C3s(const std::vector<std::vector<int>> &routes,
						const std::vector<double> &frac_routes,
						std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateR1C1s(const std::vector<std::vector<int>> &non_ele_routes,
					 const std::vector<double> &frac_none_ele_routes,
					 std::vector<std::pair<std::vector<int>, double>> &cut_info) const;

  void newGenerateR1C1s(const std::vector<std::vector<int>> &non_ele_routes,
						const std::vector<double> &frac_none_ele_routes,
						std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateHighDimensionR1Cs(const std::vector<std::vector<int>> &routes,
								 const std::vector<double> &frac_routes,
								 const std::vector<int> &cut_type,
								 std::vector<std::pair<std::vector<int>, double>> &cut_info);

#ifdef WRITE_ENUMERATION_TREES
  std::unordered_map<yzzLong, std::tuple<std::vector<int>, double, size_t>>
	  enumeration_col_idx{};// sequence, cost, index
  void writeEnuTree(BbNode *node);
  void populateEnuCols(std::vector<size_t> &col_idx, const std::vector<std::pair<size_t, double>> &col_infos);
  void writeEnuCols();
  void writeEnuCuts(BbNode *node,
					const std::vector<size_t> &col_idx,
					const std::vector<std::pair<size_t, double>> &col_infos);
#endif

#ifdef READ_ENUMERATION_TREES
  std::string tree_path{};
  std::string col_pool_path{};
  void solveEnumTree();
  void readEnumTree(BbNode *node, std::vector<std::pair<size_t, double>> &col_idx);
  void readColumnPool(BbNode *node, const std::vector<std::pair<size_t, double>> &col_idx, std::vector<double> &objs);
  void restoreModel(BbNode *&node);
  void deleteBranchCutsSpecifiedForReadEnumTree(BbNode *const node) const;
#endif

#ifdef BRANCH_FASHION_MEM_SAVING
  StackTree bbt{};
  StackTree sub_bbt{};
  std::priority_queue<double, std::vector<double>, std::greater<>> supp_queue;
  std::unordered_set<double, MyDoubleHash, MyDoubleEqual> supp_set{};
  bool checkTimeFail();
  void printInfo(BbNode *node);
  void SB_MemSave(StackTree &tree);
  void addBrCut(BbNode *node, BbNode *&node2, const std::pair<int, int> &info);
  void tackleUnsolvedNode(BbNode *&node);
  void solveModel();
  void terminateNodeMemSave(BbNode *&root_node);
#else
  QueueTree bbt;
  QueueTree sub_bbt;
#endif
#ifdef if_draw_BBT_graph
  std::unordered_map<int, std::pair<std::pair<int, int>, std::string>> BBT_nodes_name;// lchild node, rchild node
  void drawBBTgraph();
  void constructBBNodeName(BbNode *node, const std::string &terminateReason = "");
#endif
#ifdef MASTER_VALVE_ML

  void getTrainingDataPhase1(BbNode *node);

  void useModelInPhase2(BbNode *node, int num);

  void getTrainingDataInPhase2(BbNode *node);

  void useMLInGeneralFramework(BbNode *node);

  void generateModelInPhase2(BbNode *node);

  void generateModelInPhase1(BbNode *node);

#ifdef  USE_M_DYNAMICS
  void useMLVaryLPTesting(BbNode *node);
  void giveDynamicM(BbNode *node, int &num);
  void evaluateM1();
  std::vector<std::pair<int, int>> cp_branch_pair{};
  double alpha{};
  double est_m{};
  double f{};
  double opt_k{};
  bool is_use_full_k{};
#endif

#endif

  void printBrDecisions();

  static void printInfoLabeling(int iter, int num_added_col, int num_col,
								int num_row,
								double mt, double spt, double et,
								double lp, double lb_val, double ub);

  static void printInfoLabeling(int num_col, int num_row, double et, double lp, double ub);
  void printOptIntSol();

#ifdef GENERATE_COL_BY_MIP
  int num_mip_var{};
  int num_mip_edge{};
  int num_mip_constr{};
  std::unordered_map<std::pair<int, int>, int, PairHasher> map_edge_mip_idx{};
  std::unordered_map<int, std::pair<int, int>> map_mip_idx_edge{};
  void findNgSubTour(const std::vector<double> &sol, std::vector<int> &tour);
  void findBucketArcTour();
  void setObjCoeffs(BbNode *node);
  void buildMIPModel(BbNode *node);
  void cbSolve(BbNode *node);
#endif
};

bool operator==(const Rcc &lhs, const Rcc &rhs);

auto CmpLabelRCLess(const Label *l1, const Label *l2) -> bool;

std::string readInstanceFile(const std::string &file_name, int line);

std::string generateInstancePath(int argc, char *argv[]);

float sqrt_self(float x);

double pow_self(double x, int n);

void self_mkdir(const std::string &path);
#endif
