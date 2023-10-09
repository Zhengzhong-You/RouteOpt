//
// Created by Zhengzhong You on 1/20/22.
//

#ifndef CVRP_CVRP_HPP
#define CVRP_CVRP_HPP

#include "Experiment.hpp"
#include "BBNODE.hpp"
#include "MACRO.hpp"
#include "CONFIG.hpp"
#include "solver.hpp"
#include "capsep.h"
//#include "ML.hpp"
//#include "ML2.hpp"
#include "ML3.hpp"
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
//#include <glob.h>
#include <climits>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <list>
//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <limits>
//#include <unistd.h>
#include <cstdio>
#include <experimental/filesystem>

#ifdef HGS_APPLIED
#include "my_hgs.hpp"
#endif

#ifdef DELUXING_APPLIED
#include "deLuxing.hpp"
#endif

using yzzLong = std::bitset<MaxNum_Customers>;
using R1CMem = std::bitset<MaxNum_R1Cs>;

struct LABEL {
  double Sum_MainResource{};
  //start from 0 and indicates the index of the last element of sequence
  int IdxEndSeq{};
  int EndVertex{};
  double Cost{};
  //rank1_std
  R1CMem Rank1CutMem{};
  int numValidRank1Cut{};
  int *validRank1Cut{};
  //rank1_multi
  int *Rank1CutMem_multi{};
  int *validRank1Cut_multi{};
  int numValidRank1Cut_multi{};
  //
  double RC{};
  yzzLong PI{};
  bool if_extended{};
  int *Seq{};
  LABEL *PLabel{};
};

struct PtrAllR1Cs {
  double *r1c_to_pi{};
  int num_r1c_nzero{};
  double *r1c_multi_to_pi{};
  int num_r1c_multi_nzero{};

  PtrAllR1Cs(BBNODE *node, CVRP *cvrp);

  ~PtrAllR1Cs();
};

struct rank1_multi_label {
  std::vector<int> c;
  std::vector<int> w_no_c;
  int plan_idx{};
  double vio{};
  char search_dir{};
  //constructor
  rank1_multi_label(std::vector<int> c, std::vector<int> w_no_c, int plan_idx, double vio, char search_dir) :
	  c(std::move(c)),
	  w_no_c(std::move(w_no_c)),
	  plan_idx(plan_idx),
	  vio(vio),
	  search_dir(search_dir) {}
  rank1_multi_label() = default;
};

class CmpNodeVal {
 public:
  auto operator()(BBNODE *a, BBNODE *b) -> bool {
	return a->Val > b->Val;
  }
};

class CmpTupleDbl2 {
 public:
  auto operator()(std::tuple<int, int, double> a, std::tuple<int, int, double> b) -> bool {
	return std::get<2>(a) > std::get<2>(b);
  }
};

class CmpTupleDbl2Less {
 public:
  auto operator()(std::tuple<int, int, double> a, std::tuple<int, int, double> b) -> bool {
	return std::get<2>(a) < std::get<2>(b);
  }
};

struct VectorHasher {
  size_t operator()(const std::vector<int> &V) const {
	size_t hash = 0;
	for (auto &i : V) {
	  hash *= MaxNum_Customers;
	  hash += i;
	}
	return hash;
  }
};

struct Rank1_multi_pair_Hasher {
  size_t operator()(const std::pair<std::vector<int>, int> &V) const {
	size_t hash = 0;
	for (auto &i : V.first) {
	  hash *= MaxNum_Customers;
	  hash += i;
	}
	hash *= MaxNum_Customers;
	hash += V.second;
	return hash;
  }
};

using VecLabel = std::pair<std::vector<LABEL *>, int>;
using QueueTree = std::priority_queue<BBNODE *, std::vector<BBNODE *>, CmpNodeVal>;
using StackTree = std::stack<BBNODE *>;

struct Bucket {
  std::vector<int> BucketArcs;
  //first is the jump resource and the second is the destination
  std::vector<std::pair<double, int>> JumpArcs;
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
  /**
   * 仅仅是测试的代码，时间段是 jun 27，过了一周就删除
   */

//  size_t num_dominance_test{};
//  std::pair<double, size_t> succ_ratio{};
//  std::pair<double, size_t> fail_ratio{};

#ifdef NominalBranchingInEnu
  std::vector<int> BranchingColSet;
  QueueTree subBBT2;
  void NominalBranching(BBNODE *node, bool &if_suc);
#endif

  std::unordered_map<std::vector<int>,
					 std::vector<std::vector<std::vector<int>>>, VectorHasher> pureMap;

  std::unordered_map<std::vector<int>, std::vector<std::vector<int>>, VectorHasher> rank1_multi_mem_plan_map;
  int NumValidR1C_multi_InCG{};

  SOLVER Solver{};

  VecLabel **LabelArrayInForwardSense{};

  VecLabel **IfExistExtraLabelsInForwardSense{};

  std::chrono::time_point<std::chrono::high_resolution_clock> GloBeg,
	  GloEnd;
  double GloEps{};

  int NumBr{};
  int NumExploredNodes{};
  char if_IntNFeasible{};
  double MaxVio{};
  double *Demand{};

  int Rank1MemSizeLimit{};

  int BranchTimes{};
  int RealDim{};
  int NumCol{}, NumRow{};
#ifdef BranchFashion_MemSaving
  StackTree BBT;
  StackTree subBBT;
  std::priority_queue<double, std::vector<double>, std::greater<>> supp_queue;
  std::unordered_set<double, MyDoubleHash, MyDoubleEqual> supp_set;
  bool checkTimeFail();
  void printInfo(BBNODE *node);
  void SB_MemSave(StackTree &tree);
  void addBrCut(BBNODE *node, BBNODE *&node2, const std::pair<int, int> &info);
  void tackleUnsolvedNode(BBNODE *&node);
  void solveModel();
  void terminateNodeMemSave(BBNODE *&root_node);
#else
  QueueTree BBT;
  QueueTree subBBT;
#endif
#ifdef if_draw_BBT_graph
  std::unordered_map<int, std::pair<std::pair<int, int>, std::string>> BBT_nodes_name;// lchild node, rchild node
#endif
  void updateLB(double val);
  void resetEnv(BBNODE *node);
  double LPVal{}, LB{};
  int MaxNumRoute4MIP{};
  bool if_MIP_enumeration_suc{};
  double LB_transformed{}, UB{};
  double Obj4FirstCol{};
  int transformed_number{};
  int IdxNode{};
  double *X{}, *Pi{}, *Slack{}, *RC{};
  bool if_exact_labeling_CG{};
  std::vector<double> Pi4Labeling;//this pi is used for labeling and cutting, need special care!
  std::vector<double> Slack4Labeling;//this slack is used for labeling and cutting, need special care!
  int *ColPool4Pricing{};
  int *copyColPool4Pricing{};//used in heuristic find UB!
  //End index for col_pool
  size_t PoolBeg4Pricing{};
  size_t PoolBeg4copyPricing{};
  size_t PriorPoolBeg4Pricing{};
  //keep mem globally
  int *ColPool4Mem{};
  size_t PoolBeg4Mem{};
  //IPOptSol[0] records the size of the mem and the other records the sequence
  std::vector<std::vector<int>> IPOptSol;
//  int *IPOptSol{};

  std::vector<std::pair<int, int>> Branch_pair;
  std::unordered_map<std::pair<int, int>, double, PairHasher> branch_LP;
  std::vector<std::pair<int, int>> Branch_pair_from_pseudo;
  std::vector<std::pair<int, int>> Branch_pair_from_fractional;
  double SmallestRC{};
  double GapBetweenLastSmallestRCAndRCThreshold{};
  double rc_std{};
  double Cap{};
  double MaxMainResource{}, MeetPointResourceInBiDir{}, MeetPointResourceInBiDirEnu{};
  int Dim{}, K{};

  std::vector<std::vector<double>> InfoVertex;
  std::vector<std::vector<double>> CostMat4Vertex;
  std::vector<yzzLong> NGMem4Vertex;
  std::vector<yzzLong> Rank1SepHeurMem4Vertex;
  std::vector<int> SizeNGMem4Vertex;
  //Reduced cost related to edges
  std::vector<std::vector<double>> ChgCostMat4Vertex;
  std::vector<int> PriorMemSets;
  std::vector<double> lb4Vertex;
  std::vector<double> ub4Vertex;

  std::vector<std::tuple<LABEL *, LABEL *, double>> NegativeRCLabelTuple;
  std::unordered_map<yzzLong, std::tuple<LABEL *, LABEL *, double>> Map4EachNegativeRCRoute;

  int IdxGlo{};
  int SizeNGMem{};
  std::string FileName;

  LABEL *AllLabel{};
  int *AllSeq{};
  int SeqBeg{};
  int *AllValidR1Cs{};
  int *AllValidR1C_multi{};
  int *AllR1C_multi_Mem{};

  //Solver
  size_t MaxNonZeroEntry{};
  size_t *solver_beg{};
  int *solver_ind{};
  int *solver_ind2{};
  double *solver_val{};
  double *solver_obj{};

  int *const_for_branching{};
  int MaxNumEdge{500};
  double **ArcGraph{};
  //only for branching use
  double **ArcGraph_revised{};
  int NumEdge{};
  double NumForwardLabelsInEnu{}, NumBackwardLabelsInEnu{};
  //This will keep decreasing according to the first bin with label to the last bin with label
  double **RC2TillThisBinInForwardSense{};
  double **RC2BinInForwardSense{};
  int MaxLengthEleRoute{};

  //cut-self, cut-mem in (bitset and vector) //cut mem should always include cut-self
  std::vector<std::tuple<std::vector<int>, R1CMem, std::vector<int>>> Vertex2ActiveInOnePricingR1Cs;
  //cut-self and cut-mem in (bitset)
  std::vector<std::pair<std::vector<int>, R1CMem>> Vertex2AllInOneLPR1Cs;

  //cut-self (with augmented unit and denominator), no-mem and cut mem
  std::vector<std::tuple<std::vector<std::tuple<int, int, int>>, std::vector<int>, std::vector<int>>>
	  Vertex2ActiveInOnePricingR1C_multi;

  //cut-self (with augmented unit) and no-mem
  std::vector<std::pair<std::vector<std::pair<int, int>>, std::vector<int> >> Vertex2AllInOneLPR1C_multi;

  std::vector<int> R1C_multi_denominator_InCG;//valid in CG
  std::vector<int> R1C_multi_denominator_InLP;//valid in CG

  std::vector<int> CstrIndex;//used in enumeration, delete rows

  size_t Mem4Pricing{};
  size_t Mem4Mem{};
  int NumArtiVars{};

  int NumStrongArtiVars{};

  int SeqSizeArtiVars{};

  int MaxNumForwardGraphArc{};

  bool if_stopArcElimination{};

  std::unordered_map<size_t, int> TellWhichBin4ArcEliminationInForwardSense;

  std::unordered_map<std::pair<int, int>, std::vector<std::pair<LABEL *, double>>, PairHasher>
	  concatenateLabelsInBackwardCG;
  int MaxNumBackwardGraphArc{};
  std::vector<std::vector<double>> MainResourceAcrossArcsInBackwardSense;
  std::unordered_map<size_t, int> TellWhichBin4ArcEliminationInBackwardSense;
  double **RC2TillThisBinInBackwardSense{};
  double **RC2BinInBackwardSense{};
  VecLabel **LabelArrayInBackwardSense{};
  VecLabel **IfExistExtraLabelsInBackwardSense{};

  int MaxNumEnuColPool{};

  double MaxGapImprovedInR1C{TOLERANCE};

  double GapTolerance4ArcEliminationNEnumeration{};

  double CutGenTimeThresholdInPricing{};

  bool If_in_Enu_State{};

  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> Map_Edge_ColIdx_in_Enu;
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> RealImprovement_up;// sum and cnt
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> RealImprovement_down;// sum and cnt
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> HeuristicImprovement_up;
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> HeuristicImprovement_down;
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> LPTestingImprovement_up;
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> LPTestingImprovement_down;
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Pair_Val;
  std::vector<std::pair<int, int>> Branch_selected_from_model1;
  std::vector<std::pair<int, int>> Branch_selected_from_model2;
  std::unordered_map<std::pair<int, int>, int, PairHasher> BranchChoice;

#ifdef AccuracyTest
  std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> Map_Edge_SB_rank;
  void generateMap_Edge_SB_rankByExactCG(BBNODE *node);
  void write2Files(int node_dep);
  void tellAccuracy(int node_dep, int phase);
#ifdef Stage1
  std::vector<std::pair<int, int>> Branching_selected_by_model1;
  std::vector<std::pair<int, int>> Branching_selected_by_model2;
#endif
  int rank_phase1{};
  int rank_phase2{};
  double max_SB_scores{};
  std::pair<double, int> average_rank_phase1;
  std::pair<double, int> average_rank_phase2;
  std::pair<double, int> average_ratio_phase1;
  std::pair<double, int> average_ratio_phase2;
#endif

  std::vector<std::vector<double>> MainResourceAcrossArcsInForwardSense;

  int NumBucketsPerVertex{};

  bool if_force_not_regenerate_bucket_graph{};

  double guessed_UB{};

  bool if_enumeration_suc{};

  bool if_findBetterUsingHeurEnumeration_suc{};

  double StepSize{};

  double HardRollBackFactor{};

  double LastMaxTimeLabeling{};

  bool if_populate_DumpIdx4BrCol{};
  std::vector<size_t> DumpIdx4LeftBrCol;//record the size_t idx in colPool
  std::vector<size_t> DumpIdx4RightBrCol;

  //reset cut memory in cut rollback
  std::vector<std::tuple<bool, int, yzzLong>> ResetCutMem;//cut index in r1c & type & memory

  std::unordered_map<std::pair<int, int>, std::vector<std::pair<LABEL *, double>>, PairHasher>
	  concatenateLabelsInForwardCG;

  size_t LabelAssign{};
  size_t RouteInMemAssign{};
  size_t RouteInPricingAssign{};
  int AverRouteLength{};
  bool if_ArcEliminationSucceed{};
  bool if_ArcEliminationTriedButFailed{};

#ifdef test_time
  std::unordered_map<std::string, double> test_time_map;
#endif

  //==0 nothing;
  //==1 return to previous state;
  //==2 resolve with larger mem;
  //==3 cuts tail-off
  int Rollback{};

  bool ForceNotRollback{true};

  double Gap4LastFailDue2DifficultyInArcElimination{1};

  bool enumerationMode{};//0: use info provided by arc elimination; 1: standard enumeration

  double LastEnumerationFailGap{1};

  double TimeResolveLP4Iter{};

  double TimePricing4Iter{};

  bool if_ban_convertDual{};

  double OldUB{};

  bool if_fill_mem{};

  std::unordered_map<yzzLong, std::pair<std::vector<int>, double>> ColSetUsedInRootMIP;

  Eigen::RowVectorXd DualInRootMIP;

  double roundUpTolerance{};

  int CGRollBackTimes{};

  bool if_can_arc_eliminationByExactCG{};

  bool final_decision_4_arc_elimination{};

  bool final_decision_4_enumeration{};

  double GlobalGap{1};

  std::vector<std::vector<std::vector<std::vector<int>>>> record_map_rank1_combinations;

  int Count4Tolerance4tryEnumerationWhenArcEliminationFails{};

  std::unordered_map<int, std::vector<std::vector<int>>> recorded_combinations4R1C;

  //cut size, plan type[i]= (coefficient_vec, denominator, rhs )
  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, int>>> map_rank1_multiplier;

  double NumDominanceChecks{};

  std::unordered_map<yzzLong, std::unordered_set<int>> cut_record;

  void selectCuts(const std::vector<std::pair<std::vector<int>, int>> &tmp_cuts,
				  std::vector<std::pair<std::vector<int>, int>> &cuts,
				  int numCuts);


//  double skip_labels{};

  double OptGap{};

  std::unordered_map<int, std::unordered_set<int>> map_rank1_multiplier_dominance;

#ifdef MASTER_VALVE_ML
  ML3 ml;
#endif

  std::pair<double, int> Ratio_DominanceChecks_NonDominant{};//sum, int
  //enumerator, denominator, rhs
  //####################################################################################################################


  void applyRCF(BBNODE *node, int round, bool if_verbose);

  void writeRCFIP(double ub, const SOLVER &solver) const;

#ifdef if_draw_BBT_graph
  void drawBBTgraph();
  void constructBBNodeName(BBNODE *node, const std::string &terminateReason = "");
#endif

#if defined(find_missed_solution) || defined (find_lost_arcs)
  void findWhySolutionDisappear(BBNODE *node, CVRP *cvrp, std::vector<int> &data, bool &if_opt, bool if_just_check_opt_node=false);
#endif

  void check_ij_bucket(BBNODE *node, int i, int j) const;

  void deleteArcByBrC_false(BBNODE *node, std::pair<int, int> edge) const;

  void regenerateBucketGraph(BBNODE *node);

  void writeMap_Edge_ColIdx_in_Enu(BBNODE *node);

  void doSB_default(BBNODE *node, std::pair<int, int> &info);

  [[nodiscard]] double calculateGapImprovement(double nowVal, double b4Val) const;

#ifdef TestIfPhase0IsWorthy
  void do_SB_openTest(BBNODE *node, std::pair<int, int> &info);
  void do_SB_openTest_stage0(BBNODE *node);
  void do_SB_openTest_stage1(BBNODE *node);
#endif

#ifdef TestIfMLIsWorthy
  void doSB_openMLtest(BBNODE *node, std::pair<int, int> &info);
  void doSB_Stage1_openMLtest(BBNODE *node);
  void doSB_Stage2_openMLtest(BBNODE *node);
#endif

  void doSB_Stage0(BBNODE *node);

  void doSB_Stage1(BBNODE *node);

  void doSB_Stage2(BBNODE *node);

#ifdef OpenHeurCGInUseModel
  void doSB_Stage3(BBNODE *node);
#endif

  void doSB_Stage1_useModel(BBNODE *node);

  void doSB_Stage2_useModel(BBNODE *node);

  void takeTheOneWithMaxColumns(BBNODE *node);

  void constructBranchingSet(BBNODE *node);

  void writeISAcc(const std::vector<std::pair<int, int>> &tmp_br);

  void writeUseModelAcc(const std::vector<std::pair<int, int>> &tmp_br, int top_n);

  void generateModelInPhase1_combinedFashion(BBNODE *node);

  void testInitialScreening(BBNODE *node);

  void InitialScreening(BBNODE *node,
						bool if_record_source,
						bool if_fill_candidates_gap,
						int num,
						double pseudo_frac);

  void readSolFile(bool if_force);

  void reviseBranch_pair();

  void updateUB_EdgeSol();

  void simulateWriteLPPseudocost(BBNODE *node);

  std::unordered_set<std::pair<int, int>, PairHasher> EdgeSol;

  void LPTesting(BBNODE *node, int num, bool if_writeBranch_pair, bool if_record_sb_scores = false,
				 bool if_record_LP_improvement = true);

#ifdef MASTER_VALVE_ML

  void useModelInPhase1(BBNODE *node, int num, bool if_distinguish_source);// to keep compatible to the old version

  void useModelInPhase1(BBNODE *node, int num, int if_force_frac = Combined_Force_Frac);

  void getTrainingDataInPhase1(BBNODE *node);

#endif

  void randomPickLP(BBNODE *node, int num);

  void initialPickLP(BBNODE *node, int num);

  void useModelInPhase1_test_one_model(BBNODE *node, int num);

  void addSearch_crazy(int plan_idx,
					   const std::vector<int> &c,
					   const std::vector<int> &w_no_c,
					   double &new_vio,
					   int &add_j,
					   const std::vector<std::unordered_map<int, int>> &v_r_map,
					   const std::vector<double> &frac_routes,
					   std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void construct_cuts_crazy(
	  const std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio,
	  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
	  std::vector<std::pair<std::vector<int>, int>> &cuts
  );

  void removeSearch_crazy(int plan_idx,
						  const std::vector<int> &c,
						  double &new_vio,
						  int &remove_j,
						  const std::vector<std::unordered_map<int, int>> &v_r_map,
						  const std::vector<double> &frac_routes,
						  std::unordered_map<yzzLong,
											 std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void swapSearch_crazy(int plan_idx,
						const std::vector<int> &c,
						const std::vector<int> &w_no_c,
						double &new_vio,
						std::pair<int, int> &swap_i_j,
						const std::vector<std::unordered_map<int, int>> &v_r_map,
						const std::vector<double> &frac_routes,
						std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void operations_crazy(
	  rank1_multi_label &label,
	  const std::vector<std::unordered_map<int, int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  int &i,
	  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
	  std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio
  );

  void start_seed_crazy(
	  const std::vector<std::unordered_map<int, int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  const std::vector<std::pair<std::vector<int>, std::vector<int>>> &c_N_noC,
	  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
	  std::vector<rank1_multi_label> &rank1_multi_label_pool,
	  std::unordered_map<
		  yzzLong, std::vector<
			  std::pair<std::vector<int>, double>>> &map_cut_plan_vio,
	  int &num_label
  );

  void addR1C_atOnceInEnu(BBNODE *node,
						  const std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generaten4colpoolMapInEnu(BBNODE *node,
								 std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> &colmap);

  void addBrC2ColPoolInEnuBycolmap(BBNODE *node,
								   const std::pair<int, int> &edge);

  void writeSB_scores2Files(int node_dep);

  void changeModel4BetterDual(BBNODE *node);

  void CGTesting(BBNODE *node,
				 bool if_exact_CG,
				 bool if_record_product_val,
				 bool if_record_improvement,
				 bool if_force_complete = false);

  void exactCGTesting(BBNODE *node);

  void ExperimentalMode(BBNODE *node, bool if_use_initial, bool if_use_sep);

  void useModelInPhase1_sep(BBNODE *node, int num);

  void calculate_acc(const std::vector<std::pair<int, int>> &top_n,
					 const std::vector<std::pair<int, int>> &branch_pair,
					 std::unordered_map<int, std::pair<int, double>> &top_n_per);

  std::unordered_map<int, std::pair<int, double>> top_n_percentage_stage1;
  std::vector<std::unordered_map<int, std::pair<int, double>>> top_n_percentage_stage2_takeout_n;

  std::unordered_map<int, int> map_level_numCandidates;

  double Ratio_pseudo{};

  void deleteColTesting_Enu_specialize(BBNODE *node, bool if_record_product_val);

#ifdef MASTER_VALVE_ML
  void useModelInPhase2(BBNODE *node, int num);

  void getTrainingDataInPhase2(BBNODE *node);

  void veryLightHeuristicCG4ML(BBNODE *node,
							   double org_val,
							   int BeforeNumRow,
							   const std::pair<int, int> &edge,
							   bool dir,
							   int col_start,
							   int iter);

  void testuseML(BBNODE *node);

#endif

  void onlyUseStageOneModel(BBNODE *node, bool if_sep);

#ifdef  useM_dynamically
  void useMLVaryLPTesting(BBNODE *node);
  void giveDynamicM(BBNODE *node, int &num);
  void evaluateM1();
  std::vector<std::pair<int, int>> cp_Branch_pair;
  std::vector<std::pair<int, int>> tested_Branch_pair;
  double m{};
  double F{};
  double k{};
  bool if_use_full_k{};
#endif

  void useMLInGeneralFramework(BBNODE *node, bool if_sep);

  void do_SB(BBNODE *node, std::pair<int, int> &info);

  void search4ParameterInPhase1(BBNODE *node);

  void generateModelInPhase2(BBNODE *node, bool if_sep);

  void SpecializedCGTesting4Enu(BBNODE *node);

  void useSpecializedBrInEnu(BBNODE *node);

  void fill_mem_first(BBNODE *node,
					  const std::vector<std::vector<int>> &routes,
					  const std::vector<double> &frac_routes,
					  std::vector<std::pair<std::vector<int>, int>> &cuts);

  void useDefaultSB(BBNODE *node);

  void findAccInitial(BBNODE *node);

  void findAccML(BBNODE *node);

  void DeleteDumpCols(BBNODE *node);

  void generateModelInPhase1(BBNODE *node, bool if_sep);

  void generateEnuModelInStageTwo(BBNODE *node);

  void doSB_enumeration(BBNODE *node, std::pair<int, int> &info);

  void findNGMemSets(BBNODE *node, bool &if_empty);

  void deleteColByNGMem(BBNODE *node);

  void augmentNGRound(BBNODE *node);

  void chgEnuMatByCuts(BBNODE *node);

  void deleteNewAddedActiveCutsByDual_N_Mem(BBNODE *node, int oldNum);

  void addR1C_multiInEnu(BBNODE *node, const std::pair<std::vector<int>, int> &cut);

  void addR1C_multi(BBNODE *node, const std::pair<std::vector<int>, int> &cut, const std::set<int> &mem, int cut_index);

  void findMem4Rank1_multi(const std::vector<std::vector<int>> &routes,
						   const std::vector<std::unordered_map<int, int>> &v_r_map,
						   const std::pair<std::vector<int>, int> &cut_pair,
						   std::set<int> &mem,
						   bool &if_suc);

  void printCutsInfo(BBNODE *node) const;

  void findPlan4Rank1_multi(const std::vector<int> &vis, int denominator, yzzLong &mem,
							std::vector<std::set<int>> &segment, std::vector<std::vector<int>> &plan);

  void combinationUtil_addOne(const std::vector<int> &arr,
							  std::vector<int> &tmp,
							  std::vector<std::vector<int>> &data,
							  int start,
							  int end,
							  int index,
							  int r);
  void findCombs4rankCutSet2getBestMultiplier(int plan_idx,
											  const std::vector<int> &cset,
											  std::vector<int> &new_cset,//the c will be sorted!
											  int c,
											  double &new_vio,
											  const std::vector<std::vector<int>> &v_r_map,
											  const std::vector<double> &frac_routes,
											  bool &if_succeed);
  void pure_map(std::vector<std::vector<std::vector<int>>> &map, const std::vector<int> &sorted_info);

  void comb_all(int n, int r, std::vector<std::vector<int>> &data);

  void combinationUtilOuter(std::vector<std::vector<int>> &data, int n);

  void addInRank1HeurSep(int plan_idx,
						 const std::vector<int> &c,
						 const std::vector<int> &w_no_c,
						 double &new_vio,
						 int &choice_j,
						 const std::vector<std::vector<int>> &v_r_map,
						 const std::vector<double> &frac_routes);

  void removeInRank1HeurSep(int plan_idx,
							const std::vector<int> &c,
							double &new_vio,
							int &choice_j,
							const std::vector<std::vector<int>> &v_r_map,
							const std::vector<double> &frac_routes);

  void generatePermutations(std::unordered_map<int, int> &count_map,
							std::vector<int> &result,
							std::vector<std::vector<int>> &results,
							int remaining);

  void exactFindBestPermuatation4Oneplan(const std::vector<std::unordered_map<int, int>> &v_r_map,
										 const std::vector<double> &frac_routes,
										 const std::vector<int> &cut,
										 int plan_idx,
										 double &vio,
										 std::unordered_map<
											 yzzLong, std::vector<
												 std::pair<std::vector<int>, double>>> &map_cut_plan_vio);

  void swapInRank1HeurSep(int plan_idx,
						  const std::vector<int> &c,
						  const std::vector<int> &w_no_c,
						  double &new_vio,
						  std::pair<int, int> &choice_i_j,
						  const std::vector<std::vector<int>> &v_r_map,
						  const std::vector<double> &frac_routes);

  void construct_operations_aggressive(
	  rank1_multi_label &label,
	  const std::vector<std::vector<int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  int &i,
	  std::unordered_map<int,
						 std::vector<std::tuple<std::vector<int>,
												int,
												double>>> &generated_rank1_multi_pool);

  void checkR1Ctotally(BBNODE *node);

  void construct_cuts_aggressive(std::unordered_map<int,
													std::vector<std::tuple<std::vector<int>,
																		   int,
																		   double>>> &generated_rank1_multi_pool,
								 std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateR1C3s_aggressive(const std::vector<std::vector<int>> &ele_routes,
								const std::vector<std::vector<int>> &non_ele_routes,
								const std::vector<double> &frac_ele_routes,
								const std::vector<double> &frac_non_ele_routes,
								std::vector<std::pair<std::vector<int>, int>> &cuts) const;

//  void construct_v_r_vec_aggressive(
//      const std::vector<std::vector<int>> &routes,
//      std::vector<std::vector<int>> &v_r_map,
//      std::vector<std::unordered_map<int, int>> &v_r_map2,
//      std::vector<std::vector<std::vector<int>>> &c,
//      std::vector<std::vector<std::vector<int>>> &w_no_c,
//      std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map);


  void construct_v_r_vec_aggressive(const std::vector<std::vector<int>> &routes,
									std::vector<std::vector<int>> &v_r_map,
									std::vector<std::unordered_map<int, int>> &v_r_map2,
									std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> &c_n_w_no_c,
									std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map);

  void findR1C_multi_aggressive(const std::vector<std::vector<int>> &routes,
								const std::vector<double> &frac_routes,
								std::vector<std::unordered_map<int, int>> &v_r_map,
								std::vector<std::pair<std::vector<int>, int>> &cuts);

  void construct_seed_aggressive(
	  const std::vector<std::vector<int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, double>>> &generated_rank1_multi_pool,
	  std::vector<rank1_multi_label> &rank1_multi_label_pool,
	  std::unordered_map<int, std::vector<std::pair<int, int>>> &c_map,
	  std::vector<std::vector<std::pair<std::vector<int>, std::vector<int>>>> &c_n_w_no_c,
	  int &num_label);

  void generateOptimalMultiplier4R1C_multi();

  void useMLInStageOne(BBNODE *node);

  virtual void getLowerBoundofMinimumNumCars();

  void findR1C_multi(const std::vector<std::vector<int>> &routes,
					 const std::vector<double> &frac_routes,
					 const std::vector<std::vector<int>> &vbasis,
					 std::vector<std::pair<std::vector<int>, double>> &cut_info,
					 std::vector<std::tuple<std::vector<int>, int, double>> &multi_cut_info);

  void newfindR1C_multi(const std::vector<std::vector<int>> &routes,
						const std::vector<double> &frac_routes,
						std::vector<std::unordered_map<int, int>> &v_r_map,
						std::vector<std::pair<std::vector<int>, int>> &cuts);

  void crazySearch(const std::vector<std::vector<int>> &routes,
				   const std::vector<double> &frac_routes,
				   std::vector<std::unordered_map<int, int>> &v_r_map,
				   std::vector<std::pair<std::vector<int>, int>> &cuts);

  void construct_v_r_map_n_seed_crazy(const std::vector<std::vector<int>> &routes,
									  std::vector<std::unordered_map<int, int>> &v_r_map,
									  std::vector<std::pair<std::vector<int>,
															std::vector<int>>> &c_N_noC);

  void construct_v_r_vec(const std::vector<std::vector<int>> &routes,
						 std::vector<std::vector<int>> &v_r_vec,
						 std::vector<std::unordered_map<int, int>> &v_r_map) const;

  void construct_cuts(const std::vector<std::unordered_map<int, int>> &v_r_map,
					  const std::vector<std::vector<int>> &seed,
					  const std::vector<double> &frac_routes,
					  std::vector<std::pair<std::vector<int>, int>> &cuts// cut self, plan idx
  );

  void construct_cuts_by_resourceGap_only(const std::vector<std::unordered_map<int, int>> &v_r_map,
										  const std::vector<std::vector<int>> &seed,
										  const std::vector<double> &frac_routes,
										  std::vector<std::pair<std::vector<int>, int>> &cuts// cut self, plan idx
  );

  void construct_mem(BBNODE *node, const std::vector<std::vector<int>> &routes,
					 const std::vector<std::unordered_map<int, int>> &v_r_map,
					 const std::vector<std::pair<std::vector<int>, int>> &cuts,
					 std::vector<std::tuple<std::vector<int>, int, int, std::set<int>>> &full_cuts);

  template<bool if_symmetry>
  void construct_v_resourceGap(const std::vector<std::vector<int>> &routes, std::vector<double> &v_resourceGap);

  void construct_seed(const std::vector<std::vector<int>> &routes,
					  const std::vector<std::vector<int>> &v_r_vec,
					  const std::vector<double> &v_resourceGap,
					  std::vector<std::vector<int>> &seed);

  static void checkR1C_mul_look(BBNODE *node, int start);

  template<bool dir, bool if_symmetry>
  int enumerateHalfwardRoutes(BBNODE *node,
							  const double *r1c_to_pi,
							  const double *r1c_multi_to_pi,
							  std::unordered_map<yzzLong, std::tuple<LABEL *, LABEL *, double>> &Tags,
							  std::vector<LABEL *> **copy_bucket,
							  int &num_routes_now);

  void generateAllR1Cs(BBNODE *node, bool if_use_MIP);

  void NewgenerateAllR1Cs(BBNODE *node, int aggressive_level);

  void addR1C(BBNODE *node, const std::vector<int> &cut, const std::set<int> &mem, int cut_index);

  void addR1CInEnu(BBNODE *node, const std::vector<int> &cut);

  void heuristicFindBestPermuatation4Oneplan(
	  const std::vector<std::unordered_map<int, int>> &v_r_map,
	  const std::vector<double> &frac_routes,
	  const std::vector<int> &cut, int plan_idx,
	  std::vector<int> &new_cut,
	  double &vio);

  void addR1C_atOnce(BBNODE *node,
					 const std::vector<std::tuple<std::vector<int>, int, int, std::set<int>>> &full_cuts);

  void runCertainNumberOfHeurCGs(BBNODE *node, int iteration);

  void generateAllR1CsInEnu(BBNODE *node, bool if_use_MIP);

  void findBetterUsingHeurEnumeration(BBNODE *&old_node);

  void refineArcsHeur(BBNODE *node) const;

  bool if_use_heur_enumeration{};

#ifdef PrintOptimalSolution

  void printOptIntSol();

#endif

#ifdef VERBOSE

  void printBrDecisions();

  static void printInfoLabeling(int iter, int num_added_col, int num_col,
								int num_row,
								double mt, double spt, double et,
								double lp, double lb_val, double ub);

#endif

  explicit CVRP(const InstanceData &instanceData);

  virtual ~CVRP();

  void buildModel();

  bool addBrCut2Unsolved(BBNODE *node, const std::pair<int, int> &info);

  void sepHybridCuts(BBNODE *&node);

  void findMem4R1CsMode3ByEnumeration_N_MIP(yzzLong v_comb,
											std::set<int> &mem,
											const std::vector<std::vector<int>> &routes,
											bool &if_suc);

  void getMemByMIP(const std::vector<std::vector<std::vector<int>>> &array,
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
  void sepRCCs(BBNODE *&node);

  void generateRCCs(BBNODE *node);

  void randomSelection(BBNODE *node);

  void solveLPInLabeling(BBNODE *node, bool if_open_heur = true, bool if_open_exact = true, bool if_record_sol = true);

  void priceLabeling(BBNODE *node);

  void rollbackEaseWay(BBNODE *node, int old_num);

  void assignMem();

  bool solveSBModel();

  void recordOptCol(BBNODE *node, bool if_force_rewrite = false);

  int optimizeLP4OneIter(BBNODE *node, double prior_value);

  void reformulateIPSol(BBNODE *node);

  void construct_seed_mode3(const std::vector<std::vector<int>> &routes,
							const std::vector<std::vector<int>> &v_r_vec,
							const std::vector<double> &v_resourceGap,
							std::vector<std::vector<int>> &seed);

  void NewfindR1C_multi_mode3(const std::vector<std::vector<int>> &routes,
							  const std::vector<double> &frac_routes,
							  std::vector<std::unordered_map<int, int>> &v_r_map,
							  std::vector<std::pair<std::vector<int>, int>> &cuts);

  //Refers to the scope from which to start the reduction
  //In general, the reduction starts from one's own generation.
  //When the col exceeds a certain upper limit, the reduction starts from fx_dimension.
  void cleanIdxCol4Node(BBNODE *node, int beg, bool if_only_rcfixing = false);

  [[nodiscard]] double calculateDif(double tmp_val, double prior_val, bool if_chg = false) const;

  void writeCol2Mem(BBNODE *node);

  void constructMap(BBNODE *node, int beg) const;

  void getNewCstrCoeffByEdge(BBNODE *node, std::pair<int, int> edge, int *c_ind, double *c_val, int &num_nz);

  void getCoeffRCC(BBNODE *node, RCC &rcc, int *c_ind, double *c_val, int &num_nz) const;

  void convertVertex2R1CsInPricing(BBNODE *node);

  void convertVertex2R1CsInOneLP(BBNODE *node);

  //write mem in one pricing
  //if return -1, this node is infeasible
  //if return 0, this pricing is finished
  template<bool if_symmetry>
  int generateColsByBidir(BBNODE *node);

  void addCols(BBNODE *node, int &ccnt);//直接从pri一直遍历到curr吧

  int generateColsByLighterHeur(BBNODE *node);

  int generateColsByHeavierHeur(BBNODE *node);

  double calculateOptGap(BBNODE *node) const;

  void eliminateArcs(BBNODE *node);

  template<bool dir>
  void populateRC2TillThisBinNRC2Bin(BBNODE *node) const;

  void initializeLabels(BBNODE *node, int mode,
						bool if_resetLabelPoint, std::tuple<bool, int, bool> control_cleanAllPtr);
  void deleteNonActiveCutsBySlack(BBNODE *node, bool if_rcc, bool if_r1c);

  void terminateByMIP(BBNODE *node);

  void enumerateMIP(BBNODE *&node);

  void cleanAllPtr(BBNODE *node, int mode, bool if_clear_concatenate);

  void cleanColsNonEle(BBNODE *node);

  void cleanColsRCLargerThanOptGapInLP(BBNODE *node);

  void organizeColsInMem2Pricing(BBNODE *node);

  void solveLPByInspection(BBNODE *node, bool if_only_need_value,
						   bool if_heuristic, bool if_record_sol);

  void deleteBrCsNR1C1s(BBNODE *node);

//  int generateColsByInspection(BBNODE *node,
//                               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &another_mat,
//                               Eigen::RowVectorXd &another_cost,
//                               Eigen::VectorXi &another_ptr,
//                               Eigen::RowVectorXd &rc,
//                               bool if_only_need_value,
//                               bool if_heuristic);

  int generateColsByInspection(BBNODE *node,
							   bool if_only_need_value);

  void addColsByInspection(BBNODE *node, const std::vector<int> &Col_added);

  void regenerateEnuMat(BBNODE *node, BBNODE *node2, bool if_force = false);

  void generateVertex2IdxCols_N_Edge2IdxCols(BBNODE *node);

  void recoverR1CsInEnu(BBNODE *node);

  void reallocateSolverPtr(size_t num);

  [[nodiscard]] int checkSolverPtr(size_t num) const;

  void optimizeLP4OneIterInEnu(BBNODE *node);

  void recordOptColInEnu(BBNODE *node);

  void terminateNode(BBNODE *&root_node);

  void sepHybridCutsInEnu(BBNODE *node);

  void addBrCut2UnsolvedInEnu(BBNODE *node, const std::pair<int, int> &info);

  void reviseEnuColInfoByBrC(BBNODE *node, BBNODE *out_node, const BrC &bf);

  void cleanColsInLPByNewBrC(BBNODE *node, const BrC &bf);

  void generateRCCsInEnu(BBNODE *node);

  void generateR1C3sInEnu(BBNODE *node, int &num_r1c3);

  void deleteNonActiveCutsInEnu(BBNODE *node, int old_num);

  void deleteNonActiveCutsSafely(BBNODE *node, int old_num);

  void writeColsInPricingPool(BBNODE *node, int &index);

  void tellIfArcElimination(BBNODE *node);

  template<bool dir, bool if_last_half, bool if_symmetry, bool if_std_optgap, bool if_res_updated>
  void updateLabel(double res, LABEL *&ki, int i, int j, int &bj,
				   const double *r1c_to_pi,
				   const double *r1c_multi_to_pi,
				   int min_sort_b,
				   bool &if_suc);

  template<bool dir>
  void doDominance(LABEL *&ki, int j, int bj, const double *r1c_to_pi,
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
  void concatenateOneLabelWithOtherLabels(LABEL *ki, int j, int arr_bj, double tmp_rc, double tmp_mainResource,
										  const double *r1c_to_pi,
										  const double *r1c_multi_to_pi,
										  int &if_state);

  template<bool dir, bool if_symmetry>
  void eliminatebuketArcs(BBNODE *node, const double *r1c_to_pi,
						  const double *r1c_multi_to_pi,
						  int dim_sq,
						  bool *stateBetween2Buckets,
						  int *latest_bucket);

  template<bool dir, bool if_symmetry>
  void eliminatebuketArc4Depot(BBNODE *node);

  template<bool dir, bool if_symmetry>
  void concatenatePhaseInArcElimination(const double *r1c_to_pi,
										const double *r1c_multi_to_pi);

  void addPathByRC(double path_rc, LABEL *ki, LABEL *kj, int num);

  void tellIfEnumeration(BBNODE *node);

  template<bool dir>
  void populateTellWhichBin4ArcElimination();

  void getInfoEdge(BBNODE *node, bool if_br) const;

  void deleteNonActiveCutsByDual(BBNODE *node, bool if_rcc_by_slack);

  void deleteNewAddedNonActiveCutsBySlack(BBNODE *node, int olNum, bool if_keep_rcc);

  void deleteNewAddedNonActiveCutsBySlackInEnu(BBNODE *node, int olNum, bool if_delete_rcc);

  void initialBucketGraph4Node(BBNODE *node);

  //Mem section
  [[nodiscard]] int checkMaxNum_Customers() const;

  [[nodiscard]] static int checkMaxNum_R1Cs(int num_r1cs);

  [[nodiscard]] int checkCST_LIMIT() const;

  [[nodiscard]] static int checkAssignGreaterMAXINT(const size_t &value, const std::string &str);

  void reallocateLabel();

  void reallocateMemPool();

  [[nodiscard]]  int checkMemPool() const;

  void reallocatePricingPool();

  [[nodiscard]]  int checkPricingPool() const;

  void rollback2PreState(BBNODE *node, const std::vector<RCC> &old_rcc,
						 const std::vector<R1C> &old_r1c,
						 const std::vector<R1C_multi> &old_r1c_multi,
						 const std::vector<BrC> &old_brc, const std::vector<size_t> &old_col_idx);

  static int inverseLastBranchconstr(char sense, double rhs, SOLVER &solver);

  static int addBranchconstr(int numnz,
							 int *cind,
							 double *cval,
							 char sense,
							 double rhs,
							 const char *constrname, SOLVER &solver);

  void chgBranchconstr(double *val,
					   int *cind,
					   int *vind,
					   int numnz,
					   int *old_ind,
					   double *old_val,
					   char sense,
					   double rhs,
					   SOLVER &solver);

  void initialBucketGraph();

  void checkIfCutsLegal(BBNODE *node) const;

  void lateProcessing();

  virtual double transformCost(double x);

  void solveMIP(BBNODE *node, bool if_inEnu);

  void giveInitialUBByMIP(BBNODE *&node);

  void copyColGeneratedInRoot4MIP();

  void collectDualInRoot4MIP(BBNODE *node);
  [[nodiscard]] double ceil_transformed_number_related(double x) const;

  void initialLabels();

  void writeAllColsInIP4Check(BBNODE *node) const;

  template<typename T, bool dir, bool if_last_half, bool if_symmetry>
  int extendKernel4Exact(LABEL *&ki,
						 int i,
						 double res,
						 const std::vector<T> &arc,
						 const double *r1c_to_pi,
						 const double *r1c_multi_to_pi,
						 int min_sort_b);

  template<bool dir>
  void obtainjumpArcs(BBNODE *node, std::bitset<2> **bitMap) const;

  template<typename T>
  int extendKernel4Exact_Backward(LABEL *&ki,
								  int i,
								  double res,
								  const std::vector<T> &arc,
								  const double *r1c_to_pi,
								  const double *r1c_multi_to_pi);

  void runHalfBackwardLabeling(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  int concatCols_NONE_Symmetry(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  int generateColsByBidir_NONE_Symmetry(BBNODE *node);

  bool decreaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end);
  template<typename T>
  int extendKernel4ArcElimination_last_half_Backward(LABEL *&ki,
													 int i,
													 double res,
													 const std::vector<T> &arc,
													 const double *r1c_to_pi,
													 const double *r1c_multi_to_pi);

  template<bool dir>
  void checkIfDominated(LABEL *&ki, int i, int b, const double *r1c_to_pi,
						const double *r1c_multi_to_pi,
						bool &if_suc);

  template<bool dir>
  void checkIfNoLabelsLeft(int &i, int b, int &min_sorted_b, bool &if_suc);

  template<bool dir, bool if_last_half, bool if_symmetry>
  void runLabeling(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  int extendKernel4ArcElimination_inner_Backward(LABEL *&ki,
												 int i, int j, double &tmp_mainResource,
												 const double *r1c_to_pi,
												 const double *r1c_multi_to_pi);

  int lastHalfBackwardInArcElimination(BBNODE *node, const double *r1c_to_pi,
									   const double *r1c_multi_to_pi);

  int backwardConcatenateInArcElimination(const double *r1c_to_pi,
										  const double *r1c_multi_to_pi);

  void eliminateBackwardArcs(BBNODE *node, const double *r1c_to_pi,
							 const double *r1c_multi_to_pi,
							 int dim_sq,
							 bool *stateBetween2Buckets,
							 int *latest_bucket);

  void obtainBackwardJumpArcs(BBNODE *node, std::bitset<2> **bitMap) const;

  int enumerateHalfBackwardRoutes(BBNODE *node,
								  const double *r1c_to_pi,
								  const double *r1c_multi_to_pi,
								  std::vector<LABEL *> **copy_bucket);

  void setTailOffStd_N_RollBackStd() const;

  int concatenateRoutes_prior_forward_InEnumeration(BBNODE *node,
													const double *r1c_to_pi,
													const double *r1c_multi_to_pi,
													std::unordered_map<yzzLong,
																	   std::tuple<LABEL *, LABEL *, double>> &Tags,
													int &num_routes_now
  );
  template<typename T>
  void extendKernel4LightHeur(LABEL *&ki,
							  int i,
							  double res,
							  const std::vector<T> &arc,
							  const double *r1c_to_pi,
							  const double *r1c_multi_to_pi,
							  bool if_return_to_depot);
  template<typename T>
  void extendKernel4HeavierHeur(LABEL *&ki,
								int i,
								double res,
								const std::vector<T> &arc,
								const double *r1c_to_pi,
								const double *r1c_multi_to_pi,
								bool if_return_to_depot);

  int lastHalfForwardInArcElimination(BBNODE *node,
									  const double *r1c_to_pi,
									  const double *r1c_multi_to_pi);

  int forwardConcatenateInArcElimination(const double *r1c_to_pi, const double *r1c_multi_to_pi);

  void eliminatebuketArcs(BBNODE *node,
						  const double *r1c_to_pi,
						  const double *r1c_multi_to_pi,
						  int dim_sq,
						  bool *stateBetween2Buckets,
						  int *latest_bucket);

  template<typename T>
  int extendKernel4ArcElimination_last_half_Forward(LABEL *&ki,
													int i,
													double res,
													const std::vector<T> &arc,
													const double *r1c_to_pi,
													const double *r1c_multi_to_pi);
  int extendKernel4ArcElimination_inner_Forward(LABEL *&ki,
												int i, int j, double &tmp_mainResource,
												const double *r1c_to_pi,
												const double *r1c_multi_to_pi);

//  void addRCC2LP(BBNODE *node, std::unordered_map<int, std::vector<int>> &vertex_map, RowMatrixXd &mat) const;

  bool runColumnGenerationType(BBNODE *node, int mode);

//  void changeModel4DeletingSlackDual_zeroCuts(BBNODE *node);

  void findNonactiveCuts(BBNODE *node);

  void deleteNonactiveCuts(BBNODE *node, std::vector<int> &nonactive_cuts);

//  bool if_upsideDownDual_in_enumeration{};

//  std::vector<int> target_cuts;

#ifdef DEBUG_MIP_FILE
  void writeMIP(SOLVER &solver, const std::string &file_type) const;
#endif

#ifdef DEBUG_LP_FILE
  void writeLP_N_tell_if_LP_corrected(SOLVER &solver) const;
#endif

  int runLighterHeurLabeling(BBNODE *node);

  int runHeavierHeurLabeling(BBNODE *node);

  void runHalfForwardLabeling(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  bool increaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end);

  virtual void specialize_MaxLengthEleRoute();

  template<bool if_symmetry>
  int concatenateCols_prior_forward(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  void runLabeling4ArcElimination(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  void eliminateBucketArcs(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  void obtainJumpArcs(BBNODE *node) const;

  bool enumerateRoutes(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs);

  virtual void setResourceInBucketGraph();

  void generateR1C3s(const std::vector<std::vector<int>> &ele_routes,
					 const std::vector<std::vector<int>> &non_ele_routes,
					 const std::vector<double> &frac_ele_routes,
					 const std::vector<double> &frac_non_ele_routes,
					 std::vector<std::pair<std::vector<int>, double>> &cut_info) const;

  void NewgenerateR1C3s(const std::vector<std::vector<int>> &routes,
						const std::vector<double> &frac_routes,
						std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateR1C1s(const std::vector<std::vector<int>> &non_ele_routes,
					 const std::vector<double> &frac_none_ele_routes,
					 std::vector<std::pair<std::vector<int>, double>> &cut_info) const;

  void NewgenerateR1C1s(const std::vector<std::vector<int>> &non_ele_routes,
						const std::vector<double> &frac_none_ele_routes,
						std::vector<std::pair<std::vector<int>, int>> &cuts);

  void generateHighDimR1Cs(const std::vector<std::vector<int>> &routes,
						   const std::vector<double> &frac_routes,
						   const std::vector<int> &cut_type,
						   std::vector<std::pair<std::vector<int>, double>> &cut_info);

  template<typename T, bool if_enu>
  void addSelectedR1C_N_multiCuts(BBNODE *node,
								  std::vector<std::tuple<int, std::set<int>, double, int, int>> &cut_info_set,
	  //cut_idx, mem, vio, lp_index(if max, then new cut!), deducted_memSize
								  T &cut_info);
  //DEBUG
#ifdef DEBUG
  void checksol_(BBNODE *node) const;
#endif

#ifdef debugVio
  std::unordered_map<yzzLong, std::vector<double>> vio_map;
#endif

#ifdef writeEnumerationTrees
  std::unordered_map<yzzLong, std::tuple<std::vector<int>, double, size_t>> enumeration_col_idx;// sequence, cost, index
  void writeEnuTree(BBNODE *node);
  void populateEnuCols(std::vector<size_t> &col_idx, const std::vector<std::pair<size_t, double>> &col_infos);
  void writeEnuCols();
  void writeEnuCuts(BBNODE *node, const std::vector<size_t> &col_idx);
#endif

  bool if_only_read_enumerationTree{};
#ifdef readEnumerationTrees
  std::string tree_path{};
  std::string colPool_path{};
  void solveEnuTree();
  void readEnuTree(BBNODE *node, std::vector<size_t> &col_idx);
  void readColPool(BBNODE *node, const std::vector<size_t> &col_idx, std::vector<double> &objs);
  void restoreModel(BBNODE *&node);
  void deleteBrCs_spec4readEnuTree(BBNODE *const node);
#endif
};

bool operator==(const RCC &lhs, const RCC &rhs);

auto CmpLabelRCLess(const LABEL *l1, const LABEL *l2) -> bool;

std::string readInstanceFile(const std::string &file_name, int line);

std::string generateInstancePath(int argc, char *argv[]);

std::string doubleToString(double d);

float sqrt_self(float x);

double pow_self(double x, int n);

void self_mkdir(const std::string &path);
#endif
