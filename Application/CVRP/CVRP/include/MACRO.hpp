//
// Created by Zhengzhong You on 5/31/22.
//

#ifndef CVRP_MACRO_HPP
#define CVRP_MACRO_HPP

#include "FindInsByCertainSetting.hpp"
#include "Experiment.hpp"
#include "techniqueControl.hpp"

//#define writeEnumerationTrees
//#define readEnumerationTrees

#if defined(writeEnumerationTrees) || defined(readEnumerationTrees)
#define tree_folder std::string("EnuTree")
#define col_pool_folder std::string("EnuColPool")
#endif

#define branch_mode 2
/**
 * 1. default 100, 3; 15, 1
 * 2. ml (100), 15, 2; (100), 15, 2
 * 3. 3pb 100, 2, 100, 2
 */

#define do_SB_mode 1
/**
 * 1. normal
 * 2. mode 1-4 (IS/ Sep?)
 * 3. test initialS is sufficient, please set ML_state == 1, test overall Acc, please set ML_state == 3, and set testAcc
 * 4. please use if_sep instead!!! test why separating data is good (generate data 1 only) please set ML_state == 1,
 */

#define MASTER_VALVE_ML

#ifdef MASTER_VALVE_ML
#define NoDifEnu_b4_af

/**
 * testing module: in default situation, please shut all down
 */
//#define getOneStageModel 1 // 1 means get one stage model(use 2 in ML_state), 2 means read one stage model (use 3 in ML_state), please incorporate with ML_state

//#define NoInitialScreening
//#define checkSimilar 2

#define ML_state 3

//#define testAcc
/*
 * 1. get data 1
 * 2. get data 2 // also used for generating one-stage model, please change the tree level
 * 3. use model
 * 4. only use phase1 model
 * 5. use it in 3 but with \hat \theta changing dynamically
 */


#if defined(ML_state) && ML_state == 5
#define useM_dynamically
#endif

#endif

#if defined(MASTER_VALVE_ML) && (ML_state == 1 || ML_state == 2)
#define TrainingDataTreeLevel 8
#endif

#define Combined_Force_Frac 1 //please set 1 or 0 // 1 means forcing to use it like a sep mode, (not really mean it is sep mode)

//#define changeModel
//#define changeModel_mustUse
//#define debugVio
//#define checkR1C_Mul_State
//#define if_force_read_sol
//#define if_use_sol_when_br
//#define random_selection // just for testing
//#define checkR1CTotally //this checks the R1C by comparing the coefficients
//#define onlyExactCG// just for benchmarking
#ifdef onlyExactCG
#define if_tell_model_works//right now only for two models
#ifdef if_tell_model_works
#define writeAccOut2File "acc"
#endif
#define if_tell_LPTesting_works
#ifdef if_tell_LPTesting_works
#define writeAccOut2File "acc2"
#endif
#if defined(if_tell_model_works) && defined(if_tell_LPTesting_works)
#error "if_tell_model_works and if_tell_LPTesting_works cannot be defined at the same time"
#endif
#endif
#if defined(MASTER_VALVE_ML) && defined(random_selection)
#error "random_selection and MASTER_VALVE_ML cannot be defined at the same time"
#endif

//#define test_time
//#define find_missed_solution
#ifdef find_missed_solution
#define traverse_required_nodes// only calculate the optimal nodes
//always incorporate with indicate_search_tree to identify the branch choice!
#endif
//#define check_if_cuts_are_legal
//#define indicate_search_tree
//#define find_lost_arcs// should indicate i and bin
//#define check_lower_bound
//#define check_enumeration_pool_unsatisfied_cols_by_br
//#define if_draw_BBT_graph
//#define STD_BRANCHING
//#define if_constraint_on_small_ins
//#define draw_residual_bucket_graph
#define writeSolOut2File "sol"

//#define SYMMETRY_PROHIBIT
//#define SOLVER_VRPTW

#define CAPACITY_AS_MAIN_RESOURCE

//#define Use_heuristic_UB

//#ifdef Use_heuristic_UB
//#define Resolve_Ins_with_Optimal_Val (not well implemented yet!)
//#endif

//#define DETAILED_EXACT_PRINT_INFO
//#define DETAILED_RANK1_PRINT_INFO
//#define HGS_APPLIED
//#define DELUXING_APPLIED
#ifdef DELUXING_APPLIED
#define DELUXING_ROUND 10
//#define writeIP "IP"
#endif
#define DUAL_TOLERANCE 1e-6

#define NewWay2AddCuts
#define PRINT_LABELING_STEP_SIZE 20
#define LARGEFLOAT 3e30
#define MaxInt 2000000000
#define VERBOSE
#define ReadNoLine (-1)
#define GlobalTimeLimit 1000000000
#define NUM_THREAD_LP 1
#define PrintOptimalSolution
#define COL_KEPT_FRAC 0.67
#define TOLERANCE 1e-6
#define CutVioTolerance 1e-6
#define ten_FOLAT_TOLERANCE 1e-4
#define RC_TOLERANCE (-1e-6)
#define MIP_TOLERANCE 1e-4
/** note for my vrptw instance, the gap should be adjust to 1e-6
 */
#ifdef SOLVER_VRPTW
#define MIPGap 1e-9
#else
#define MIPGap 1e-8
#endif
#define MAXNOOFCUTS 100
#define LIMIT_NEW_ADDED_COLS 10000
#define LPCol_FINAL_LIMIT 10000
#define InitialNumCol4MapInRootMIP 100000
#define MAXTIMELIMT4MIP 1000000000 //no status check for now
#define MAXCUTOFF 1000000000
#define MAXTIMELIMT4R1CMIP 1 //no status check for now
#define MAXMATRIXCOLS4MULTIPLICATION 65536
#define FIND_MEM_USE_ENUMERATION_OR_MIP 1000
#define Initial_rank1_multi_label_pool_size 2048
#define MaxNumColsBasic 1000000
#define TimeLimit4MIPFindMem 0.2

//#define DEBUG_LP_FILE
//#define DEBUG_MIP_FILE

//#define DEBUG
//#define debugMLInputData


//section control all mem-related hyperparameters

//these are strictly satisfied!
#ifdef SOLVER_VRPTW
#define MaxNum_Customers 512 //exit(-1)
//print out the MaxNum_Customers
#else
#define MaxNum_Customers 600 //exit(-1)
#endif
#define MaxNum_R1Cs 1024//exit(-2)
#define MaxNum_R1C_multi 1024//exit(-3)
#define CST_LIMIT 3000//exit(-3)
#define FACTOR_NUM_LABEL 100
//#define BranchFashion_MemSaving //save the memory

//these can be changed in size!
//max_int ~ 2.1e9

#define LABEL_ASSIGN 8000000
#define MAX_ROUTE_MEM 10000000
#define MAX_ROUTE_PRICING 10000000
#define LABEL_LIMIT_PER_BIN 128

#define MAX_NG_SIZE 8
#define PATH_LENGTH_TOLERANCE_FOR_NG_AUGMENTATION 5
#define MAX_OPTCOLUMN_SELECT 100

#define AUGMENT_NG 1

#define SMALL_PHASE_SEPARATION "--------------------------------------------------------------------------------------------------------------\n"
#define BIG_PHASE_SEPARATION "##############################################################################################################\n"
#define MID_PHASE_SEPARATION "**************************************************************************************************************\n"

#define TERMINATED_MESSAGE_PROMISING_VEHICLES "terminate reason: node is guaranteed to be unpromising when solving lp!\n"
#define TERMINATED_MESSAGE_PROMISING_UPDATE_UB "terminate reason: node is guaranteed to be unpromising after finding UB!\n"
#define TERMINATED_MESSAGE_SEP_RCC "terminate reason: node is guaranteed to be unpromising when separating RCCs!\n"
#define TERMINATED_MESSAGE_BRANCHING "terminate reason: node is guaranteed to be unpromising when two-side branching terminated!\n"

//indicate where might be the problem
#define MemWarning "Notice there is a mem-warning sign!\n"

#define safe_solver(call){ int solver_err = (call);\
if (solver_err != 0) {                             \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call+"\n ERROR CODE= " +std::to_string(solver_err)); \
}}                                                 \

#define safe_Hyperparameter(call){ int Hyperparameter_err = (call);\
if (Hyperparameter_err != 0) {                             \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call+"\n ERROR CODE= " +std::to_string(Hyperparameter_err)); \
}}                                                                 \

#ifdef INCORPORATE_ML
//#define GENERATE_ML_EXACT_DATA//generate exact data // need train_EXACT folder
//#define ONLY_USE_MODEL//use model
#endif

#include <bitset>
using yzzLong = std::bitset<MaxNum_Customers>;

struct PairHasher {
  size_t operator()(const std::pair<int, int> &V) const {
	return V.first * MaxNum_Customers + V.second;
  }
};

#if !defined(NoAdjustmentMIPTimeLimit) && defined (SolveMIPNotFixAll)
#error "AdjustmentMIPTimeLimit is not compatible with SolveMIPNotFixAll For now!"
#endif

#endif //CVRP_MACRO_HPP
