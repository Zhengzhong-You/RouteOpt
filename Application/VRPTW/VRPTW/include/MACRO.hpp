
#ifndef CVRP_MACRO_HPP
#define CVRP_MACRO_HPP

#include <cstddef>
#include <utility>
#include "techniqueControl.hpp"

#define MAX_NUM_CUSTOMERS 512
#define DUAL_TOLERANCE 1e-6
#define PRINT_LABELING_STEP_SIZE 20
#define LARGE_FLOAT 3e30
#define MAX_INT 2000000000
#define VERBOSE_MODE
#define READ_NO_LINE (-1)
#define GLOBAL_TIME_LIMIT 1000000000
#define NUM_THREADS_LP 1
#define COL_KEEP_FRAC 0.67
#define TOLERANCE 1e-6
#define CUT_VIO_TOLERANCE 1e-6
#define RC_TOLERANCE (-1e-6)
#define MIP_GAP_TOLERANCE 1e-9
#define MAX_NUM_OF_CUTS 100
#define LIMIT_NEW_ADDED_COLUMNS 10000
#define LP_COL_FINAL_LIMIT 10000
#define MAX_TIME_LIMIT_FOR_MIP 1000000000 //no status check for now
#define FIND_MEM_USE_ENUMERATION_OR_MIP 1000 //when to use MIP to find the least memory
#define INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE 2048
#define MAX_NUM_R1CS 512//exit(-2)
#define MAX_NUM_R1C_MULTI 512//exit(-3)
#define CST_LIMIT 3000//exit(-3)
#define FACTOR_NUM_LABEL 100//assign labels in the bucket graph
#define LABEL_ASSIGN 800000
#define MAX_ROUTE_MEMORY 100000
#define MAX_ROUTE_PRICING 1600000
#define LABEL_LIMIT_PER_BIN 128
#define WRITE_SOL_OUT_TO_FILE "sol"

#ifdef FIND_MEM_USE_ENUMERATION_OR_MIP
#define TIME_LIMIT_FOR_MIP_FIND_MEM 0.2
#endif

#ifdef WRITE_ENUMERATION_TREES
#define BRANCH_FASHION_MEM_SAVING //save the memory
#endif

#define SMALL_PHASE_SEPARATION "--------------------------------------------------------------------------------------------------------------\n"
#define BIG_PHASE_SEPARATION "##############################################################################################################\n"
#define MID_PHASE_SEPARATION "**************************************************************************************************************\n"

#define TERMINATED_MESSAGE_PROMISING_VEHICLES "terminate reason: node is guaranteed to be unpromising when solving lp!\n"
#define TERMINATED_MESSAGE_PROMISING_UPDATE_UB "terminate reason: node is guaranteed to be unpromising after finding ub!\n"
#define TERMINATED_MESSAGE_SEP_RCC "terminate reason: node is guaranteed to be unpromising when separating rccs!\n"

#define safe_solver(call){ int solver_err = (call);\
if (solver_err != 0) {                             \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call+"\n ERROR CODE= " +std::to_string(solver_err)); \
}}                                                 \

#define safe_Hyperparameter(call){ int Hyperparameter_err = (call);\
if (Hyperparameter_err != 0) {                             \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call+"\n ERROR CODE= " +std::to_string(Hyperparameter_err)); \
}}                                                                 \

struct PairHasher {
  size_t operator()(const std::pair<int, int> &V) const {
	return V.first * MAX_NUM_CUSTOMERS + V.second;
  }
};

#if defined(WRITE_ENUMERATION_TREES) || defined(READ_ENUMERATION_TREES)
#define TREE_FOLDER std::string("EnuTree")
#define COL_POOL_FOLDER std::string("EnuColPool")
#endif

#if defined(MASTER_VALVE_ML) && (ML_STATE == 1 || ML_STATE == 2)
#define TRAINING_DATA_TREE_LEVEL 8
#endif

#ifdef SOLVER_VRPTW
#define SYMMETRY_PROHIBIT
#else
#define CAPACITY_AS_MAIN_RESOURCE
#endif

#if defined(WRITE_ENUMERATION_TREES) && defined(READ_ENUMERATION_TREES)
#error "WRITE_ENUMERATION_TREES and READ_ENUMERATION_TREES cannot be defined at the same time"
#endif

#if defined(ML_STATE) && ML_STATE == 4
#define USE_M_DYNAMICS
#endif

#endif //CVRP_MACRO_HPP
