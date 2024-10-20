
#ifndef CVRP_MACRO_HPP
#define CVRP_MACRO_HPP

#include <cstddef>
#include <utility>
#include "techniqueControl.hpp"
#include <bitset>
#include <unordered_map>
#include <vector>

#define RESOURCE_FACTOR 100 //times 100 to make it integer
#define MAX_NUM_REGENERATE_BUCKET 4
#define MAX_NUM_CUSTOMERS 400
#define DUAL_TOLERANCE 1e-6
#define PRINT_LABELING_STEP_SIZE 30
#define LARGE_FLOAT 3e30
#define MAX_INT 2000000000
#define READ_NO_LINE (-1)
#define GLOBAL_TIME_LIMIT 1000000000
#define NUM_THREADS_LP 1
#define COL_KEEP_FRAC 0.67
#define TOLERANCE 1e-6
#define CUT_VIO_TOLERANCE 1e-6
#define RC_TOLERANCE (-1e-6)
#define MIP_GAP_TOLERANCE 1e-9
#define MAX_NUM_OF_CUTS 100
#define LP_COL_FINAL_LIMIT 10000
#define MAX_TIME_LIMIT_FOR_MIP 1000000000 //no status check for now
#define FIND_MEM_USE_ENUMERATION_OR_MIP 1000 //when to use MIP to find the least memory
#define INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE 2048
#define MAX_NUM_R1CS_IN_PRICING 1024
#define MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX 128
#define FACTOR_NUM_LABEL 100//assign labels in the bucket graph
#define LABEL_ASSIGN 8000000
#define MAX_ROUTE_PRICING 6400000
#define LABEL_LIMIT_PER_BIN 128
//#define CUT_ANGLE 0.9
#define WRITE_SOL_OUT_TO_FILE "sol"
#define ROLL_BACK_INIT 0
#define ROLL_BACK_LONG_TIME 1
#define ROLL_BACK_LACK_MEMORY 2
#define ROLL_BACK_TAIL_OFF 3
#define SPARSITY_ESTIMATOR 0.05
#define FRESH_MARK (-1)

using res_int = int;//cannot be unsigned type
using yzzLong = std::bitset<MAX_NUM_CUSTOMERS>;
using R1CINDEX = std::bitset<MAX_NUM_R1CS_IN_PRICING>;

#ifdef FIND_MEM_USE_ENUMERATION_OR_MIP
#define TIME_LIMIT_FOR_MIP_FIND_MEM 0.2
#endif

#ifdef WRITE_ENUMERATION_TREES
#define BRANCH_FASHION_MEM_SAVING //save the memory
#endif

#define SMALL_PHASE_SEPARATION "---------------------------------------------------------------------------------------------------------------\n"
#define BIG_PHASE_SEPARATION "##################################################################################################################\n"
#define MID_PHASE_SEPARATION "******************************************************************************************************************\n"

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
#define VRPTW_DISTANCE_TOLERANCE 1e-6
#if SOLVER_VRPTW == 2
#define CAPACITY_AS_MAIN_RESOURCE
#endif
#if SOLVER_VRPTW == 2 || SOLVER_VRPTW == 3
#define USE_TWO_RESOURCE
#endif
#else //CVRP
#define CAPACITY_AS_MAIN_RESOURCE
#endif

#if defined(WRITE_ENUMERATION_TREES) && defined(READ_ENUMERATION_TREES)
#error "WRITE_ENUMERATION_TREES and READ_ENUMERATION_TREES cannot be defined at the same time"
#endif

#if (defined(ML_STATE) && ML_STATE == 4) || (BRANCH_CANDIDATES_TYPE == 4)
#define USE_M_DYNAMICS
#endif

#define RANK1_INVALID (-1)

#ifdef EXTRA_ARC_MEMORY
using ArcMap = std::unordered_map<std::pair<int, int>, int, PairHasher>;
#endif

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

#if LIMITED_MEMORY_TYPE == 2
#ifdef SOLVER_VRPTW
#define FIXED_MEET_POINT_FOR_RESOURCE
#endif
#endif

#ifdef USE_M_DYNAMICS
#define MAX_M_COEFF 0.8
#define MIN_M_COEFF 0.2
#endif

#if BRANCH_CANDIDATES_TYPE == 4 && defined(MASTER_VALVE_ML)
#error "BRANCH_CANDIDATES_TYPE 4 and MASTER_VALVE_ML cannot be defined at the same time"
#endif
/**
 * hidden techniques
 */

//#define MEMORY_SELECTION

#if LIMITED_MEMORY_TYPE == 2
#define PRICING_HARD
//#define EXTRA_ARC_MEMORY
#endif

//#define CHECK_PRICING_LABELS

#endif //CVRP_MACRO_HPP
