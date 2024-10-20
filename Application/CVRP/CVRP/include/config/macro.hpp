#ifndef CVRP_MACRO_HPP
#define CVRP_MACRO_HPP

#include <cstddef>
#include <utility>
#include "technique_control.hpp"
#include <bitset>
#include <unordered_map>
#include <vector>
#include "setting2.hpp"

#define RESOURCE_FACTOR 100 //times 100 to make it integer
#define MAX_NUM_REGENERATE_BUCKET 4
#define MAX_NUM_CUSTOMERS 400
#define DUAL_TOLERANCE 1e-6
#define PRINT_LABELING_STEP_SIZE 30 //30
#define LARGE_FLOAT 3e30
#define READ_NO_LINE (-1)
#define GLOBAL_TIME_LIMIT 1000000000
#define NUM_THREADS_LP 1
#define COL_KEEP_FRAC 0.67
#define TOLERANCE 1e-6
#define CUT_VIO_TOLERANCE 1e-6
#define RC_TOLERANCE (-1e-4)
#define MIP_GAP_TOLERANCE 1e-9
#define MAX_NUM_OF_CUTS 100
#define LP_COL_FINAL_LIMIT 10000
#define MAX_TIME_LIMIT_FOR_MIP 1000000000 //no status check for now

#define FACTOR_NUM_LABEL 100//assign labels in the bucket graph
#define LABEL_ASSIGN 8000000
#define MAX_ROUTE_PRICING 6400000
#define LABEL_LIMIT_PER_BIN 128
#define WRITE_SOL_OUT_TO_FILE "sol"
#define VRPTW_DISTANCE_TOLERANCE 1e-6

using res_int = int; //cannot be unsigned type
using yzzLong = std::bitset<MAX_NUM_CUSTOMERS>;



#define SMALL_PHASE_SEPARATION "---------------------------------------------------------------------------------------------------------------\n"
#define BIG_PHASE_SEPARATION "##################################################################################################################\n"
#define MID_PHASE_SEPARATION "******************************************************************************************************************\n"

#define TERMINATED_MESSAGE_PROMISING_VEHICLES "terminate reason: node is guaranteed to be unpromising when solving lp!\n"
#define TERMINATED_MESSAGE_PROMISING_UPDATE_UB "terminate reason: node is guaranteed to be unpromising after finding ub!\n"
#define TERMINATED_MESSAGE_SEP_RCC "terminate reason: node is guaranteed to be unpromising when separating rccs!\n"

#define safe_solver(call){ int solver_err = (call);\
if (solver_err != 0) {                             \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call+"\n ERROR CODE= " +std::to_string(solver_err)); \
}}
#define safe_Hyperparameter(call){ int Hyperparameter_err = (call);\
if (Hyperparameter_err != 0) {                             \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call+"\n ERROR CODE= " +std::to_string(Hyperparameter_err)); \
}}
#define safe_eigen(func) \
    try { \
        func; \
    } catch (const std::exception &e) { \
        std::ostringstream oss; \
        oss << "Error in file " << __FILE__ << " at line " << __LINE__ << ": " << e.what(); \
        throw std::runtime_error(oss.str()); \
    }

struct PairHasher {
    size_t operator()(const std::pair<int, int> &V) const {
        return V.first * MAX_NUM_CUSTOMERS + V.second;
    }
};

#define TREE_FOLDER std::string("EnuTree")
#define COL_POOL_FOLDER std::string("EnuColPool")

#if defined(SOLVER_VRPTW) && SOLVER_VRPTW == 0
#undef SOLVER_VRPTW
#endif

#ifdef SOLVER_VRPTW
#if SOLVER_VRPTW == 0
#undef SOLVER_VRPTW
#endif
#define SYMMETRY_PROHIBIT
#if SOLVER_VRPTW == 2
#define CAPACITY_AS_MAIN_RESOURCE
#endif
#if SOLVER_VRPTW == 2 || SOLVER_VRPTW == 3
#define USE_TWO_RESOURCE
#endif
#else //CVRP
#define CAPACITY_AS_MAIN_RESOURCE
#endif

#ifdef SOLVER_RVRPSTW
#define SYMMETRY_PROHIBIT
#endif

#if defined(WRITE_ENUMERATION_TREES) && defined(READ_ENUMERATION_TREES)
#error "WRITE_ENUMERATION_TREES and READ_ENUMERATION_TREES cannot be defined at the same time"
#endif

#define RANK1_INVALID (-1)

struct Rank1MultiPairHasher {
    size_t operator()(const std::pair<std::vector<int>, int> &V) const {
        size_t hash = 0;
        for (auto &i: V.first) {
            hash *= MAX_NUM_CUSTOMERS;
            hash += i;
        }
        hash *= MAX_NUM_CUSTOMERS;
        hash += V.second;
        return hash;
    }
};

#define MAX_M_COEFF 0.8
#define MIN_M_COEFF 0.2

/**
 * hidden techniques
 */


//#define CHECK_PRICING_LABELS

#ifdef DUAL_SMOOTHING
#define dual_smoothing_call(...) __VA_ARGS__;
#else
#define dual_smoothing_call(method);
#endif

#if SETTING == 2
#define setting_2_call(method) Setting2::method;
#else
#define setting_2_call(method);
#endif

#ifdef SYMMETRY_PROHIBIT
#define symmetry_prohibit_call(...) __VA_ARGS__;
#else
#define symmetry_prohibit_call(...);
#endif

#ifdef SOLVER_VRPTW
#define solver_vrptw_call(method) method;
#else
#define solver_vrptw_call(method);
#endif

#define ml_call(state, ...) \
    if (state) { \
           __VA_ARGS__;\
        };

#ifdef USE_M_DYNAMICS
#define dynamic_call(...) __VA_ARGS__;
#else
#define dynamic_call(...)
#endif

#if  VERBOSE_MODE == 1
#define verbose_call(...) __VA_ARGS__;
#else
#define verbose_call(...)
#endif

#ifdef WRITE_ENUMERATION_TREES
#define write_enumeration_trees_call(...) __VA_ARGS__;
#else
#define write_enumeration_trees_call(...)
#endif

#ifdef FASTER_DUAL
#define faster_dual_call(...) __VA_ARGS__;
#else
#define faster_dual_call(...)
#endif

#ifdef COMBINE_FASTER_SMOOTHING
#define combine_faster_call(...) __VA_ARGS__;
#else
#define combine_faster_call(...)
#endif

#define ML_GET_DATA_1 1
#define ML_GET_DATA_2 2
#define ML_USE_MODEL 3
#define ML_USE_MODEL_1 4

#define BARRIER_SELECTION
#ifdef BARRIER_SELECTION
#define barrier_call(...) __VA_ARGS__
#else
#define barrier_call(...)
#endif

#ifdef HEURISTIC
#define heuristic_call(...) __VA_ARGS__;
#else
#define heuristic_call(...)
#endif

#ifdef READ_ENUMERATION_TREES
#define read_enumeration_trees_call(...) __VA_ARGS__;
#else
#define read_enumeration_trees_call(...)
#endif


#ifdef DYNAMICS_TYPE
#define easy_dynamic_call(...) __VA_ARGS__;
#else
#define easy_dynamic_call(...)
#endif

#ifdef READ_DUAL_SOLUTION
#define read_dual_sol_call(...) __VA_ARGS__;
#else
#define read_dual_sol_call(...)
#endif

#endif //CVRP_MACRO_HPP
