/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_L2B_MACRO_HPP
#define ROUTE_OPT_L2B_MACRO_HPP
#include <string>

namespace RouteOpt::Application::CVRP {
    constexpr int L2B_PHASE0 = 100;
    constexpr int MAX_NUM_CLUSTER = 10;
    constexpr int MAX_NUM_RUN_OPTIMAL_K = 10;
    constexpr int ML_RANDOM_SEED = 42;
    constexpr int L2B_SIMULATE_NUM = 30;
    constexpr int MAX_TREE_LEVEL = 5;
    constexpr int MAX_R_SECOND_STAGE = 2; //scale the prediction value to 0-2!
    constexpr int MAX_RANK_DEVIATION = (MAX_R_SECOND_STAGE + 2);
    constexpr int BEST_PRE = (-1);
    constexpr int TRUST_PRE_BAR = (MAX_R_SECOND_STAGE - 1);
    constexpr double THRESHOLD_DELTA = 1.;

    constexpr std::string_view TRAIN_FOLDER_LP = "train_lp";
    constexpr std::string_view TRAIN_FOLDER_EXACT = "train_exact";

#define PseudoMark std::string("pseudo-")
#define SAFE_XGBOOST_CALL(call) {  int Xerr = (call);\
if (Xerr != 0) { \
throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call + ":" + XGBGetLastError()); \
}\
}
// #define DEBUG_ML_INPUT_DATA
#ifdef DEBUG_ML_INPUT_DATA
#define DEBUG_INPUT_DATA_CALL(...) __VA_ARGS__;
#else
#define DEBUG_INPUT_DATA_CALL(...)
#endif

#define L2B_VERBOSE

#ifdef L2B_VERBOSE
#define   L2B_VERBOSE_EXEC(...) __VA_ARGS__;
#else
    L2B_VERBOSE_EXEC();
#endif
}

#endif // ROUTE_OPT_L2B_MACRO_HPP
