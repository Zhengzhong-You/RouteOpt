//
// Created by You, Zhengzhong on 5/19/24.
//

#ifndef INCLUDE_HEURISTIC_HEURISTIC_HPP_
#define INCLUDE_HEURISTIC_HEURISTIC_HPP_

#include "cvrp.hpp"
#include "template_functors.hpp"

#define HEURISTIC_ENUMERATION_ON //
#define HEURISTIC_DELUXING_ON //

#define GUESSED_OPT_GAP_PER 0.001 // overly large value will lead to a fail heuristic
#define NUM_COL_HEURISTIC_MIP 10000
#define HEURISTIC_MIP_THREADS_NUM 1

#define HEURISTIC_DELUXING_ROUND 20
#define HEURISTIC_DELUXING_BETA1 1000
#define HEURISTIC_DELUXING_BETA2 100000
#define HEURISTIC_DELUXING_TIME_LIMIT 100000
#define HEURISTIC_DELUXING_VERBOSE 0

#ifdef HEURISTIC_DELUXING_ON
#include "deluxing.hpp"
#define heuristic_deluxing_call(...) __VA_ARGS__;
#else
#define heuristic_deluxing_call(...)
#endif

#ifdef HEURISTIC_ENUMERATION_ON
#define heuristic_enumeration_call(...) __VA_ARGS__;
#else
#define heuristic_enumeration_call(...)
#endif

class Heuristic {
public:
    static double fake_ub;
    static bool if_fail;
    static double heuristic_mip_time;
    static std::vector<SequenceInfo> heuristic_col_info;

    static void cleanHeuristicMIPCols();

    static void addHeuristicMIPCols(CVRP *cvrp, BbNode *node);

    static void heuristicMIP(CVRP *cvrp, BbNode *&node);

    static void enumerateHeuristicMIP(CVRP *cvrp, BbNode *node);

    template<bool if_symmetry>
    static void completeHeuristic(CVRP *cvrp, BbNode *node);

    static void heuristicApplyRCF(Solver &local_solver,
                                  int round,
                                  bool if_verbose,
                                  double ub,
                                  int real_dim,
                                  std::vector<int> &col_map);

    static void rmNonCapCols(const double *demand, double cap);

    static void rmNonEleCols();
};

#endif //INCLUDE_HEURISTIC_HEURISTIC_HPP_
