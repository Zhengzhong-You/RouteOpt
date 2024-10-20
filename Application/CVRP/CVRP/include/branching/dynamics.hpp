//
// Created by You, Zhengzhong on 5/5/24.
//

#ifndef INCLUDE_BRANCHING_DYNAMICS_HPP_
#define INCLUDE_BRANCHING_DYNAMICS_HPP_

#include <vector>
#include <unordered_map>
#include <tuple>
#include <utility>
#include "macro.hpp"

#define RATIO_HEURISTIC_LP 20
#define MAX_HEURISTIC_ALPHA 0.8
#define ABS_HEURISTIC_NO_MORE 4
#define MIP_TIME_ENU_TIME_RATIO 8
#define ENUMERATION_TIME_CONCERN_LIMIT 10
class CVRP;
class BbNode;

class Dynamics {
public:
    static CVRP *cvrp;
    static BbNode *node;
    static double alpha;
    static double est_m;
    static double f;
    static double opt_k;
    static double solve_cg_time;
    static bool is_use_more_k;
    static std::vector<std::pair<std::pair<double, int>, std::pair<double, int> > > r_star_depth;

    static double solve_a_node_time;
    static double t_for_one_lp_b4;
    static double c_b4;
    static int n_nodes_b4;

    static double t_for_one_lp_enu;
    static double c_enu;
    static int n_nodes_enu;


    static double alpha_heuristic;
    static double est_m_heuristic;
    static double t_for_one_heuristic_b4;
    static double t_for_one_heuristic_enu;

    static double r_best;

    int k_node{0};

    int &getKNode() { return k_node; }

    static void init(CVRP *cvrp);

    static void updateNode(BbNode *node);

    static void updateState(double new_value, double &old_value, int n);

    static void updateStateWithWeights(double new_value, double &old_value, int n);

    static void updateStateAverage(double new_value, double &old_value, int n);

    static void calculateRStar(double lift);

    static void calculateF(double eps, double node_value = 0);

    static void evaluateM1();

    static void giveDynamicK(int &num, bool if_lp = true);

    static void getAverageT4LPNHeuristic(double average_t, bool if_lp = true);

    static void copyDynamicData4EachNode(const Dynamics &d, Dynamics &d2);
};

#endif //INCLUDE_BRANCHING_DYNAMICS_HPP_
