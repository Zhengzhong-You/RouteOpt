//
// Created by You, Zhengzhong on 4/26/24.
//

#ifndef INCLUDE_BRANCHING_BRANCHING_HPP_
#define INCLUDE_BRANCHING_BRANCHING_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"

#define HEURISTIC_LIGHT_TESTING_MAX_COLUMN_RATIO 1.5// over 1.5* 10000, the testing will be stopped
#define HEURISTIC_HEAVY_TESTING_MAX_COLUMN_RATIO 2//over 2 * 10000, the testing will be stopped
#define ADJUSTMENT_SCORE_RATIO 3

#define OUTPUT_OPT 1
#define OUTPUT_TIME_LIMIT 2
#define OUTPUT_INIT 3
#define OUTPUT_BIG_SEP 4
#define OUTPUT_COLUMN 5
#define OUTPUT_DYNAMIC_SEARCH 6
#define OUTPUT_NODE_INFO 7
#define OUTPUT_BUILD_MODEL 8
#define OUTPUT_SUBTREE 9

struct BranchingHistory {
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> exact_improvement_up{}; // sum and cnt
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> exact_improvement_down{}; // sum and cnt
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> heuristic_improvement_up{};
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> heuristic_improvement_down{};
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> lp_testing_improvement_up{};
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> lp_testing_improvement_down{};
    std::unordered_map<std::pair<int, int>, int, PairHasher> branch_choice{};
    std::vector<std::pair<std::pair<double, int>, std::pair<double, int> > > increase_depth; //left, cnt, right, cnt
};

struct CurrentBranchingInfo {
    std::vector<std::pair<std::pair<int, int>, double> > branch_pair_val{};
    std::vector<std::pair<int, int> > branch_selected_from_model1{};
    std::vector<std::pair<int, int> > branch_selected_from_model2{};
    std::vector<std::pair<int, int> > branch_pair{};
    std::unordered_map<std::pair<int, int>, double, PairHasher> branch_lp{};
    std::vector<std::pair<int, int> > branch_pair_from_pseudo{};
    std::vector<std::pair<int, int> > branch_pair_from_fractional{};
};

struct EdgeScoreInfo {
    std::pair<int, int> edge{};
    double dif1{}, dif2{}; //left and right
    double ratio{1}; //max/min
    bool right_max{};
};

class BaseBranching {
public:
    static CVRP *cvrp;
    //
    static BbNode *node;
    static double lb_transformed, ub, lb;
    static double sub_lb_transformed, sub_lb;
    static int num_explored_nodes, num_br, nd_rmn;
    static BranchingHistory branching_history;
    static CurrentBranchingInfo current_branching_info;
    static bool if_terminate_tree;
    static bool if_test_cg;
    static std::chrono::time_point<std::chrono::high_resolution_clock> glo_beg,
            glo_end;
    static double glo_eps;
    static double global_gap;

    static void reviseExtremeUnbalancedScore(
        std::vector<EdgeScoreInfo> &edge_info,
        bool if_record_LP_improvement
    );


    static void init(CVRP *pr_cvrp);

    static void solveNode();

    static void tellIfTerminateTree();

    static void prepare();

    static void recordBranchingHistory(const Brc &brc, double old_val);

    static double tellBranchingValue();

    static void freeNode();

    static void printInfo(int mode);

    static void resetEnv();

    template<typename SubTreeType>
    static void prepareEnuTree(SubTreeType &sub_bbt);

    static void updateLowerBound();

    static void initialScreen(bool if_record_source,
                              bool if_fill_candidates_gap,
                              int num,
                              double pseudo_frac);

    static double calculateDifference(double tmp_val, double prior_val, bool if_consider_ub = false);

    static void testLP(int num, bool if_writeBranch_pair, bool if_record_sb_scores = false,
                       bool if_record_LP_improvement = true);

    static void testCG(bool if_exact_CG, bool if_record_product_val,
                       bool if_record_improvement, bool if_force_complete = false);

    static void controlBrSelection(std::pair<int, int> &info);

    static void useDefaultSB();
};

template<typename SubTreeType>
void BaseBranching::prepareEnuTree(SubTreeType &sub_bbt) {
    cvrp->prepareSolveEnuTree(node);

    sub_bbt.push(node);
    --num_explored_nodes;
}

#endif //INCLUDE_BRANCHING_BRANCHING_HPP_
