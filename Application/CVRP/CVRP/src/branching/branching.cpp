//
// Created by You, Zhengzhong on 4/26/24.
//

#include "branching.hpp"

#if SOLUTION_TYPE == 1
#include "best_bound_first_branching.hpp"
#elif SOLUTION_TYPE == 2
#include "depth_first_branching.hpp"
#endif

#ifdef WRITE_ENUMERATION_TREES
#include "write_enumeration_tree.hpp"
#endif

#ifdef HEURISTIC
#include "heuristic.hpp"
#endif

#ifdef READ_ENUMERATION_TREES
#include "read_enumeration_tree.hpp"
#endif

#ifdef FASTER_DUAL
#include "faster_dual.hpp"
#endif

#ifdef DUAL_SMOOTHING
#include "dual_smoothing.hpp"
#endif

#ifdef COMBINE_FASTER_SMOOTHING
#include "combine_faster_smoothing.hpp"
#endif

#include "write_node_out.hpp"

#ifdef MASTER_VALVE_ML
#include "train.hpp"
#include "predict.hpp"
#endif

#include "read_node_in.hpp"
#include "robust_control.hpp"

CVRP *BaseBranching::cvrp{};
BbNode *BaseBranching::node{};
double BaseBranching::lb{};
double BaseBranching::lb_transformed{};
double BaseBranching::sub_lb{};
double BaseBranching::sub_lb_transformed{};
double BaseBranching::ub{};
int BaseBranching::num_explored_nodes{};
int BaseBranching::num_br{};
int BaseBranching::nd_rmn{};
BranchingHistory BaseBranching::branching_history{};
CurrentBranchingInfo BaseBranching::current_branching_info{};
bool BaseBranching::if_terminate_tree{};
bool BaseBranching::if_test_cg{};
double BaseBranching::glo_eps{};
std::chrono::high_resolution_clock::time_point BaseBranching::glo_beg{};
std::chrono::high_resolution_clock::time_point BaseBranching::glo_end{};
double BaseBranching::global_gap{};


using namespace std;
using namespace std::chrono;

void getGeomean(double &old_v, int &n, double new_v);

void BaseBranching::tellIfTerminateTree() {
    if (cvrp->getIfInEnuState()) {
        if (sub_lb_transformed >= ub) {
            sub_lb_transformed = ub;
            printInfo(OUTPUT_SUBTREE);
            freeNode();
            if_terminate_tree = true;
            goto QUIT;
        }
    } else {
        if (lb_transformed >= ub) {
            lb_transformed = ub;
            printInfo(OUTPUT_OPT);
            freeNode();
            if_terminate_tree = true;
            goto QUIT;
        }
    }

    if (glo_eps > GLOBAL_TIME_LIMIT) {
        printInfo(OUTPUT_TIME_LIMIT);
        freeNode();
        if_terminate_tree = true;
    }
QUIT:;
}

void BaseBranching::init(CVRP *pr_cvrp) {
    cvrp = pr_cvrp;
}

void BaseBranching::prepare() {
    heuristic_call(Heuristic::cleanHeuristicMIPCols();)

    cvrp->initialProcessing();
    glo_beg = std::chrono::high_resolution_clock::now();

#ifdef READ_ENUMERATION_TREES
  ReadEnumerationTree::init(cvrp);
  ReadEnumerationTree::getPath();
  ReadEnumerationTree::restoreModel();
#elif defined(READ_NODE_IN)
    ReadNodeIn::init(cvrp);
    ReadNodeIn::recoverNodeInfo();
#else
    cvrp->buildModel();
#endif

    glo_end = chrono::high_resolution_clock::now();
    glo_eps = duration<double>(glo_end - glo_beg).count();

    printInfo(OUTPUT_BUILD_MODEL);

#ifdef WRITE_ENUMERATION_TREES
  WriteEnumerationTree::init(cvrp);
#endif

    RobustControl::init(cvrp);
    faster_dual_call(FasterDual::init(cvrp);)
    dual_smoothing_call(DualSmoothing::init(cvrp))
    combine_faster_call(CombineFasterSmoothing::init(cvrp))
    ml_call(MASTER_VALVE_ML, MachineLearning::init(cvrp))
    ml_call(MASTER_VALVE_ML, MachineLearning::calculatePrerequisites())
    ml_call(MASTER_VALVE_ML == ML_GET_DATA_1 || MASTER_VALVE_ML == ML_GET_DATA_2, GetTrainingData::initOutputPath())
    ml_call(MASTER_VALVE_ML == ML_GET_DATA_2, Predict::loadModel(1))
    ml_call(
        MASTER_VALVE_ML == ML_USE_MODEL || MASTER_VALVE_ML == ML_USE_MODEL_1,
        (Predict::loadModel(1), Predict::loadModel(2)))
    dynamic_call(Dynamics::init(cvrp))
    write_node_out_call(WriteNodeOut::init(cvrp))
}


void BaseBranching::solveNode() {
    dynamic_call(auto beg_c = std::chrono::high_resolution_clock::now())
    double old_val = node->obtainParentNodeValue();
    ml_call((MASTER_VALVE_ML == ML_GET_DATA_1 || MASTER_VALVE_ML == ML_GET_DATA_2),
            if (node->getTreeLevel() > MAX_TREE_LEVEL) {
            verbose_call(cout << "The tree level is too deep for training data collection!" << endl)
            freeNode();
            goto QUIT;
            })


    cvrp->startSolveNode(node);


    if (node->getIfTerminated()) {
        freeNode();
        goto QUIT;
    }

    if (node->getTreeLevel())recordBranchingHistory(node->getBrCs().back(), old_val);

    if (!cvrp->getIfInEnuState()) updateLowerBound();
    cvrp->postSolveNode(node);
    if (!node) goto QUIT;


    cvrp->enterCuttingPhase(node);
    if (!node) goto QUIT;

    verbose_call(printInfo(OUTPUT_COLUMN))

    dynamic_call({
        auto end_c = std::chrono::high_resolution_clock::now();
        auto eps = duration<double>(end_c - beg_c).count();
        Dynamics::solve_a_node_time = eps;
        verbose_call(printInfo(OUTPUT_DYNAMIC_SEARCH))
        if (node->getTreeLevel()) {
        Dynamics::calculateRStar(node->getCurrentNodeVal() - old_val);
        }
        })
    printInfo(OUTPUT_BIG_SEP);
QUIT:
    cvrp->finishSolveNode();
}

void BaseBranching::freeNode() {
    delete node;
    node = nullptr;
}

void BaseBranching::printInfo(int mode) {
    switch (mode) {
        case OUTPUT_OPT: cout << "Pruned! Optimality has been proven!  lb=ub= " << ub << endl;
            cout << SMALL_PHASE_SEPARATION;
            break;
        case OUTPUT_TIME_LIMIT: cout << SMALL_PHASE_SEPARATION;
            cout << "solution process shuts down due to reaching lm. gt.= " << GLOBAL_TIME_LIMIT << endl;
            cout << "lb= " << lb << "  ub= " << ub << "  gap(ub-lb/ub)= "
                    << (ub - lb) / (ub) * 100 << "%  gt= "
                    << glo_eps << " nd= "

                    << num_explored_nodes << "  br= " << num_br << endl;
            break;
        case OUTPUT_BIG_SEP: cout << BIG_PHASE_SEPARATION;
            break;
        case OUTPUT_COLUMN: cout << BIG_PHASE_SEPARATION;
            cout << "Run column reduction in memory... ncol= " << cvrp->getNumCol() << " nrow= " << cvrp->getNumRow()
                    << endl;
            cout << BIG_PHASE_SEPARATION;
            break;
        case OUTPUT_DYNAMIC_SEARCH: verbose_call(
                cout << "c_b4= " << Dynamics::c_b4 << " c_enu= " << Dynamics::c_enu << " t_b4= " << Dynamics::
                t_for_one_lp_b4
                << " t_enu= " << Dynamics::t_for_one_lp_enu << " heu_b4= " << Dynamics::t_for_one_heuristic_b4
                << " heu_enu= " << Dynamics::t_for_one_heuristic_enu<<endl;
            )
            break;
        case OUTPUT_NODE_INFO: cout << BIG_PHASE_SEPARATION;
            cout << "nd_ind= " << node->getNodeIdx() << "  nd_col= " << cvrp->getNumCol() << "  nd_val= " << node->
                    getCurrentNodeVal() <<
                    "  nd_dep= "
                    << node->getTreeLevel() << "  et= " << glo_eps << "  lb= " << lb << "  ub= "
                    << ub << "  nd_rmn= " << nd_rmn << endl;
            node->printBrCInfo();
            break;
        case OUTPUT_BUILD_MODEL:
            cout << "Build model...   ub= " << ub << "  et= " << glo_eps << endl;
            cout << BIG_PHASE_SEPARATION;
            break;
        case OUTPUT_SUBTREE:
            cout << "The subtree has been terminated for tree lb= " << node->getCurrentNodeVal() << " but ub= "
                    << ub << endl;
        default: break;
    }
}

void BaseBranching::resetEnv() {
    safe_solver(node->getSolver().updateModel())
    safe_solver(node->getSolver().getNumRow(&cvrp->getNumRow()))
    safe_solver(node->getSolver().getNumCol(&cvrp->getNumCol()))
    cvrp->getVCutMapLP(node);
#if VERBOSE_MODE == 1
    printInfo(OUTPUT_INIT);
#endif
    printInfo(OUTPUT_NODE_INFO);

    RobustControl::updateNode(node)
            faster_dual_call(FasterDual::updateNode(node))dual_smoothing_call(DualSmoothing::updateNode(node))
    combine_faster_call(CombineFasterSmoothing::updateNode(node))
    ml_call(MASTER_VALVE_ML, MachineLearning::updateNode(node))
    dynamic_call(Dynamics::updateNode(node))
    write_node_out_call(WriteNodeOut::updateNode(node))
}

void BaseBranching::recordBranchingHistory(const Brc &brc, double old_val) {
    double dif = max(node->getCurrentNodeVal() - old_val, TOLERANCE);
    auto &exact_improvement_up = branching_history.exact_improvement_up;
    auto &exact_improvement_down = branching_history.exact_improvement_down;
    auto &increase_depth = branching_history.increase_depth;
    if (increase_depth.size() <= node->getTreeLevel()) increase_depth.resize(node->getTreeLevel() + 1);
    auto &r_star = increase_depth[node->getTreeLevel()];
    auto dir = brc.br_dir;
    auto &recordings = dir ? r_star.second.second : r_star.first.second;
    auto &increase = dir ? r_star.second.first : r_star.first.first;
    if (node->getIfJustEnterEnu()) {
        node->getNodeBrValueImproved() = 0;
        return;
    }
    if (node->getCurrentNodeVal() - old_val < RC_TOLERANCE * old_val) {
        cout << "node->getCurrentNodeVal()=" << node->getCurrentNodeVal() << " old_val=" << old_val << endl;
        throw runtime_error("branching leads to a decreasing in lb, plz check the code!");
    }
    if (brc.br_dir) {
        exact_improvement_up[brc.edge].first += dif;
        ++exact_improvement_up[brc.edge].second;
    } else {
        exact_improvement_down[brc.edge].first += dif;
        ++exact_improvement_down[brc.edge].second;
    }
    getGeomean(increase, recordings, dif);
    // verbose_call(cout << "Branching history: " << endl;
    //     for (int i = 0; i < increase_depth.size(); ++i) {
    //     cout << "Tree level: " << i << " increase: " << increase_depth[i].first.first << " "
    //     << increase_depth[i].second.first << endl;
    //     cout << "left recoding: " << increase_depth[i].first.second << " right recording: "
    //     << increase_depth[i].second.second << endl;
    //     })
    node->getNodeBrValueImproved() = tellBranchingValue();
}

void BaseBranching::updateLowerBound() {
#if SOLUTION_TYPE == 1
    BestBoundFirstBranching::updateLowerBound();
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::updateLowerBound();
#endif
}

double BaseBranching::tellBranchingValue() {
    if (branching_history.increase_depth.size() <= node->getTreeLevel())
        throw runtime_error(
            "the tree level is too deep!");
    auto &left_increase = branching_history.increase_depth[node->getTreeLevel()].first;
    auto &right_increase = branching_history.increase_depth[node->getTreeLevel()].second;
    double r_star;
    if (left_increase.second == 0 || right_increase.second == 0) {
        r_star = numeric_limits<float>::max();
        for (auto &[fst, snd]: branching_history.increase_depth) {
            if (fst.second == 0 || snd.second == 0) continue;
            r_star = min(r_star, sqrt(fst.first * snd.first));
        }
        if (r_star == numeric_limits<float>::max()) {
            const double val = (left_increase.second == 0) ? right_increase.first : left_increase.first;
            r_star = val * INITIAL_CUTTING_BRANCHING_RATIO;
            write_node_out_call(r_star/=INITIAL_CUTTING_BRANCHING_RATIO) //reverse the effect of the ratio
        }
    } else {
        r_star = sqrt(left_increase.first * right_increase.first);
    }
    return r_star * R_DISCOUNT;
}

void BaseBranching::reviseExtremeUnbalancedScore(
    std::vector<EdgeScoreInfo> &edge_info,
    bool if_record_LP_improvement
) {
    double geo_ratio = 0;

    for (auto &edge: edge_info) {
        edge.right_max = (edge.dif2 >= edge.dif1);
        if (edge.right_max) {
            edge.ratio = edge.dif2 / edge.dif1;
        } else {
            edge.ratio = edge.dif1 / edge.dif2;
        }
        geo_ratio += log(edge.ratio);
    }
    geo_ratio = std::exp(geo_ratio / (double) edge_info.size()) * ADJUSTMENT_SCORE_RATIO;
    double max_dif = 0;
    for (auto &edge: edge_info) {
        if (edge.ratio > geo_ratio) continue;
        if (edge.right_max) {
            max_dif = max(max_dif, edge.dif2);
        } else {
            max_dif = max(max_dif, edge.dif1);
        }
    }

    for (auto &edge: edge_info) {
        if (edge.ratio <= geo_ratio) continue;
        if (edge.right_max) {
            edge.dif2 = min(max_dif, edge.dif2);
            cout << "edge: " << edge.edge.first << "-" << edge.edge.second << " dif2: " << edge.dif2 << endl;
        } else {
            edge.dif1 = min(max_dif, edge.dif1);
            cout << "edge: " << edge.edge.first << "-" << edge.edge.second << " dif1: " << edge.dif1 << endl;
        }
    }

    auto &branch_lp = current_branching_info.branch_lp;
    branch_lp.clear();
    for (auto &edge: edge_info) {
        branch_lp[edge.edge] = edge.dif1 * edge.dif2;
    }

    if (if_record_LP_improvement) {
        auto &lp_testing_improvement_down = branching_history.lp_testing_improvement_down;
        auto &lp_testing_improvement_up = branching_history.lp_testing_improvement_up;
        for (auto &edge: edge_info) {
            lp_testing_improvement_down[edge.edge].first += edge.dif1;
            ++lp_testing_improvement_down[edge.edge].second;
            lp_testing_improvement_up[edge.edge].first += edge.dif2;
            ++lp_testing_improvement_up[edge.edge].second;
        }
    }
}

void getGeomean(double &old_v, int &n, double new_v) {
    if (new_v < TOLERANCE) return;
    if (n == 0) old_v = new_v;
    old_v = double(n) / (n + 1) * log(old_v) + 1. / (n + 1) * log(new_v);
    old_v = exp(old_v);
    ++n;
}
