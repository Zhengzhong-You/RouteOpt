#ifndef INCLUDE_ML_HPP_
#define INCLUDE_ML_HPP_

#include "branching_node.hpp"
#include "cvrp.hpp"
#include "branching.hpp"


#if MASTER_VALVE_ML==0
using bst_ulong = unsigned long;
// Dummy type definitions to mimic XGBoost types.
typedef void *BoosterHandle; // Dummy handle for boosters.
typedef void *DMatrixHandle; // Dummy handle for DMatrix.

// Dummy functions for XGBoost-related calls.
inline int XGBoosterCreate(const void *cache, int len, BoosterHandle *out) { return 0; }
inline int XGBoosterLoadModel(BoosterHandle handle, const char *fname) { return 0; }
inline int XGBoosterFree(BoosterHandle handle) { return 0; }

inline int XGBoosterPredict(BoosterHandle handle, DMatrixHandle dmat, int option_mask, unsigned ntree_limit,
                            int training, bst_ulong *out_len, const float **out_result) { return 0; }

inline int XGDMatrixCreateFromMat(const float *data, bst_ulong nrow, bst_ulong ncol, float missing,
                                  DMatrixHandle *out) { return 0; }

inline int XGDMatrixFree(DMatrixHandle handle) { return 0; }
inline const char *XGBGetLastError() { return "dummy error"; }
#else
#include <xgboost/c_api.h>
#endif


//#define DEBUG_ML_INPUT_DATA
#define MAX_NUM_CLUSTER 10
#define MAX_NUM_RUN_OPTIMAL_K 10
#define ML_RANDOM_SEED 42

struct TmpEdgeRelatedData {
    double sb_scores{};
    std::vector<std::pair<std::string, double> > basic_features{};
    std::vector<std::pair<std::string, double> > extra_features_edge0{};
    std::vector<std::pair<std::string, double> > extra_features_edge1{};
    std::vector<std::pair<std::string, double> >
    resolving_lp_features{}; // 1 is the left branch, 3 is the right branch, and will be used!
};

struct LongEdgeRelatedData {
    std::pair<double, int> aver_edge_lp{};
    std::pair<double, int> aver_exact_lp_discrepancy_down{}; //1- exact/lp, number of times
    std::pair<double, int> aver_exact_lp_discrepancy_up{}; //1- exact/lp, number of times
};

struct DualRC {
    std::pair<int, int> edge{};
    double dual1{}, rc1{};
    double dual2{}, rc2{};
};

class MachineLearning {
public:
    static std::unordered_map<std::pair<int, int>, std::pair<double, double>, PairHasher> edge_lp_change;

    static void collectOneSideEdgeFeatures();

    static void updateOneSideLPChange(const std::pair<int, int> &edge, double obj_change, bool if_left);

public:
    static CVRP *cvrp;
    static BbNode *node;
    static double max_edge_cost;
    static double max_mid_point_edge_cord_2_depot;
    static double cluster_coeff;
    static double depot_2_center;
    static std::pair<int, double> average_route_length;
    static std::vector<std::vector<std::pair<double, double> > > mid_point_edge_cord;
    static std::vector<std::vector<double> > mid_point_edge_cord_2_depot;
    static std::vector<std::vector<std::vector<int> > > node_density_in_std_dis_vec_form;
    static std::vector<std::vector<double> > edge_2_other_convert_dis;
    static std::unordered_map<std::pair<int, int>, TmpEdgeRelatedData, PairHasher> edge_tmp_info;
    static std::unordered_map<std::pair<int, int>, LongEdgeRelatedData, PairHasher> edge_long_info;
    static std::vector<double> is_in_solution;
    static std::unordered_map<std::pair<int, int>, double, PairHasher> edge_val;

    static void init(CVRP *pr_cvrp);

    static void updateNode(BbNode *pr_node);

    static void calculatePrerequisites();

    static void recordEdgeLongInfo(int i, int j);

    static void recordDiscrepancyLongInfo(const std::pair<int, int> &edge, double cg_change, bool if_left);

    static void cleanLastData();

    static void getFeatureDataPhase1();

    static void getFeatureDataPhase2();

    static void calculateAverageRouteLength();

    static void calculateClusteringCoefficient();

    static void calculateDisDepot2Center();

    static void collectStaticFeatures();

    static void collectEdgeRelatedFeatures(double org_val);

    static void collectVariableRelatedFeatures(
        const std::pair<int, int> &edge,
        const int *solver_ind,
        int BeforeNumRow,
        int numnz,
        double org_val);

    static void collectResolvingDualRC(DualRC &dual_rc, const std::pair<int, int> &edge,
                                       int BeforeNumRow, bool if_left);

    static void collectResolvingFeatures(const std::vector<EdgeScoreInfo> &edge_info,
                                         const std::vector<DualRC> &dual_rc);

    static void collectScore();

    static void findDiscrepancyResolvingFeatures(const std::pair<int, int> &edge, bool if_left);

    static void debugInputData(const std::pair<std::string, double> &fs);

    static void printFeatures();

    static void chooseBestNCandidate(int num);
};

#define PseudoMark std::string("pseudo-")

#define safe_xgboost(call) {  int Xerr = (call);\
if (Xerr != 0) { \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call + ":" + XGBGetLastError()); \
}\
}

#ifdef DEBUG_ML_INPUT_DATA
#define debug_input_data_call(...) __VA_ARGS__;
#else
#define debug_input_data_call(...)
#endif

#endif //INCLUDE_ML_HPP_
