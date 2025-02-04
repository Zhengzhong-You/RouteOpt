#ifndef CVRP_CVRP_HPP
#define CVRP_CVRP_HPP

#include "branching_node.hpp"
#include "macro.hpp"
#include "config.hpp"
#include "solver_interface/solver.hpp"
#include "capsep.h"
#include "cnstrmgr.h"
#include "instance_data.hpp"
#include <cstdio>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <tuple>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <random>
#include <algorithm>
#include <climits>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <list>
#include <Eigen/Sparse>
#include <limits>
#include <cstdio>
#include <filesystem>
#include "label.hpp"

#ifdef HGS_APPLIED
#include "my_hgs.hpp"
#endif

#ifdef DELUXING_APPLIED
#include "deluxing.hpp"
#endif

using Label =
#ifdef SOLVER_RVRPSTW
	RVRPSTWLabel;
#else
CVRPLabel;
#endif


class CmpNodeValue {
public:
    auto operator()(BbNode *a, BbNode *b) -> bool {
        return a->value > b->value;
    }
};

struct VectorHash {
    size_t operator()(const std::vector<int> &V) const {
        size_t hash = 0;
        for (auto &i: V) {
            hash *= MAX_NUM_CUSTOMERS;
            hash += i;
        }
        return hash;
    }
};

using VecLabel = std::pair<std::vector<Label *>, int>;
using ListLabel = std::list<Label *>;
using QueueTree = std::priority_queue<BbNode *, std::vector<BbNode *>, CmpNodeValue>;
using StackTree = std::stack<BbNode *>;

struct Bucket {
    std::vector<int> bucket_arcs{};
    std::vector<std::pair<double, int> > jump_arcs{};
    int i{};
};

struct MyDoubleHash {
    size_t operator()(double d) const {
        return std::hash<int>()((int) std::round(d / TOLERANCE));
    }
};

struct MyDoubleEqual {
    bool operator()(double a, double b) const {
        return std::abs(a - b) < TOLERANCE;
    }
};

class CVRP {
public:
    yzzLong can_leave_depot_forward{};
    yzzLong can_leave_depot_backward{};

    std::unordered_map<std::vector<int>,
        std::vector<std::vector<std::vector<int> > >, VectorHash> pureMap;

    std::unordered_map<std::vector<int>, std::vector<std::vector<int> >, VectorHash> rank1_multi_mem_plan_map;

    Solver solver{};

    ListLabel **label_array_in_forward_sense{};

    VecLabel **if_exist_extra_labels_in_forward_sense{};

    char if_int_n_feasible{};
    double max_vio{};
    double *demand{};

    int real_dim{};
    int num_col{}, num_row{};

    double lp_val{};
    int MaxNumRoute4Mip{};
    bool if_mip_enumeration_suc{};
    int transformed_number{};
    int idx_node{};
    bool if_exact_labeling_cg{};
    std::vector<double> pi4_labeling{}; //this pi is used for labeling and cutting, need special care!
    int *col_pool4_pricing{};
    int *copy_col_pool4_pricing{}; //used in heuristic find ub!
    size_t pool_beg4_pricing{};
    std::vector<std::vector<int> > ip_opt_sol{};

    double rc_std{};
    double cap{};
    double meet_point_resource_in_bi_dir{}, meet_point_resource_in_bi_dir_enu{};
    ResTuple resource{};
    int dim{}, num_vehicle{}, max_num_vehicle{};
#if  SOLVER_VRPTW == 1
    bool if_force_keep_rcc{};
#endif

    std::vector<std::vector<double> > info_vertex{};
    std::vector<std::vector<double> > cost_mat4_vertex{};
    std::vector<yzzLong> ng_mem4_vertex{};
    std::vector<std::vector<int> > ng_mem4_vertex_sort_dist{};
    std::vector<std::vector<double> > chg_cost_mat4_vertex{};
    std::vector<int> prior_mem_sets{};
    std::vector<ResTuple> lb4_vertex{};
    std::vector<ResTuple> ub4_vertex{};

    std::vector<std::tuple<Label *, Label *, double> > negative_rc_label_tuple{};

    int idx_glo{};
    std::string file_name;

    Label *all_label{};
    size_t label_int_space_len{};
    int *label_int_space{};

    double **arc_graph{};
    double **arc_graph_revised{};
    int num_edge{};
    double num_forward_labels_in_enu{}, num_backward_labels_in_enu{};
    double **rc2_till_this_bin_in_forward_sense{};
    double **rc2_bin_in_forward_sense{};

    std::vector<double> rank1_dual{}; //nonzero

    bool if_exact_cg_finished{};
    double max_enumeration_success_gap{};
    double min_enumeration_fail_gap{1};
    std::pair<double, double> max_bucket_arc_suc_enumeration{}; //forward and backward
    std::pair<double, double> min_bucket_arc_fail_enumeration{1e9, 1e9};
    bool if_force_enumeration_suc{};

    std::vector<int> cstr_index{}; //used in enumeration, delete rows

    size_t mem4_pricing{};


    int max_num_forward_graph_arc{};

    bool if_stop_arc_elimination{};

    std::unordered_map<size_t, int> tell_which_bin4_arc_elimination_in_forward_sense{};

    std::unordered_map<std::pair<int, int>, std::vector<std::pair<Label *, ResTuple> >, PairHasher>
    concatenate_labels_in_backward_cg{};
    int max_num_backward_graph_arc{};
    std::vector<std::vector<ResTuple> > resource_across_arcs_in_backward_sense{};
    std::unordered_map<size_t, int> tell_which_bin4_arc_elimination_in_backward_sense{};
    double **rc2_till_this_bin_in_backward_sense{};
    double **rc2_bin_in_backward_sense{};
    ListLabel **label_array_in_backward_sense{};
    VecLabel **if_exist_extra_labels_in_backward_sense{};

    int max_num_enu_col_pool{};

    double gap_tolerance4_arc_elimination_n_enumeration{};

    bool if_in_enu_state{};

    std::vector<std::vector<ResTuple> > resource_across_arcs_in_forward_sense{};

    int num_buckets_per_vertex{};

    bool if_force_not_regenerate_bucket_graph{};

    double guessed_ub{};

    bool if_enumeration_suc{};

    res_int step_size{};

    double max_labeling_time_for_node_so_far{};
    double max_labeling_time_last_round_for_node_so_far{};

    std::vector<size_t> dump_idx4_left_br_col{}; //record the size_t idx in colPool
    std::vector<size_t> dump_idx4_right_br_col{};

    // std::vector<std::tuple<int, std::vector<int>, std::vector<std::pair<std::vector<int>, int> > > >
    // reset_cut_mem{}; //cut index in r1c  & memory

    std::unordered_map<std::pair<int, int>, std::vector<std::pair<Label *, ResTuple> >, PairHasher>
    concatenate_labels_in_forward_cg{};

    size_t label_assign{};
    size_t route_in_pricing_assign{};
    double aver_route_length{};
    bool if_arc_elimination_succeed{};
    bool if_arc_elimination_tried_but_failed{};

    bool if_roll_back{};
    bool if_tail_off{};
    bool if_short_memory{};

    bool force_not_rollback{true};

    std::pair<double, int> success_enumeration_gap{};

    double time_resolve_lp4_iter{};

    double time_pricing4_iter{};

    double old_ub{};


    double round_up_tolerance{};

    bool if_can_arc_elimination_by_exact_cg{};

    bool final_decision_4_arc_elimination{};

    bool final_decision_4_enumeration{};


    std::unordered_map<int, std::vector<std::vector<int> > > recorded_combinations4_r1c{};


    double num_dominance_checks{};

    std::unordered_map<yzzLong, std::unordered_set<int> > cut_record{};


    double opt_gap{};

    std::pair<double, int> ratio_dominance_checks_non_dominant{}; //sum, int

    sparseRowMatrixXd row_basic_matrix;

    bool if_reset_label_ptr{};
    bool if_clean_all_ptr_list{};
    bool if_clean_concatenate{};

    bool if_allow_change_col{true};

    std::vector<SequenceInfo> new_cols{};
    std::vector<R1CUseStates> cg_v_cut_map{};
    std::vector<std::vector<std::vector<int> > > cg_v_v_use_states{}; //0 remember ;-1 forget; >0 add
    std::vector<int> cg_r1c_denominator{};

    std::vector<int> lp_r1c_map; //map the cut idx to row idx
    std::vector<std::pair<R1CINDEX, std::vector<int> > > lp_v_cut_map{}; //only record > 0 and sparse rep
    std::vector<std::vector<std::vector<int> > > lp_v_v_use_states{}; //0 remember ;-1 forget; >0 add
    std::vector<int> lp_r1c_denominator{};


    double arc_elimination_time{};
    double gap_improved_4_arc_elimination_n_enumeration{};
    std::vector<double> optimal_dual_vector{};
    Solver rollback_solver{};
    std::vector<SequenceInfo> rollback_cols{};

#ifdef CHECK_PRICING_LABELS
    std::pair<double, double> inner_bin_len{};
    std::pair<double, double> outer_bin_len{};
    std::pair<double, double> outer_bin_but_keep_len{};
    std::map<std::vector<int>, double> seq_rc{};
#endif

    void adjustEnumerationStdBucketArcs(BbNode *node, bool if_suc);

    void getRank1DualsInCG(BbNode *node, const std::vector<double> &pi_vector);

    void getVCutMapLP(BbNode *node);

    void giveMemInNode(BbNode *node,
                       const std::vector<R1c> &full_cuts,
                       std::vector<int> &idx
    ) const;

    void resizePoolWarning(size_t &pricing_warning);

    void applyRCF(BbNode *node, int round, bool if_verbose);

    void checkBucketIJ(BbNode *node, int i, int j) const;

    void deleteArcByFalseBranchConstraint(Bucket **buckets, const std::pair<int, int> &edge) const;

    void regenerateGraphBucket(BbNode *node);


    [[nodiscard]] double calculateGapImprovement(double nowVal, double b4Val) const;

    void readSolutionFile(bool if_force);

    void setSolverEnv();

    [[nodiscard]] double getGapStdTryEnumeration() const;

    void adjustEnumerationStdGap(double local_gap, bool if_suc);

    void separateNAddRank1Cuts(BbNode *node);

    void addR1CAtOnceInEnum(BbNode *node, const std::vector<R1c> &new_cuts);

    void generateN4ColpoolMapInEnum(BbNode *node,
                                    std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> &colmap);

    void addBranchConstraint2ColPoolInEnumByColMap(BbNode *node,
                                                   const std::pair<int, int> &edge) const;

    void findNGMemorySets(BbNode *node, bool &if_empty);

    void deleteColumnByNGMemory(BbNode *node, int start, bool if_full_mem = false);

    void augmentNonGuillotineRound(BbNode *node);

    void changeEnumMatByCuts(BbNode *node);

    void printCutsInformation(BbNode *node) const;

    virtual void getLowerBoundofMinimumNumberCars();


    template<bool dir, bool if_symmetry>
    int enumerateHalfwardRoutes(BbNode *node,
                                std::unordered_map<yzzLong, std::tuple<Label *, Label *, double> > &Tags,
                                std::vector<Label *> **copy_bucket,
                                int &num_routes_now);


    void addR1C(BbNode *node, const std::vector<int> &cut, const std::set<int> &mem, int cut_index);

    void addR1CInEnum(BbNode *node, const std::vector<int> &cut);


    void addR1CAtOnce(BbNode *node,
                      const std::vector<int> &idx);


    explicit CVRP(const InstanceData &instanceData);

    virtual ~CVRP();

    virtual void buildModel();

    void addBranchCutToUnsolved(BbNode *node, const std::pair<int, int> &info);

    void separateHybridCuts(BbNode *&node);

    void separateRCCs(BbNode *&node);

    void generateRCCs(BbNode *node);

    void solveLPInLabeling(BbNode *node,
                           bool if_open_heur = true,
                           bool if_open_exact = true,
                           bool if_record_sol = true);

    virtual void priceLabeling(BbNode *node, const std::vector<double> &pi_vector);

    void priceRCC(BbNode *node, const std::vector<double> &pi_vector);

    void priceBRC(BbNode *node, const std::vector<double> &pi_vector);

    virtual bool tellIfInt(BbNode *node, const std::vector<double> &X); //node is very essential

    void rollbackEasyWay(BbNode *node, int old_num);

    void assignMemory();

    void recordOptimalColumn(BbNode *node, bool if_force_rewrite = false);

    int optimizeLPForOneIteration(BbNode *node, double prior_value, int lp_method = SOLVER_PRIMAL_SIMPLEX);

    void cleanIndexColForNode(BbNode *node, bool if_only_rcfixing);

    void constructMap(BbNode *node, int beg) const;

    void getNewConstraintCoefficientByEdge(BbNode *node,
                                           const std::pair<int, int> &edge,
                                           std::vector<int> &ind,
                                           std::vector<double> &val);

    void getCoefficientRCC(BbNode *node, Rcc &rcc, std::vector<int> &ind, std::vector<double> &val);

    template<bool if_symmetry>
    int generateColsByBidir(BbNode *node);

    void addColumns(BbNode *node, int &ccnt, bool if_check_rc = true);

    template<bool if_symmetry>
    int generateColumnsByLighterHeuristic(BbNode *node);

    template<bool if_symmetry>
    int generateColumnsByHeavierHeuristic(BbNode *node);

    double calculateOptimalGap(BbNode *node) const;

    virtual void eliminateArcs(BbNode *node);

    template<bool dir>
    void populateRC2TillThisBinNRC2Bin();

    void initializeLabels(BbNode *node, int mode,
                          bool if_resetLabelPoint, std::tuple<bool, int, bool> control_cleanAllPtr);

    virtual void initializeExtraLabels(Label *label, bool if_forward) {
    };

    void terminateByMIP(BbNode *node);

    virtual void enumerateMIP(BbNode *&node);

    void cleanAllPointers(BbNode *node, int mode, bool if_clear_concatenate);


    void solveLPByInspection(BbNode *node, bool if_only_need_value,
                             bool if_heuristic, bool if_record_sol);

    void deleteBranchCutsAndR1C1s(BbNode *node);

    int generateColumnsByInspection(BbNode *node,
                                    bool if_only_need_value);

    void addColumnsByInspection(BbNode *node, const std::vector<int> &Col_added);

    void regenerateEnumMat(BbNode *node, BbNode *node2, bool if_force = false);

    void generateVertex2IndexColsAndEdge2IndexCols(BbNode *node);

    void recoverR1CsInEnum(BbNode *node);

    void recoverRCCsInEnum(BbNode *node);

    void buildRCCInEnuMatrix(BbNode *node,
                             std::vector<Eigen::Triplet<double> > &triplets,
                             int old_num = 0) const;

    void buildAllR1CInEnuMatrix(BbNode *node,
                                std::vector<Eigen::Triplet<double> > &triplets, int oldNum = 0);

    void createBasicMatrix(BbNode *node) const;

    void rmColByBranchInEnuMatrix(BbNode *node, std::vector<bool> &deleted_columns_in_enumeration_pool,
                                  bool if_test_not_use,
                                  const std::vector<Brc> &brcs) const;

    void optimizeLPForOneIterationInEnum(BbNode *node);

    void addBranchCutToUnsolvedInEnu(BbNode *node, const std::pair<int, int> &info);

    void reviseEnumColInfoByBrC(BbNode *node, BbNode *out_node, const Brc &bf);

    void cleanColumnsInLPByNewBranchCut(BbNode *node, const Brc &bf);

    void generateRCCsInEnum(BbNode *node);

    void deleteNonActiveCutsSafely(BbNode *node, int old_num, bool &if_give_up_this_cutting_round);

    void writeColumnsInPricingPool();

    void determineIfArcElimination(BbNode *node);

    template<bool dir, bool if_last_half, bool if_complete, bool if_symmetry, bool if_std_optgap, bool if_res_updated,
        int heuristic_level>
    void updateLabel(BbNode *node, const ResTuple &res, Label *ki, int i, int j, int &bj,
                     bool &if_suc);

    template<bool dir, int heuristic_level>
    void doDominance(Label *ki, int j, int bj, bool &if_suc);

    template<typename T, bool dir, bool if_symmetry>
    void concatenateTestKernelInArcElimination(int i,
                                               int b,
                                               const std::vector<T> &arc,
                                               int dim_sq,
                                               bool *stateBetween2Buckets,
                                               int *latest_bucket);

    template<bool dir, bool if_symmetry, bool if_std_optgap>
    void concatenateOneLabelWithOtherLabels(Label *ki, int j, int arr_bj, double tmp_rc, const ResTuple &tmp_res,
                                            int &if_state);

    template<char type>
    [[nodiscard]] bool tellResTupleRelations(const ResTuple &res1, const ResTuple &res2) const;

    template<bool dir, bool if_symmetry>
    void eliminateBucketArcs(BbNode *node,
                             int dim_sq,
                             bool *stateBetween2Buckets,
                             int *latest_bucket);

    template<bool dir, bool if_symmetry>
    void eliminateBuketArc4Depot(BbNode *node);

    template<bool dir, bool if_symmetry>
    void concatenatePhaseInArcElimination(BbNode *node);

    void addPathByRC(double path_rc, Label *ki, Label *kj, int num);

    virtual void determineIfEnumeration(BbNode *node);

    template<bool dir>
    void populateTellWhichBin4ArcElimination();

    void getEdgeInfo(BbNode *node, bool if_br) const;

    void deleteNonActiveCutsByDual(BbNode *node, bool if_rcc_by_slack);

    void deleteNewAddedNonActiveCutsBySlack(BbNode *node, int olNum, bool if_keep_rcc);

    void initializeBucketGraphForNode(BbNode *node);

    [[nodiscard]] int checkMaximumNumberCustomers() const;

    [[nodiscard]] static int checkMaximumNumberR1Cs(int num_r1cs);

    [[nodiscard]] static int checkIfAssignGreaterMAXINT(const size_t &value, const std::string &str);

    void reallocateLabel();

    virtual void assignMemory4Label() {
    };

    void reallocatePricingPool(size_t num = 0);

    [[nodiscard]] int checkPricingPool() const;

    void rollbackToPreviousState(BbNode *node, const std::vector<Rcc> &old_rcc,
                                 const std::vector<R1c> &old_r1c,
                                 const std::vector<Brc> &old_brc, const std::vector<size_t> &old_col_idx);

    static int inverseLastBranchConstraint(char sense, double rhs, Solver &solver);

    static int addBranchConstraint(std::vector<int> &cind,
                                   std::vector<double> &cval,
                                   char sense,
                                   double rhs,
                                   const char *constrname,
                                   Solver &local_solver);

    void changeBranchConstraint(const std::vector<int> &vind,
                                const std::vector<double> &vval,
                                char sense,
                                double rhs,
                                int row_idx,
                                Solver &local_solver);

    void initializeBucketGraph();

    void checkIfCutsLegal(BbNode *node) const;

    void initialProcessing();

    virtual double transformCost(double x);

    void solveMIP(BbNode *node, bool if_inEnu);

    [[nodiscard]] double ceilTransformedNumberRelated(double x) const;

    void initializeLabels();

    template<typename T, bool dir, bool if_last_half, bool if_complete, bool if_symmetry, int heuristic_level>
    int extendKernel4Exact(BbNode *node, Label *ki,
                           int i,
                           ResTuple res,
                           const std::vector<T> &arc);

    template<bool dir>
    void obtainjumpArcs(BbNode *node, std::bitset<2> **bitMap) const;

    template<typename T>
    int extendKernel4Exact_Backward(Label *&ki,
                                    int i,
                                    double res,
                                    const std::vector<T> &arc);

    bool decreaseMainResourceConsumption(const ResTuple &nowResource,
                                         ResTuple &newResource,
                                         int start,
                                         int end);

    template<bool dir, int heuristic_level>
    void checkIfDominated(Label *&ki, int i, int b,
                          bool &if_suc);

    template<bool dir, bool if_last_half, bool if_complete, bool if_symmetry, int heuristic_level>
    void runLabeling(BbNode *node, bool if_default = true);

    void obtainBackwardJumpArcs(BbNode *node, std::bitset<2> **bitMap) const;

    int enumerateHalfBackwardRoutes(BbNode *node,
                                    std::vector<Label *> **copy_bucket);

    void assignInitialLabelingMemory() const;

    int concatenateRoutesPriorForwardInEnumeration(BbNode *node,
                                                   std::unordered_map<yzzLong,
                                                       std::tuple<Label *, Label *, double> > &Tags,
                                                   int &num_routes_now
    );

    template<typename T, bool dir>
    void extendKernel4LightHeur(BbNode *node, Label *&ki,
                                int i,
                                ResTuple res,
                                const std::vector<T> &arc,
                                bool if_return_to_depot);

    template<typename T, bool dir>
    void extendKernel4HeavierHeur(BbNode *node, Label *&ki,
                                  int i,
                                  ResTuple res,
                                  const std::vector<T> &arc,
                                  bool if_return_to_depot);

    int lastHalfForwardInArcElimination(BbNode *node);

    int forwardConcatenateInArcElimination();

    void eliminateBucketArcs(BbNode *node,
                             int dim_sq,
                             bool *stateBetween2Buckets,
                             int *latest_bucket);

    template<typename T>
    int extendKernel4ArcElimination_last_half_Forward(Label *&ki,
                                                      int i,
                                                      double res,
                                                      const std::vector<T> &arc);

    bool runColumnGenerationType(BbNode *node, int mode);

    void findNonActiveCuts(BbNode *node, bool if_optimal_dual = true);

    void deleteNonactiveCuts(BbNode *node, std::vector<int> &nonactive_cuts);

    void runHalfForwardLabeling(BbNode *node);

    bool increaseMainResourceConsumption(const ResTuple &nowMainResource,
                                         ResTuple &newMainResource,
                                         int start,
                                         int end);

    template<bool if_symmetry>
    int concatenateCols_prior_forward(BbNode *node);

    void runLabelingForArcElimination(BbNode *node);

    void eliminateBucketArcs(BbNode *node);

    void obtainJumpArcs(BbNode *node) const;

    virtual bool enumerateRoutes(BbNode *node);

    virtual void setResourceInBucketGraph();

    virtual void checkSolutionFeasibleByCapacity(bool &feasible) {
    };

    virtual void popArcGraph(BbNode *node) {
    };

    virtual void cleanColumnsCapInfeasible(BbNode *node) {
    };

    void tellIfEnterMIP(BbNode *node);

    void printBrDecisions();

    static void printInfoLabeling(int iter, int num_added_col, int num_col,
                                  int num_row,
                                  double mt, double spt, double et,
                                  double lp, double lb_val, double ub);

    void printOptIntSol();

    /**
     * new supply functions
     */

    void rmLPCols(BbNode *node, const std::vector<int> &col_idx);

    void updateEdgeColMap(BbNode *node, bool if_br);

    void updateIPOptSol(BbNode *node, const std::vector<double> &X);

    void getFullCoeff(const std::vector<std::pair<std::vector<int>, int> > &cuts,
                      const std::vector<std::vector<int> > &routes,
                      Eigen::MatrixXd &dense_coeff);

    template<typename MatrixType>
    void getLimitedR1CCoeffs(const std::vector<SequenceInfo> &seq_info, MatrixType &mat);

    void getLimitedR1CPre(BbNode *node, const std::vector<int> &idx);

    void addLimitedMemoryR1CsNodeBased(BbNode *node,
                                       const std::vector<R1c> &full_cuts);


    void updateR1CStates(double &rc, R1CPricingStat &out_states,
                         const R1CPricingStat &in_states,
                         int from, int to);

    virtual void updateExtraStates(Label *p_label, bool if_forward) {
    };

    virtual void lighterHeuristic(BbNode *node, int &num);

    virtual void heavierHeuristic(BbNode *node, int &num);

    virtual void exactLabeling(BbNode *node, int &num);

    template<bool dir>
    bool doRCTermDominance(Label *ki, Label *kj);

    virtual bool doExtraDominance(double &gap, Label *ki, Label *kj, bool if_forward) { return true; }

    bool doR1CDominance(double &gap,
                        const R1CPricingStat &out_states,
                        const R1CPricingStat &in_states);

    bool concatenateR1CStates(double &rc, double req, R1CPricingStat &out_states,
                              const R1CPricingStat &in_states, int out, int in);

    void getCoefficientExtendR1C(std::vector<int> &states,
                                 std::vector<int> &sparse_rep,
                                 std::unordered_map<int, int> &cnt,
                                 int &valid_sparse_num,
                                 int from,
                                 int to
    );


    void writeIntoNode(BbNode *node, std::vector<R1c> &full_cuts, std::vector<int> &idx) const;

    void addLimitedMemoryR1Cs(BbNode *node, std::vector<R1c> &full_cuts);

    std::vector<double> revised_rank1_dual;

    std::vector<int> cg_cut_map; //map back the cg cut to lp cut!

    template<bool dir>
    void getTopologicalOrder4OneBin(BbNode *node, int b);

    void getTopologicalOrder(BbNode *node);

    template<bool dir>
    void sortLabelsInBinByRC(int i, int b);

    template<bool dir, int heuristic_level>
    bool dominanceCore(Label *ki, Label *kj);

    template<bool dir>
    void updateEnumerationLabel(Label *ki, int i, int j, int &bj);

    template<bool dir>
    bool dominanceCoreInEnumeration(Label *ki, Label *kj);

    template<bool dir>
    void doDominanceEnumerationLabel(Label *ki, int i, int j, int bj, bool &if_suc);

    template<bool dir, bool if_symmetry>
    void extendKernel4Enumeration(BbNode *node, int i, int b, Label *ki, std::vector<Label *> **copy_bucket,
                                  std::unordered_map<yzzLong, std::tuple<Label *, Label *, double> > &Tags,
                                  int &num_routes_now, int &status);

    virtual bool checkFeasibleExtensionInUpdateLabelByBackwardCompletionBounds(Label *ki, int j,
                                                                               double which_rc) { return true; };

    virtual void tryCallHeuristicB4Branching(BbNode *&node) {
    };

    virtual void getNewRhs4ChangingDual(BbNode *node, std::vector<double> &new_model_rhs);

    virtual bool stopChangingDualCondition(BbNode *node);

    void calculateColumnCoefficientsB4Enumeration(BbNode *node, sparseColMatrixXd &mat,
                                                  Eigen::RowVectorXd &cost,
                                                  std::unordered_map<std::pair<int, int>,
                                                      std::vector<int>,
                                                      PairHasher> &edge_map);

    void prepareSolveEnuTree(BbNode *node);

    bool &getIfInEnuState();

    void startSolveNode(BbNode *node);

    void postSolveNode(BbNode *&node);

    void cuttingAfterEnu(BbNode *&node);

    void enterCuttingPhase(BbNode *&node);

    void finishSolveNode();

    int &getNumCol();

    int &getNumRow();

    double &getMaxEnumerationSuccessGap();

    void preprocess4TestCG(BbNode *node);

    void postprocess4TestCG(BbNode *node);

    bool &getIfAllowChangeCol();
};

bool operator==(const Rcc &lhs, const Rcc &rhs);

auto CmpLabelRCLess(const Label *l1, const Label *l2) -> bool;

res_int roundAndConvertResLong(double value);

std::string readInstanceFile(const std::string &file_name, int line);

std::string generateInstancePath(int argc, char *argv[]);

double pow_self(double x, int n);

void self_mkdir(const std::string &path);
#endif
