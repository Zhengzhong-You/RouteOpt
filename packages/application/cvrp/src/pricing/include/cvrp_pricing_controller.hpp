/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_CVRP_PRICING_CONTROLLER_HPP
#define ROUTE_OPT_CVRP_PRICING_CONTROLLER_HPP
#include <list>
#include "cvrp_macro.hpp"
#include "label.hpp"
#include "rank1_rc_controller.hpp"
#include "rank1_data_shared.hpp"
#include "bucket.hpp"
#include "pricing_macro.hpp"
#include "cuts_definition.hpp"


namespace RouteOpt::Application::CVRP {
    using VecLabel = std::pair<std::vector<Label *>, int>;
    using ListLabel = std::list<Label *>;

    struct StatCount {
        std::pair<double, double> average{};

        double getAverage() const {
            if (average.second == 0) return 0;
            return average.first / average.second;
        }

        void updateAverage(double val) {
            average.first += val;
            ++average.second;
        }
    };

    class CVRP_Pricing {
    public:
        CVRP_Pricing(const int &dim,
                     const int &max_num_vehicle,
                     double cap,
                     const std::vector<double> &demand,
                     double H,
                     const std::vector<double> &earliest_time,
                     const std::vector<double> &latest_time,
                     const std::vector<double> &service_time,
                     std::vector<std::vector<double> > &cost_mat4_vertex,
                     Rank1Cuts::RCGetter::Rank1RCController &rank1_rc_controller)
            : dim(dim),
              max_num_vehicle_ref(max_num_vehicle),
              cost_mat4_vertex_ref(cost_mat4_vertex),
              rank1_rc_controller_ref(rank1_rc_controller) {
            chg_cost_mat4_vertex.resize(dim, std::vector<double>(dim));
            getNG();
            setCapResourceInBucketGraph<!IF_SYMMETRY_PROHIBIT>(cap, demand);
            setTWResourceInBucketGraph(H, earliest_time, latest_time, service_time, cost_mat4_vertex);
        }

        void initLabelingMemory() {
            assignMemory();
            initializeBucketGraph<!IF_SYMMETRY_PROHIBIT>();
            initializeLabels<!IF_SYMMETRY_PROHIBIT>();
        }


        void updatePtr(Bucket **all_forward_buckets, Bucket **all_backward_buckets,
                       std::vector<std::vector<std::vector<int> > > *topological_order_forward_ptr,
                       std::vector<std::vector<std::vector<int> > > *topological_order_backward_ptr);


        template<bool if_symmetry>
        void initializeBucketGraphForNode(
            Bucket **&all_forward_buckets, Bucket **&all_backward_buckets,
            int &num_forward_bucket_arcs, int &num_backward_bucket_arcs);


        void initializeOnceB4WholePricing();

        template<bool if_symmetry>
        int generateColumnsByExact(double time_limit);

        template<bool if_symmetry>
        void adjustResourceMeetPointInPricing();

        void setTerminateMarker(double val, double ub, bool &if_terminate);

        template<bool if_symmetry, PRICING_LEVEL pricing_level>
        int generateColumnsByHeuristic();

        //getters

        [[nodiscard]] bool getIfCompleteCG() {
            return if_exact_cg_finished;
        }

        const auto &getStepSize() const {
            return step_size;
        }


        const auto &getNewCols() const {
            return new_cols;
        }

        const auto &getNegativeLabelTuple() const {
            return negative_rc_label_tuple;
        }

        const auto &getNGMem() const {
            return ng_mem4_vertex;
        }

        const auto &getNumBucketPerVertex() const {
            return num_buckets_per_vertex;
        }

        const auto &getResourceForward() const {
            return resource_across_arcs_in_forward_sense;
        }

        const auto &getResourceBackward() const {
            return resource_across_arcs_in_backward_sense;
        }


        auto getColumnPoolPtr() const {
            return col_pool4_pricing;
        }

        const auto &getExistLabelsInForward() const {
            return if_exist_extra_labels_in_forward_sense;
        }

        const auto &getExistLabelsInBackward() const {
            return if_exist_extra_labels_in_backward_sense;
        }

        const auto &getMaxEnumerationSuccessGap() const {
            return max_enumeration_success_gap;
        }

        const auto &getSuccessEnumerationGap() const {
            return success_enumeration_gap;
        }

        const auto &getMinEnumerationFailGap() const {
            return min_enumeration_fail_gap;
        }

        const auto &getMaxGap2TryEnumeration() const {
            return max_gap2_try_enumeration_enumeration;
        }

        const auto &getStabGamma() const {
            return stab_gamma;
        }

        const auto &getStabDelta() const {
            return stab_delta;
        }

        const auto &getIncumbentDualSolution() const {
            return incumbent_dual_solution;
        }

        const auto &getTimeLP() const {
            return time_lp;
        }

        const auto &getTimePricing() const {
            return time_pricing;
        }

        //refers

        auto &refSeqRCMap() {
            return seq_rc;
        }

        auto &refNG() {
            return ng_mem4_vertex;
        }

        auto &refExistLabelsInForward() {
            return if_exist_extra_labels_in_forward_sense;
        }

        auto &refExistLabelsInBackward() {
            return if_exist_extra_labels_in_backward_sense;
        }

        auto &refNumBucketPerVertex() {
            return num_buckets_per_vertex;
        }

        auto &refStepSize() {
            return step_size;
        }

        auto &refColumnPoolPtr() {
            return col_pool4_pricing;
        }

        auto &refMem4Pricing() {
            return mem4_pricing;
        }

        auto &refPoolBeg4Pricing() {
            return pool_beg4_pricing;
        }

        auto &refMaxEnumerationSuccessGap() {
            return max_enumeration_success_gap;
        }

        auto &refSuccessEnumerationGap() {
            return success_enumeration_gap;
        }

        auto &refMinEnumerationFailGap() {
            return min_enumeration_fail_gap;
        }

        auto &refMaxGap2TryEnumeration() {
            return max_gap2_try_enumeration_enumeration;
        }

        auto &refNewCols() {
            return new_cols;
        }

        auto &refStabGamma() {
            return stab_gamma;
        }

        auto &refStabDelta() {
            return stab_delta;
        }

        auto &refIncumbentDualSolution() {
            return incumbent_dual_solution;
        }

        auto &refTimeLP() {
            return time_lp;
        }

        auto &refTimePricing() {
            return time_pricing;
        }

        auto &refNumColGeneratedUB() {
            return num_col_generated_ub;
        }

        auto &refAverageRouteLength() {
            return aver_route_length;
        }

        //

        void priceConstraints(const std::vector<Rcc> &rccs,
                              const std::vector<R1c> &r1cs,
                              const std::vector<Brc> &brcs,
                              const std::vector<double> &pi_vector);

        template<bool if_symmetry>
        void considerRegenerateBucketGraph(
            Bucket **&all_forward_buckets, Bucket **&all_backward_buckets,
            int &num_forward_bucket_arcs, int &num_forward_jump_arcs,
            int &num_backward_bucket_arcs, int &num_backward_jump_arcs);


        template<bool dir>
        void populateTellWhichBin4ArcElimination();

        template<bool if_symmetry>
        void getTopologicalOrder();

        template<bool if_symmetry>
        void eliminateArcs(
            const std::vector<Rcc> &rccs,
            const std::vector<R1c> &r1cs,
            const std::vector<Brc> &brcs,
            const std::vector<double> &optimal_dual_vector,
            double ub,
            double opt_gap,
            double &last_gap,
            int &num_forward_bucket_arcs,
            int &num_backward_bucket_arcs,
            int &num_forward_jump_arcs,
            int &num_backward_jump_arcs
        );

        template<bool if_symmetry>
        bool enumerateMIP(const std::vector<Rcc> &rccs,
                          const std::vector<R1c> &r1cs,
                          const std::vector<Brc> &brcs,
                          const std::vector<double> &optimal_dual_vector,
                          double ub,
                          double opt_gap,
                          int num_forward_bucket_arcs,
                          int num_backward_bucket_arcs,
                          bool &if_in_enu_state,
                          RowVectorXT &index_columns_in_enumeration_column_pool,
                          RowVectorXd &cost_for_columns_in_enumeration_column_pool,
                          int &valid_size,
                          bool if_fix_meet_point,
                          const std::vector<double> &optional_demand_testifier,
                          double optional_cap_testifier);


        double getAverageRouteLength();

        double getSmallestRC();

        CVRP_Pricing() = delete;

        ~CVRP_Pricing() {
            freeMemory();
        }

    private:
        //private
        int dim{};
        int initial_ng_size{};
        Label *all_label{};
        bool if_exact_labeling_cg{};
        ListLabel **label_array_in_forward_sense{};
        VecLabel **if_exist_extra_labels_in_forward_sense{};
        int *col_pool4_pricing{};
        int *copy_col_pool4_pricing{}; //used in heuristic std::find ub!
        size_t pool_beg4_pricing{};
        std::vector<std::vector<double> > chg_cost_mat4_vertex{};
        std::unordered_map<int, double> adjust_brc_dual4_single_route{}; //customer, dual
        std::vector<std::vector<Resource> > resource_across_arcs_in_forward_sense{};
        std::vector<std::vector<Resource> > resource_across_arcs_in_backward_sense{};
        Resource resource{};
        double meet_point_resource_in_bi_dir{};
        std::vector<Resource> lb4_vertex{};
        std::vector<Resource> ub4_vertex{};
        double **rc2_till_this_bin_in_forward_sense{};
        double **rc2_bin_in_forward_sense{};

        double **rc2_till_this_bin_in_backward_sense{};
        double **rc2_bin_in_backward_sense{};

        ListLabel **label_array_in_backward_sense{};
        VecLabel **if_exist_extra_labels_in_backward_sense{};

        res_int step_size{};

        int num_buckets_per_vertex{InitialNumBuckets};
        std::vector<routeOptLong> ng_mem4_vertex{};
        routeOptLong can_leave_depot_forward{};
        routeOptLong can_leave_depot_backward{};

        size_t label_assign{};

        size_t mem4_pricing{};

        StatCount aver_route_length{};
        std::unordered_map<size_t, int> tell_which_bin4_arc_elimination_in_forward_sense{},
                tell_which_bin4_arc_elimination_in_backward_sense{};
        size_t max_num_forward_graph_arc{}, max_num_backward_graph_arc{};
        double NumExistedLabels{};
        double NumExistedLabel_back{};

        //ptrs

        Bucket **all_forward_buckets{};
        Bucket **all_backward_buckets{};
        std::vector<std::vector<std::vector<int> > > *topological_order_forward_ptr{};
        std::vector<std::vector<std::vector<int> > > *topological_order_backward_ptr{};


        //refers
        std::reference_wrapper<const int> max_num_vehicle_ref;
        std::reference_wrapper<std::vector<std::vector<double> > > cost_mat4_vertex_ref;
        std::reference_wrapper<Rank1Cuts::RCGetter::Rank1RCController> rank1_rc_controller_ref;

        //initial each cg
        int idx_glo{};
        double rc_std{};
        double opt_gap{};
        std::unordered_map<std::pair<int, int>, std::vector<std::pair<Label *, Resource> >, PairHasher>
        concatenate_labels_in_forward_cg{};

        std::unordered_map<std::pair<int, int>, std::vector<std::pair<Label *, Resource> >, PairHasher>
        concatenate_labels_in_backward_cg{};
        std::vector<std::tuple<Label *, Label *, double> > negative_rc_label_tuple{};
        double num_dominance_checks{};
        bool if_short_memory{};
        bool if_exact_cg_finished{};
        bool if_exact_labeling_finished{};
        bool if_arc_elimination_succeed{};
        double arc_elimination_time{HardTimeThresholdInArcEliminationLastHalf};
        int num_col_generated_ub{};
        std::vector<SequenceInfo> new_cols{};
        //initial once in whole pricing
        std::pair<double, int> ratio_dominance_checks_non_dominant{}; //sum, int

        //arc elimination
        double gap_improved_4_arc_elimination_n_enumeration{GapImproved4ArcEliminationNEnumeration};
        double gap_tolerance4_arc_elimination_n_enumeration{InitGapTolerance4ArcEliminationNEnumeration};
        bool if_stop_arc_elimination{};
        double old_ub{std::numeric_limits<float>::max()};

        //enumeration
        double num_forward_labels_in_enu{}, num_backward_labels_in_enu{};
        double max_gap2_try_enumeration_enumeration{InitialMaxGap2TryEnumeration};
        double meet_point_resource_in_bi_dir_enu{};
        double max_enumeration_success_gap{};
        double min_enumeration_fail_gap{1.};
        bool if_force_enumeration_suc{};
        int max_label_in_enumeration{};
        int max_route_in_enumeration{};
        std::vector<int> label_per_bin_in_enumeration{};
        res_int step_size_label_check_in_enumeration{};

        double stab_gamma{};
        double stab_delta{};
        StatCount time_lp{};
        StatCount time_pricing{};
        std::pair<std::vector<double>, double> incumbent_dual_solution{}; // dual sol & lb

        std::pair<double, int> success_enumeration_gap{};
        std::pair<double, double> max_bucket_arc_suc_enumeration{}; //forward and backward
        std::pair<double, double> min_bucket_arc_fail_enumeration{
            std::numeric_limits<float>::max(), std::numeric_limits<float>::max()
        };

        //monitor pricing data
        std::pair<double, double> inner_bin_len{};
        std::pair<double, double> outer_bin_len{};
        std::pair<double, double> outer_bin_but_keep_len{};
        std::map<std::vector<int>, double> seq_rc{};

        //enumeration
        template<bool if_symmetry>
        bool determineIfEnumeration(double ub, double opt_gap, int num_forward_bucket_arcs,
                                    int num_backward_bucket_arcs);

        template<bool dir>
        dominanceCoreInEnumeration_STATE dominanceCoreInEnumeration(Label *ki, Label *kj);

        template<bool dir>
        void doDominanceEnumerationLabel(int j, int bj, bool &if_suc);

        void updateEnumerationLabel(Label *ki, int i, int j, int &bj);

        template<bool dir, bool if_symmetry>
        void extendKernel4Enumeration(int i, int b, Label *ki,
                                      std::vector<Label *> **copy_bucket,
                                      std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &
                                      Tags,
                                      int &num_routes_now, ENUMERATION_STATE &status);

        template<bool dir, bool if_symmetry>
        ENUMERATION_STATE enumerateHalfwardRoutes(
            std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &Tags,
            std::vector<Label *> **copy_bucket,
            int &num_routes_now);

        template<bool if_symmetry>
        void checkGroupInner(
            ListLabel &label_arr,
            double &tmp_rc, Resource &tmp_Resource,
            Label *ki, int i, int j,
            std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &Tags,
            int &num_routes_now, ENUMERATION_STATE &status
        );

        template<bool if_symmetry>
        ENUMERATION_STATE concatenateRoutesPriorForwardInEnumeration(
            std::unordered_map<routeOptLong, std::tuple<Label *, Label *, double> > &Tags,
            int &num_routes_now);

        template<bool if_symmetry>
        bool enumerateRoutes(
            bool &if_in_enu_state,
            ENUMERATION_STATE &status,
            RowVectorXT &index_columns_in_enumeration_column_pool,
            RowVectorXd &cost_for_columns_in_enumeration_column_pool,
            int &valid_size,
            const std::vector<double> &optional_demand_testifier,
            double &optional_cap_testifier);

        //end

        template<bool if_symmetry>
        void obtainJumpArcs(int &num_forward_jump_arcs, int &num_backward_jump_arcs);

        template<bool dir>
        void obtainJumpArcs(std::bitset<2> **bitMap,
                            int &num_forward_jump_arcs,
                            int &num_backward_jump_arcs);

        template<bool if_symmetry>
        void eliminateBucketArcs(int &num_forward_bucket_arcs, int &num_backward_bucket_arcs);

        template<bool dir, bool if_symmetry>
        void eliminateBucketArcs(
            int dim_sq,
            bool *stateBetween2Buckets,
            int *latest_bucket,
            int &num_forward_bucket_arcs,
            int &num_backward_bucket_arcs);

        template<bool dir, bool if_symmetry>
        void eliminateBuketArc4Depot();

        template<typename T, bool dir, bool if_symmetry>
        void concatenateTestKernelInArcElimination(int i,
                                                   int b,
                                                   const std::vector<T> &arc,
                                                   int dim_sq,
                                                   bool *stateBetween2Buckets,
                                                   int *latest_bucket);

        template<bool if_symmetry>
        void runLabelingForArcElimination();

        bool determineIfArcElimination(double ub, double opt_gap, double &last_gap);

        template<bool dir, bool if_symmetry>
        void concatenatePhaseInArcElimination();

        void freeMemory() const;

        template<bool if_symmetry>
        void initializeBucketGraph();

        void assignMemory();

        void getNG();

        template<bool if_symmetry>
        void setCapResourceInBucketGraph(double cap, const std::vector<double> &demand);

        void setTWResourceInBucketGraph(double H, const std::vector<double> &e,
                                        const std::vector<double> &l, const std::vector<double> &s,
                                        const std::vector<std::vector<double> > &cost_mat4_vertex);

        template<bool dir>
        void getTopologicalOrder4OneBin(int b);

        template<bool if_symmetry>
        void updateDominanceStatics();

        template<bool if_symmetry>
        void regenerateGraphBucket(
            Bucket **&all_forward_buckets, Bucket **&all_backward_buckets,
            int &num_forward_bucket_arcs, int &num_forward_jump_arcs,
            int &num_backward_bucket_arcs, int &num_backward_jump_arcs);

        void pricePartitioning(const std::vector<double> &pi_vector);

        void priceBRC(const std::vector<Brc> &Brcs, const std::vector<double> &pi_vector);

        int checkPricingPool() const;

        bool increaseMainResourceConsumption(const Resource &nowMainResource,
                                             Resource &newMainResource,
                                             int start,
                                             int end);

        bool decreaseMainResourceConsumption(const Resource &nowMainResource,
                                             Resource &newMainResource,
                                             int start,
                                             int end);

        void addPathByRC(double path_rc, Label *ki, Label *kj, int num);

        template<char type>
        bool tellResTupleRelations(const Resource &res1, const Resource &res2) const;

        template<bool dir, bool if_last_half, bool if_complete, bool if_symmetry, bool if_std_optgap, bool
            if_res_updated,
            PRICING_LEVEL pricing_level>
        void updateLabel(const Resource &res, Label *ki, int i, int j, int &bj,
                         bool &if_suc);

        template<bool dir, bool if_symmetry, bool if_std_optgap>
        void concatenateOneLabelWithOtherLabels(Label *ki, int j, int arr_bj, double tmp_rc,
                                                const Resource &tmp_res,
                                                int &if_state);


        template<bool dir>
        bool doRCTermDominance(Label *ki, Label *kj);


        template<bool dir, PRICING_LEVEL pricing_level>
        bool dominanceCore(Label *ki, Label *kj);

        template<bool dir, PRICING_LEVEL pricing_level>
        void doDominance(Label *ki, int j, int bj, bool &if_suc);


        template<bool dir, PRICING_LEVEL pricing_level>
        void checkIfDominated(Label *&ki, int i, int b,
                              bool &if_suc);

        template<typename T, bool dir, bool if_last_half, bool if_complete, bool if_symmetry, PRICING_LEVEL
            heuristic_level>
        PRICING_STATE extendKernel4Exact(Label *ki,
                                         int i,
                                         Resource res,
                                         const std::vector<T> &arc);

        template<bool dir>
        void sortLabelsInBinByRC(int i, int b);

        template<bool dir, bool if_last_half, bool if_complete, bool if_symmetry, PRICING_LEVEL heuristic_level>
        void runLabeling(double time_limit);

        template<bool dir, bool if_symmetry, bool if_reset_label_point, bool if_clear_all, bool if_clear_concatenate>
        void initializeLabels();

        template<bool dir, bool if_clear_concatenate>
        void cleanAllPointers();

        template<bool dir>
        void populateRC2TillThisBinNRC2Bin();

        void writeColumnsInPricingPool();

        void reallocatePricingPool(size_t num = 0);

        template<bool if_symmetry>
        int concatenateCols_prior_forward();

        void reallocateLabel();

        template<bool if_symmetry>
        void initializeLabels();

        void inspectColumns();

        void resizePoolWarning(size_t &pricing_warning);
    };
}

#include "set_resource_in_bucket_graph.hpp"
#include "regenerate_bucket_graph.hpp"
#include "get_topological_order.hpp"
#include "preprocess_arc_elimination.hpp"
#include "pricing_functors.hpp"
#include "run_labeling.hpp"
#include "heuristic_labeling.hpp"
#include "exact.hpp"
#include "write_columns_from_pricing.hpp"
#include "arc_elimination.hpp"
#include "enumeration.hpp"
#endif // ROUTE_OPT_CVRP_PRICING_CONTROLLER_HPP
