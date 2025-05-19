/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_NODE_HPP
#define ROUTE_OPT_NODE_HPP
#include <vector>
#include <deque>
#include "cvrp_macro.hpp"
#include "solver.hpp"
#include "route_opt_macro.hpp"
#include "bucket.hpp"
#include "cuts_definition.hpp"
#include "bbt_controller.hpp"
#include "pricing_macro.hpp"
#include "rank1_coefficient_controller.hpp"
#include "get_rcc_coefficient.hpp"

namespace RouteOpt::Application::CVRP {
    class BbNode {
    public:
        //user defined functions
        static bool defineBetterNode(BbNode *a, BbNode *b) {
            if (a->if_in_enu_state != b->if_in_enu_state)
                return a->if_in_enu_state > b->if_in_enu_state;

            if (!equalFloat(a->value, b->value))
                return a->value < b->value;

            return a->idx < b->idx;
        }

        static auto &refNodeBrIncreaseVal(BbNode *node) {
            return node->br_value_improved;
        }

        static int getNodeIdx(BbNode *node) {
            return node->idx;
        }

        static double getNodeValue(BbNode *node) {
            return node->value;
        }

        static std::tuple<std::pair<int, int>, bool, int> getLastBrc(BbNode *node) {
            if (node->brcs.empty())
                return {std::make_pair(0, 0), false, 0};
            return {node->brcs.back().edge, node->brcs.back().br_dir, static_cast<int>(node->brcs.size())};
        }


        static bool getNodeIfTerminate(BbNode *node) {
            return node->if_terminate;
        }

        static bool getNodeIfInEnumState(BbNode *node) {
            return node->if_in_enu_state;
        }

        static void setDim(int dim) {
            BbNode::dim = dim;
        }

        //user defined functions end

        static void buildModel(int num_vehicle, int dim, Solver *solver, BbNode *node);

        static void regenerateEnumMat(BbNode *node, BbNode *node2, bool if_force, std::vector<double> &duals);

        static std::unordered_map<std::pair<int, int>, double, PairHasher> obtainSolEdgeMap(BbNode *node);


        //refers
        Solver &refSolver() {
            return solver;
        }

        auto &refCols() {
            return cols;
        }

        auto &refRCCs() {
            return rccs;
        }

        auto &refR1Cs() {
            return r1cs;
        }

        auto &refBrCs() {
            return brcs;
        }

        auto &refIndexColPool() {
            return index_columns_in_enumeration_column_pool;
        }

        auto &refCostColPool() {
            return cost_for_columns_in_enumeration_column_pool;
        }

        //


        BbNode() = default;

        ~BbNode();

        // //getters
        [[nodiscard]] int getIdx() const {
            return idx;
        }


        [[nodiscard]] double getValue() const {
            return value;
        }

        [[nodiscard]] bool getIfTerminate() const {
            return if_terminate;
        }

        [[nodiscard]] const std::vector<SequenceInfo> &getCols() const {
            return cols;
        }

        [[nodiscard]] const auto &getBrCs() const {
            return brcs;
        }

        [[nodiscard]] const auto &getR1Cs() const {
            return r1cs;
        }

        [[nodiscard]] const auto &getRCCs() const {
            return rccs;
        }

        [[nodiscard]] const auto &getIfRootNode() const {
            return if_root_node;
        }

        [[nodiscard]] auto getIfInEnumState() const {
            return if_in_enu_state;
        }

        [[nodiscard]] const auto &getIndexColPool() const {
            return index_columns_in_enumeration_column_pool;
        }

        [[nodiscard]] const auto &getCostColPool() const {
            return cost_for_columns_in_enumeration_column_pool;
        }

        [[nodiscard]] const auto &getMatrixColPool() const {
            return matrix_in_enumeration;
        }

        [[nodiscard]] auto getBrValueImproved() const {
            return br_value_improved;
        }

        [[nodiscard]] auto getTreeSize() const {
            return static_cast<int>(brcs.size());
        }

        [[nodiscard]] const auto &getAllForwardBuckets() const {
            return all_forward_buckets;
        }

        [[nodiscard]] const auto &getAllBackwardBuckets() const {
            return all_backward_buckets;
        }

        [[nodiscard]] const auto &getDeletedColumnsInEnumerationColumnPool() const {
            return deleted_columns_in_enumeration_pool;
        }

        [[nodiscard]] auto getColPoolSize() const {
            return static_cast<int>(index_columns_in_enumeration_column_pool.size());
        }


        //refers

        auto &refIfInEnumState() {
            return if_in_enu_state;
        }


        bool &refIfTerminate() {
            return if_terminate;
        }

        auto &refValue() {
            return value;
        }

        auto &refNumForwardBucketArcs() {
            return num_forward_bucket_arcs;
        }

        auto &refNumForwardJumpArcs() {
            return num_forward_jump_arcs;
        }

        auto &refNumBackwardBucketArcs() {
            return num_backward_bucket_arcs;
        }

        auto &refNumBackwardJumpArcs() {
            return num_backward_jump_arcs;
        }

        auto &refLastGap() {
            return last_gap;
        }

        auto &refAllForwardBuckets() {
            return all_forward_buckets;
        }

        auto &refAllBackwardBuckets() {
            return all_backward_buckets;
        }

        auto &refValidSize() {
            return valid_size;
        }

        auto &refIndexColumnsInEnumerationColumnPool() {
            return index_columns_in_enumeration_column_pool;
        }

        auto &refCostForColumnsInEnumerationColumnPool() {
            return cost_for_columns_in_enumeration_column_pool;
        }

        auto &refDeletedColumnsInEnumerationColumnPool() {
            return deleted_columns_in_enumeration_pool;
        }

        auto &refMatrixColPool() {
            return matrix_in_enumeration;
        }

        auto &refIdx() {
            return idx;
        }

        auto &refBasicMatrix() {
            return basic_matrix;
        }

        auto &refRowBasicMatrix() {
            return row_basic_matrix;
        }


        //ptrs

        auto ptrTopologicalOrderForward() {
            return &topological_order_forward;
        }

        auto ptrTopologicalOrderBackward() {
            return &topological_order_backward;
        }

        void optimizeLPForOneIteration(double &prior_value, bool if_allow_delete_col,
                                       int lp_method = SOLVER_PRIMAL_SIMPLEX);

        [[nodiscard]] double calculateOptimalGap(double ub) const {
            return ub - value + round_up_tolerance;
        }

        void preprocessEnumeration(
            const int *col_pool4_pricing,
            Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter,
            const std::vector<routeOptLong> &ng_mem4_vertex, const std::vector<double> &optional_demand,
            double optional_cap);

        void cleanIndexColForNode(double ratio = COL_KEEP_FRAC);

        bool validateCuts(int old_num, bool if_care_lb_improvement, double cutting_branching_ratio);

        void mergeR1Cs();

        void changeEnumMatByCuts(const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter);

        int generateColumnsByInspection(const std::vector<double> &pi, std::vector<int> &col_added, RowVectorXd &rc);

        void cleanIndexColForNode(double ub, std::vector<double> &optimal_dual_vector);

        [[nodiscard]] bool tellIf4MIP() const {
            if (!if_in_enu_state || if_terminate) return false;
            int num_col;
            solver.getNumCol(&num_col);
            return valid_size + num_col <= MaxNumRoute4Mip;
        }

        void solveMIP(double ub);

        void findNonActiveCuts(std::vector<double> &optional_optimal_dual, std::vector<int> &cstr_index);

        int countActiveCuts(const std::vector<double> &optimal_dual) const;

        void deleteNonactiveCuts(std::vector<int> &nonactive_cuts, std::vector<double> &optimal_dual_vector,
                                 std::vector<int> &cstr_index);

        void obtainBrcCoefficient(const std::pair<int, int> &edge, std::vector<int> &ind,
                                  std::vector<double> &val);


        void addBranchConstraint2ColPoolInEnumByColMap(const std::pair<int, int> &edge, const int *col_pool4_pricing);

        void rmLPCols(const std::vector<int> &col_idx);

        void obtainColIdxNotAllowByEdge(const std::pair<int, int> &edge, std::vector<int> &col_idx);

        void clearEdgeMap() {
            edge_col_map.clear();
        }

        BbNode(BbNode *node,
               bool if_symmetry,
               int num_buckets_per_vertex,
               const Brc &bf);

        BbNode(BbNode *node, const Brc &bf);

        void rmColByBranchInEnuMatrix(
            std::vector<bool> &deleted_columns_in_enumeration_pool,
            bool if_not_use_by_impose_br,
            const std::vector<Brc> &brcs,
            const int *col_pool4_pricing);

        void deleteColumnByNGMemory(int start, const std::vector<routeOptLong> &ng_mem4_vertex, bool if_full_mem);

        void generateVertex2IndexColsAndEdge2IndexCols(const int *col_pool4_pricing,
                                                       const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter
                                                       &rank1_coefficient_getter,
                                                       bool if_rm_cols = true);


        void getNewIdx() {
            idx = ++node_idx_counter;
        }

    private:
        static int dim;
        static int node_idx_counter;
        static sparseRowMatrixXd row_basic_matrix;

        double last_gap{1.};
        double br_value_improved{};
        bool if_root_node{};
        bool if_terminate{};
        int idx{};
        bool if_in_enu_state{};
        std::vector<bool> deleted_columns_in_enumeration_pool{};
        int valid_size{};
        double value{};
        Solver solver{};
        std::vector<SequenceInfo> cols{};
        Bucket **all_forward_buckets{}, **all_backward_buckets{};
        std::vector<std::vector<std::vector<int> > > topological_order_forward{}, topological_order_backward{};

        RowVectorXT index_columns_in_enumeration_column_pool{};
        RowVectorXd cost_for_columns_in_enumeration_column_pool{};
        sparseColMatrixXd basic_matrix{};
        std::deque<sparseRowMatrixXd> matrix_in_enumeration{};
        std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int> >, PairHasher> edge_col_map{};
        //edge, col idx, cnt


        int num_forward_bucket_arcs{}, num_forward_jump_arcs{}, num_backward_bucket_arcs{}, num_backward_jump_arcs{};
        std::vector<Rcc> rccs{};
        std::vector<R1c> r1cs{};
        std::vector<Brc> brcs{};

        void deleteCuts(std::vector<int> &deleted_cstrs, std::vector<int> &local_cstr_index);


        void deleteBranchCutsAndR1C1s();


        void createBasicMatrix();

        void cleanColumnsCapInfeasible(const std::vector<double> &demand, double cap);

        void obtainBrcMap();

        void buildRCCInEnuMatrix(
            std::vector<Eigen::Triplet<double> > &triplets, int old_num = 0);

        void buildAllR1CInEnuMatrix(
            std::vector<Eigen::Triplet<double> > &triplets,
            const Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter,
            int old_num = 0
        );
    };
}


#endif // ROUTE_OPT_NODE_HPP
