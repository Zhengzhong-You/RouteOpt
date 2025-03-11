/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_L2B_CONTROLLER_HPP
#define ROUTE_OPT_L2B_CONTROLLER_HPP
#include <vector>
#include <unordered_map>
#include "solver.hpp"
#include "branching_macro.hpp"
#include "label.hpp"
#include "node.hpp"

/*
 * include this hpp will result a mandatory use of BrCType=std::pair<int, int> and Hasher=PairHasher
 */
namespace RouteOpt::Application::CVRP {
    namespace L2BDetail {
        inline std::vector<int> getBrConstraintNonzeroIdx(BbNode *node, const std::pair<int, int> &edge) {
            std::vector<int> ind;
            std::vector<double> val;
            node->obtainBrcCoefficient(edge, ind, val);
            return ind;
        }
    }

    struct TmpEdgeRelatedData {
        double sb_scores{};
        std::vector<std::pair<std::string, double> > basic_features{};
        std::vector<std::pair<std::string, double> > extra_features_edge0{};
        std::vector<std::pair<std::string, double> > extra_features_edge1{};
        std::vector<std::pair<std::string, double> >
        resolving_lp_features{}; // 1 is the std::left branch, 3 is the right branch, and will be used!
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

    struct CustomerLocation {
        double x{};
        double y{};
    };

    template<typename Node, typename BrCType, typename Hasher>
    class Learning2Branch {
    public:
        //needs to clean
        std::unordered_map<BrCType, TmpEdgeRelatedData, Hasher> edge_tmp_info{};
        std::vector<BrCType> branch_pair_from_pseudo{};
        std::vector<BrCType> branch_pair_from_fractional{};
        std::vector<DualRC> dual_rc{};

        Learning2Branch(int dim,
                        const std::vector<std::vector<double> > &cost_mat4_vertex,
                        const std::vector<std::vector<Resource> > &
                        resource_across_arcs_in_forward_sense,
                        const std::vector<std::vector<Resource> > &
                        resource_across_arcs_in_backward_sense,
                        const std::function<std::vector<int> (Node *, const BrCType &)> &
                        getBrConstraintNonzeroIdx,
                        const std::vector<std::vector<double> > &
                        infor_vector, bool if_init);


        void getFeatureDataPhase1(
            const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            Node *node,
            int tree_level,
            const std::vector<double> &lp_solution,
            const std::vector<int> &route_length
        );

        void getFeatureDataPhase2(
            Node *node,
            Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            Branching::CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &branching_testing,
            double current_local_gap
        );

        void distinguishBranchPairSource(
            const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            const Branching::BranchingHistory<BrCType, Hasher> &branching_history);

        void collectOneSideEdgeFeatures(
            const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            double current_local_gap
        );

        void recordEdgeLongInfo(
            const std::unordered_map<BrCType, double, Hasher> &edge_map);

        void collectResolvingDualRC(
            const Solver &solver,
            const BrCType &candidate,
            int BeforeNumRow,
            bool if_left);

        void recordDiscrepancyLongInfo(const BrCType &candidate, double cg_change,
                                       bool if_left);

        void cleanLastData();

        //getters
        const auto &getEdgeLPChange() const {
            return edge_lp_change;
        }


        Learning2Branch() = delete;

        ~Learning2Branch() = default;

    private:
        //refers
        int dim{};
        std::reference_wrapper<const std::vector<std::vector<double> >> cost_mat4_vertex_ref;
        std::reference_wrapper<const std::vector<std::vector<Resource> >> resource_across_arcs_in_forward_sense_ref;
        std::reference_wrapper<const std::vector<std::vector<Resource> >> resource_across_arcs_in_backward_sense_ref;


        //private
        double max_edge_cost{};
        double max_mid_point_edge_cord_2_depot{};
        double cluster_coeff{};
        double depot_2_center{};
        std::pair<int, double> average_route_length{};
        std::vector<std::vector<std::pair<double, double> > > mid_point_edge_cord{};
        std::vector<std::vector<double> > mid_point_edge_cord_2_depot{};
        std::vector<std::vector<std::vector<int> > > node_density_in_std_dis_vec_form{};
        std::vector<std::vector<double> > edge_2_other_convert_dis{};
        std::unordered_map<BrCType, LongEdgeRelatedData, Hasher> edge_long_info{};
        std::vector<double> is_in_solution{};
        std::unordered_map<BrCType, std::pair<double, double>, Hasher> edge_lp_change{};

        //function ptrs
        std::function<std::vector<int> (Node *, const BrCType &)> getBrConstraintNonzeroIdx{};

        void getMaxEdgeCost();

        void getMidPointEdgeCord(const std::vector<CustomerLocation> &customer_location); //0 must be depot

        void getEdge2OtherConvertDis();

        void calculateClusteringCoefficient(const std::vector<CustomerLocation> &customer_location);

        void calculateDisDepot2Center(const std::vector<CustomerLocation> &customer_location);

        void collectVariableRelatedFeatures(
            const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            const Branching::BranchingHistory<BrCType, Hasher> &branching_history,
            const BrCType &candidate,
            const std::vector<int> &nonzero_idx);

        void collectEdgeRelatedFeatures(
            const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            int tree_level
        );

        void collectStaticFeatures(
            const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
            const std::vector<double> &lp_solution,
            const std::vector<int> &route_length);

        void collectResolvingFeatures(
            const std::vector<Branching::CandidateScoreInfo<BrCType> > &edge_info);

        void updateOneSideLPChange(const BrCType &candidate, double obj_change, bool if_left);

        void findDiscrepancyResolvingFeatures(const BrCType &candidate, bool if_left);
    };
}

#include "helper_l2b.hpp"
#include "t_l2b_phase1.hpp"
#include "t_l2b_phase2.hpp"
#include "t_l2b.hpp"
#include  "t_l2b_k_means.hpp"
#include "t_l2b_update.hpp"

#endif // ROUTE_OPT_L2B_CONTROLLER_HPP
