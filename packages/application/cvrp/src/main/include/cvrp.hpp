/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_CVRP_HPP
#define ROUTE_OPT_CVRP_HPP
#include <chrono>

#include "cvrp_pricing_controller.hpp"
#include "node.hpp"
#include "read_data_controller.hpp"
#include "hgs_controller.hpp"
#include "initial_controller.hpp"
#include "rank1_separation_controller.hpp"
#include "rank1_coefficient_controller.hpp"
#include "add_column_controller.hpp"
#include "l2b_train.hpp"
#include "rcc_coefficient_controller.hpp"
#include "rcc_separation_controller.hpp"
#include "rcc_rc_controller.hpp"
#include "two_stage_controller.hpp"

namespace RouteOpt::Application::CVRP {
    using CandidateSelectionFuncType = std::function<std::pair<int, int>(
        BbNode *,
        Branching::BranchingHistory<std::pair<int, int>, PairHasher> &,
        Branching::BranchingDataShared<std::pair<int, int>, PairHasher> &,
        Branching::CandidateSelector::BranchingTesting<BbNode, std::pair<int, int>, PairHasher> &
    )>;

    using OutNodeFuncType = std::function<void(
        BbNode *,
        const Branching::BranchingHistory<std::pair<int, int>, PairHasher> &,
        const Branching::BKF::BKFDataShared &)>;

    using InNodeFuncType = std::function<void(
        BbNode *,
        Branching::BranchingHistory<std::pair<int, int>, PairHasher> &,
        Branching::BKF::BKFDataShared &)>;

    class CVRPSolver {
    public:
        CVRPSolver(int argc, char *argv[]): read_data_controller(argc, argv, f_name,
                                                                 ins_name, tree_path, num_vehicle, cap, dim,
                                                                 info_vertex, ub),
                                            hgs_controller(IF_HGS_HEURISTIC, f_name, dim, ub, ip_opt_sol),
                                            initial_controller(
                                                dim,
                                                solver, demand,
                                                H, earliest_time, latest_time, service_time,
                                                cost_mat4_vertex,
                                                info_vertex,
                                                [this] {
                                                    this->getLowerBoundofMinimumNumberCars();
                                                }),
                                            rank1_cuts_data_shared(dim),
                                            rank1_coefficient_getter(rank1_cuts_data_shared),
                                            rank1_rc_controller(rank1_cuts_data_shared),
                                            rank1_separation_controller(
                                                rank1_cuts_data_shared, MAX_ROW_RANK1,
                                                MAX_NUM_R1C3_PER_ROUND, MAX_NUM_R1C_PER_ROUND, solver,
                                                cost_mat4_vertex),
                                            pricing_controller(
                                                dim, max_num_vehicle, cap, demand,
                                                H, earliest_time, latest_time, service_time,
                                                cost_mat4_vertex,
                                                rank1_rc_controller),
                                            add_column_controller(
                                                dim, pricing_controller.getNewCols(), cost_mat4_vertex,
                                                pricing_controller.getNegativeLabelTuple(),
                                                rank1_coefficient_getter,
                                                pricing_controller.refSeqRCMap(),
                                                pricing_controller.refColumnPoolPtr()),
                                            l2b_controller(dim, cost_mat4_vertex,
                                                           pricing_controller.getResourceForward(),
                                                           pricing_controller.getResourceBackward(),
                                                           L2BDetail::getBrConstraintNonzeroIdx,
                                                           info_vertex, ml_type != ML_TYPE::ML_NO_USE),
                                            l2b_predict_controller(l2b_controller,
                                                                   [this](BbNode *node, const std::pair<int, int> &edge,
                                                                          double &dif1, bool dir) {
                                                                       this->processOneSideCGTesting<false>(
                                                                           node, edge, dif1,
                                                                           dir);
                                                                   }, ml_type == ML_TYPE::ML_GET_DATA_2 || ml_type ==
                                                                      ML_TYPE::ML_USE_MODEL,
                                                                   ml_type == ML_TYPE::ML_USE_MODEL),
                                            l2b_train_controller(l2b_controller, ins_name,
                                                                 ml_type == ML_TYPE::ML_GET_DATA_1 || ml_type ==
                                                                 ML_TYPE::ML_GET_DATA_2) {
            std::cout << "instance name= " << ins_name << std::endl;
        }

        void callCutting(BbNode *node);

        void imposeBranching(BbNode *node, const std::pair<int, int> &brc, std::vector<BbNode *> &children);

        void processLPTesting(BbNode *node, const std::pair<int, int> &edge, double &dif1, double &dif2);

        template<bool if_exact>
        void processCGTesting(BbNode *node, const std::pair<int, int> &edge, double &dif1, double &dif2);

        template<bool if_exact>
        void processOneSideCGTesting(BbNode *node, const std::pair<int, int> &edge, double &dif1,
                                     bool dir);

        std::pair<int, int> callMLCandidateSelection(BbNode *node,
                                                     Branching::BranchingHistory<std::pair<int, int>,
                                                         PairHasher> &history,
                                                     Branching::BranchingDataShared<std::pair<int, int>,
                                                         PairHasher> &data_shared,
                                                     Branching::CandidateSelector::BranchingTesting<BbNode,
                                                         std::pair<int, int>, PairHasher> &tester);

        void callWriteNodeOut(BbNode *node,
                              const Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                              const Branching::BKF::BKFDataShared &bkf_data_shared) const;

        void callReadNodeIn(BbNode *node,
                            Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                            Branching::BKF::BKFDataShared &bkf_data_shared);


        virtual void getLowerBoundofMinimumNumberCars();

        CVRPSolver() = delete;

        virtual ~CVRPSolver() {
            solver.freeEnv();
        }

        //getters
        auto getDim() const {
            return dim;
        }

        auto getNumVehicle() const {
            return num_vehicle;
        }

        auto getCap() const {
            return cap;
        }

        auto &getDemand() const {
            return demand;
        }

        auto &getDisMat() const {
            return cost_mat4_vertex;
        }

        //refers; allow for modification
        Solver &refSolver() {
            return solver;
        }

        auto &refUB() {
            return ub;
        }

        void terminateByMIP(BbNode *node);

        void printOptSol(std::ostream &os, int num_nodes, double lower_bound);

        void callPricingAtBeg(BbNode *node);

        void solveLPInLabeling(BbNode *node, bool if_open_heur, bool if_open_exact,
                               bool if_update_node_val,
                               bool if_consider_regenerate_bucket_graph, bool if_possible_terminate_early,
                               bool if_fix_row, bool if_fix_meet_point, bool if_allow_delete_col, double time_limit);

        void solveLPByInspection(BbNode *node, bool if_update_column_pool, bool if_allow_delete_col);

        void augmentNGRound(BbNode *node, std::vector<routeOptLong> &ng_mem4_vertex);

        void initSolver();

    private:
        Solver solver{};
        int num_vehicle{}, max_num_vehicle{};
        int dim{};
        std::vector<double> demand{};
        std::vector<std::vector<double> > info_vertex{};
        double cap{};
        double H{};
        std::vector<double> earliest_time{};
        std::vector<double> latest_time{};
        std::vector<double> service_time{};
        //read from file
        std::string f_name{};
        std::string ins_name{};
        std::string tree_path{};

        //refer
        double ub{};

        //processed
        std::vector<std::vector<int> > ip_opt_sol{};
        std::vector<std::vector<double> > cost_mat4_vertex{};

        //cvrp private controller
        std::vector<double> optimal_dual_vector{};


        //controllers
        CVRP_ReadDataController read_data_controller;
        HGSController hgs_controller;
        InitialController initial_controller;


        //r1c
        Rank1Cuts::Rank1CutsDataShared rank1_cuts_data_shared;
        Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter rank1_coefficient_getter;
        Rank1Cuts::RCGetter::Rank1RCController rank1_rc_controller;
        Rank1Cuts::Separation::Rank1SeparationController rank1_separation_controller;


        CVRP_Pricing pricing_controller;
        CVRP_AddColumnController add_column_controller;

        //ML
        Learning2Branch<BbNode, std::pair<int, int>, PairHasher> l2b_controller; // must put before train_controller
        Predict<BbNode, std::pair<int, int>, PairHasher> l2b_predict_controller;
        GetTrainingData<BbNode, std::pair<int, int>, PairHasher> l2b_train_controller;


        void callEnumeration(BbNode *node);

        void setEnv(BbNode *node);

        void updateIntegerSolution(double val, const std::vector<double> &X,
                                   const std::vector<SequenceInfo> &cols, bool &if_integer, bool &if_feasible);

        void runColumnGenerationType(BbNode *node, PRICING_LEVEL cg_mode, double time_limit,
                                     bool if_possible_terminate_early,
                                     bool if_fix_row, bool if_fix_meet_point, bool if_allow_delete_col,
                                     bool if_last_cg_type, bool if_stabilization);

        virtual void checkSolutionFeasibility(const std::vector<double> &X,
                                              const std::vector<SequenceInfo> &cols,
                                              bool &feasible) {
            feasible = true;
        };

        virtual void addFeasibilityCuts(int &num_row,
                                        const std::vector<double> &x,
                                        const std::vector<SequenceInfo> &cols,
                                        std::vector<Rcc> &existing_rccs,
                                        Solver &solver) {
        };

        void callInspection(BbNode *node, double &time_4_pure_pricing);

        void callLabeling(BbNode *node, double labeling_time_limit, double &time_4_pure_pricing);

        void callPricing(BbNode *node, double labeling_time_limit, double &time_4_pure_pricing);
    };

    namespace ConstructInitialColumnsDetail {
        void constructInitialColumns(const CVRPSolver &cvrp_solver, BbNode *node);
    }
}

#include "call_branching.hpp"
#endif // ROUTE_OPT_CVRP_HPP
