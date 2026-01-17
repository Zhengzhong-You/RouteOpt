/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_TWO_STAGE_CONTROLLER_HPP
#define ROUTE_OPT_TWO_STAGE_CONTROLLER_HPP
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif
#include "cvrp_pricing_controller.hpp"
#include "label.hpp"

namespace RouteOpt::Application::CVRP {
    class TwoStageController {
    public:
        template<bool if_symmetry>
        static void writeNodeOut(
            const std::string &ins_name,
            const Solver &node_solver,
            const std::vector<SequenceInfo> &cols,
            const std::vector<Rcc> &rccs,
            const std::vector<R1c> &r1cs,
            const std::vector<Brc> &brcs,
            double val, int idx, bool if_enu,
            int dim, int num_buckets_per_vertex, res_int step_size,
            VecLabel *const *if_exist_extra_labels_in_forward_sense,
            Bucket *const *all_forward_buckets,
            VecLabel *const *if_exist_extra_labels_in_backward_sense,
            Bucket *const *all_backward_buckets,
            const std::vector<bool> &deleted_columns_in_enumeration_pool,
            const RowVectorXT &index_columns_in_enumeration_column_pool,
            const RowVectorXd &cost_for_columns_in_enumeration_column_pool,
            const int *col_pool4_pricing,
            double max_enumeration_success_gap,
            const std::pair<double, int> &success_enumeration_gap,
            double min_enumeration_fail_gap,
            double max_gap2try_enumeration,
            const Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
            const Branching::BKF::BKFDataShared &bkf_data_shared);


        //node in

        template<bool if_symmetry>
        static void readNodeIn(const std::string &node_file_name,
                               std::vector<double> &cost,
                               std::vector<SequenceInfo> &cols,
                               std::vector<Rcc> &rccs,
                               std::vector<R1c> &r1cs,
                               std::vector<Brc> &brcs,
                               double &val, int &idx, bool &if_enu,
                               int dim, int &num_buckets_per_vertex, res_int &step_size,
                               VecLabel **&if_exist_extra_labels_in_forward_sense,
                               Bucket **&all_forward_buckets,
                               VecLabel **&if_exist_extra_labels_in_backward_sense,
                               Bucket **&all_backward_buckets,
                               std::vector<bool> &deleted_columns_in_enumeration_pool,
                               RowVectorXT &index_columns_in_enumeration_column_pool,
                               RowVectorXd &cost_for_columns_in_enumeration_column_pool,
                               int *&col_pool4_pricing,
                               size_t &mem4_pricing,
                               size_t &pool_beg4_pricing,
                               double &max_enumeration_success_gap,
                               std::pair<double, int> &success_enumeration_gap,
                               double &min_enumeration_fail_gap,
                               double &max_gap2try_enumeration,
                               Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                               Branching::BKF::BKFDataShared &bkf_data_shared);


        //static update ub
        static void updateUB(const std::string &ub_name, double &ub);

        static void deleteInFile(const std::string &node_file_name);

        TwoStageController() = default;

        ~TwoStageController() = default;
    };
}


#include "t_node_out.hpp"
#include "t_node_in.hpp"

#endif // ROUTE_OPT_TWO_STAGE_CONTROLLER_HPP
