/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_ADD_COLUMN_CONTROLLER_HPP
#define ROUTE_OPT_ADD_COLUMN_CONTROLLER_HPP
#include <vector>
#include <deque>
#include "route_opt_macro.hpp"
#include "label.hpp"
#include "rank1_coefficient_controller.hpp"
#include "cuts_definition.hpp"
#include "cvrp_macro.hpp"

namespace RouteOpt::Application::CVRP {
    class CVRP_AddColumnController {
    public:
        CVRP_AddColumnController(
            const int &dim,
            const std::vector<SequenceInfo> &new_cols,
            const std::vector<std::vector<double> > &cost_mat4_vertex,
            const std::vector<std::tuple<Label *, Label *, double> > &negative_rc_label_tuple,
            Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter &rank1_coefficient_getter,
            std::map<std::vector<int>, double> &seq_rc_map,
            int *&col_pool4_pricing)
            : dim(dim),
              cost_mat4_vertex_ref(cost_mat4_vertex),
              negative_rc_label_tuple_ref(negative_rc_label_tuple),
              new_cols_ref(new_cols),
              rank1_coefficient_getter_ref(rank1_coefficient_getter),
              seq_rc_map_ref(seq_rc_map),
              col_pool4_pricing_ref(col_pool4_pricing) {
        }


        void updatePtr(
            std::vector<SequenceInfo> &existing_cols,
            Solver *solver,
            const std::vector<Rcc> &rccs,
            const std::vector<R1c> &r1cs,
            const std::vector<Brc> &brcs,
            const RowVectorXT &idx_col_pool,
            const RowVectorXd &cost_col_pool,
            const std::deque<sparseRowMatrixXd> &matrix_col_pool,
            const int *col_pool4_pricing) {
            this->existing_cols_ptr = &existing_cols;
            this->solver_ptr = solver;
            this->rccs_ptr = &rccs;
            this->r1cs_ptr = &r1cs;
            this->brcs_ptr = &brcs;
            this->idx_col_pool_ptr = &idx_col_pool;
            this->cost_col_pool_ptr = &cost_col_pool;
            this->matrix_col_pool_ptr = &matrix_col_pool;
        }

        void addColumns(int &ccnt, const std::vector<double> &pi4_labeling, bool if_check_rc, bool if_enu = false);

        void addColumnsByInspection(const std::vector<int> &Col_added);

    private:
        int num_row{};

        int dim{};
        //ref
        std::reference_wrapper<const std::vector<std::vector<double> >> cost_mat4_vertex_ref;
        std::reference_wrapper<const std::vector<std::tuple<Label *, Label *, double> >> negative_rc_label_tuple_ref;
        std::reference_wrapper<const std::vector<SequenceInfo>> new_cols_ref;
        //cannot be const
        std::reference_wrapper<Rank1Cuts::CoefficientGetter::Rank1CoefficientGetter> rank1_coefficient_getter_ref;
        std::reference_wrapper<std::map<std::vector<int>, double> > seq_rc_map_ref;
        std::reference_wrapper<int *> col_pool4_pricing_ref;


        //ptr
        const std::vector<Rcc> *rccs_ptr{};
        const std::vector<R1c> *r1cs_ptr{};
        const std::vector<Brc> *brcs_ptr{};
        const RowVectorXT *idx_col_pool_ptr{};
        const RowVectorXd *cost_col_pool_ptr{};
        const std::deque<sparseRowMatrixXd> *matrix_col_pool_ptr{};
        //cannot be const
        Solver *solver_ptr{};
        std::vector<SequenceInfo> *existing_cols_ptr{};

        void calculateColumnCoefficientsB4Enumeration(bool if_enu,
                                                      sparseColMatrixXd &mat,
                                                      Eigen::RowVectorXd &cost,
                                                      std::unordered_map<std::pair<int, int>,
                                                          std::vector<int>
                                                          , PairHasher> &
                                                      edge_map);
    };
}

#endif // ROUTE_OPT_ADD_COLUMN_CONTROLLER_HPP
