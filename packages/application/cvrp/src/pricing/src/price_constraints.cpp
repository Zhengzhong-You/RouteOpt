/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cuts_definition.hpp"
#include "cvrp_pricing_controller.hpp"
#include "rcc_rc_controller.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRP_Pricing::priceConstraints(const std::vector<Rcc> &rccs,
                                        const std::vector<R1c> &r1cs,
                                        const std::vector<Brc> &brcs,
                                        const std::vector<double> &pi_vector) {
        pricePartitioning(pi_vector);
        RCCs::RCGetter::RCCRCController::priceRCC(rccs, pi_vector, chg_cost_mat4_vertex);
        priceBRC(brcs, pi_vector);
        //price rank1 cuts;
        rank1_rc_controller_ref.get().getRank1DualsInCG(r1cs, pi_vector);
    }

    void CVRP_Pricing::pricePartitioning(const std::vector<double> &pi_vector) {
        auto &cm = cost_mat4_vertex_ref.get();
        auto real_dim = dim - 1;
        for (int i = 1; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                chg_cost_mat4_vertex[i][j] = cm[i][j] - 0.5 * (pi_vector[i - 1] + pi_vector[j - 1]);
            }
        }
        for (int i = 1; i < dim; ++i) {
            chg_cost_mat4_vertex[0][i] =
                    cm[0][i] - 0.5 * (pi_vector[i - 1] + pi_vector[real_dim]);
        }
        for (int i = 1; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                chg_cost_mat4_vertex[j][i] = chg_cost_mat4_vertex[i][j];
            }
        }
        for (int i = 1; i < dim; ++i) {
            chg_cost_mat4_vertex[i][0] = chg_cost_mat4_vertex[0][i];
        }
    }

    void CVRP_Pricing::priceBRC(const std::vector<Brc> &brcs, const std::vector<double> &pi_vector) {
        adjust_brc_dual4_single_route.clear();
        for (auto &brc: brcs) {
            if (!brc.br_dir) {
                chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] = std::numeric_limits<float>::max();
                chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] =
                        std::numeric_limits<float>::max(); //do not use double since the number will overflow
            } else {
                if (brc.edge.first == 0) adjust_brc_dual4_single_route[brc.edge.second] = pi_vector[brc.idx_brc];
                chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] -= pi_vector[brc.idx_brc];
                chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] -= pi_vector[brc.idx_brc];
            }
        }
    }
}
