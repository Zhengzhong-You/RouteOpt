/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <route_opt_macro.hpp>
#include <vector>
#include "rcc_macro.hpp"
#include "rcc_rc_controller.hpp"

namespace RouteOpt::RCCs::RCGetter {
    void RCCRCController::priceRCC(const std::vector<Rcc> &rccs, const std::vector<double> &pi_vector,
                                   std::vector<std::vector<double> > &chg_cost_mat4_vertex) {
        double rc;
        for (auto &rcc: rccs) {
            if (rcc.form_rcc == static_cast<int>(RCCForm::RCC_FORM_1)) {
                auto &info = rcc.info_rcc_customer;
                rc = pi_vector[rcc.idx_rcc];
                for (auto i = info.begin(); i != info.end(); ++i) {
                    auto j = i;
                    ++j;
                    for (; j != info.end(); ++j) {
                        chg_cost_mat4_vertex[*i][*j] -= rc;
                        chg_cost_mat4_vertex[*j][*i] -= rc;
                    }
                }
            } else if (rcc.form_rcc == static_cast<int>(RCCForm::RCC_FORM_2)) {
                auto &outside_customer_info = rcc.info_rcc_outside_customer;
                auto &customer_info = rcc.info_rcc_customer;
                rc = pi_vector[rcc.idx_rcc];
                for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
                    auto j = i;
                    ++j;
                    for (; j != outside_customer_info.end(); ++j) {
                        chg_cost_mat4_vertex[*i][*j] -= rc;
                        chg_cost_mat4_vertex[*j][*i] -= rc;
                    }
                }
                double half_rc = 0.5 * rc;
                for (auto it: outside_customer_info) {
                    chg_cost_mat4_vertex[0][it] -= half_rc;
                    chg_cost_mat4_vertex[it][0] -= half_rc;
                }
                for (auto it: customer_info) {
                    chg_cost_mat4_vertex[0][it] += half_rc;
                    chg_cost_mat4_vertex[it][0] += half_rc;
                }
            } else
                THROW_RUNTIME_ERROR("rcc form does not match")
        }
    }
}
