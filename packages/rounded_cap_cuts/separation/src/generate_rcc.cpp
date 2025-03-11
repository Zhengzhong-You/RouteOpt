/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */


#include <vector>
#include <algorithm>
#include "route_opt_macro.hpp"
#include "rcc_macro.hpp"
#include "rcc_separation_macro.hpp"
#include "rcc_separation_controller.hpp"
#include "capsep.h"
#include "cnstrmgr.h"
#include "rank1_macro.hpp"


namespace RouteOpt::RCCs::Separation {
    void getMap(const std::vector<double> &sol_x,
                const std::vector<SequenceInfo> &sols,
                int dim,
                std::vector<int> &edge_tail, std::vector<int> &edge_head,
                std::vector<double> &edge_value) {
        std::unordered_map<std::pair<int, int>, double, PairHasher> coeff_map;
        for (int i = 0; i < sol_x.size(); ++i) {
            double val = sol_x[i];
            const auto &c = sols[i].col_seq;
            int b4 = 0;
            for (int j: c) {
                if (b4 < j) {
                    coeff_map[{b4, j}] += val;
                } else {
                    coeff_map[{j, b4}] += val;
                }
                b4 = j;
            }
            coeff_map[{0, b4}] += val;
        }
        edge_head.resize(1);
        edge_tail.resize(1);
        edge_value.resize(1);

        for (const auto &pr: coeff_map) {
            if (pr.first.first == 0) {
                edge_tail.emplace_back(dim);
            } else {
                edge_tail.emplace_back(pr.first.first);
            }
            edge_head.emplace_back(pr.first.second);
            edge_value.emplace_back(pr.second);
        }
    }

    void RCCSeparationController::generateRCCs(int dim,
                                               double cap,
                                               const std::vector<double> &demand,
                                               bool if_keep_rcc,
                                               bool if_force_first_form,
                                               bool if_strengthen_rcc,
                                               const std::vector<double> &sol_x,
                                               const std::vector<SequenceInfo> &sols,
                                               const std::vector<Rcc> &existing_rccs,
                                               std::vector<Rcc> &new_rccs) {
        CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
        CMGR_CreateCMgr(&MyCutsCMP, MAX_NUM_OF_CUTS);
        CMGR_CreateCMgr(&MyOldCutsCMP, MAX_NUM_OF_CUTS);
        std::vector<int> solver_ind;
        std::vector<double> solver_val;

        std::vector<int> edge_tail;
        std::vector<int> edge_head;
        std::vector<double> edge_value;

        getMap(sol_x, sols, dim, edge_tail, edge_head, edge_value);

        char if_int_n_feasible;
        double max_vio = MAX_VIO;

        CAPSEP_SeparateCapCuts(dim - 1, demand.data(), cap,
                               static_cast<int>(edge_tail.size()) - 1,
                               edge_tail.data(),
                               edge_head.data(), edge_value.data(),
                               MyOldCutsCMP,
                               MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                               &if_int_n_feasible, &max_vio, MyCutsCMP);


        if (!MyCutsCMP->Size) goto QUIT;

        new_rccs.clear();
        for (int i = 0; i < MyCutsCMP->Size; ++i) {
            Rcc rcc;
            auto &tmp_customerInfo = rcc.info_rcc_customer;
            for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
                tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
            }

            if (if_force_first_form || tmp_customerInfo.size() <= dim / 2 || if_strengthen_rcc) {
                rcc.form_rcc = if_strengthen_rcc
                                   ? static_cast<int>(RCCForm::RCC_FORM_3)
                                   : static_cast<int>(RCCForm::RCC_FORM_1);
                rcc.rhs = MyCutsCMP->CPL[i]->RHS;
                if (if_strengthen_rcc) rcc.rhs = static_cast<double>(tmp_customerInfo.size()) - rcc.rhs;
            } else {
                rcc.form_rcc = static_cast<int>(RCCForm::RCC_FORM_2);
                rcc.rhs = dim - 1 - static_cast<double>(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
                auto &tmp_NoCustomerInfo = rcc.info_rcc_outside_customer;
                std::vector<bool> tmp(dim, false);
                for (int j: tmp_customerInfo) tmp[j] = true;
                for (int j = 1; j < dim; ++j) {
                    if (!tmp[j]) {
                        tmp_NoCustomerInfo.emplace_back(j);
                    }
                }
            }

            if (std::find(existing_rccs.begin(), existing_rccs.end(), rcc) != existing_rccs.end()) {
                continue;
            }

            if (if_keep_rcc) rcc.if_keep = true;

            new_rccs.emplace_back(rcc);
        }
    QUIT:
        for (int i = 0; i < MyCutsCMP->Size; ++i) {
            CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
        }
        MyCutsCMP->Size = 0;

        CMGR_FreeMemCMgr(&MyOldCutsCMP);
        CMGR_FreeMemCMgr(&MyCutsCMP);
    }
}
