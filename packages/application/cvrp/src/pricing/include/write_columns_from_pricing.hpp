/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_WRITE_COLUMNS_FROM_PRICING_HPP
#define ROUTE_OPT_WRITE_COLUMNS_FROM_PRICING_HPP
#include "cvrp_pricing_controller.hpp"
#include "pricing_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<bool if_symmetry>
    void CVRP_Pricing::initializeLabels() {
        all_label[0].end_vertex = 0;
        all_label[0].p_label = nullptr;
        Rank1Cuts::RCGetter::Rank1RCController::initializeLabel(all_label[0].r1c);
        for (int i = 1; i < dim; ++i) {
            all_label[i].pi.set(i);
            all_label[i].cost = cost_mat4_vertex_ref.get()[0][i];
            all_label[i].end_vertex = i;
            all_label[i].p_label = all_label;
            increaseMainResourceConsumption({}, all_label[i].res, 0, i);
        }
        if constexpr (!if_symmetry) {
            int max_num = 2 * dim - 1;
            for (int i = dim; i < max_num; ++i) {
                int point = i - dim + 1;
                all_label[i].pi.set(point);
                all_label[i].cost = cost_mat4_vertex_ref.get()[0][point];
                all_label[i].end_vertex = point;
                all_label[i].p_label = all_label;
                decreaseMainResourceConsumption(resource, all_label[i].res, 0, point);
            }
        }
    }

    inline void CVRP_Pricing::reallocateLabel() {
        try {
            delete[] all_label;
            label_assign *= 2;
            all_label = new Label[label_assign];
            rank1_rc_controller_ref.get().assignLabelMem(all_label, label_assign, &Label::r1c);
        } catch (const std::bad_alloc &) {
            THROW_RUNTIME_ERROR("Label assigned=" + std::to_string(label_assign) +
                " failed in reallocateLabel. Not enough memory.");
        }
    }

    inline int CVRP_Pricing::checkPricingPool() const {
        if (pool_beg4_pricing >= static_cast<size_t>(FracMemTolerance * static_cast<double>(mem4_pricing))) {
            return true;
        }
        return false;
    }

    inline void CVRP_Pricing::reallocatePricingPool(size_t num) {
        try {
            if (num == 0) {
                mem4_pricing *= 2;
            } else {
                if (num < mem4_pricing) return;
                mem4_pricing = static_cast<size_t>(2 * static_cast<double>(num));
            }
            auto tmp = col_pool4_pricing;
            col_pool4_pricing = new int[mem4_pricing];
            std::copy(tmp, tmp + pool_beg4_pricing, col_pool4_pricing);
            delete[] tmp;
        } catch (const std::bad_alloc &) {
            THROW_RUNTIME_ERROR("Memory allocation failed during reallocatePricingPool. Not enough memory.");
        }
    }


    inline void cleanNegativeTuple(std::vector<std::tuple<Label *, Label *, double> > &negative_rc_label_tuple,
                                   int num) {
        if (negative_rc_label_tuple.empty()) return;
        std::sort(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(),
                  [](const std::tuple<Label *, Label *, double> &a, const std::tuple<Label *, Label *, double> &b) {
                      return std::get<2>(a) < std::get<2>(b);
                  });
        int cnt = 0;
        for (auto it = negative_rc_label_tuple.begin() + 1; it != negative_rc_label_tuple.end();) {
            if (std::abs(std::get<2>(*it) - std::get<2>(*(it - 1))) < TOLERANCE) {
                it = negative_rc_label_tuple.erase(it);
            } else {
                ++it;
                ++cnt;
                if (cnt >= num) break;
            }
        }
        if (negative_rc_label_tuple.size() > num)
            negative_rc_label_tuple.resize(num);
    }

    inline void CVRP_Pricing::writeColumnsInPricingPool() {
        Label *p;
        if (checkPricingPool()) reallocatePricingPool();

        cleanNegativeTuple(negative_rc_label_tuple, num_col_generated_ub);
        new_cols.clear();
        new_cols.resize(negative_rc_label_tuple.size());
        int col_idx = 0;

        for (auto &i: negative_rc_label_tuple) {
            auto &ki = std::get<0>(i);
            auto &kj = std::get<1>(i);
            auto &col = new_cols[col_idx].col_seq;

            col.reserve(dim);

            p = ki;
            while (p && p->end_vertex) {
                col.emplace_back(p->end_vertex);

                p = p->p_label;
            }
            std::reverse(col.begin(), col.end());

            new_cols[col_idx].forward_concatenate_pos = static_cast<int>(col.size()) - 1;

            if (kj) {
                p = kj;
                while (p && p->end_vertex) {
                    col.emplace_back(p->end_vertex);

                    p = p->p_label;
                }
            }
            if constexpr (CHECK_PRICING_LABELS)seq_rc[col] = std::get<2>(i);

            ++col_idx;
        }
        if constexpr (INSPECT_COLUMN_FEASIBILITY) {
            PRINT_DEBUG("inspect columns");
            inspectColumns();
        }
    }

    inline double CVRP_Pricing::getSmallestRC() {
        return negative_rc_label_tuple.empty() ? 0 : std::get<2>(negative_rc_label_tuple[0]);
    }

    inline void CVRP_Pricing::setTerminateMarker(double val, double ub, bool &if_terminate) {
        if (if_exact_labeling_cg && if_exact_labeling_finished) {
            if (ceilTransformedNumberRelated(getSmallestRC() * max_num_vehicle_ref + val + RC_TOLERANCE) + TOLERANCE
                >= ub) {
                if_terminate = true;
            }
        }
    }

    inline void CVRP_Pricing::inspectColumns() {
        for (const auto &col: new_cols) {
            auto &col_seq = col.col_seq;
            auto res = Resource{};
            int b4 = 0;
            for (int j: col_seq) {
                if (!increaseMainResourceConsumption(res, res, b4, j)) {
                    THROW_RUNTIME_ERROR("infeasible solution in extension");
                }
                b4 = j;
            }
            if (res.resources[0] > resource.resources[0]) {
                THROW_RUNTIME_ERROR("infeasible solution");
            }
        }
    }
}

#endif // ROUTE_OPT_WRITE_COLUMNS_FROM_PRICING_HPP
