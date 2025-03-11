/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_HELPER_L2B_HPP
#define ROUTE_OPT_T_HELPER_L2B_HPP
#include <iostream>
#include <unordered_map>
#include <cmath>
#include "route_opt_macro.hpp"
#include "l2b_macro.hpp"
#include "solver.hpp"

namespace RouteOpt::Application::CVRP::L2BDetail {
    inline void checkOptimalStatus(const Solver &solver) {
        int status;
        SAFE_SOLVER(solver.getStatus(&status))
        if (status != SOLVER_OPTIMAL && status != SOLVER_SUBOPTIMAL) {
            if (status == SOLVER_NUMERIC) {
                int crossover;
                SAFE_SOLVER(solver.getCrossOver(&crossover))
                if (crossover == SOLVER_CROSSOVER_DOWN) {
                    SAFE_SOLVER(solver.setEnvCrossOver(SOLVER_CROSSOVER_DEFAULT))
                }
                SAFE_SOLVER(solver.reoptimize(SOLVER_BARRIER))
                SAFE_SOLVER(solver.getStatus(&status))
                if (crossover == SOLVER_CROSSOVER_DOWN) {
                    SAFE_SOLVER(solver.setEnvCrossOver(SOLVER_CROSSOVER_DOWN))
                }
                if (status != SOLVER_OPTIMAL && status != SOLVER_SUBOPTIMAL) {
                    std::cout << "WARNING: status is not optimal or suboptimal again: " << status << std::endl;
                }
            } else {
                throw std::runtime_error(
                    "status is not optimal or suboptimal or numeric issue:" + std::to_string(status));
            }
        }
    }

    inline void debugInputData(const std::pair<std::string, double> &fs) {
        if (std::isnan(fs.second)) {
            PRINT_WARNING("nan in data: " + fs.first + " " + std::to_string(fs.second));
            std::cout << fs.first << " is nan" << std::endl;
        } else if (std::isinf(fs.second)) {
            PRINT_WARNING("inf in data: " + fs.first + " " + std::to_string(fs.second));
            std::cout << fs.first << " is inf" << std::endl;
        }
    }

    inline void printAllInputData(const std::pair<std::string, double> &fs) {
        std::cout << "feature: " << fs.first << " value: " << fs.second << std::endl;
    }


    inline void checkMLInputData(int num_row, float *data, int num_features) {
        PRINT_DEBUG("check in MLInputData");
        size_t tmp_len = num_row * num_features;
        for (int i = 0; i < tmp_len; i++) {
            if (std::isnan(data[i]))
                THROW_RUNTIME_ERROR("nan in data");
            if (std::isinf(data[i]))
                THROW_RUNTIME_ERROR("inf in data");
        }
    }


    template<typename BrCType, typename Hasher>
    void printFeatures(const std::unordered_map<BrCType, TmpEdgeRelatedData, Hasher> &edge_tmp_info) {
        for (auto &tmp_info: edge_tmp_info) {
            int cnt = 0;
            for (auto &feature: tmp_info.second.basic_features) {
                std::cout << " " << cnt << ":" << feature.first << " " << feature.second << std::endl;
                ++cnt;
            }
            for (auto &feature: tmp_info.second.resolving_lp_features) {
                std::cout << " " << cnt << ":" << feature.first << " " << feature.second << std::endl;
                ++cnt;
            }
            std::cout << BIG_PHASE_SEPARATION;
        }
        auto &edge = *edge_tmp_info.begin();
        int cnt = 0;
        std::vector<int> f_set;
        for (auto &feature: edge.second.basic_features) {
            if (feature.first.find(PseudoMark) != std::string::npos) {
                std::cout << " " << cnt << ":" << feature.first << " " << feature.second << std::endl;
                f_set.emplace_back(cnt);
            }
            ++cnt;
        }
        std::cout << "f_set: " << std::endl;
        std::cout << "[";
        for (auto &tmp: f_set) {
            std::cout << "," << tmp;
        }
        std::cout << "]";
    }

    template<typename BrCType, typename Hasher>
    void addEdgeInfo(
        std::vector<std::pair<std::string, double> > &edge_info,
        const std::unordered_map<BrCType, std::pair<double, int>, Hasher> &history,
        const std::pair<int, int> &edge,
        const std::string &indicator_key,
        const std::string &change_key) {
        auto it = history.find(edge);
        if (it != history.end()) {
            const auto &[first, second] = it->second; // structured binding for clarity
            edge_info.emplace_back(indicator_key, 1);
            edge_info.emplace_back(change_key, first / second);
        } else {
            edge_info.emplace_back(indicator_key, false);
            edge_info.emplace_back(change_key, 0);
        }
    }
}

#endif // ROUTE_OPT_T_HELPER_L2B_HPP
