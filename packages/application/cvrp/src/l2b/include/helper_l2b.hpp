/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HELPER_L2B_HPP
#define ROUTE_OPT_HELPER_L2B_HPP
#include <unordered_map>
#include "solver.hpp"


namespace RouteOpt::Application::CVRP::L2BDetail {
    void printAllInputData(const std::pair<std::string, double> &fs);

    void checkMLInputData(int num_row, float *data, int num_features);

    void checkOptimalStatus(const Solver &solver);

    void debugInputData(const std::pair<std::string, double> &fs);

    template<typename BrCType, typename Hasher>
    void printFeatures(const std::unordered_map<BrCType, TmpEdgeRelatedData, Hasher> &edge_tmp_info);

    template<typename BrCType, typename Hasher>
    void addEdgeInfo(
        std::vector<std::pair<std::string, double> > &edge_info,
        const std::unordered_map<BrCType, std::pair<double, int>, Hasher> &history,
        const std::pair<int, int> &edge,
        const std::string &indicator_key,
        const std::string &change_key);
}

#include "t_helper_l2b.hpp"
#endif // ROUTE_OPT_HELPER_L2B_HPP
