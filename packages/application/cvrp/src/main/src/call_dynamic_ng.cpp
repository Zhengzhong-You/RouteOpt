/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <set>
#include <chrono>
#include <boost/math/tools/config.hpp>

#include "cvrp.hpp"

namespace RouteOpt::Application::CVRP {
    namespace DynamicNGDetail {
        void findNGMemorySets(BbNode *node,
                              int dim,
                              std::vector<routeOptLong> &ng_mem4_vertex,
                              bool &if_empty) {
            std::unordered_map<int, std::vector<std::pair<int, routeOptLong> > > map_size_cycle;
            std::vector<size_t> tmp(dim);

            int num_col;
            SAFE_SOLVER(node->refSolver().getNumCol(&num_col))
            std::vector<double> x(num_col);
            SAFE_SOLVER(node->refSolver().getX(0, num_col, x.data()))

            std::vector<std::pair<int, double> > opt_col_idx;

            for (int i = 0; i < num_col; ++i) {
                if (x[i] < SOL_NG_X_TOLERANCE) continue;
                opt_col_idx.emplace_back(i, x[i]);
            }

            std::sort(opt_col_idx.begin(), opt_col_idx.end(),
                      [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
                          return a.second > b.second;
                      });

            int re_cols = 0;
            for (auto &pr: opt_col_idx) {
                std::fill(tmp.begin(), tmp.end(), -1);
                bool if_re_col = false;
                auto &route = node->getCols()[pr.first].col_seq;
                for (int i = 0; i < route.size(); ++i) {
                    int curr_node = route[i];
                    if (tmp[curr_node] != -1) {
                        if_re_col = true;
                        auto length = static_cast<int>(i - tmp[curr_node] - 1);
                        map_size_cycle[length].emplace_back(curr_node, 0);
                        auto &mem = map_size_cycle[length].back().second;
                        for (auto k = tmp[curr_node] + 1; k < i; ++k) {
                            mem.set(route[k]);
                        }
                    }
                    tmp[curr_node] = i; //not in else!
                }
                if (if_re_col) {
                    ++re_cols;
                    if (re_cols > MaxNumColsInNGAug) break;
                }
            }
            std::vector<std::pair<int, std::vector<std::pair<int, routeOptLong> > > > size_cycle(
                map_size_cycle.begin(), map_size_cycle.end());
            std::sort(size_cycle.begin(), size_cycle.end(), [](
                  const std::pair<int, std::vector<std::pair<int, routeOptLong> > > &a,
                  const std::pair<int, std::vector<std::pair<int, routeOptLong> > > &b) {
                          return a.first < b.first;
                      });
            if (size_cycle.empty()) {
                PRINT_REMIND("no small cycles are found!");
                if_empty = true;
                return;
            }
            if_empty = false;
            for (auto &i: size_cycle[0].second) {
                int re = i.first;
                for (int j = 1; j < dim; ++j) {
                    if (i.second[j]) {
                        auto &mem = ng_mem4_vertex[j];
                        mem.set(re);
                        if (mem.count() == MAX_NG_SIZE) break;
                    }
                }
            }
            for (int i = 1; i < size_cycle.size(); ++i) {
                if (size_cycle[i].first > CYCLE_SIZE) break;
                for (auto &j: size_cycle[i].second) {
                    int re = j.first;
                    for (int k = 1; k < dim; ++k) {
                        if (j.second[k]) {
                            auto &mem = ng_mem4_vertex[k];
                            mem.set(re);
                            if (mem.count() == MAX_NG_SIZE) break;
                        }
                    }
                }
            }

            node->deleteColumnByNGMemory(1, ng_mem4_vertex, false);
        }

        void reduceNGSet(BbNode *node, int dim, std::vector<routeOptLong> &ng_mem4_vertex) {
            auto sol_edge_map = BbNode::obtainSolEdgeMap(node);
            for (int i = 1; i < dim; ++i) {
                std::set<int> tmp;
                std::set<int> tmp2;
                for (int j = 1; j < dim; ++j) {
                    if (i == j) continue;
                    auto pr = i < j ? std::make_pair(i, j) : std::make_pair(j, i);
                    if (sol_edge_map[pr] > TOLERANCE) {
                        tmp.emplace(j);
                    }
                    if (ng_mem4_vertex[i][j]) tmp2.emplace(j);
                }
                std::set<int> tmp3;
                set_intersection(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end(), inserter(tmp3, tmp3.begin()));
                ng_mem4_vertex[i] = 0;
                ng_mem4_vertex[i].set(i);
                for (auto n: tmp3) {
                    ng_mem4_vertex[i].set(n);
                }
            }
        }

        void dssrRound(BbNode *node, CVRPSolver *cvrp, std::vector<routeOptLong> &ng_mem4_vertex, bool &if_empty,
                       double &max_time) {
            int rd = 0;
            double old_val = node->getValue();
            while (true) {
                std::cout << SMALL_PHASE_SEPARATION;
                std::cout << "DSSR: " << rd++ << std::endl;
                findNGMemorySets(node, cvrp->getDim(), ng_mem4_vertex, if_empty);
                if (if_empty) break;
                auto eps = TimeSetter::measure([&]() {
                    cvrp->solveLPInLabeling(node, true, true, true,
                                            false, false,
                                            false, false, true,
                                            std::numeric_limits<float>::max());
                });
                max_time = std::max(max_time, eps);
                if (node->getValue() + TOLERANCE > old_val) break;
            }
        }

        void augmentNG(BbNode *node, CVRPSolver *cvrp, double std_time, std::vector<routeOptLong> &ng_mem4_vertex) {
            int rd = 0;
            auto old_val = node->getValue();
            auto improved = TOLERANCE;
            while (true) {
                std::cout << "NG Augmentation: " << rd++ << std::endl;
                bool if_empty;
                auto old_mem = ng_mem4_vertex;
                findNGMemorySets(node, cvrp->getDim(), ng_mem4_vertex, if_empty);
                if (if_empty) break;

                auto eps = TimeSetter::measure([&]() {
                    cvrp->solveLPInLabeling(node, true, true, true,
                                            false, false,
                                            false, false, true,
                                            NGAugTimeHardThresholdFactor * std_time);
                });

                if (eps > NGAugTimeHardThresholdFactor * std_time) {
                    std::cout << "eps= " << eps << " std_time= " << std_time << " NGAugTimeHardThresholdFactor= " <<
                            NGAugTimeHardThresholdFactor << std::endl;
                    PRINT_REMIND("NG Augmentation: too time consuming, rollback!");
                    ng_mem4_vertex = old_mem;
                    node->deleteColumnByNGMemory(1, ng_mem4_vertex, false);
                    cvrp->solveLPInLabeling(node, true, true, true,
                                            false, false,
                                            false, false, true,
                                            std::numeric_limits<float>::max());
                    return;
                }
                if (eps > NGAugTimeSoftThresholdFactor * std_time) {
                    break;
                }
                double tmp = std::abs(node->getValue() - old_val) / node->getValue();
                if (tmp > improved) improved = tmp;
                else if (tmp < improved * NGAugTailOff) break;
                else old_val = node->getValue();
            }
        }
    }

    void CVRPSolver::augmentNGRound(BbNode *node, std::vector<routeOptLong> &ng_mem4_vertex) {
        if (node->getIfTerminate()) return;
        std::vector<routeOptLong> old_mem = ng_mem4_vertex;
        DynamicNGDetail::reduceNGSet(node, dim, ng_mem4_vertex);
        bool if_empty;
        double max_time = 0.;
        DynamicNGDetail::dssrRound(node, this, ng_mem4_vertex, if_empty, max_time);
        if (if_empty) return;
        DynamicNGDetail::augmentNG(node, this, max_time, ng_mem4_vertex);
    }
}
