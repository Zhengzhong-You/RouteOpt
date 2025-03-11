/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_GET_TOPOLOGICAL_ORDER_HPP
#define ROUTE_OPT_GET_TOPOLOGICAL_ORDER_HPP
#include <stack>
#include <queue>

namespace RouteOpt::Application::CVRP {
    namespace detail {
        inline void fillOrder(const std::unordered_map<int, std::vector<int> > &graph,
                              int node,
                              std::unordered_set<int> &visited,
                              std::stack<int> &finishStack) {
            std::vector<bool> onStack(graph.size() + 1, false); //since we start from 1
            std::stack<int> dfsStack;
            dfsStack.push(node);

            while (!dfsStack.empty()) {
                int curr = dfsStack.top();
                if (!visited.count(curr)) {
                    visited.insert(curr);
                    onStack[curr] = true;
                    if (graph.find(curr) != graph.end()) {
                        for (int neighbor: graph.at(curr)) {
                            if (!visited.count(neighbor) && !onStack[neighbor]) {
                                dfsStack.push(neighbor);
                            }
                        }
                    }
                } else {
                    dfsStack.pop();
                    if (onStack[curr]) {
                        finishStack.push(curr);
                        onStack[curr] = false;
                    }
                }
            }
        }

        inline void getTranspose(const std::unordered_map<int, std::vector<int> > &graph,
                                 std::unordered_map<int, std::vector<int> > &transposeGraph) {
            transposeGraph.clear();
            for (auto &pr: graph) {
                for (int neighbor: pr.second) {
                    transposeGraph[neighbor].push_back(pr.first);
                }
            }
        }

        inline void dfsTransposed(const std::unordered_map<int, std::vector<int> > &graph,
                                  int node,
                                  std::unordered_set<int> &visited,
                                  std::vector<int> &component) {
            std::stack<int> dfsStack;
            dfsStack.push(node);

            while (!dfsStack.empty()) {
                int curr = dfsStack.top();
                dfsStack.pop();
                if (!visited.count(curr)) {
                    visited.insert(curr);
                    component.push_back(curr);
                    if (graph.count(curr)) {
                        for (int neighbor: graph.at(curr)) {
                            if (!visited.count(neighbor)) {
                                dfsStack.push(neighbor);
                            }
                        }
                    }
                }
            }
        }

        inline void kosaraju(const std::unordered_map<int, std::vector<int> > &graph,
                             std::vector<std::vector<int> > &scc) {
            std::stack<int> orderStack;
            std::unordered_set<int> visited;
            scc.clear();

            for (auto &pr: graph) {
                if (!visited.count(pr.first)) {
                    fillOrder(graph, pr.first, visited, orderStack);
                }
            }

            std::unordered_map<int, std::vector<int> > transposedGraph;
            getTranspose(graph, transposedGraph);

            visited.clear();
            while (!orderStack.empty()) {
                int curr = orderStack.top();
                orderStack.pop();
                if (!visited.count(curr)) {
                    std::vector<int> component;
                    dfsTransposed(transposedGraph, curr, visited, component);
                    scc.emplace_back(component);
                }
            }
        }

        inline void buildCondensedGraph(const std::vector<std::vector<int> > &scc,
                                        const std::unordered_map<int, std::vector<int> > &graph,
                                        std::unordered_map<int, std::vector<int> > &condensedGraph) {
            std::unordered_map<int, int> nodeToComponent;
            int index = 0;
            for (const auto &component: scc) {
                for (int node: component) {
                    nodeToComponent[node] = index;
                }
                index++;
            }

            condensedGraph.clear();
            for (int i = 0; i < scc.size(); ++i) {
                std::unordered_set<int> seen;
                for (int node: scc[i]) {
                    if (graph.find(node) != graph.end()) {
                        for (int neighbor: graph.at(node)) {
                            int neighborComponent = nodeToComponent[neighbor];
                            if (neighborComponent != i && seen.find(neighborComponent) == seen.end()) {
                                condensedGraph[i].push_back(neighborComponent);
                                seen.insert(neighborComponent);
                            }
                        }
                    }
                }
            }
        }

        inline void topologicalSort(const std::unordered_map<int, std::vector<int> > &graph, std::vector<int> &order) {
            std::unordered_map<int, int> inDegree;
            for (const auto &pr: graph) {
                for (int neighbor: pr.second) {
                    inDegree[neighbor]++;
                }
            }

            std::queue<int> q;
            for (const auto &pr: graph) {
                if (inDegree[pr.first] == 0) {
                    q.push(pr.first);
                }
            }

            order.clear();
            while (!q.empty()) {
                int node = q.front();
                q.pop();
                order.push_back(node);
                if (graph.find(node) == graph.end())continue;
                for (int neighbor: graph.at(node)) {
                    if (--inDegree[neighbor] == 0) {
                        q.push(neighbor);
                    }
                }
            }
        }
    }


    template<bool dir>
    void CVRP_Pricing::getTopologicalOrder4OneBin(int b) {
        auto &order = dir ? topological_order_forward_ptr->at(b) : topological_order_backward_ptr->at(b);
        order.clear();
        Resource f_r{};
        dir ? f_r.resources[0] = {b * step_size} : f_r.resources[0] = {(b + 1) * step_size - 1};
        res_int max_res;
        dir ? max_res = (b + 1) * step_size : max_res = b * step_size - 1;
        std::unordered_map<int, std::vector<int> > graph;
        for (int i = 1; i < dim; ++i) {
            routeOptLong connect = 0;
            for (auto &j: dir
                              ? all_forward_buckets[i][b].bucket_arcs
                              : all_backward_buckets[i][b].bucket_arcs
            )
                connect.set(j);
            graph[i] = std::vector<int>();
            for (int j = 1; j < dim; ++j) {
                Resource tmp_res{};
                if (!connect.test(j)) continue;
                if (i == j) continue;
                dir
                    ? increaseMainResourceConsumption(f_r, tmp_res, i, j)
                    : decreaseMainResourceConsumption(f_r, tmp_res, i, j);
                if (dir ? tmp_res.resources[0] >= max_res : tmp_res.resources[0] <= max_res) continue;
                graph[i].emplace_back(j);
            }
        }
        std::vector<std::vector<int> > scc;
        detail::kosaraju(graph, scc);
        std::unordered_map<int, std::vector<int> > condensedGraph;
        detail::buildCondensedGraph(scc, graph, condensedGraph);
        std::vector<int> topoOrder;
        detail::topologicalSort(condensedGraph, topoOrder);
        std::unordered_set<int> tmpSet;
        tmpSet.reserve(topoOrder.size());
        order.resize(scc.size());
        for (int i = 0; i < topoOrder.size(); ++i) {
            int index = topoOrder[i];
            order[i] = scc[index];
            tmpSet.emplace(index);
        }
        std::vector<std::vector<int> > tmp;
        tmp.reserve(scc.size());
        for (int i = 0; i < scc.size(); ++i) {
            if (tmpSet.find(i) == tmpSet.end()) {
                tmp.emplace_back(scc[i]);
            }
        }
        std::sort(tmp.begin(), tmp.end(), [](const std::vector<int> &a, const std::vector<int> &b) {
            return a.size() > b.size();
        });

        auto index = topoOrder.size();
        for (auto &vec: tmp)order[index++] = vec;
    }

    template<bool if_symmetry>
    void CVRP_Pricing::getTopologicalOrder() {
        topological_order_forward_ptr->resize(num_buckets_per_vertex);
        for (int b = 0; b < num_buckets_per_vertex; ++b) {
            getTopologicalOrder4OneBin<true>(b);
        }
        if constexpr (!if_symmetry) {
            topological_order_backward_ptr->resize(num_buckets_per_vertex);
            for (int b = 0; b < num_buckets_per_vertex; ++b) {
                getTopologicalOrder4OneBin<false>(b);
            }
        }
    }
}

#endif // ROUTE_OPT_GET_TOPOLOGICAL_ORDER_HPP
