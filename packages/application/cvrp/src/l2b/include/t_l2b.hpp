/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_HPP
#define ROUTE_OPT_T_L2B_HPP
#include <numeric>
#include <iostream>
#include "route_opt_macro.hpp"
#include "l2b_controller.hpp"
#include "l2b_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::getMaxEdgeCost() {
        max_edge_cost = 0;
        auto real_dim = dim - 1;
        for (int i = 0; i < real_dim; ++i) {
            double pre_cost = *std::max_element(cost_mat4_vertex_ref.get()[i].begin() + i + 1,
                                                cost_mat4_vertex_ref.get()[i].end());
            max_edge_cost = std::max(max_edge_cost, pre_cost);
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::getMidPointEdgeCord(
        const std::vector<CustomerLocation> &customer_location) {
        mid_point_edge_cord.resize(dim, std::vector<std::pair<double, double> >(dim));
        for (int i = 0; i < dim; ++i) {
            for (int j = i; j < dim; ++j) {
                mid_point_edge_cord[i][j].first = (customer_location[i].x + customer_location[j].x) / 2;
                mid_point_edge_cord[i][j].second = (customer_location[i].y + customer_location[j].y) / 2;
                mid_point_edge_cord[j][i] = mid_point_edge_cord[i][j];
            }
        }
        mid_point_edge_cord_2_depot.resize(dim, std::vector<double>(dim));
        for (int i = 0; i < dim; ++i) {
            for (int j = i; j < dim; ++j) {
                mid_point_edge_cord_2_depot[i][j] =
                        std::sqrt(float(
                            (mid_point_edge_cord[i][j].first - customer_location[0].x)
                            * (mid_point_edge_cord[i][j].first - customer_location[0].x) +
                            (mid_point_edge_cord[i][j].second - customer_location[0].y)
                            * (mid_point_edge_cord[i][j].second - customer_location[0].y)));
                mid_point_edge_cord_2_depot[j][i] = mid_point_edge_cord_2_depot[i][j];
            }
        }
        max_mid_point_edge_cord_2_depot = 0;
        auto real_dim = dim - 1;
        for (int i = 0; i < real_dim; ++i) {
            double max_dis = *std::max_element(mid_point_edge_cord_2_depot[i].begin() + i + 1,
                                               mid_point_edge_cord_2_depot[i].end());
            max_mid_point_edge_cord_2_depot = std::max(max_mid_point_edge_cord_2_depot, max_dis);
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::getEdge2OtherConvertDis() {
        double geo_dis = std::exp(std::accumulate(cost_mat4_vertex_ref.get()[0].begin() + 1,
                                                  cost_mat4_vertex_ref.get()[0].end(),
                                                  0.0,
                                                  [](double a, double b) {
                                                      return a + std::log(b);
                                                  }) / (dim - 1));
        std::vector<routeOptLong> v_neighbor(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (cost_mat4_vertex_ref.get()[i][j] <= geo_dis) {
                    v_neighbor[i].set(j);
                }
            }
        }
        std::vector<std::vector<routeOptLong> > density_std_dis(dim, std::vector<routeOptLong>(dim));
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                density_std_dis[i][j] = v_neighbor[i] & v_neighbor[j];
                density_std_dis[j][i] = density_std_dis[i][j];
            }
        }
        node_density_in_std_dis_vec_form.resize(dim, std::vector<std::vector<int> >(dim));
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                for (int k = 0; k < dim; ++k) {
                    if (density_std_dis[i][j].test(k)) {
                        node_density_in_std_dis_vec_form[i][j].emplace_back(k);
                    }
                }
                node_density_in_std_dis_vec_form[j][i] = node_density_in_std_dis_vec_form[i][j];
            }
        }

        edge_2_other_convert_dis.resize(dim, std::vector<double>(dim));
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                double aver_dis = 0;
                for (int k = 0; k < node_density_in_std_dis_vec_form[i][j].size(); ++k) {
                    for (int l = k + 1; l < node_density_in_std_dis_vec_form[i][j].size(); ++l) {
                        double dif_x =
                                mid_point_edge_cord[node_density_in_std_dis_vec_form[i][j][k]][
                                    node_density_in_std_dis_vec_form[i][j][l]].first
                                - mid_point_edge_cord[i][j].first;
                        double dif_y =
                                mid_point_edge_cord[node_density_in_std_dis_vec_form[i][j][k]][
                                    node_density_in_std_dis_vec_form[i][j][l]].second
                                - mid_point_edge_cord[i][j].second;
                        aver_dis += std::sqrt(static_cast<float>(dif_x * dif_x + dif_y * dif_y));
                    }
                }
                if (node_density_in_std_dis_vec_form[i][j].empty())
                    aver_dis = 0;
                else
                    aver_dis /=
                            (double) std::pow(node_density_in_std_dis_vec_form[i][j].size(), 2) / 2;
                edge_2_other_convert_dis[i][j] = aver_dis / geo_dis;
                edge_2_other_convert_dis[j][i] = edge_2_other_convert_dis[i][j];
            }
        }
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::calculateDisDepot2Center(
        const std::vector<CustomerLocation> &customer_location) {
        double center_x = 0;
        double center_y = 0;
        for (int i = 1; i < dim; ++i) {
            center_x += customer_location[i].x;
            center_y += customer_location[i].y;
        }
        center_x /= dim - 1;
        center_y /= dim - 1;
        depot_2_center = std::sqrt(
            static_cast<float>((center_x - customer_location[0].x) * (center_x - customer_location[0].x) +
                               (center_y - customer_location[0].y) * (center_y - customer_location[0].y)));
        depot_2_center /= max_edge_cost;
    }

    template<typename Node, typename BrCType, typename Hasher>
    Learning2Branch<Node, BrCType, Hasher>::Learning2Branch(int dim,
                                                            const std::vector<std::vector<double> > &
                                                            cost_mat4_vertex,
                                                            const std::vector<std::vector<Resource> > &
                                                            resource_across_arcs_in_forward_sense,
                                                            const std::vector<std::vector<Resource> > &
                                                            resource_across_arcs_in_backward_sense,
                                                            const std::function<std::vector<int> (Node *,
                                                                const BrCType &)> &getBrConstraintNonzeroIdx,
                                                            const std::vector<std::vector<double> > &
                                                            infor_vector, bool if_init): dim(dim),
        cost_mat4_vertex_ref(cost_mat4_vertex),
        resource_across_arcs_in_forward_sense_ref(
            resource_across_arcs_in_forward_sense),
        resource_across_arcs_in_backward_sense_ref(
            resource_across_arcs_in_backward_sense),
        getBrConstraintNonzeroIdx(
            getBrConstraintNonzeroIdx) {
        if (!if_init) return;
        std::vector<CustomerLocation> customer_location(dim);
        for (int i = 0; i < dim; ++i) {
            customer_location[i].x = infor_vector[i][1];
            customer_location[i].y = infor_vector[i][2];
        }
        getMaxEdgeCost();
        getMidPointEdgeCord(customer_location);
        getEdge2OtherConvertDis();
        calculateClusteringCoefficient(customer_location);
        calculateDisDepot2Center(customer_location);
    }

    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::distinguishBranchPairSource(
        const Branching::BranchingDataShared<BrCType, Hasher> &branching_data_shared,
        const Branching::BranchingHistory<BrCType, Hasher> &branching_history) {
        const auto &branch_pair = branching_data_shared.getBranchPair();
        branch_pair_from_pseudo.clear();
        branch_pair_from_fractional.clear();
        for (auto &edge: branch_pair) {
            if (branching_history.isRecordedCandidate(edge)) {
                branch_pair_from_pseudo.emplace_back(edge);
            } else {
                branch_pair_from_fractional.emplace_back(edge);
            }
        }
    }
}

#endif // ROUTE_OPT_T_L2B_HPP
