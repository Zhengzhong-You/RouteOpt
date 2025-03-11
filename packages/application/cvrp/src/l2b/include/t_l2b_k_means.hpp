/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_L2B_K_MEANS_HPP
#define ROUTE_OPT_T_L2B_K_MEANS_HPP
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <numeric>
#include "l2b_macro.hpp"
#include "l2b_controller.hpp"

namespace RouteOpt::Application::CVRP {
    namespace L2BDetail {
        inline double calculateDistance(const CustomerLocation &a, const CustomerLocation &b) {
            return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
        }

        inline void kMeansPlusPlus(const std::vector<CustomerLocation> &data, int k,
                                   std::vector<CustomerLocation> &centroids,
                                   int seed) {
            auto n = static_cast<int>(data.size());
            centroids.resize(k);

            std::mt19937 gen(seed);
            std::uniform_int_distribution<> dis(0, n - 1);

            centroids[0] = data[dis(gen)];

            std::vector<double> minDistances(n, std::numeric_limits<double>::max());

            for (int i = 1; i < k; ++i) {
                double totalDist = 0.0;

                for (int j = 0; j < n; ++j) {
                    double dist = calculateDistance(data[j], centroids[i - 1]);
                    if (dist < minDistances[j]) {
                        minDistances[j] = dist;
                    }
                    totalDist += minDistances[j];
                }

                double target = std::uniform_real_distribution<>(0, totalDist)(gen);
                for (int j = 0; j < n; ++j) {
                    if ((target -= minDistances[j]) <= 0) {
                        centroids[i] = data[j];
                        break;
                    }
                }
            }
        }

        inline void kMeans(const std::vector<CustomerLocation> &data,
                           int k,
                           std::vector<int> &labels,
                           std::vector<CustomerLocation> &centroids,
                           int seed) {
            auto n = static_cast<int>(data.size());
            labels.assign(n, -1);

            kMeansPlusPlus(data, k, centroids, seed);

            bool changed = true;
            while (changed) {
                changed = false;

                for (int i = 0; i < n; ++i) {
                    double minDist = std::numeric_limits<double>::max();
                    int bestCluster = -1;

                    for (int j = 0; j < k; ++j) {
                        double dist = calculateDistance(data[i], centroids[j]);
                        if (dist < minDist) {
                            minDist = dist;
                            bestCluster = j;
                        }
                    }

                    if (labels[i] != bestCluster) {
                        labels[i] = bestCluster;
                        changed = true;
                    }
                }

                std::vector<CustomerLocation> newCentroids(k);
                std::vector<int> counts(k, 0);

                for (int i = 0; i < n; ++i) {
                    int cluster = labels[i];
                    newCentroids[cluster].x += data[i].x;
                    newCentroids[cluster].y += data[i].y;
                    counts[cluster]++;
                }

                for (int j = 0; j < k; ++j) {
                    if (counts[j] > 0) {
                        newCentroids[j].x /= counts[j];
                        newCentroids[j].y /= counts[j];
                    }
                }
                centroids = newCentroids;
            }
        }

        inline double calculateWCSS(const std::vector<CustomerLocation> &data,
                                    const std::vector<int> &labels,
                                    const std::vector<CustomerLocation> &centroids) {
            double wcss = 0.0;
            for (int i = 0; i < data.size(); ++i) {
                int cluster = labels[i];
                wcss += std::pow(calculateDistance(data[i], centroids[cluster]), 2);
            }
            return wcss;
        }

        inline int findOptimalCluster(const std::vector<CustomerLocation> &data, int maxK, int seed) {
            std::vector<double> wcssValues;

            for (int k = 1; k <= maxK; ++k) {
                std::vector<int> labels;
                std::vector<CustomerLocation> centroids;
                kMeans(data, k, labels, centroids, seed);
                double wcss = calculateWCSS(data, labels, centroids);
                wcssValues.emplace_back(wcss);
            }

            int optimalCluster = 1;
            double maxDelta = 0.0;

            for (int k = 1; k < static_cast<int>(wcssValues.size()) - 1; ++k) {
                double delta = wcssValues[k - 1] - wcssValues[k];
                if (delta > maxDelta) {
                    maxDelta = delta;
                    optimalCluster = k + 1;
                }
            }

            return optimalCluster;
        }
    }


    template<typename Node, typename BrCType, typename Hasher>
    void Learning2Branch<Node, BrCType, Hasher>::calculateClusteringCoefficient(
        const std::vector<CustomerLocation> &customer_location) {
        std::vector<double> k_vec(MAX_NUM_RUN_OPTIMAL_K);
        for (int i = 0; i < MAX_NUM_RUN_OPTIMAL_K; ++i) {
            k_vec[i] = L2BDetail::findOptimalCluster(customer_location, MAX_NUM_CLUSTER, ML_RANDOM_SEED);
        }
        auto optimalCluster = static_cast<int>(std::round(
            std::accumulate(k_vec.begin(), k_vec.end(), 0.0) / MAX_NUM_RUN_OPTIMAL_K));

        std::vector<int> labels;
        std::vector<CustomerLocation> centroids;
        L2BDetail::kMeans(customer_location, optimalCluster, labels, centroids, ML_RANDOM_SEED);

        std::unordered_map<int, std::vector<int> > cluster_map;
        for (int i = 0; i < customer_location.size(); ++i) {
            cluster_map[labels[i]].emplace_back(i);
        }

        double sum_dis_in_cluster = 0;
        int cnt = 0;
        for (auto &pr: cluster_map) {
            auto &vec = pr.second;
            double dis_in = 0;
            if (vec.size() <= 1) continue;
            ++cnt;
            for (int i = 0; i < vec.size(); ++i) {
                for (int j = i + 1; j < vec.size(); ++j) {
                    dis_in += L2BDetail::calculateDistance(customer_location[vec[i]], customer_location[vec[j]]);
                }
            }
            dis_in /= static_cast<double>(vec.size() * static_cast<int>(vec.size() - 1)) / 2;
            sum_dis_in_cluster += dis_in;
        }
        sum_dis_in_cluster /= cnt * max_edge_cost;

        cluster_coeff = 1 / sum_dis_in_cluster;
    }
}

#endif // ROUTE_OPT_T_L2B_K_MEANS_HPP
