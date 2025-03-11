/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_SET_RESOURCE_IN_BUCKET_GRAPH_HPP
#define ROUTE_OPT_SET_RESOURCE_IN_BUCKET_GRAPH_HPP
#include "cvrp_pricing_controller.hpp"


namespace RouteOpt::Application::CVRP {
    template<bool if_symmetry>
    void CVRP_Pricing::initializeBucketGraphForNode(
        Bucket **&all_forward_buckets, Bucket **&all_backward_buckets,
        int &num_forward_bucket_arcs, int &num_backward_bucket_arcs) {
        max_num_forward_graph_arc = num_buckets_per_vertex * (dim - 2) * (dim - 1);
        if (all_forward_buckets == nullptr) {
            all_forward_buckets = new Bucket *[dim];
            for (int i = 0; i < dim; ++i) {
                all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
            }
            std::vector<int> bucketArc(dim - 2);
            for (int i = 1; i < dim; ++i) {
                int tmp = 0;
                for (int j = 1; j < i; ++j) bucketArc[tmp++] = j;
                for (int j = i + 1; j < dim; ++j) bucketArc[tmp++] = j;
                for (int b = 0; b < num_buckets_per_vertex; ++b) {
                    auto &bucket = all_forward_buckets[i][b];
                    bucket.bucket_arcs = bucketArc;
                    bucket.i = i;
                }
            }
            auto &bucket = all_forward_buckets[0][0];
            bucket.bucket_arcs.resize(dim - 1);
            int cnt = 0;
            for (int i = 1; i < dim; ++i) {
                Resource res{};
                if (!increaseMainResourceConsumption({}, res, 0, i)) continue;
                bucket.bucket_arcs[cnt++] = i;
            }
            bucket.bucket_arcs.resize(cnt);
            this->all_forward_buckets = all_forward_buckets;
            num_forward_bucket_arcs = static_cast<int>(max_num_forward_graph_arc);
        } else {
            num_forward_bucket_arcs = 0;
            for (int i = 1; i < dim; ++i) {
                for (int b = 0; b < num_buckets_per_vertex; ++b) {
                    num_forward_bucket_arcs += static_cast<int>(all_forward_buckets[i][b].bucket_arcs.size());
                }
            }
        }


        if constexpr (!if_symmetry) {
            max_num_backward_graph_arc = num_buckets_per_vertex * (dim - 2) * (dim - 1);
            if (all_backward_buckets == nullptr) {
                all_backward_buckets = new Bucket *[dim];
                for (int i = 0; i < dim; ++i) {
                    all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
                }
                for (int i = 1; i < dim; ++i) {
                    for (int b = 0; b < num_buckets_per_vertex; ++b) {
                        all_backward_buckets[i][b] = all_forward_buckets[i][0];
                    }
                }
                auto &bucket_back = all_backward_buckets[0][0];
                int cnt = 0;
                bucket_back.bucket_arcs.resize(dim - 1);
                for (int i = 1; i < dim; ++i) {
                    Resource res{};
                    if (!decreaseMainResourceConsumption(resource, res, 0, i)) continue;
                    bucket_back.bucket_arcs[cnt++] = i;
                }
                bucket_back.bucket_arcs.resize(cnt);
                this->all_backward_buckets = all_backward_buckets;
                num_backward_bucket_arcs = static_cast<int>(max_num_backward_graph_arc);
            } else {
                num_backward_bucket_arcs = 0;
                for (int i = 1; i < dim; ++i) {
                    for (int b = 0; b < num_buckets_per_vertex; ++b) {
                        num_backward_bucket_arcs += static_cast<int>(all_backward_buckets[i][b].bucket_arcs.size());
                    }
                }
            }
        }

        getTopologicalOrder<if_symmetry>();
    }

    template<bool if_symmetry>
    void CVRP_Pricing::initializeBucketGraph() {
        if (step_size == 0) {
            auto tmp = static_cast<int>(static_cast<double>(resource.resources[0]) / num_buckets_per_vertex /
                                        std::pow(2, MAX_NUM_REGENERATE_BUCKET));
            step_size = static_cast<res_int>(std::max(tmp * std::pow(2, MAX_NUM_REGENERATE_BUCKET), 1.));
            //at least be 1!
            num_buckets_per_vertex = static_cast<int>(std::floor(
                                         resource.resources[0] / step_size)) + 1;
        }


        if (label_array_in_forward_sense == nullptr) {
            label_array_in_forward_sense = new ListLabel *[dim];
            for (int i = 0; i < dim; ++i) {
                label_array_in_forward_sense[i] = new ListLabel[num_buckets_per_vertex];
            }
        }

        if (if_exist_extra_labels_in_forward_sense == nullptr) {
            if_exist_extra_labels_in_forward_sense = new VecLabel *[dim];
            for (int i = 0; i < dim; ++i) {
                if_exist_extra_labels_in_forward_sense[i] = new VecLabel[num_buckets_per_vertex];
                for (int j = 0; j < num_buckets_per_vertex; ++j) {
                    if_exist_extra_labels_in_forward_sense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
                }
            }
        }

        if (rc2_till_this_bin_in_forward_sense == nullptr) {
            rc2_till_this_bin_in_forward_sense = new double *[dim];
            for (int i = 0; i < dim; ++i) {
                rc2_till_this_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
            }
        }

        if (rc2_bin_in_forward_sense == nullptr) {
            rc2_bin_in_forward_sense = new double *[dim];
            for (int i = 0; i < dim; ++i) {
                rc2_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
            }
        }


        if constexpr (!if_symmetry) {
            if (label_array_in_backward_sense == nullptr) {
                label_array_in_backward_sense = new ListLabel *[dim];
                for (int i = 0; i < dim; ++i) {
                    label_array_in_backward_sense[i] = new ListLabel[num_buckets_per_vertex];
                }
            }

            if (if_exist_extra_labels_in_backward_sense == nullptr) {
                if_exist_extra_labels_in_backward_sense = new VecLabel *[dim];
                for (int i = 0; i < dim; ++i) {
                    if_exist_extra_labels_in_backward_sense[i] = new VecLabel[num_buckets_per_vertex];
                    for (int j = 0; j < num_buckets_per_vertex; ++j) {
                        if_exist_extra_labels_in_backward_sense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
                    }
                }
            }

            if (rc2_till_this_bin_in_backward_sense == nullptr) {
                rc2_till_this_bin_in_backward_sense = new double *[dim];
                for (int i = 0; i < dim; ++i) {
                    rc2_till_this_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
                }
            }

            if (rc2_bin_in_backward_sense == nullptr) {
                rc2_bin_in_backward_sense = new double *[dim];
                for (int i = 0; i < dim; ++i) {
                    rc2_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
                }
            }
        }
    }

    template<bool if_symmetry>
    void CVRP_Pricing::setCapResourceInBucketGraph(double cap, const std::vector<double> &demand) {
        resource.resources[CapResourceIdx] = roundAndConvertResLong(cap);
        if constexpr (CapResourceIdx == 0)
            meet_point_resource_in_bi_dir =
                    static_cast<double>(resource.resources[CapResourceIdx]) / 2;
        if (resource_across_arcs_in_forward_sense.empty()) {
            resource_across_arcs_in_forward_sense.resize(dim);
            for (auto &vertex: resource_across_arcs_in_forward_sense) vertex.resize(dim);
        }
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if constexpr (!if_symmetry) {
                    resource_across_arcs_in_forward_sense[i][j].resources[CapResourceIdx] =
                            roundAndConvertResLong(demand[i]);
                } else {
                    resource_across_arcs_in_forward_sense[i][j].resources[CapResourceIdx] =
                            roundAndConvertResLong((demand[i] + demand[j]) / 2);
                }
            }
        }
        if constexpr (!if_symmetry) {
            if (resource_across_arcs_in_backward_sense.empty()) {
                resource_across_arcs_in_backward_sense.resize(dim);
                for (auto &vertex: resource_across_arcs_in_backward_sense) vertex.resize(dim);
            }
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    resource_across_arcs_in_backward_sense[i][j].resources[CapResourceIdx] = roundAndConvertResLong(
                        demand[j]);
                }
            }
        }
        if (lb4_vertex.empty()) lb4_vertex.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lb4_vertex[i].resources[CapResourceIdx] = 0;
        }
        if (ub4_vertex.empty()) ub4_vertex.resize(dim);
        for (int i = 0; i < dim; ++i) {
            ub4_vertex[i].resources[CapResourceIdx] = roundAndConvertResLong(cap);
        }
    }

    inline void CVRP_Pricing::setTWResourceInBucketGraph(double H, const std::vector<double> &e,
                                                         const std::vector<double> &l, const std::vector<double> &s,
                                                         const std::vector<std::vector<double> > &cost_mat4_vertex) {
        if constexpr (TWResourceIdx >= NUM_RESOURCE) {
            return;
        }
        resource.resources[TWResourceIdx] = roundAndConvertResLong(H);
        if constexpr (TWResourceIdx == 0)
            meet_point_resource_in_bi_dir =
                    static_cast<double>(resource.resources[TWResourceIdx]) / 2;
        if (resource_across_arcs_in_forward_sense.empty()) {
            resource_across_arcs_in_forward_sense.resize(dim);
            for (auto &vertex: resource_across_arcs_in_forward_sense) vertex.resize(dim);
        }
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                resource_across_arcs_in_forward_sense[i][j].resources[TWResourceIdx] =
                        roundAndConvertResLong(s[i] + cost_mat4_vertex[i][j]);
            }
        }

        if (resource_across_arcs_in_backward_sense.empty()) {
            resource_across_arcs_in_backward_sense.resize(dim);
            for (auto &vertex: resource_across_arcs_in_backward_sense) vertex.resize(dim);
        }
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                resource_across_arcs_in_backward_sense[i][j].resources[TWResourceIdx] =
                        roundAndConvertResLong(s[j] + cost_mat4_vertex[i][j]);
            }
        }
        if (lb4_vertex.empty()) lb4_vertex.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lb4_vertex[i].resources[TWResourceIdx] = roundAndConvertResLong(e[i]);
        }
        if (ub4_vertex.empty()) ub4_vertex.resize(dim);
        for (int i = 0; i < dim; ++i) {
            ub4_vertex[i].resources[TWResourceIdx] = roundAndConvertResLong(l[i]);
        }
    }
}

#endif // ROUTE_OPT_SET_RESOURCE_IN_BUCKET_GRAPH_HPP
