/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_REGENERATE_BUCKET_GRAPH_HPP
#define ROUTE_OPT_REGENERATE_BUCKET_GRAPH_HPP
#include <stack>
#include "cvrp_pricing_controller.hpp"
#include "pricing_macro.hpp"

namespace RouteOpt::Application::CVRP {
    template<bool if_symmetry>
    void CVRP_Pricing::considerRegenerateBucketGraph(
        Bucket **&all_forward_buckets, Bucket **&all_backward_buckets,
        int &num_forward_bucket_arcs, int &num_forward_jump_arcs,
        int &num_backward_bucket_arcs, int &num_backward_jump_arcs) {
        double aver_ratio = ratio_dominance_checks_non_dominant.first / ratio_dominance_checks_non_dominant.
                            second;
        if (aver_ratio > BucketResizeFactorRatioDominanceChecksNonDominant) {
            double numArcs = num_forward_bucket_arcs + num_forward_jump_arcs;
            numArcs += num_backward_bucket_arcs + num_backward_jump_arcs;
            if (numArcs / dim < BucketResizeFactorNumBucketArcsPerVertex) {
                regenerateGraphBucket<if_symmetry>(all_forward_buckets, all_backward_buckets, num_forward_bucket_arcs,
                                                   num_forward_jump_arcs,
                                                   num_backward_bucket_arcs, num_backward_jump_arcs);
            }
        }
    }

    template<bool if_symmetry>
    void CVRP_Pricing::regenerateGraphBucket(
        Bucket **&all_forward_buckets, Bucket **&all_backward_buckets,
        int &num_forward_bucket_arcs, int &num_forward_jump_arcs,
        int &num_backward_bucket_arcs, int &num_backward_jump_arcs) {
        if (step_size / 2 < 1)
            return;
        if (step_size % 2) {
            std::cout << "regenerateGraphBucket is banned since step_size is odd!" << std::endl;
            return;
        }

        if_stop_arc_elimination = true;
        step_size /= 2;
        num_buckets_per_vertex *= 2;
        max_num_forward_graph_arc *= 2;
        num_forward_bucket_arcs *= 2;
        num_forward_jump_arcs *= 2;
        if constexpr (!if_symmetry) {
            num_backward_bucket_arcs *= 2;
            num_backward_jump_arcs *= 2;
            max_num_backward_graph_arc *= 2;
        }

        auto reallocate2DDouble = [this](double **&arr, int newSize) {
            for (int i = 0; i < dim; ++i) {
                delete[] arr[i];
                arr[i] = new double[newSize];
            }
        };

        auto reallocate2DListLabel = [this](ListLabel **&arr, int newSize) {
            for (int i = 0; i < dim; ++i) {
                delete[] arr[i];
                arr[i] = new ListLabel[newSize];
            }
        };

        auto regenerateExtraLabels = [this](VecLabel **oldExtra, int newSize) -> VecLabel **{
            VecLabel **newExtra = new VecLabel *[dim];
            for (int i = 0; i < dim; ++i) {
                newExtra[i] = new VecLabel[newSize];
                for (int j = 0; j < newSize; ++j)
                    newExtra[i][j].first.resize(oldExtra[i][j / 2].first.size() / 2);
            }
            return newExtra;
        };

        auto regenerateBuckets = [this](Bucket **oldBuckets, int newSize) -> Bucket **{
            Bucket **newBuckets = new Bucket *[dim];
            for (int i = 0; i < dim; ++i) {
                newBuckets[i] = new Bucket[newSize];
                for (int j = 0; j < newSize; ++j)
                    newBuckets[i][j] = oldBuckets[i][j / 2];
            }
            return newBuckets;
        };

        VecLabel **new_IfExistExtraLabelsInForwardSense = regenerateExtraLabels(
            if_exist_extra_labels_in_forward_sense, num_buckets_per_vertex);
        reallocate2DDouble(rc2_till_this_bin_in_forward_sense, num_buckets_per_vertex);
        reallocate2DDouble(rc2_bin_in_forward_sense, num_buckets_per_vertex);
        reallocate2DListLabel(label_array_in_forward_sense, num_buckets_per_vertex);
        for (int i = 0; i < dim; ++i)
            delete[] if_exist_extra_labels_in_forward_sense[i];
        delete[] if_exist_extra_labels_in_forward_sense;
        if_exist_extra_labels_in_forward_sense = new_IfExistExtraLabelsInForwardSense;

        if constexpr (!if_symmetry) {
            VecLabel **new_IfExistExtraLabelsInBackwardSense = regenerateExtraLabels(
                if_exist_extra_labels_in_backward_sense, num_buckets_per_vertex);
            reallocate2DDouble(rc2_till_this_bin_in_backward_sense, num_buckets_per_vertex);
            reallocate2DDouble(rc2_bin_in_backward_sense, num_buckets_per_vertex);
            reallocate2DListLabel(label_array_in_backward_sense, num_buckets_per_vertex);
            for (int i = 0; i < dim; ++i)
                delete[] if_exist_extra_labels_in_backward_sense[i];
            delete[] if_exist_extra_labels_in_backward_sense;
            if_exist_extra_labels_in_backward_sense = new_IfExistExtraLabelsInBackwardSense;
        }

        Bucket **new_AllForwardBuckets = regenerateBuckets(all_forward_buckets, num_buckets_per_vertex);
        for (int i = 0; i < dim; ++i)
            delete[] all_forward_buckets[i];
        delete[] all_forward_buckets;
        all_forward_buckets = new_AllForwardBuckets;

        if constexpr (!if_symmetry) {
            Bucket **new_AllBackwardBuckets = regenerateBuckets(all_backward_buckets, num_buckets_per_vertex);
            for (int i = 0; i < dim; ++i)
                delete[] all_backward_buckets[i];
            delete[] all_backward_buckets;
            all_backward_buckets = new_AllBackwardBuckets;
        }

        this->all_forward_buckets = all_forward_buckets;
        this->all_backward_buckets = all_backward_buckets;

        tell_which_bin4_arc_elimination_in_forward_sense.clear();
        populateTellWhichBin4ArcElimination<true>();
        if constexpr (!if_symmetry) {
            tell_which_bin4_arc_elimination_in_backward_sense.clear();
            populateTellWhichBin4ArcElimination<false>();
        }

        getTopologicalOrder<if_symmetry>();
        std::cout << SMALL_PHASE_SEPARATION;
        PRINT_REMIND("new generated bucket graph: #buckets per vertex: " + std::to_string(num_buckets_per_vertex));
        std::cout << SMALL_PHASE_SEPARATION;
    }
}

#endif // ROUTE_OPT_REGENERATE_BUCKET_GRAPH_HPP
