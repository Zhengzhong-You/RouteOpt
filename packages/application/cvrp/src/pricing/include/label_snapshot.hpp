/*
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_LABEL_SNAPSHOT_HPP
#define ROUTE_OPT_LABEL_SNAPSHOT_HPP
#include <algorithm>
#include <list>
#include <span>
#include <vector>
#include "label.hpp"

namespace RouteOpt::Application::CVRP {
    // Preserve bucket order without copying labels; their reduced costs must remain unchanged.
    struct LabelSnapshot {
        std::vector<Label *> labels;
        std::vector<size_t> offsets;
        int num_buckets{};

        LabelSnapshot() = default;

        LabelSnapshot(std::list<Label *> **buckets, int dim, int num_buckets)
            : offsets(static_cast<size_t>(dim) * num_buckets + 1), num_buckets(num_buckets) {
            size_t index = 0;
            for (int i = 0; i < dim; ++i) {
                for (int b = 0; b < num_buckets; ++b, ++index) {
                    offsets[index + 1] = offsets[index] + buckets[i][b].size();
                }
            }
            labels.resize(offsets.back());
            index = 0;
            for (int i = 0; i < dim; ++i) {
                for (int b = 0; b < num_buckets; ++b, ++index) {
                    std::copy(buckets[i][b].begin(), buckets[i][b].end(), labels.begin() + offsets[index]);
                }
            }
        }

        [[nodiscard]] std::span<Label * const> bucket(int i, int b) const {
            const auto index = static_cast<size_t>(i) * num_buckets + b;
            return std::span<Label * const>(labels).subspan(offsets[index], offsets[index + 1] - offsets[index]);
        }
    };
}

#endif // ROUTE_OPT_LABEL_SNAPSHOT_HPP
