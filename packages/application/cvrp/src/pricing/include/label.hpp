/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_LABEL_HPP
#define ROUTE_OPT_LABEL_HPP
#include "route_opt_macro.hpp"
#include "pricing_macro.hpp"
#include "cvrp_macro.hpp"
#include "rank1_rc_controller.hpp"

namespace RouteOpt::Application::CVRP {
    using res_int = int; //cannot be unsigned type
    /**
     * Resource tuple for label states
     * N: total number of resources, M: resources compared in dominance
     */
    template<int N, int M> //all resource, resource that is compared in dominance rule;
    struct ResTuple {
        std::array<res_int, N> resources;

        ResTuple operator+(const ResTuple &other) const {
            ResTuple result;
            for (int i = 0; i < N; ++i) {
                result.resources[i] = this->resources[i] + other.resources[i];
            }
            return result;
        }

        ResTuple operator-(const ResTuple &other) const {
            ResTuple result;
            for (int i = 0; i < N; ++i) {
                result.resources[i] = this->resources[i] - other.resources[i];
            }
            return result;
        }

        bool operator <=(const ResTuple &other) const {
            /**
         * first resource is checked at last since it is sorted to some extent
         */
            for (int i = 1; i < M; ++i) {
                if (resources[i] > other.resources[i]) {
                    return false;
                }
            }
            return resources[0] <= other.resources[0];
        }

        bool operator==(const ResTuple &other) const {
            /**
         * first resource is checked at last since it is sorted to some extent
         */
            for (int i = 1; i < M; ++i) {
                if (resources[i] != other.resources[i]) {
                    return false;
                }
            }
            return resources[0] == other.resources[0];
        }

        void takeLarger(const ResTuple &other) {
            adjust(other, std::less<>());
        }

        void takeSmaller(const ResTuple &other) {
            adjust(other, std::greater<>());
        }

    private:
        template<typename Compare>
        void adjust(const ResTuple &standard, Compare comp) {
            for (int i = 0; i < N; ++i) {
                if (comp(resources[i], standard.resources[i])) {
                    resources[i] = standard.resources[i];
                }
            }
        }
    };


    /**
     * Base label structure for dynamic programming states
     */
    template<typename DerivedLabel>
    struct LabelBase {
        bool is_extended{};
        int end_vertex{};
        double rc{};
        double cost{};
        routeOptLong pi{};
        Rank1Cuts::RCGetter::R1CPricingStat r1c{};
        DerivedLabel *p_label{};
        ResTuple<NUM_RESOURCE, NUM_RESOURCE_CMP_IN_DOMINANCE> res{};

        LabelBase() = default;

        virtual ~LabelBase() = default;
    };

    struct CVRPLabel : public LabelBase<CVRPLabel> {
    };

    using Label = CVRPLabel;
    using Resource = ResTuple<
        NUM_RESOURCE,
        NUM_RESOURCE_CMP_IN_DOMINANCE>;


    inline res_int roundAndConvertResLong(double value) {
        double rounded = std::round(value * RESOURCE_FACTOR) / RESOURCE_FACTOR;
        if (std::abs(rounded - value) > TOLERANCE)
            THROW_RUNTIME_ERROR("RESOURCE_FACTOR is too small");
        rounded *= RESOURCE_FACTOR;
        if (rounded > std::numeric_limits<res_int>::max()) {
            throw std::overflow_error("Value exceeds the range of res_int");
        }
        return static_cast<res_int>(rounded);
    }
}

#endif // ROUTE_OPT_LABEL_HPP
