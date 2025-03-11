/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_RANK1_MEMORY_FINDER_HPP
#define ROUTE_OPT_RANK1_MEMORY_FINDER_HPP
#include <vector>
#include "solver.hpp"
#include "rank1_data_shared.hpp"


namespace RouteOpt::Rank1Cuts::Separation {
    class MemGenerator {
    public:
        MemGenerator(const Rank1CutsDataShared &rank1CutsDataShared, DataShared &sharedData) : rank1CutsDataShared_ref(
                std::ref(rank1CutsDataShared)),
            sharedData_ref(std::ref(sharedData)) {
        }

        void fillMemory() const;

        void findMemory4Cuts();

        //setters
        void setSolver(Solver &value) { solver.getEnv(&value); }
        void setPricingHardLevel(PRICING_HARD_LEVEL value) { pricing_hard_level = value; }


        //getters
        const Solver &getSolver() const { return solver; }
        auto getPricingHardLevel() const { return pricing_hard_level; }

        void cleanData() {
            rank1_multi_mem_plan_map.clear();
            pricing_hard_level = PRICING_HARD_LEVEL::EASY;
        }

        MemGenerator() = delete;

        ~MemGenerator() = default;

    private:
        //need to keep
        Solver solver{};

        //need to clean
        PRICING_HARD_LEVEL pricing_hard_level{};
        std::unordered_map<std::vector<int>, std::vector<std::vector<int> >, VectorHashInRank1>
        rank1_multi_mem_plan_map{};

        //ref
        const std::reference_wrapper<const Rank1CutsDataShared> rank1CutsDataShared_ref;
        std::reference_wrapper<DataShared> sharedData_ref;
    };
}

#endif // ROUTE_OPT_RANK1_MEMORY_FINDER_HPP
