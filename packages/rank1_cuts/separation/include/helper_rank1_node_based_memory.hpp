/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HELPER_RANK1_NODE_BASED_MEMORY_HPP
#define ROUTE_OPT_HELPER_RANK1_NODE_BASED_MEMORY_HPP
#include "rank1_data_shared.hpp"
#include "solver.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    struct other_ {
        int beg{};
        int tor{};
        std::vector<int> left_c{};
        std::vector<int> mem_c{};

        other_(int beg, int tor, std::vector<int> left_c,
               std::vector<int> mem_c) : beg(beg), tor(tor),
                                         left_c(std::move(left_c)),
                                         mem_c(std::move(mem_c)) {
        }

        other_() = default;
    };

    void constructMemoryVertexBased(const Rank1CutsDataShared &rank1CutsDataShared,
                                    DataShared &sharedData,
                                    PRICING_HARD_LEVEL pricing_hard_level, Solver &solver,
                                    std::unordered_map<std::vector<int>,
                                        std::vector<std::vector<int> >,
                                        VectorHashInRank1> &
                                    rank1_multi_mem_plan_map);

    void findMem(const std::vector<std::vector<std::vector<int> > > &arr,
                 const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                 std::unordered_set<int> &mem
    );

    void findPlan4R1CMulti(const Rank1CutsDataShared &rank1CutsDataShared,
                           std::unordered_map<std::vector<int>,
                               std::vector<std::vector<int> >, VectorHashInRank1> &
                           rank1_multi_mem_plan_map,
                           const std::vector<int> &vis, int denominator, cutLong &mem,
                           std::vector<std::unordered_set<int> > &segment,
                           std::vector<std::vector<int> > &plan);

    void combinationUtil(const std::vector<int> &arr,
                         std::vector<int> &tmp,
                         std::vector<std::vector<int> > &data,
                         int start,
                         int end,
                         int index,
                         int r);

    void combinations(const std::vector<std::vector<std::vector<int> > > &arr,
                      const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                      int i,
                      const std::vector<int> &accum,
                      const std::unordered_set<int> &mem,
                      int &record_min,
                      std::unordered_set<int> &new_mem
    );

    void getMemoryByMIP(const Rank1CutsDataShared &rank1CutsDataShared,
                        Solver &solver,
                        const std::vector<std::vector<std::vector<int> > > &arr,
                        const std::vector<std::vector<std::unordered_set<int> > >
                        &vec_segment,
                        std::unordered_set<int> &mem, bool &if_suc);
} // namespace RouteOpt::Rank1Cuts::Separation


#endif // ROUTE_OPT_HELPER_RANK1_NODE_BASED_MEMORY_HPP
