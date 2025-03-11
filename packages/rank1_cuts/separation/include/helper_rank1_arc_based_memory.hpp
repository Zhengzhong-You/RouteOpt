/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HELPER_RANK1_ARC_BASED_MEMORY_HPP
#define ROUTE_OPT_HELPER_RANK1_ARC_BASED_MEMORY_HPP
#include "rank1_separation_macro.hpp"
#include "rank1_data_shared.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    struct State {
        int coeff{};
        int state{};
        int end_segment{};
        std::bitset<MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN> bit{};
    };

    struct Arcs {
        std::vector<std::unordered_set<std::pair<int, int>, PairHasher>
        >
        arc_plan;
    };

    void constructMemoryArcBased(const Rank1CutsDataShared &rank1CutsDataShared, DataShared &sharedData);

    void findLeastMemoryArcBased(const Rank1CutsDataShared &rank1CutsDataShared, const DataShared &sharedData,
                                 const sparseRowMatrixXI &sol_matrix, R1c &cut, bool &if_suc);


    void reduceArcs(std::vector<Arcs> &all_arcs,
                    std::unordered_set<std::pair<int, int>, PairHasher> &existing_arcs,
                    bool &if_suc);

    void getLeastMemory(std::vector<Arcs> &all_arcs,
                        std::unordered_set<std::pair<int, int>, PairHasher> &existing_arcs,
                        bool &if_suc);

    void findLeastPlans2MakeCoeffRight(const std::vector<int> &vertex_states,
                                       int denominator,
                                       std::vector<std::vector<int> > &plans,
                                       bool &if_suc);

    void getVertexStates(const std::vector<int> &sequence,
                         int forward_pos,
                         const std::unordered_map<int, int> &mp_,
                         std::vector<int> &vertex_states,
                         std::vector<std::unordered_set<std::pair<int, int>,
                             PairHasher> > &arcs);
} // namespace RouteOpt::Rank1Cuts::Separation


#endif // ROUTE_OPT_HELPER_RANK1_ARC_BASED_MEMORY_HPP
