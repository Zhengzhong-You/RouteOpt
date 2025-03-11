/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_RANK1_CUTS_SELECTOR_HPP
#define ROUTE_OPT_RANK1_CUTS_SELECTOR_HPP
#include "rank1_data_shared.hpp"


namespace RouteOpt::Rank1Cuts::Separation {
    class CutSelector {
    public:
        CutSelector(const Rank1CutsDataShared &rank1CutsDataShared, DataShared &sharedData) : rank1CutsDataShared_ref(
                std::ref(rank1CutsDataShared)),
            sharedData_ref(std::ref(sharedData)) {
        }

        void selectR1CsByVioNMemory();

        void cleanData() {
        };

        CutSelector() = delete;

        ~CutSelector() = default;

    private:
        //ref
        const std::reference_wrapper<const Rank1CutsDataShared> rank1CutsDataShared_ref;
        std::reference_wrapper<DataShared> sharedData_ref;
    };
}


#endif // ROUTE_OPT_RANK1_CUTS_SELECTOR_HPP
