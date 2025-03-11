/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <bitset>
#include <iostream>
#include "rank1_macro.hpp"

namespace RouteOpt::Rank1Cuts {
    void Rank1CutsDataShared::generateOptimalMultiplier() {
        map_rank1_multiplier[MAX_RANK_ROW + 1].resize(7, {});
        for (int i = 1; i <= MAX_RANK_ROW; ++i) {
            auto &it = map_rank1_multiplier[i];
            it.resize(7, {});
            if (i % 2) {//plan 0
                std::get<0>(it[0]) = std::vector<int>(i, 1);
                std::get<1>(it[0]) = 2;
                std::get<2>(it[0]) = static_cast<int>(i / 2);
            }
            if ((i - 2) % 3 == 0 && i >= 5) {//plan 1
                std::get<0>(it[1]) = std::vector<int>(i, 1);
                std::get<1>(it[1]) = 3;
                std::get<2>(it[1]) = static_cast<int>(i / 3);
            }
            if (i >= 5) {
                auto &tmp = std::get<0>(it[2]);//plan 2
                tmp.resize(i, 1);
                tmp[0] = i - 3;
                tmp[1] = i - 3;
                tmp[2] = 2;
                std::get<1>(it[2]) = i - 2;
                std::get<2>(it[2]) = 2;
                auto &tmp2 = std::get<0>(it[3]);//plan 3
                tmp2.resize(i, 1);
                tmp2[0] = i - 2;
                tmp2[1] = i - 2;
                tmp2[2] = 2;
                tmp2[3] = 2;
                std::get<1>(it[3]) = i - 1;
                std::get<2>(it[3]) = 2;
                auto &tmp3 = std::get<0>(it[4]);//plan 4
                tmp3.resize(i, 1);
                tmp3[0] = i - 3;
                tmp3[1] = 2;
                std::get<1>(it[4]) = i - 1;
                std::get<2>(it[4]) = 1;
                auto &tmp4 = std::get<0>(it[6]);//plan 6
                tmp4.resize(i, 1);
                tmp4[0] = i - 2;
                tmp4[1] = 2;
                tmp4[2] = 2;
                std::get<1>(it[6]) = i;
                std::get<2>(it[6]) = 1;
            }
            if (i >= 4) {//plan 5
                auto &tmp5 = std::get<0>(it[5]);
                tmp5.resize(i, 1);
                tmp5[0] = i - 2;
                std::get<1>(it[5]) = i - 1;
                std::get<2>(it[5]) = 1;
            }
        }
    }
}
