//
// Created by You, Zhengzhong on 4/26/24.
//

#ifndef INCLUDE_BRANCHING_LOWER_BOUND_PRIORITIZE_BRANCHING_HPP_
#define INCLUDE_BRANCHING_LOWER_BOUND_PRIORITIZE_BRANCHING_HPP_


#include "branching.hpp"

class BestBoundFirstBranching : public BaseBranching {
public:
    static QueueTree bbt;
    static QueueTree sub_bbt;

    static void updateLowerBound();

    static void solve(QueueTree &tree);
};

#endif //INCLUDE_BRANCHING_LOWER_BOUND_PRIORITIZE_BRANCHING_HPP_
