/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "node.hpp"

namespace RouteOpt::Application::CVRP {
    void BbNode::optimizeLPForOneIteration(double &prior_value, bool if_allow_delete_col, int lp_method) {
        int num_col;
        SAFE_SOLVER(solver.reoptimize(lp_method))
        SAFE_SOLVER(solver.getNumCol(&num_col))

        double val;
        SAFE_SOLVER(solver.getObjVal(&val))

        if (!if_in_enu_state && if_allow_delete_col && num_col > LP_COL_FINAL_LIMIT) {
            double tol = std::max(TOLERANCE * val, TOLERANCE);
            if (std::abs(prior_value - val) > tol) {
                cleanIndexColForNode();
            }
        }
        prior_value = val;
    }
}
