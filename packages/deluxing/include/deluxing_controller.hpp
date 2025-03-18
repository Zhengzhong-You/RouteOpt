/*
 * Copyright (c) 2025 Yu Yang & Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file deluxing.hpp
 * @brief Header file for the DeLuxingController class.
 *
 * This header defines the DeLuxingController class, which is used to perform the DeLuxing algorithm
 * to reduce the number of enumerated columns.
 */

#ifndef ROUTE_OPT_DELUXING_HPP
#define ROUTE_OPT_DELUXING_HPP


#include <vector>
#include "solver.hpp"

namespace RouteOpt::DeLuxing {
    /**
     * @class DeLuxingController
     * @brief Provides functionality for performing the DeLuxing procedure to remove
     *        unnecessary variables in solver-based optimization tasks.
     */
    class DeLuxingController {
    public:
        /**
         * @brief Performs the DeLuxing procedure on the given solver.
         *
         * @param solver      A reference to the solver object representing the current model state.
         * @param UB          Upper bound on the objective value (often from a known feasible solution).
         * @param NClust      Number of clusters for grouping variables.
         * @param beta1       Threshold controlling the initial aggressiveness of variable removal.
         * @param beta2       Threshold controlling deeper iteration for variable removal.
         * @param Idxdel      A reference to a vector that will store the indices of removed variables.
         * @param Timelimit   The maximum time allowed (in seconds) for the DeLuxing procedure.
         * @param Tolerance   A small value to handle floating-point comparisons when fixing variables.
         * @param Verbose     If true, prints additional logs about the DeLuxing process.
         */
        static void deLuxing(Solver &solver,
                             double UB,
                             int NClust,
                             int beta1,
                             int beta2,
                             std::vector<int> &Idxdel,
                             double Timelimit,
                             double Tolerance,
                             bool Verbose);
    };
}


#endif
