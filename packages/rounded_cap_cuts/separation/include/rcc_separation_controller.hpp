/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rcc_separation_controller.hpp
 * @brief Controller for Rounded Cap Cuts (RCCs) separation in RouteOpt.
 *
 * This header defines the RCCSeparationController class, which provides a static
 * method to generate RCCs based on instance dimensions, capacity, demand, and solution
 * information. The controller supports various options including keeping existing cuts,
 * enforcing the first form, and strengthening the generated cuts.
 */

#ifndef ROUTE_OPT_RCC_SEPARATION_CONTROLLER_HPP
#define ROUTE_OPT_RCC_SEPARATION_CONTROLLER_HPP

namespace RouteOpt::RCCs::Separation {
    /**
     * @brief Controller for the RCC separation process.
     *
     * The RCCSeparationController class is responsible for calling CVRPSEP package.
     * It provides a static method to collect new RCCs given the problem data and solution information.
     */
    class RCCSeparationController {
    public:
        /**
         * @brief Collects new RCCs based on the provided problem and solution data.
         *
         * This static method collects new Rounded Cap Cuts (RCCs) using the following parameters:
         * - The problem dimension (dim) which determines the size of the problem.
         * - A capacity parameter (cap) used in the cut generation logic.
         * - A demand vector (demand) representing the demand at each customer.
         * - A flag (if_keep_rcc) setting this flag to newly generated RCCs.
         * - A flag (if_strengthen_rcc) that indicates whether to apply strengthening form of RCC.
         * - A solution vector (sol_x) containing fractional values from the LP relaxation.
         * - A vector of SequenceInfo objects (sols) representing solution routes or sequences.
         * - A vector of existing RCCs (existing_rccs) that may be used to avoid duplicates.
         * - An output vector (new_rccs) to store the newly generated RCCs.
         *
         * @param dim Problem dimension (number of variables or customers).
         * @param cap Capacity constraint used in generating cuts.
         * @param demand Vector of demands for each customer.
         * @param if_keep_rcc Setter of flag of newly generated RCCs.
         * @param if_strengthen_rcc Flag to determine if RCCs should be strengthened form, i.e., form 3.
         * @param sol_x Vector containing fractional solution values.
         * @param sols Vector of SequenceInfo objects representing solution sequences.
         * @param existing_rccs Vector of already generated RCCs to check against duplicates.
         * @param new_rccs [in,out] Output vector for storing the new RCCs.
         */
        static void generateRCCs(int dim,
                                 double cap,
                                 const std::vector<double> &demand,
                                 bool if_keep_rcc,
                                 bool if_strengthen_rcc,
                                 const std::vector<double> &sol_x,
                                 const std::vector<SequenceInfo> &sols,
                                 const std::vector<Rcc> &existing_rccs,
                                 std::vector<Rcc> &new_rccs);

        // Delete the default constructor to enforce usage of the static method.
        RCCSeparationController() = delete;

        // Default destructor.
        ~RCCSeparationController() = default;
    };
} // namespace RouteOpt::RCCs::Separation

#endif // ROUTE_OPT_RCC_SEPARATION_CONTROLLER_HPP
