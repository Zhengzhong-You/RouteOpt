/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rcc_rc_controller.hpp
 * @brief Controller for RCC pricing and RC updates in RouteOpt.
 *
 * This header defines the RCCRCController class, which provides a static method to price
 * Rounded Cap Cuts (RCCs) based on dual values. The controller calculates changes in the cost matrix
 * for customers using the provided RCCs and dual prices.
 */

#ifndef ROUTE_OPT_RCC_RC_CONTROLLER_HPP
#define ROUTE_OPT_RCC_RC_CONTROLLER_HPP

namespace RouteOpt::RCCs::RCGetter {
    /**
     * @brief Controller for RCC pricing.
     *
     * The RCCRCController class is responsible for pricing RCCs by using dual values obtained from the LP relaxation.
     * Its static method updates the change in cost matrix for customers based on the RCCs and the dual price vector.
     */
    class RCCRCController {
    public:
        /**
         * @brief Prices RCCs and updates the change cost matrix for customers.
         *
         * This static method performs pricing for Rounded Cap Cuts (RCCs) by computing
         * the change in cost for each vertex based on the provided dual prices. It takes as input
         * a vector of RCCs, a dual price vector (pi_vector), and outputs an updated change cost matrix.
         *
         * @param rccs Vector of RCCs used in the pricing process.
         * @param pi_vector Vector of dual prices from the LP relaxation.
         * @param chg_cost_mat4_vertex [in,out] Matrix of cost changes for customers, which is updated based on RCC pricing.
         */
        static void priceRCC(const std::vector<Rcc> &rccs,
                             const std::vector<double> &pi_vector,
                             std::vector<std::vector<double> > &chg_cost_mat4_vertex);

        // Default constructor.
        RCCRCController() = default;

        // Default destructor.
        ~RCCRCController() = default;
    };
} // namespace RouteOpt::RCCs::RCGetter

#endif // ROUTE_OPT_RCC_RC_CONTROLLER_HPP
