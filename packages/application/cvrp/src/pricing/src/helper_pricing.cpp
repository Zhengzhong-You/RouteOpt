/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    double CVRP_Pricing::getAverageRouteLength() {
        if (aver_route_length.getAverage() == 0) return EXPECTED_AVER_ROUTE_LENGTH;
        return aver_route_length.getAverage();
    }

    void CVRP_Pricing::resizePoolWarning(size_t &pricing_warning) {
        if (pool_beg4_pricing >= pricing_warning) {
            std::cout << SMALL_PHASE_SEPARATION;
            PRINT_REMIND("pricing pool is almost full!");
            std::cout << "pool_beg4_pricing=" << pool_beg4_pricing << std::endl;
            std::cout << "mem4_pricing=" << mem4_pricing << std::endl;
            std::cout << "we reallocate the pricing pool!" << std::endl;
            reallocatePricingPool();
            pricing_warning = static_cast<size_t>(PricingWarningThreshold * static_cast<double>(mem4_pricing));
            std::cout << "the new mem4_pricing=" << mem4_pricing << std::endl;
            std::cout << SMALL_PHASE_SEPARATION;
        }
    }
}
