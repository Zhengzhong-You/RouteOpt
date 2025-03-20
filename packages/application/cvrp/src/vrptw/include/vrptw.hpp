/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_VRPTW_HPP
#define ROUTE_OPT_VRPTW_HPP
#include "cvrp.hpp"

namespace RouteOpt::Application::CVRP {
    class VRPTW : public CVRPSolver {
    public:
        VRPTW(int argc, char *argv[]): CVRPSolver(argc, argv) {
            std::cout << "VRPTW MODE= " << (vrptw_type == VRPTW_TYPE::SINGLE_RESOURCE
                                                ? "SINGLE_RESOURCE"
                                                : vrptw_type == VRPTW_TYPE::CAP_MAIN_RESOURCE
                                                      ? "CAP_MAIN_RESOURCE"
                                                      : "TW_MAIN_RESOURCE") << std::endl;
        }

        ~VRPTW() override = default;

    private:
        void checkSolutionFeasibility(const std::vector<double> &X,
                                      const std::vector<SequenceInfo> &cols,
                                      bool &feasible) override;

        void addFeasibilityCuts(int &num_row,
                                const std::vector<double> &x,
                                const std::vector<SequenceInfo> &cols,
                                std::vector<Rcc> &existing_rccs,
                                Solver &solver) override;
    };
}


#endif // ROUTE_OPT_VRPTW_HPP
