/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_HGS_CONTROLLER_HPP
#define ROUTE_OPT_HGS_CONTROLLER_HPP
#include <string>
#include <vector>

namespace RouteOpt::Application::CVRP {
    class HGSController {
    public:
        HGSController(bool if_run, std::string &f_name, int &dim, double &ub,
                      std::vector<std::vector<int> > &ip_opt_sol): f_name(f_name),
                                                                   dim(dim), ub_ref(ub),
                                                                   ip_opt_sol_ref(ip_opt_sol) {
            if (if_run) runHGS();
        }

        HGSController() = delete;

        ~HGSController() = default;

    private:
        int dim{};
        std::string f_name{};

        std::reference_wrapper<double> ub_ref;
        std::reference_wrapper<std::vector<std::vector<int> > > ip_opt_sol_ref;

        void runHGS();
    };
}

#endif // ROUTE_OPT_HGS_CONTROLLER_HPP
