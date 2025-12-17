/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_INITIAL_CONTROLLER_HPP
#define ROUTE_OPT_INITIAL_CONTROLLER_HPP
#include <functional>


#include "node.hpp"

namespace RouteOpt::Application::CVRP {
    /**
     * Controller for initial solution and preprocessing
     * Sets up initial LP model and root node
     */
    class InitialController {
    public:
        InitialController(int &dim,
                          Solver &solver,
                          std::vector<double> &demand,
                          double &H,
                          std::vector<double> &earliest_time,
                          std::vector<double> &latest_time,
                          std::vector<double> &service_time,
                          std::vector<std::vector<double> > &cost_mat4_vertex,
                          std::vector<std::vector<double> > &info_vertex,
                          const std::function<void()> &getLowerBoundFn)
            : dim(dim),
              demand_ref(demand),
              cost_mat4_vertex_ref(cost_mat4_vertex),
              info_vertex_ref(info_vertex),
              getLowerBoundofMinimumNumberCarsRef(getLowerBoundFn),
              solver_ref(solver) {
            initialProcessing();
            if (info_vertex[0].size() >= 7) {
                earliest_time.resize(dim);
                latest_time.resize(dim);
                service_time.resize(dim);
                for (int i = 0; i < dim; ++i) {
                    earliest_time[i] = info_vertex[i][4];
                    latest_time[i] = info_vertex[i][5];
                    service_time[i] = info_vertex[i][6];
                }
                H = info_vertex[0][5];
            }
        }

        InitialController() = delete;

        ~InitialController() = default;

    private:
        int dim{};
        const std::reference_wrapper<std::vector<double> > demand_ref;
        const std::reference_wrapper<std::vector<std::vector<double> > > cost_mat4_vertex_ref;
        const std::reference_wrapper<std::vector<std::vector<double> > > info_vertex_ref;
        std::function<void()> getLowerBoundofMinimumNumberCarsRef; //std::copy here but refer to the original function

        std::reference_wrapper<Solver> solver_ref;

        void setSolverEnv() const;

        void initialProcessing();
    };
}

#endif // ROUTE_OPT_INITIAL_CONTROLLER_HPP
