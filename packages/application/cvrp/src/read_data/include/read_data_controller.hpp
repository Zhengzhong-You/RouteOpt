/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_READ_DATA_CONTROLLER_HPP
#define ROUTE_OPT_READ_DATA_CONTROLLER_HPP
#include <string>
#include <fstream>
#include "read_data_macro.hpp"
#include "cvrp_macro.hpp"

namespace RouteOpt::Application::CVRP {
    class CVRP_ReadDataController {
    public:
        CVRP_ReadDataController(int argc, char *argv[],
                                std::string &f_name, std::string &ins_name, std::string &tree_path, int &num_vehicle,
                                double &cap, int &dim,
                                std::vector<std::vector<double> > &info_vertex, double &ub): f_name_ref(f_name),
            ins_name_ref(ins_name),
            tree_path_ref(tree_path), num_vehicle_ref(num_vehicle),
            cap_ref(cap),
            dim_ref(dim),
            info_vertex_ref(info_vertex), ub_ref(ub) {
            getData(argc, argv);
        }


        CVRP_ReadDataController() = delete;

        ~CVRP_ReadDataController() = default;

    private:
        std::ifstream file{};
        //ref
        std::reference_wrapper<std::string> f_name_ref;
        std::reference_wrapper<std::string> ins_name_ref;
        std::reference_wrapper<std::string> tree_path_ref;
        std::reference_wrapper<int> num_vehicle_ref;
        std::reference_wrapper<double> cap_ref;
        std::reference_wrapper<int> dim_ref;
        std::reference_wrapper<std::vector<std::vector<double> > > info_vertex_ref;
        std::reference_wrapper<double> ub_ref;


        void getData(int argc, char *argv[]) {
            ub_ref.get() = INITIAL_UB;
            generateInstancePath(argc, argv);
            takeDataFromFile();
        }

        void takeDataFromFile();

        void extractName();

        void parseVRPTWOrRVRPSTW();

        void parseVRPTWSecondType();

        void parseCVRPSolver();

        void readInstanceFile(int line);

        void generateInstancePath(int argc, char *argv[]);
    };
}

#endif // ROUTE_OPT_READ_DATA_CONTROLLER_HPP
