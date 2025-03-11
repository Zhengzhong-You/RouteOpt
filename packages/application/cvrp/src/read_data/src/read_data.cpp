/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "read_data_controller.hpp"
#include "read_data_macro.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRP_ReadDataController::generateInstancePath(int argc, char *argv[]) {
        int n = READ_NO_LINE;
        bool if_type1 = false;
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-d" && i + 1 < argc) {
                f_name_ref.get() = argv[++i];
                if_type1 = true;
            } else if (arg == "-n" && i + 1 < argc) {
                std::stringstream convert(argv[++i]);
                if (!(convert >> n)) {
                    printf("Invalid number: %s\n", argv[i]);
                    continue;
                }
                if_type1 = true;
            } else if (arg == "-u" && i + 1 < argc) {
                std::stringstream convert(argv[++i]);
                if (!(convert >> ub_ref.get())) {
                    printf("Invalid number: %s\n", argv[i]);
                }
            } else if (arg == "-b" && i + 1 < argc) {
                std::stringstream convert(argv[++i]);
                if (!(convert >> tree_path_ref.get())) {
                    //share the tree path
                    printf("Invalid std::string: %s\n", argv[i]);
                } else std::cout << tree_path_ref.get() << std::endl;
            }
        }

        if (if_type1) {
            readInstanceFile(n);
            return;
        }
        if (argc > 1) {
            f_name_ref.get() = argv[1];
            return;
        }
        throw std::invalid_argument("Not enough arguments");
    }

    void splitString(const std::string &str, const std::string &splits,
                     std::vector<std::string> &res) {
        if (str.empty()) return;
        auto strs = str + splits;
        size_t pos = strs.find(splits);
        auto step = static_cast<int>(splits.size());
        while (pos != std::string::npos) {
            std::string tmp = strs.substr(0, pos);
            res.push_back(tmp);
            strs = strs.substr(pos + step, strs.size());
            pos = strs.find(splits);
        }
    }

    void CVRP_ReadDataController::readInstanceFile(int line) {
        if (line == READ_NO_LINE)
            THROW_RUNTIME_ERROR("READ_NO_LINE error");
        std::fstream file(f_name_ref.get());
        if (!file.is_open()) {
            THROW_RUNTIME_ERROR("File not found: " + f_name_ref.get());
        }
        file.seekg(std::ios::beg);
        for (int i = 0; i < line; ++i) {
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        std::string line_str;
        getline(file, line_str);
        std::vector<std::string> strList;
        splitString(line_str, " ", strList);
        if (strList.size() == 2) {
            ub_ref.get() = strtod(strList[1].c_str(), nullptr);
        } else
            THROW_RUNTIME_ERROR("Read error");
        file.close();
        f_name_ref.get() = strList[0];
    }

    void CVRP_ReadDataController::takeDataFromFile() {
        file.open(f_name_ref.get());
        if (!file.is_open()) {
            THROW_RUNTIME_ERROR("File not found: " + f_name_ref.get());
        }

        extractName();

        if constexpr (app_type == APPLICATION_TYPE::VRPTW) {
            parseVRPTWSecondType();
            if (dim_ref.get() == 0) parseVRPTWOrRVRPSTW();
        } else {
            parseCVRPSolver();
        }

        file.close();
    }

    void CVRP_ReadDataController::extractName() {
        bool if_pop = false;
        std::string tmp_name;
        ins_name_ref.get().clear();
        for (auto i = f_name_ref.get().end() - 2; i >= f_name_ref.get().begin(); --i) {
            if (*(i + 1) == '.') {
                if_pop = true;
            }
            if (*i == '/' || *i == '\\') break;
            if (if_pop) {
                ins_name_ref.get().push_back(*i);
            }
            tmp_name.push_back(*i);
        }
        if (!if_pop) ins_name_ref.get() = tmp_name;
        std::reverse(ins_name_ref.get().begin(), ins_name_ref.get().end());
    }

    void CVRP_ReadDataController::parseVRPTWOrRVRPSTW() {
        std::string tmp_string;
        while (true) {
            if (tmp_string.find("NUMBER") != std::string::npos) {
                std::getline(file, tmp_string);
                std::stringstream temp_ss(tmp_string);
                temp_ss >> num_vehicle_ref.get();
                temp_ss >> cap_ref.get();
                break;
            }
            std::getline(file, tmp_string);
        }

        while (true) {
            std::getline(file, tmp_string);
            if (tmp_string.find("CUST NO.") != std::string::npos) {
                std::getline(file, tmp_string);
                break;
            }
        }

        int num = MAX_NUM_CUSTOMERS, tmp_j = 0;
        auto &info_vertex = info_vertex_ref.get();
        info_vertex.resize(num, std::vector<double>(7));
        while (std::getline(file, tmp_string) && tmp_j < num) {
            if (tmp_string.empty()) break;
            std::stringstream temp_ss(tmp_string);
            temp_ss >> info_vertex[tmp_j][0]; // index
            temp_ss >> info_vertex[tmp_j][1]; // x
            temp_ss >> info_vertex[tmp_j][2]; // y
            temp_ss >> info_vertex[tmp_j][3]; // demand
            temp_ss >> info_vertex[tmp_j][4]; // ready time
            temp_ss >> info_vertex[tmp_j][5]; // due date
            temp_ss >> info_vertex[tmp_j][6]; // service time
            ++tmp_j;
        }
        PRINT_DEBUG("");
        cap_ref.get() = 300;
        dim_ref.get() = tmp_j;
        info_vertex.resize(dim_ref.get());
    }

    bool checkForTravelTimeMatrix(std::ifstream &file) {
        std::string line;
        while (std::getline(file, line)) {
            if (line.find("Travel Time Matrix") != std::string::npos) {
                return true;
            }
        }
        return false;
    }

    void CVRP_ReadDataController::parseVRPTWSecondType() {
        file.clear();
        file.seekg(0, std::ios::beg);

        if (!checkForTravelTimeMatrix(file)) {
            file.clear();
            file.seekg(0, std::ios::beg);
            return;
        }

        file.clear();
        file.seekg(0, std::ios::beg);

        std::string tmp_string;

        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            if (line.find("Dimension:") != std::string::npos) {
                ss.ignore(11); // Ignore the label part
                ss >> dim_ref.get();
            } else if (line.find("Capacity:") != std::string::npos) {
                ss.ignore(10); // Ignore the label part
                ss >> cap_ref.get();
            } else if (line.find("ub:") != std::string::npos) {
                ss.ignore(4); // Ignore the label part
                ss >> ub_ref.get();
            } else if (line.find("Vertex Information:") != std::string::npos) {
                break; // Stop after reading all required headers
            }
        }

        std::getline(file, line); // Skip the header line with ID, X, Y, etc.
        auto &info_vertex = info_vertex_ref.get();
        info_vertex.resize(dim_ref.get(), std::vector<double>(7));
        int index = 0;
        while (std::getline(file, line) && index < dim_ref.get()) {
            if (line.empty()) continue; // Skip empty lines
            std::stringstream ss(line);
            for (int j = 0; j < 7; ++j) {
                ss >> info_vertex[index][j];
            }
            ++index;
        }
    }

    void CVRP_ReadDataController::parseCVRPSolver() {
        std::string line, name, tmp_string, tmp_string_2;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string key;
            iss >> key;
            if (key == "NAME") {
                iss >> key; // Discard the ":"
                iss >> name;
                if (name.find('k') != std::string::npos) {
                    size_t k_pos = name.find_last_of('k');
                    num_vehicle_ref.get() = std::stoi(name.substr(k_pos + 1));
                } else {
                    num_vehicle_ref.get() = -1;
                }
            } else if (key == "DIMENSION") {
                iss >> key; // Discard the ":"
                iss >> dim_ref.get();
                break;
            }
        }

        if (num_vehicle_ref.get() == -1) {
            num_vehicle_ref.get() = dim_ref.get() - 1;
        }

        std::getline(file, tmp_string);
        while (true) {
            if (tmp_string.find("CAPACITY") != std::string::npos) {
                std::stringstream temp_ss(tmp_string);
                temp_ss >> tmp_string_2;
                temp_ss >> tmp_string_2;
                temp_ss >> cap_ref.get();
                break;
            }
            std::getline(file, tmp_string);
        }

        int tmp_j = 0;
        auto &info_vertex = info_vertex_ref.get();
        info_vertex.resize(dim_ref.get(), std::vector<double>(4));
        std::getline(file, tmp_string);
        while (true) {
            std::getline(file, tmp_string);
            if (tmp_string.find("DEMAND_SECTION") != std::string::npos || tmp_j == dim_ref.get()) {
                break;
            }
            std::stringstream temp_ss(tmp_string);
            temp_ss >> info_vertex[tmp_j][0];
            temp_ss >> info_vertex[tmp_j][1];
            temp_ss >> info_vertex[tmp_j][2];
            ++tmp_j;
        }

        while (true) {
            if (tmp_string.find("DEMAND_SECTION") != std::string::npos) {
                break;
            }
            std::getline(file, tmp_string);
        }

        int tmp_int;
        tmp_j = 0;
        while (true) {
            std::getline(file, tmp_string);
            if (tmp_string.find("DEPOT_SECTION") != std::string::npos || tmp_j == dim_ref.get()) {
                break;
            }
            std::stringstream temp_ss(tmp_string);
            temp_ss >> tmp_int;
            temp_ss >> info_vertex[tmp_j][3];
            ++tmp_j;
        }
        info_vertex.resize(dim_ref.get());
    }
}
