/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "read_data_controller.hpp"
#include "read_data_macro.hpp"
#include "route_opt_macro.hpp"

namespace RouteOpt::Application::CVRP {
    namespace {
        std::string trimCopy(const std::string &input) {
            const auto first = input.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) {
                return {};
            }
            const auto last = input.find_last_not_of(" \t\r\n");
            return input.substr(first, last - first + 1);
        }

        std::string toUpperAscii(std::string input) {
            std::transform(input.begin(), input.end(), input.begin(), [](const unsigned char ch) {
                return static_cast<char>(std::toupper(ch));
            });
            return input;
        }

        bool startsWithKey(const std::string &line, const std::string &key) {
            const auto trimmed = trimCopy(line);
            const auto upper_line = toUpperAscii(trimmed);
            const auto upper_key = toUpperAscii(key);
            if (upper_line.size() < upper_key.size() || upper_line.compare(0, upper_key.size(), upper_key) != 0) {
                return false;
            }
            if (upper_line.size() == upper_key.size()) {
                return true;
            }
            const unsigned char next = static_cast<unsigned char>(upper_line[upper_key.size()]);
            return std::isspace(next) || next == ':';
        }

        bool isSectionLine(const std::string &line, const std::string &section_name) {
            return toUpperAscii(trimCopy(line)) == toUpperAscii(section_name);
        }

        std::string extractHeaderValue(const std::string &line) {
            const auto colon_pos = line.find(':');
            if (colon_pos != std::string::npos) {
                return trimCopy(line.substr(colon_pos + 1));
            }

            std::istringstream iss(line);
            std::string key;
            iss >> key;
            std::string remainder;
            std::getline(iss, remainder);
            return trimCopy(remainder);
        }

        std::size_t findSectionIndex(const std::vector<std::string> &lines, const std::string &section_name) {
            for (std::size_t i = 0; i < lines.size(); ++i) {
                if (isSectionLine(lines[i], section_name)) {
                    return i;
                }
            }
            return lines.size();
        }

    }

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
        file.clear();
        file.seekg(0, std::ios::beg);

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(file, line)) {
            lines.push_back(line);
        }

        std::string name;
        bool has_capacity = false;
        bool has_vehicles = false;
        num_vehicle_ref.get() = -1;

        for (const auto &raw_line: lines) {
            if (startsWithKey(raw_line, "NAME")) {
                name = extractHeaderValue(raw_line);
                if (!has_vehicles) {
                    const auto k_pos = name.find_last_of('k');
                    if (k_pos != std::string::npos && k_pos + 1 < name.size()) {
                        num_vehicle_ref.get() = std::stoi(name.substr(k_pos + 1));
                    }
                }
            } else if (startsWithKey(raw_line, "VEHICLES")) {
                std::stringstream ss(extractHeaderValue(raw_line));
                if (!(ss >> num_vehicle_ref.get())) {
                    THROW_RUNTIME_ERROR("Invalid VEHICLES header: " + raw_line);
                }
                has_vehicles = true;
            } else if (startsWithKey(raw_line, "DIMENSION")) {
                std::stringstream ss(extractHeaderValue(raw_line));
                if (!(ss >> dim_ref.get())) {
                    THROW_RUNTIME_ERROR("Invalid DIMENSION header: " + raw_line);
                }
            } else if (startsWithKey(raw_line, "CAPACITY")) {
                std::stringstream ss(extractHeaderValue(raw_line));
                if (!(ss >> cap_ref.get())) {
                    THROW_RUNTIME_ERROR("Invalid CAPACITY header: " + raw_line);
                }
                has_capacity = true;
            }
        }

        if (dim_ref.get() <= 0) {
            THROW_RUNTIME_ERROR("Missing or invalid DIMENSION in CVRP instance: " + f_name_ref.get());
        }
        if (!has_capacity) {
            THROW_RUNTIME_ERROR("Missing CAPACITY in CVRP instance: " + f_name_ref.get());
        }
        if (num_vehicle_ref.get() == -1) {
            num_vehicle_ref.get() = dim_ref.get() - 1;
        }

        const auto node_coord_idx = findSectionIndex(lines, "NODE_COORD_SECTION");
        const auto demand_idx = findSectionIndex(lines, "DEMAND_SECTION");
        const auto depot_idx = findSectionIndex(lines, "DEPOT_SECTION");

        if (node_coord_idx == lines.size()) {
            THROW_RUNTIME_ERROR("Missing NODE_COORD_SECTION in CVRP instance: " + f_name_ref.get());
        }
        if (demand_idx == lines.size()) {
            THROW_RUNTIME_ERROR("Missing DEMAND_SECTION in CVRP instance: " + f_name_ref.get());
        }
        if (depot_idx == lines.size()) {
            THROW_RUNTIME_ERROR("Missing DEPOT_SECTION in CVRP instance: " + f_name_ref.get());
        }
        if (!(node_coord_idx < demand_idx && demand_idx < depot_idx)) {
            THROW_RUNTIME_ERROR("Unexpected CVRP section order in instance: " + f_name_ref.get());
        }

        auto &info_vertex = info_vertex_ref.get();
        info_vertex.assign(dim_ref.get(), std::vector<double>(4, 0.0));
        std::vector<bool> coord_seen(dim_ref.get(), false);
        std::vector<bool> demand_seen(dim_ref.get(), false);

        // Parse TSPLIB sections by node id so optional headers and line ordering do not shift data.
        for (std::size_t i = node_coord_idx + 1; i < demand_idx; ++i) {
            const auto trimmed = trimCopy(lines[i]);
            if (trimmed.empty() || isSectionLine(trimmed, "EOF")) {
                continue;
            }

            std::stringstream ss(trimmed);
            int node_id = 0;
            double x = 0.0;
            double y = 0.0;
            if (!(ss >> node_id >> x >> y)) {
                THROW_RUNTIME_ERROR("Malformed NODE_COORD_SECTION entry: " + lines[i]);
            }
            if (node_id < 1 || node_id > dim_ref.get()) {
                THROW_RUNTIME_ERROR("Coordinate node id out of range: " + lines[i]);
            }
            if (coord_seen[node_id - 1]) {
                THROW_RUNTIME_ERROR("Duplicate coordinate entry for node " + std::to_string(node_id));
            }

            info_vertex[node_id - 1][0] = node_id;
            info_vertex[node_id - 1][1] = x;
            info_vertex[node_id - 1][2] = y;
            coord_seen[node_id - 1] = true;
        }

        for (int node_id = 1; node_id <= dim_ref.get(); ++node_id) {
            if (!coord_seen[node_id - 1]) {
                THROW_RUNTIME_ERROR("Missing coordinate entry for node " + std::to_string(node_id));
            }
        }

        for (std::size_t i = demand_idx + 1; i < depot_idx; ++i) {
            const auto trimmed = trimCopy(lines[i]);
            if (trimmed.empty() || isSectionLine(trimmed, "EOF")) {
                continue;
            }

            std::stringstream ss(trimmed);
            int node_id = 0;
            double demand = 0.0;
            if (!(ss >> node_id >> demand)) {
                THROW_RUNTIME_ERROR("Malformed DEMAND_SECTION entry: " + lines[i]);
            }
            if (node_id < 1 || node_id > dim_ref.get()) {
                THROW_RUNTIME_ERROR("Demand node id out of range: " + lines[i]);
            }
            if (demand_seen[node_id - 1]) {
                THROW_RUNTIME_ERROR("Duplicate demand entry for node " + std::to_string(node_id));
            }

            info_vertex[node_id - 1][3] = demand;
            demand_seen[node_id - 1] = true;
        }

        for (int node_id = 1; node_id <= dim_ref.get(); ++node_id) {
            if (!demand_seen[node_id - 1]) {
                THROW_RUNTIME_ERROR("Missing demand entry for node " + std::to_string(node_id));
            }
        }

        std::vector<int> depot_nodes;
        for (std::size_t i = depot_idx + 1; i < lines.size(); ++i) {
            const auto trimmed = trimCopy(lines[i]);
            if (trimmed.empty() || isSectionLine(trimmed, "EOF")) {
                continue;
            }

            std::stringstream ss(trimmed);
            int depot_id = 0;
            if (!(ss >> depot_id)) {
                THROW_RUNTIME_ERROR("Malformed DEPOT_SECTION entry: " + lines[i]);
            }
            if (depot_id == -1) {
                break;
            }
            depot_nodes.push_back(depot_id);
        }

        if (depot_nodes.empty()) {
            THROW_RUNTIME_ERROR("DEPOT_SECTION does not contain a depot node");
        }
        if (depot_nodes.size() != 1) {
            THROW_RUNTIME_ERROR("Only single-depot CVRP instances are supported");
        }
        if (depot_nodes.front() != 1) {
            THROW_RUNTIME_ERROR("Unsupported depot id " + std::to_string(depot_nodes.front()) +
                                "; RouteOpt expects the depot to be node 1");
        }
    }
}
