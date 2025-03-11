/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_RANK1_DATA_SHARED_HPP
#define ROUTE_OPT_RANK1_DATA_SHARED_HPP
#include <vector>
#include <unordered_set>
#include <unordered_map>


#include "rank1_separation_macro.hpp"
#include "rank1_macro.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    class DataShared {
    public:
        explicit DataShared(const std::vector<std::vector<double> > &cost_mat4_vertex) : cost_mat4_vertex_ref(
            std::ref(cost_mat4_vertex)) {
        }


        //clean data
        void cleanData() {
            limited_memory_type = MemoryType::NODE_MEMORY;
            sol.clear();
            old_cuts.clear();
            cuts.clear();
            cut_record.clear();
            v_r_map.clear();
        }


        // Setters
        void setLimitedMemoryType(MemoryType value) {
            limited_memory_type = value;
        }

        void setSol(const std::vector<RouteInfo> &value) {
            sol = value;
        }

        void setOldCuts(const std::vector<R1c> &value) { old_cuts = value; }
        void setVRMap(const std::vector<std::unordered_map<int, int> > &value) { v_r_map = value; }

        // Getters
        const std::vector<std::vector<double> > &getCostMat4Vertex() const { return cost_mat4_vertex_ref.get(); }


        MemoryType getLimitedMemoryType() const { return limited_memory_type; }
        const std::vector<RouteInfo> &getSol() const { return sol; }
        const std::vector<R1c> &getOldCuts() const { return old_cuts; }


        auto &refCuts() { return cuts; }
        auto &refCutRecord() { return cut_record; }

        const std::vector<std::unordered_map<int, int> > &getVRMap() const { return v_r_map; }

        DataShared() = delete;

        ~DataShared() = default;

    private:
        // need to keep
        const std::reference_wrapper<const std::vector<std::vector<double> >> cost_mat4_vertex_ref;

        //need to clean
        MemoryType limited_memory_type{MemoryType::NODE_MEMORY};
        std::vector<RouteInfo> sol{}; //
        std::vector<R1c> old_cuts{}; // only need for limited memory mode
        std::vector<R1c> cuts{}; // current generated cuts
        std::vector<std::unordered_map<int, int> > v_r_map{};
        std::unordered_map<std::vector<int>, std::unordered_set<int>, VectorHashInRank1> cut_record{};
    };
}

#endif // ROUTE_OPT_RANK1_DATA_SHARED_HPP
