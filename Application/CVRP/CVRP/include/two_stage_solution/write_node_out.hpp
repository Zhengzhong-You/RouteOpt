//
// Created by Ricky You on 11/20/24.
//

#ifndef INCLUDE_TWO_STAGE_SOLUTION_WRITE_NODE_OUT_HPP
#define INCLUDE_TWO_STAGE_SOLUTION_WRITE_NODE_OUT_HPP


#include "cvrp.hpp"
#include "branching_node.hpp"

#define NODE_FOLDER "nodes"
#define NODE_FILE_SUFFIX ".node"
#define NODE_FOLDER_NUM_FILE_LIMIT 360
#define MEMORY_USAGE_LIMIT 12// 16
#define SMALL_MEMORY_USE 32
#define LARGE_MEMORY_USE 64

class WriteNodeOut {
public:
    static CVRP *cvrp;
    static BbNode *node;
    static int node_counter;
    static int recommended_memory_usage;
    // static std::unordered_map<yzzLong, std::tuple<std::vector<int>, double, size_t> > enumeration_col_idx;

    static void init(CVRP *cvrp);

    static void updateNode(BbNode *node);

    ///write out node and set nullptr, this is a different node than the one in the class
    static void writeNodeOut(BbNode *&node, int recommended_memory_usage = SMALL_MEMORY_USE);

    // static void writeEnuCols();

    // static void writeCutsOut(const std::vector<std::pair<size_t, double> > &col_idx_cost);
    // static void writeCutsOut(const std::vector<std::pair<size_t, double> > &col_idx_cost);


    // static void populateEnuCols(std::vector<std::pair<size_t, double> > &col_idx_cost,
    //                             const std::vector<SequenceInfo> &lp_col_info,
    //                             const std::vector<std::pair<size_t, double> > &enu_col_idx_cost);
};


#ifdef WRITE_NODE_OUT
#define write_node_out_call(...) __VA_ARGS__;
#else
#define write_node_out_call(...)
#endif


#endif //INCLUDE_TWO_STAGE_SOLUTION_WRITE_NODE_OUT_HPP
