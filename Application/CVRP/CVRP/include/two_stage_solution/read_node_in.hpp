//
// Created by Ricky You on 11/21/24.
//

#ifndef INCLUDE_TWO_STAGE_SOLUTION_READ_NODE_IN_HPP
#define INCLUDE_TWO_STAGE_SOLUTION_READ_NODE_IN_HPP

#include "write_node_out.hpp"

class ReadNodeIn {
public:
    static CVRP *cvrp;
    static std::string original_file_name;
    // static std::string node_name;
    // static BbNode *node;
    // static int node_counter;
    // static std::unordered_map<yzzLong, std::tuple<std::vector<int>, double, size_t> > enumeration_col_idx;

    static void init(CVRP *cvrp);

    // static void updateNode(BbNode *node);

    static void recoverNodeInfo();

    static void rmNodeFile();

    static void tryUpdateUB();
};


#ifdef READ_NODE_IN
#define read_node_in_call(...) __VA_ARGS__;
#else
#define read_node_in_call(...)
#endif


#endif //INCLUDE_TWO_STAGE_SOLUTION_READ_NODE_IN_HPP
