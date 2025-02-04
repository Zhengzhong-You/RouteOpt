//
// Created by Ricky You on 11/20/24.
//

#include "write_node_out.hpp"
#include "robust_control.hpp"
#include <branching.hpp>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <sys/resource.h>

CVRP *WriteNodeOut::cvrp{nullptr};
BbNode *WriteNodeOut::node{nullptr};
int WriteNodeOut::node_counter{0};
int WriteNodeOut::recommended_memory_usage{SMALL_MEMORY_USE};

//here open the file stream and write the information out;

//write ng out
//write bucket graph information
//write cut information
//all in this small function!

std::ofstream outFile;
CVRP *local_cvrp;
BbNode *local_node;
int local_node_counter;

// void writeNgOut();

void initFileStream();

void closeFileStream();

void writeLPOut();

bool checkMemory(size_t memory_usage_limit);

void writeBucketGraphOut();

void writeCutOut();

void writeColSequence();

void writeBranchInfo();

void writeRobustControl();

void writeEnumerationTryGap();

void WriteNodeOut::init(CVRP *cvrp) {
    WriteNodeOut::cvrp = cvrp;
}


void WriteNodeOut::updateNode(BbNode *node) {
    WriteNodeOut::node = node;
}

bool checkMemory(const size_t memory_usage_limit) {
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    if (const auto memoryUsageKB = usage.ru_maxrss; memoryUsageKB > (memory_usage_limit * 1048576)) {
        std::cout << "Memory usage is " << memoryUsageKB / 1048576 << "GB, exceeds the limit of " <<
                memory_usage_limit
                << "GB." << std::endl;
        return true;
    }
    return false;
}


void WriteNodeOut::writeNodeOut(BbNode *&node, int recommended_memory_usage) {
    //check if the number of files is over the NODE_FOLDER_NUM_FILE_LIMIT
    self_mkdir(NODE_FOLDER);
    // int count;
    // if (checkMemory(MEMORY_USAGE_LIMIT)) goto WRITE_OUT;
    // count = 0;
    // // 1. get the number of files that ends with .node, in the NODE_FOLDER
    // for (const auto &entry: std::filesystem::directory_iterator(NODE_FOLDER)) {
    //     if (entry.is_regular_file() && entry.path().extension() == NODE_FILE_SUFFIX) {
    //         ++count;
    //     }
    // }
    // if (count >= NODE_FOLDER_NUM_FILE_LIMIT) return;
    // WRITE_OUT:
    ++node_counter;
    local_cvrp = cvrp;
    local_node = node;
    local_node_counter = node_counter;
    WriteNodeOut::recommended_memory_usage = recommended_memory_usage;
    initFileStream();

    writeLPOut();
    /// the following sequence cannot be changed!
    writeBucketGraphOut();
    writeColSequence();
    writeCutOut();
    writeBranchInfo();
    writeRobustControl();
    writeEnumerationTryGap();

    closeFileStream();
    delete node;
    node = nullptr;
}

void writeLPOut() {
    const auto file_name = std::string(NODE_FOLDER) + "/" + local_cvrp->file_name + "_" +
                           std::to_string(local_node_counter) + "_G" + std::to_string(
                               WriteNodeOut::recommended_memory_usage) + ".mps";
    safe_solver(local_node->getSolver().write(file_name.c_str()))
}


void writeBucketGraphOut() {
    if (!outFile.is_open())
        throw std::runtime_error("File stream is not open.");

    const auto dim = local_cvrp->dim;
    const auto num_buckets_per_vertex = local_cvrp->num_buckets_per_vertex;
    const auto step_size = local_cvrp->step_size;
    const auto &if_exist_extra_labels_in_forward_sense = local_cvrp->if_exist_extra_labels_in_forward_sense;
    const auto &all_forward_buckets = local_node->all_forward_buckets;

#ifdef SYMMETRY_PROHIBIT
    const auto &if_exist_extra_labels_in_backward_sense = local_cvrp->if_exist_extra_labels_in_backward_sense;
    const auto &all_backward_buckets = local_node->all_backward_buckets;
#endif
    outFile.write(reinterpret_cast<const char *>(&local_node->value), sizeof(local_node->value));
    outFile.write(reinterpret_cast<const char *>(&local_node->index), sizeof(local_node->index));
    outFile.write(reinterpret_cast<const char *>(&step_size), sizeof(step_size));
    outFile.write(reinterpret_cast<const char *>(&num_buckets_per_vertex), sizeof(num_buckets_per_vertex));
    // writing num_buckets_per_vertex
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < num_buckets_per_vertex; ++j) {
            auto size = static_cast<int>(if_exist_extra_labels_in_forward_sense[i][j].first.size());
            outFile.write(reinterpret_cast<const char *>(&size), sizeof(size));
            // writing size of label_array_in_forward_sense[i][j]

            const auto &[bucket_arcs, jump_arcs, s] = all_forward_buckets[i][j];
            auto arc_count = static_cast<int>(bucket_arcs.size());
            outFile.write(reinterpret_cast<const char *>(&arc_count), sizeof(arc_count));
            for (int arc: bucket_arcs) {
                outFile.write(reinterpret_cast<const char *>(&arc), sizeof(arc)); // writing arc
            }

            auto jump_count = static_cast<int>(jump_arcs.size());
            outFile.write(reinterpret_cast<const char *>(&jump_count), sizeof(jump_count)); // writing jump_count
            for (const auto &[fst, snd]: jump_arcs) {
                outFile.write(reinterpret_cast<const char *>(&fst), sizeof(fst)); // writing t node
                outFile.write(reinterpret_cast<const char *>(&snd), sizeof(snd)); // writing q resource
            }
#ifdef SYMMETRY_PROHIBIT
            auto size_backward = static_cast<int>(if_exist_extra_labels_in_backward_sense[i][j].first.size());
            outFile.write(reinterpret_cast<const char *>(&size_backward), sizeof(size_backward));
            // writing size of label_array_in_backward_sense[i][j]

            const auto &[bucket_arcs_backward, jump_arcs_backward, s_backward] = all_backward_buckets[i][j];
            auto arc_count_backward = static_cast<int>(bucket_arcs_backward.size());
            outFile.write(reinterpret_cast<const char *>(&arc_count_backward), sizeof(arc_count_backward));
            for (int arc: bucket_arcs_backward) {
                outFile.write(reinterpret_cast<const char *>(&arc), sizeof(arc)); // writing arc
            }

            auto jump_count_backward = static_cast<int>(jump_arcs_backward.size());
            outFile.write(reinterpret_cast<const char *>(&jump_count_backward), sizeof(jump_count_backward)); // writing jump_count
            for (const auto &[fst, snd]: jump_arcs_backward) {
                outFile.write(reinterpret_cast<const char *>(&fst), sizeof(fst)); // writing t node
                outFile.write(reinterpret_cast<const char *>(&snd), sizeof(snd)); // writing q resource
            }
#endif
        }
    }
}

void writeColSequence() {
    if (!outFile.is_open())
        throw std::runtime_error("File stream is not open.");
    const auto &cols = local_node->getCols();

    // Write the number of columns first
    const auto num_cols = static_cast<int>(cols.size());
    outFile.write(reinterpret_cast<const char *>(&num_cols), sizeof(num_cols));

    for (const auto &col: cols) {
        // Write each column sequence vector
        const auto num_seq = static_cast<int>(col.col_seq.size());
        outFile.write(reinterpret_cast<const char *>(&num_seq), sizeof(num_seq));
        for (const auto &i: col.col_seq) {
            outFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }

        // Write each main resource vector
        const auto num_res = static_cast<int>(col.main_res.size());
        outFile.write(reinterpret_cast<const char *>(&num_res), sizeof(num_res));
        for (const auto &i: col.main_res) {
            outFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }

        // Write the forward concatenate position
        outFile.write(reinterpret_cast<const char *>(&col.forward_concatenate_pos),
                      sizeof(col.forward_concatenate_pos));
    }
}

void writeCutOut() {
    if (!outFile.is_open())
        throw std::runtime_error("File stream is not open.");

    const auto &rccs = local_node->rccs;
    const auto &r1cs = local_node->r1cs;
    const auto &brcs = local_node->getBrCs();

    // Write rccs data
    const auto num_rccs = static_cast<int>(rccs.size());
    outFile.write(reinterpret_cast<const char *>(&num_rccs), sizeof(num_rccs));
    for (const auto &rcc: rccs) {
        outFile.write(reinterpret_cast<const char *>(&rcc.idx_rcc), sizeof(rcc.idx_rcc));
        outFile.write(reinterpret_cast<const char *>(&rcc.rhs), sizeof(rcc.rhs));
        outFile.write(reinterpret_cast<const char *>(&rcc.form_rcc), sizeof(rcc.form_rcc));
#if  SOLVER_VRPTW ==1
        outFile.write(reinterpret_cast<const char*>(&rcc.if_keep), sizeof(rcc.if_keep));
#endif
        // Customer info
        auto customer_count = static_cast<int>(rcc.info_rcc_customer.size());
        outFile.write(reinterpret_cast<const char *>(&customer_count), sizeof(customer_count));
        for (const auto &i: rcc.info_rcc_customer) {
            outFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }
        // Outside customer info
        auto outside_count = static_cast<int>(rcc.info_rcc_outside_customer.size());
        outFile.write(reinterpret_cast<const char *>(&outside_count), sizeof(outside_count));
        for (const auto &i: rcc.info_rcc_outside_customer) {
            outFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }
    }

    // Write r1cs data
    const auto num_r1cs = static_cast<int>(r1cs.size());
    outFile.write(reinterpret_cast<const char *>(&num_r1cs), sizeof(num_r1cs));
    for (const auto &r1c: r1cs) {
        outFile.write(reinterpret_cast<const char *>(&r1c.idx_r1c), sizeof(r1c.idx_r1c));
        outFile.write(reinterpret_cast<const char *>(&r1c.rhs), sizeof(r1c.rhs));
        // Info r1c
        auto info_count = static_cast<int>(r1c.info_r1c.first.size());
        outFile.write(reinterpret_cast<const char *>(&info_count), sizeof(info_count));
        for (const auto &i: r1c.info_r1c.first) {
            outFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }
        outFile.write(reinterpret_cast<const char *>(&r1c.info_r1c.second), sizeof(r1c.info_r1c.second));
        // Arc memory
        auto arc_mem_count = static_cast<int>(r1c.arc_mem.size());
        outFile.write(reinterpret_cast<const char *>(&arc_mem_count), sizeof(arc_mem_count));
        for (const auto &[arc, mem]: r1c.arc_mem) {
            auto arc_size = static_cast<int>(arc.size());
            outFile.write(reinterpret_cast<const char *>(&arc_size), sizeof(arc_size));
            for (const auto &i: arc) {
                outFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
            }
            outFile.write(reinterpret_cast<const char *>(&mem), sizeof(mem));
        }
    }

    // Write br_cs data
    const auto num_brcs = static_cast<int>(brcs.size());
    outFile.write(reinterpret_cast<const char *>(&num_brcs), sizeof(num_brcs));
    for (const auto &brc: brcs) {
        outFile.write(reinterpret_cast<const char *>(&brc.idx_brc), sizeof(brc.idx_brc));
        outFile.write(reinterpret_cast<const char *>(&brc.edge.first), sizeof(brc.edge.first));
        outFile.write(reinterpret_cast<const char *>(&brc.edge.second), sizeof(brc.edge.second));
        outFile.write(reinterpret_cast<const char *>(&brc.br_dir), sizeof(brc.br_dir));
    }
}


template<typename Map>
void writeMap(const Map &mapData) {
    const auto mapSize = static_cast<int>(mapData.size());
    outFile.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
    for (const auto &pair: mapData) {
        outFile.write(reinterpret_cast<const char *>(&pair.first.first), sizeof(pair.first.first));
        outFile.write(reinterpret_cast<const char *>(&pair.first.second), sizeof(pair.first.second));
        outFile.write(reinterpret_cast<const char *>(&pair.second.first), sizeof(pair.second.first));
        outFile.write(reinterpret_cast<const char *>(&pair.second.second), sizeof(pair.second.second));
    }
}

void writeBranchInfo() {
    // Write r_star_depth data using the vector helper
    writeMap(Dynamics::r_star_depth);
    //print r_star_depth
    // for (const auto &r_star: Dynamics::r_star_depth) {
    //     std::cout << "r_star1: " << r_star.first.first << " " << r_star.first.second << std::endl;
    //     std::cout << "r_star2: " << r_star.second.first << " " << r_star.second.second << std::endl;
    // }
    // Write exact_improvement_up using the map helper
    writeMap(BaseBranching::branching_history.exact_improvement_up);
    writeMap(BaseBranching::branching_history.exact_improvement_down);
    writeMap(BaseBranching::branching_history.heuristic_improvement_up);
    writeMap(BaseBranching::branching_history.heuristic_improvement_down);
    writeMap(BaseBranching::branching_history.lp_testing_improvement_up);
    writeMap(BaseBranching::branching_history.lp_testing_improvement_down);
    writeMap(BaseBranching::branching_history.increase_depth);
}

void writeRobustControl() {
    if (!outFile.is_open())
        throw std::runtime_error("File stream is not open.");
    outFile.write(reinterpret_cast<const char *>(&RobustControl::rank1_cuts_mode),
                  sizeof(RobustControl::rank1_cuts_mode));
    outFile.write(reinterpret_cast<const char *>(&RobustControl::if_fix_resource_point),
                  sizeof(RobustControl::if_fix_resource_point));
    outFile.write(reinterpret_cast<const char *>(&RobustControl::pricing_hard_level),
                  sizeof(RobustControl::pricing_hard_level));
    outFile.write(reinterpret_cast<const char *>(&RobustControl::if_ever_roll_back),
                  sizeof(RobustControl::if_ever_roll_back));
}

void writeEnumerationTryGap() {
    if (!outFile.is_open())
        throw std::runtime_error("File stream is not open.");
    outFile.write(reinterpret_cast<const char *>(&local_cvrp->max_enumeration_success_gap),
                  sizeof(local_cvrp->max_enumeration_success_gap));
    outFile.write(reinterpret_cast<const char *>(&local_cvrp->success_enumeration_gap.first),
                  sizeof(local_cvrp->success_enumeration_gap.first));
    outFile.write(reinterpret_cast<const char *>(&local_cvrp->success_enumeration_gap.second),
                  sizeof(local_cvrp->success_enumeration_gap.second));
    outFile.write(reinterpret_cast<const char *>(&local_cvrp->min_enumeration_fail_gap),
                  sizeof(local_cvrp->min_enumeration_fail_gap));
    outFile.write(reinterpret_cast<const char *>(&Config::MaxGap2TryEnumeration),
                  sizeof(Config::MaxGap2TryEnumeration));
}


void initFileStream() {
    const auto filename = std::string(NODE_FOLDER) + "/" + local_cvrp->file_name + "_" +
                          std::to_string(local_node_counter) + "_G" + std::to_string(
                              WriteNodeOut::recommended_memory_usage) + NODE_FILE_SUFFIX;
    outFile.open(filename, std::ios::binary | std::ios::out); // Open in out mode
    if (!outFile.is_open()) {
        throw std::runtime_error("Could not open file stream.");
    }
}

void closeFileStream() {
    outFile.flush();
    if (outFile.is_open()) {
        outFile.close();
    }
}
