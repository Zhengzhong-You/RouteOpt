//
// Created by Ricky You on 11/21/24.
//


#include "read_node_in.hpp"
#include "branching.hpp"
#include "robust_control.hpp"
#include "best_bound_first_branching.hpp"
#include "depth_first_branching.hpp"
#include <fstream>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/file.h>

CVRP *ReadNodeIn::cvrp{nullptr};
std::string ReadNodeIn::original_file_name{};
std::string file_name{};
// BbNode *ReadNodeIn::node{nullptr};
std::ifstream inFile{};
BbNode *node{nullptr};

void initInFileStream();

void closeInFileStream();

void readNodeMPS();

void readBucketGraphIn();

void readColSequenceIn();

void readCutIn();

void readBranchInfo();

void readRobustControl();

void readEnumerationTryGap();

void printInfoNode();

void deleteFile(const std::string &fileName);

void ReadNodeIn::init(CVRP *cvrp) {
    ReadNodeIn::cvrp = cvrp;
    original_file_name = cvrp->file_name;
    // X-n336-k84_3_G16; remove anything after G, and record the number after G
    auto tmp_path = Config::tree_path;
    auto pos = tmp_path.find("_G");
    if (pos != std::string::npos) {
        cvrp->file_name = tmp_path.substr(0, pos);
    }
    file_name = Config::tree_path;
}

void ReadNodeIn::recoverNodeInfo() {
    node = new BbNode(500, cvrp);
    // node = node;
    initInFileStream();
    readNodeMPS();
    readBucketGraphIn();
    readColSequenceIn();
    readCutIn();
    readBranchInfo();
    readRobustControl();
    readEnumerationTryGap();

    closeInFileStream();

    BaseBranching::ub = Config::ub;
    cvrp->old_ub = BaseBranching::ub;

#if SOLUTION_TYPE == 1
    BestBoundFirstBranching::bbt.push(node);
#elif SOLUTION_TYPE == 2
    DepthFirstBranching::addNodeIn(DepthFirstBranching::bbt, node);
#endif
    BaseBranching::lb = node->getCurrentNodeVal();
    BaseBranching::lb_transformed = cvrp->ceilTransformedNumberRelated(BaseBranching::lb);

    // printInfoNode();
}


void ReadNodeIn::rmNodeFile() {
    auto fileName = std::string(NODE_FOLDER) + "/" + file_name +
                    NODE_FILE_SUFFIX;

    deleteFile(fileName);
    fileName = std::string(NODE_FOLDER) + "/" + file_name + ".mps";
    deleteFile(fileName);
}

void ReadNodeIn::tryUpdateUB() {
    // Remove all '_' from the file_name to create the shared file name
    std::string fileName = original_file_name;
    fileName += ".ub";

    // Try to open the file in read mode
    std::ifstream inFile(fileName);
    double recordedUB = std::numeric_limits<double>::max(); // Assume UB is a double
    auto &ub = BaseBranching::ub;

    if (inFile.is_open() && inFile >> recordedUB) {
        inFile.close();

        if (recordedUB < ub - TOLERANCE) {
            ub = recordedUB;
            std::cout << "Updated UB from file: " << fileName << " to " << ub << "\n";
        } else if (recordedUB > ub + TOLERANCE) {
            // Lock and write to the file only if necessary
            int fd = open(fileName.c_str(), O_WRONLY);
            if (fd == -1) {
                throw std::runtime_error("Failed to open file for writing: " + fileName);
            }
            if (flock(fd, LOCK_EX) == -1) {
                // Exclusive lock
                close(fd);
                throw std::runtime_error("Failed to lock file: " + fileName);
            }

            std::ofstream outFile(fileName, std::ios::out);
            if (!outFile.is_open()) {
                flock(fd, LOCK_UN);
                close(fd);
                throw std::runtime_error("Failed to open file for writing: " + fileName);
            }
            outFile << ub;
            outFile.close();

            // Unlock the file
            flock(fd, LOCK_UN);
            close(fd);
            std::cout << "Updated UB in file: " << fileName << " to " << ub << "\n";
        }
    } else {
        std::cout << "File does not exist or failed to read UB. Creating a new file with current UB.\n";
        std::ofstream outFile(fileName, std::ios::out);
        if (!outFile.is_open()) {
            throw std::runtime_error("Failed to create or open file for writing: " + fileName);
        }
        outFile << ub;
        outFile.close();
    }
}


void readNodeMPS() {
    // node
    const auto &cvrp = ReadNodeIn::cvrp;
    const auto mps_name = std::string(NODE_FOLDER) + "/" + file_name +
                          ".mps";
    cvrp->setSolverEnv();
    node->getSolver().getEnv(&cvrp->solver);
    cvrp->rollback_solver.getEnv(&cvrp->solver);
    safe_solver(node->solver.readModel(mps_name.c_str()))
}

void readBucketGraphIn() {
    if (!inFile.is_open())
        throw std::runtime_error("File stream is not open.");
    const auto &cvrp = ReadNodeIn::cvrp;
    const auto dim = cvrp->dim;
    const auto real_dim = cvrp->real_dim;
    auto &num_buckets_per_vertex = cvrp->num_buckets_per_vertex;
    auto &step_size = cvrp->step_size;
    auto &label_array_in_forward_sense = cvrp->label_array_in_forward_sense;
    auto &if_exist_extra_labels_in_forward_sense = cvrp->if_exist_extra_labels_in_forward_sense;
    auto &rc2_till_this_bin_in_forward_sense = cvrp->rc2_till_this_bin_in_forward_sense;
    auto &rc2_bin_in_forward_sense = cvrp->rc2_bin_in_forward_sense;
    auto &all_forward_buckets = node->all_forward_buckets;
    auto &num_forward_bucket_arcs = node->num_forward_bucket_arcs;
    num_forward_bucket_arcs = 0;
#ifdef SYMMETRY_PROHIBIT
    auto &label_array_in_backward_sense = cvrp->label_array_in_backward_sense;
    auto &if_exist_extra_labels_in_backward_sense = cvrp->if_exist_extra_labels_in_backward_sense;
    auto &rc2_till_this_bin_in_backward_sense = cvrp->rc2_till_this_bin_in_backward_sense;
    auto &rc2_bin_in_backward_sense = cvrp->rc2_bin_in_backward_sense;
    auto &all_backward_buckets = node->all_backward_buckets;
    auto &num_backward_bucket_arcs = node->num_backward_bucket_arcs;
    num_backward_bucket_arcs=0;
#endif
    inFile.read(reinterpret_cast<char *>(&node->value), sizeof(node->value));
    inFile.read(reinterpret_cast<char *>(&node->index), sizeof(node->index));
    inFile.read(reinterpret_cast<char *>(&step_size), sizeof(step_size));
    inFile.read(reinterpret_cast<char *>(&num_buckets_per_vertex), sizeof(num_buckets_per_vertex));

    label_array_in_forward_sense = new ListLabel *[dim];
    if_exist_extra_labels_in_forward_sense = new VecLabel *[dim];
    rc2_till_this_bin_in_forward_sense = new double *[dim];
    rc2_bin_in_forward_sense = new double *[dim];
    all_forward_buckets = new Bucket *[dim];

    for (int i = 0; i < dim; ++i) {
        label_array_in_forward_sense[i] = new ListLabel[num_buckets_per_vertex];
        if_exist_extra_labels_in_forward_sense[i] = new VecLabel[num_buckets_per_vertex];
        rc2_till_this_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
        rc2_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
        all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
    }
    cvrp->max_num_forward_graph_arc = num_buckets_per_vertex * (real_dim - 1) * real_dim;

#ifdef SYMMETRY_PROHIBIT
    label_array_in_backward_sense = new ListLabel *[dim];
    if_exist_extra_labels_in_backward_sense = new VecLabel *[dim];
    rc2_till_this_bin_in_backward_sense = new double *[dim];
    rc2_bin_in_backward_sense = new double *[dim];
    all_backward_buckets = new Bucket *[dim];

    for (int i = 0; i < dim; ++i) {
        label_array_in_backward_sense[i] = new ListLabel[num_buckets_per_vertex];
        if_exist_extra_labels_in_backward_sense[i] = new VecLabel[num_buckets_per_vertex];
        rc2_till_this_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
        rc2_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
        all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
    }
    cvrp->max_num_backward_graph_arc = num_buckets_per_vertex * (real_dim - 1) * real_dim;
#endif

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < num_buckets_per_vertex; ++j) {
            int size;
            inFile.read(reinterpret_cast<char *>(&size), sizeof(size));
            // label_array_in_forward_sense[i][j].resize(size);
            if_exist_extra_labels_in_forward_sense[i][j].first.resize(size);

            ///
            int arc_count;
            inFile.read(reinterpret_cast<char *>(&arc_count), sizeof(arc_count));
            all_forward_buckets[i][j].bucket_arcs.resize(arc_count);
            for (int k = 0; k < arc_count; ++k) {
                inFile.read(reinterpret_cast<char *>(&all_forward_buckets[i][j].bucket_arcs[k]),
                            sizeof(int));
            }
            num_forward_bucket_arcs += arc_count;

            ///
            int jump_count;
            inFile.read(reinterpret_cast<char *>(&jump_count), sizeof(jump_count));
            all_forward_buckets[i][j].jump_arcs.resize(jump_count);
            for (int k = 0; k < jump_count; ++k) {
                double first;
                int second;
                inFile.read(reinterpret_cast<char *>(&first), sizeof(first));
                inFile.read(reinterpret_cast<char *>(&second), sizeof(second));
                all_forward_buckets[i][j].jump_arcs[k] = {first, second};
            }
#ifdef SYMMETRY_PROHIBIT
            int size_backward;
            inFile.read(reinterpret_cast<char *>(&size_backward), sizeof(size_backward));
            // label_array_in_backward_sense[i][j].resize(size_backward);
            if_exist_extra_labels_in_backward_sense[i][j].first.resize(size_backward);

            int arc_count_backward;
            inFile.read(reinterpret_cast<char *>(&arc_count_backward), sizeof(arc_count_backward));
            all_backward_buckets[i][j].bucket_arcs.resize(arc_count_backward);
            for (int k = 0; k < arc_count_backward; ++k) {
                inFile.read(reinterpret_cast<char *>(&all_backward_buckets[i][j].bucket_arcs[k]), sizeof(int));
            }
            num_backward_bucket_arcs += arc_count_backward;

            int jump_count_backward;
            inFile.read(reinterpret_cast<char *>(&jump_count_backward), sizeof(jump_count_backward));
            all_backward_buckets[i][j].jump_arcs.resize(jump_count_backward);
            for (int k = 0; k < jump_count_backward; ++k) {
                double first;
                int second;
                inFile.read(reinterpret_cast<char *>(&first), sizeof(first));
                inFile.read(reinterpret_cast<char *>(&second), sizeof(second));
                all_backward_buckets[i][j].jump_arcs[k] = {first, second};
            }
#endif
        }
    }
    cvrp->getTopologicalOrder(node);
}

void readColSequenceIn() {
    if (!inFile.is_open())
        throw std::runtime_error("File stream is not open.");
    int num_cols;
    inFile.read(reinterpret_cast<char *>(&num_cols), sizeof(num_cols));
    auto &cols = node->getCols();
    cols.resize(num_cols);
    for (auto &[col_seq, main_res, forward_concatenate_pos]: cols) {
        int num_seq;
        inFile.read(reinterpret_cast<char *>(&num_seq), sizeof(num_seq));
        col_seq.resize(num_seq);
        for (auto &i: col_seq) {
            inFile.read(reinterpret_cast<char *>(&i), sizeof(i));
        }

        int num_res;
        inFile.read(reinterpret_cast<char *>(&num_res), sizeof(num_res));
        main_res.resize(num_res);
        for (auto &i: main_res) {
            inFile.read(reinterpret_cast<char *>(&i), sizeof(i));
        }
        inFile.read(reinterpret_cast<char *>(&forward_concatenate_pos), sizeof(forward_concatenate_pos));
    }
}

void readCutIn() {
    if (!inFile.is_open())
        throw std::runtime_error("File stream is not open.");
    // Read rccs data
    int num_rccs;
    inFile.read(reinterpret_cast<char *>(&num_rccs), sizeof(num_rccs));
    node->rccs.resize(num_rccs);
    for (auto &rcc: node->rccs) {
        inFile.read(reinterpret_cast<char *>(&rcc.idx_rcc), sizeof(rcc.idx_rcc));
        inFile.read(reinterpret_cast<char *>(&rcc.rhs), sizeof(rcc.rhs));
        inFile.read(reinterpret_cast<char *>(&rcc.form_rcc), sizeof(rcc.form_rcc));
#if SOLVER_VRPTW == 1
        inFile.read(reinterpret_cast<char*>(&rcc.if_keep), sizeof(rcc.if_keep));
#endif
        int customer_count;
        inFile.read(reinterpret_cast<char *>(&customer_count), sizeof(customer_count));
        rcc.info_rcc_customer.resize(customer_count);
        for (int &i: rcc.info_rcc_customer) {
            inFile.read(reinterpret_cast<char *>(&i), sizeof(i));
        }
        int outside_count;
        inFile.read(reinterpret_cast<char *>(&outside_count), sizeof(outside_count));
        rcc.info_rcc_outside_customer.resize(outside_count);

        for (int &i: rcc.info_rcc_outside_customer) {
            inFile.read(reinterpret_cast<char *>(&i), sizeof(i));
        }
    }

    // Read r1cs data
    int num_r1cs;
    inFile.read(reinterpret_cast<char *>(&num_r1cs), sizeof(num_r1cs));
    node->r1cs.resize(num_r1cs);
    for (auto &r1c: node->r1cs) {
        inFile.read(reinterpret_cast<char *>(&r1c.idx_r1c), sizeof(r1c.idx_r1c));
        inFile.read(reinterpret_cast<char *>(&r1c.rhs), sizeof(r1c.rhs));

        int info_count;
        inFile.read(reinterpret_cast<char *>(&info_count), sizeof(info_count));
        r1c.info_r1c.first.resize(info_count);
        for (int &i: r1c.info_r1c.first) {
            inFile.read(reinterpret_cast<char *>(&i), sizeof(i));
        }
        inFile.read(reinterpret_cast<char *>(&r1c.info_r1c.second), sizeof(r1c.info_r1c.second));

        int arc_mem_count;
        inFile.read(reinterpret_cast<char *>(&arc_mem_count), sizeof(arc_mem_count));
        r1c.arc_mem.resize(arc_mem_count);
        for (auto &[arc, mem]: r1c.arc_mem) {
            int arc_size;
            inFile.read(reinterpret_cast<char *>(&arc_size), sizeof(arc_size));
            arc.resize(arc_size);
            for (int &i: arc) {
                inFile.read(reinterpret_cast<char *>(&i), sizeof(i));
            }
            inFile.read(reinterpret_cast<char *>(&mem), sizeof(mem));
        }
    }
    // Read brcs data
    int num_brcs;
    inFile.read(reinterpret_cast<char *>(&num_brcs), sizeof(num_brcs));
    node->brcs.resize(num_brcs);
    for (auto &brc: node->brcs) {
        inFile.read(reinterpret_cast<char *>(&brc.idx_brc), sizeof(brc.idx_brc));
        inFile.read(reinterpret_cast<char *>(&brc.edge.first), sizeof(brc.edge.first));
        inFile.read(reinterpret_cast<char *>(&brc.edge.second), sizeof(brc.edge.second));
        inFile.read(reinterpret_cast<char *>(&brc.br_dir), sizeof(brc.br_dir));
    }
    node->tree_level = static_cast<int>(node->brcs.size());
}

template<typename T1, typename T2>
void readPair(std::pair<T1, T2> &p) {
    inFile.read(reinterpret_cast<char *>(&p.first), sizeof(p.first));
    inFile.read(reinterpret_cast<char *>(&p.second), sizeof(p.second));
}

template<typename Map, bool if_map = true>
void readMap(Map &mapData) {
    int mapSize;
    inFile.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    if constexpr (!if_map) mapData.resize(mapSize);
    for (int i = 0; i < mapSize; ++i) {
        if constexpr (!if_map) {
            auto &key = mapData[i].first;
            auto &value = mapData[i].second;
            readPair(key);
            readPair(value);
        } else {
            typename Map::key_type key;
            typename Map::mapped_type value;
            readPair(key);
            readPair(value);
            mapData.insert(std::make_pair(key, value));
        }
    }
}

void readBranchInfo() {
    if (!inFile.is_open())
        throw std::runtime_error("File stream is not open.");
    readMap<decltype(Dynamics::r_star_depth), false>(Dynamics::r_star_depth);
    readMap(BaseBranching::branching_history.exact_improvement_up);
    readMap(BaseBranching::branching_history.exact_improvement_down);
    readMap(BaseBranching::branching_history.heuristic_improvement_up);
    readMap(BaseBranching::branching_history.heuristic_improvement_down);
    readMap(BaseBranching::branching_history.lp_testing_improvement_up);
    readMap(BaseBranching::branching_history.lp_testing_improvement_down);
    readMap<decltype(BaseBranching::branching_history.increase_depth), false>(
        BaseBranching::branching_history.increase_depth);
}

void readRobustControl() {
    if (!inFile.is_open())
        throw std::runtime_error("File stream is not open.");
    inFile.read(reinterpret_cast<char *>(&RobustControl::rank1_cuts_mode),
                sizeof(RobustControl::rank1_cuts_mode));
    inFile.read(reinterpret_cast<char *>(&RobustControl::if_fix_resource_point),
                sizeof(RobustControl::if_fix_resource_point));
    inFile.read(reinterpret_cast<char *>(&RobustControl::pricing_hard_level),
                sizeof(RobustControl::pricing_hard_level));
    inFile.read(reinterpret_cast<char *>(&RobustControl::if_ever_roll_back),
                sizeof(RobustControl::if_ever_roll_back));
}

void readEnumerationTryGap() {
    if (!inFile.is_open())
        throw std::runtime_error("File stream is not open.");
    const auto &local_cvrp = ReadNodeIn::cvrp;
    inFile.read(reinterpret_cast<char *>(&local_cvrp->max_enumeration_success_gap),
                sizeof(local_cvrp->max_enumeration_success_gap));
    inFile.read(reinterpret_cast<char *>(&local_cvrp->success_enumeration_gap.first),
                sizeof(local_cvrp->success_enumeration_gap.first));
    inFile.read(reinterpret_cast<char *>(&local_cvrp->success_enumeration_gap.second),
                sizeof(local_cvrp->success_enumeration_gap.second));
    inFile.read(reinterpret_cast<char *>(&local_cvrp->min_enumeration_fail_gap),
                sizeof(local_cvrp->min_enumeration_fail_gap));
    inFile.read(reinterpret_cast<char *>(&Config::MaxGap2TryEnumeration),
                sizeof(Config::MaxGap2TryEnumeration));
}

void initInFileStream() {
    const auto filename = std::string(NODE_FOLDER) + "/" + file_name +
                          NODE_FILE_SUFFIX;
    inFile.open(filename, std::ios::binary);
    if (!inFile.is_open())
        throw std::runtime_error("Unable to open file.");
}

void deleteFile(const std::string &fileName) {
    if (std::remove(fileName.c_str()) == 0) {
        std::cout << "File deleted successfully: " << fileName << std::endl;
    }
}


void closeInFileStream() {
    if (inFile.is_open()) {
        inFile.close();
    }
}

void printInfoNode() {
    using namespace std;

    // Basic node information
    cout << "Node value: " << node->value << endl;
    cout << "Node index: " << node->index << endl;
    cout << "Step size: " << ReadNodeIn::cvrp->step_size << endl;
    cout << "Number of buckets per vertex: " << ReadNodeIn::cvrp->num_buckets_per_vertex << endl;

    // Forward graph data
    cout << "Total forward bucket arcs: " << node->num_forward_bucket_arcs << endl;
    for (int i = 0; i < ReadNodeIn::cvrp->dim; ++i) {
        for (int j = 0; j < ReadNodeIn::cvrp->num_buckets_per_vertex; ++j) {
            cout << "if_exist_extra_labels_in_forward_sense: " << ReadNodeIn::cvrp->
                    if_exist_extra_labels_in_forward_sense[i][j].second << endl;
            cout << "Forward Bucket " << i << ", " << j << ":" << endl;

            cout << "  Arcs: ";
            for (int arc: node->all_forward_buckets[i][j].bucket_arcs) {
                cout << arc << " ";
            }
            cout << endl;
            cout << "  Jumps: ";
            for (const auto &jump: node->all_forward_buckets[i][j].jump_arcs) {
                cout << "(" << jump.first << ", " << jump.second << ") ";
            }
            cout << endl;
        }
    }

    // If symmetry is prohibited, print backward graph data
#ifdef SYMMETRY_PROHIBIT
    cout << "Total backward bucket arcs: " << node->num_backward_bucket_arcs << endl;
    for (int i = 0; i < ReadNodeIn::cvrp->dim; ++i) {
        for (int j = 0; j < ReadNodeIn::cvrp->num_buckets_per_vertex; ++j) {
            cout << "Backward Bucket " << i << ", " << j << ":" << endl;
            cout << "  Arcs: ";
            for (int arc : node->all_backward_buckets[i][j].bucket_arcs) {
                cout << arc << " ";
            }
            cout << endl;
            cout << "  Jumps: ";
            for (const auto& jump : node->all_backward_buckets[i][j].jump_arcs) {
                cout << "(" << jump.first << ", " << jump.second << ") ";
            }
            cout << endl;
        }
    }
#endif

    // Print Column Sequences
    auto &cols = node->getCols();
    cout << "Column Sequences:" << endl;
    for (const auto &col: cols) {
        cout << "  Sequence: ";
        for (int idx: col.col_seq) {
            cout << idx << " ";
        }
        cout << endl;
        cout << "  Main Resources: ";
        for (auto res: col.main_res) {
            cout << res << " ";
        }
        cout << endl;
        cout << "  Concat Position: " << col.forward_concatenate_pos << endl;
    }

    // Print Cut Information
    cout << "RCCs:" << endl;
    for (const auto &rcc: node->rccs) {
        cout << "  RCC index: " << rcc.idx_rcc << ", RHS: " << rcc.rhs << ", Form: " << rcc.form_rcc << endl;
    }

    cout << "R1Cs:" << endl;
    for (const auto &r1c: node->r1cs) {
        cout << "  R1C index: " << r1c.idx_r1c << ", RHS: " << r1c.rhs << endl;
        cout << "  Info: ";
        for (int info: r1c.info_r1c.first) {
            cout << info << " ";
        }
        cout << endl;
        cout << "  Arc Memory: ";
        for (const auto &arc_mem: r1c.arc_mem) {
            cout << "Arcs: ";
            for (int arc: arc_mem.first) {
                cout << arc << " ";
            }
            cout << " to: " << arc_mem.second;
        }
        cout << endl;
    }

    // Print BRCS Information if any
    cout << "BRCS:" << endl;
    for (const auto &brc: node->brcs) {
        cout << "  BRC index: " << brc.idx_brc << ", Edge: (" << brc.edge.first << ", " << brc.edge.second <<
                "), Direction: " << brc.br_dir << endl;
    }


    // Print branching information
    cout << "Branching Information:" << endl;
    cout << "  R Star Depth:" << endl;
    for (const auto &r_star: Dynamics::r_star_depth) {
        cout << "    " << r_star.first.first << ": " << r_star.first.second << ", " << r_star.second.first << ": " <<
                r_star.second.second << endl;
    }

    exit(0);
}
