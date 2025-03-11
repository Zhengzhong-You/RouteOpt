/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_NODE_OUT_HPP
#define ROUTE_OPT_T_NODE_OUT_HPP
#include <fstream>
#include <iostream>
#include <zlib.h>

#include "two_stage_macro.hpp"
#include "cvrp_pricing_controller.hpp"
#include "label.hpp"

namespace RouteOpt::Application::CVRP {
    namespace OutNodeNameSpace {
        inline int local_node_counter{};
        inline std::ofstream outFile{};

        inline void reportMemoryUsage() {
            rusage usage{};
            getrusage(RUSAGE_SELF, &usage);
            std::cout << "memory usage= " << usage.ru_maxrss / 1048576 << "GB." << std::endl;
        }

        inline void initFileStream(const std::string &file_name, int memory_usage_limit) {
            mkDir(NODE_FOLDER.data());
            reportMemoryUsage();

            std::string base_name;
            auto pos = file_name.find("_G");
            if (pos != std::string::npos) {
                base_name = file_name.substr(0, pos);
            } else {
                base_name = file_name;
            }
            const auto filename = std::string(NODE_FOLDER) + "/" + base_name + "_" +
                                  std::to_string(local_node_counter) + "_G" +
                                  std::to_string(memory_usage_limit) + std::string(NODE_FILE_SUFFIX);


            outFile.open(filename, std::ios::binary | std::ios::out);
            // Open in out mode, _0 means not in calculation yet
            if (!outFile.is_open()) {
                THROW_RUNTIME_ERROR("cannot open file stream.");
            }
            ++local_node_counter;
        }

        inline void closeFileStream() {
            outFile.flush();
            if (outFile.is_open()) {
                outFile.close();
            }
        }

        inline void writeCompressedData(std::ostream &out, const std::string &rawData) {
            // Get the size of the uncompressed data.
            size_t srcLen = rawData.size();

            // Determine the maximum compressed size.
            size_t destLen = compressBound(srcLen);
            std::vector<Bytef> compressedData(destLen);

            // Compress the data using maximum compression.
            int ret = compress2(compressedData.data(), &destLen,
                                reinterpret_cast<const Bytef *>(rawData.data()),
                                srcLen, Z_BEST_COMPRESSION);
            if (ret != Z_OK) {
                THROW_RUNTIME_ERROR("compression failed with error code: " + std::to_string(ret));
            }

            // Write header information (uncompressed size and compressed size)
            out.write(reinterpret_cast<const char *>(&srcLen), sizeof(srcLen));
            out.write(reinterpret_cast<const char *>(&destLen), sizeof(destLen));

            // Write the compressed data.
            out.write(reinterpret_cast<const char *>(compressedData.data()), static_cast<std::streamsize>(destLen));
        }

        inline void writeColSequence(const Solver &node_solver, const std::vector<SequenceInfo> &cols) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            // 1. Serialize the column sequence data into an in-memory buffer.
            std::ostringstream oss(std::ios::binary);

            // Write the number of columns first
            const int num_cols = static_cast<int>(cols.size());
            oss.write(reinterpret_cast<const char *>(&num_cols), sizeof(num_cols));

            std::vector<double> obj(num_cols);
            SAFE_SOLVER(node_solver.getObj(0, num_cols, obj.data()))

            // Write each column's data
            for (int i = 0; i < num_cols; ++i) {
                const auto &col = cols[i];
                const int num_seq = static_cast<int>(col.col_seq.size());
                oss.write(reinterpret_cast<const char *>(&num_seq), sizeof(num_seq));
                for (const auto &j: col.col_seq) {
                    oss.write(reinterpret_cast<const char *>(&j), sizeof(j));
                }
                // Write the forward concatenate position
                oss.write(reinterpret_cast<const char *>(&col.forward_concatenate_pos),
                          sizeof(col.forward_concatenate_pos));
                oss.write(reinterpret_cast<const char *>(&obj[i]), sizeof(obj[i]));
            }
            writeCompressedData(outFile, oss.str());
        }

        inline void writeCutOut(const std::vector<Rcc> &rccs,
                                const std::vector<R1c> &r1cs,
                                const std::vector<Brc> &brcs) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            // Use an ostringstream to serialize all the cut-out data.
            std::ostringstream oss(std::ios::binary);

            // Write Rccs data.
            const int num_rccs = static_cast<int>(rccs.size());
            oss.write(reinterpret_cast<const char *>(&num_rccs), sizeof(num_rccs));
            for (const auto &rcc: rccs) {
                oss.write(reinterpret_cast<const char *>(&rcc.idx_rcc), sizeof(rcc.idx_rcc));
                oss.write(reinterpret_cast<const char *>(&rcc.rhs), sizeof(rcc.rhs));
                oss.write(reinterpret_cast<const char *>(&rcc.form_rcc), sizeof(rcc.form_rcc));
                oss.write(reinterpret_cast<const char *>(&rcc.if_keep), sizeof(rcc.if_keep));

                int customer_count = static_cast<int>(rcc.info_rcc_customer.size());
                oss.write(reinterpret_cast<const char *>(&customer_count), sizeof(customer_count));
                for (const auto &i: rcc.info_rcc_customer) {
                    oss.write(reinterpret_cast<const char *>(&i), sizeof(i));
                }

                int outside_count = static_cast<int>(rcc.info_rcc_outside_customer.size());
                oss.write(reinterpret_cast<const char *>(&outside_count), sizeof(outside_count));
                for (const auto &i: rcc.info_rcc_outside_customer) {
                    oss.write(reinterpret_cast<const char *>(&i), sizeof(i));
                }
            }

            // Write R1cs data.
            const int num_r1cs = static_cast<int>(r1cs.size());
            oss.write(reinterpret_cast<const char *>(&num_r1cs), sizeof(num_r1cs));
            for (const auto &r1c: r1cs) {
                oss.write(reinterpret_cast<const char *>(&r1c.idx_r1c), sizeof(r1c.idx_r1c));
                oss.write(reinterpret_cast<const char *>(&r1c.rhs), sizeof(r1c.rhs));

                int info_count = static_cast<int>(r1c.info_r1c.first.size());
                oss.write(reinterpret_cast<const char *>(&info_count), sizeof(info_count));
                for (const auto &i: r1c.info_r1c.first) {
                    oss.write(reinterpret_cast<const char *>(&i), sizeof(i));
                }
                oss.write(reinterpret_cast<const char *>(&r1c.info_r1c.second), sizeof(r1c.info_r1c.second));

                int arc_mem_count = static_cast<int>(r1c.arc_mem.size());
                oss.write(reinterpret_cast<const char *>(&arc_mem_count), sizeof(arc_mem_count));
                for (const auto &pair: r1c.arc_mem) {
                    const auto &arc = pair.first;
                    const auto &mem = pair.second;
                    int arc_size = static_cast<int>(arc.size());
                    oss.write(reinterpret_cast<const char *>(&arc_size), sizeof(arc_size));
                    for (const auto &i: arc) {
                        oss.write(reinterpret_cast<const char *>(&i), sizeof(i));
                    }
                    oss.write(reinterpret_cast<const char *>(&mem), sizeof(mem));
                }
            }

            // Write Brcs data.
            const int num_brcs = static_cast<int>(brcs.size());
            oss.write(reinterpret_cast<const char *>(&num_brcs), sizeof(num_brcs));
            for (const auto &brc: brcs) {
                oss.write(reinterpret_cast<const char *>(&brc.idx_brc), sizeof(brc.idx_brc));
                oss.write(reinterpret_cast<const char *>(&brc.edge.first), sizeof(brc.edge.first));
                oss.write(reinterpret_cast<const char *>(&brc.edge.second), sizeof(brc.edge.second));
                oss.write(reinterpret_cast<const char *>(&brc.br_dir), sizeof(brc.br_dir));
            }
            // Compress the data and write it to outFile using the general utility.
            writeCompressedData(outFile, oss.str());
        }

        inline void writeEnuState(bool if_enu) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::ostringstream oss(std::ios::binary);
            unsigned char u = if_enu ? 1 : 0;
            oss.write(reinterpret_cast<const char *>(&u), sizeof(u));
            writeCompressedData(outFile, oss.str());
        }

        inline void writeNodeInfo(double val, int idx) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::ostringstream oss(std::ios::binary);
            oss.write(reinterpret_cast<const char *>(&val), sizeof(val));
            oss.write(reinterpret_cast<const char *>(&idx), sizeof(idx));
            writeCompressedData(outFile, oss.str());
        }

        template<bool if_symmetry>
        void writeBucketGraphOut(
            int dim, int num_buckets_per_vertex, res_int step_size,
            const VecLabel *const*if_exist_extra_labels_in_forward_sense,
            const Bucket *const*all_forward_buckets,
            const VecLabel *const*if_exist_extra_labels_in_backward_sense,
            const Bucket *const*all_backward_buckets) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::ostringstream oss(std::ios::binary);

            oss.write(reinterpret_cast<const char *>(&step_size), sizeof(step_size));
            oss.write(reinterpret_cast<const char *>(&num_buckets_per_vertex), sizeof(num_buckets_per_vertex));

            // Write forward bucket graph data.
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < num_buckets_per_vertex; ++j) {
                    // Write the size of the extra forward label array.
                    int fwdLabelSize = static_cast<int>(if_exist_extra_labels_in_forward_sense[i][j].first.size());
                    oss.write(reinterpret_cast<const char *>(&fwdLabelSize), sizeof(fwdLabelSize));

                    // Retrieve bucket information.
                    // Using structured binding (C++17) for clarity.
                    const auto &[bucket_arcs, jump_arcs, s] = all_forward_buckets[i][j];

                    int arc_count = static_cast<int>(bucket_arcs.size());
                    oss.write(reinterpret_cast<const char *>(&arc_count), sizeof(arc_count));
                    for (int arc: bucket_arcs) {
                        oss.write(reinterpret_cast<const char *>(&arc), sizeof(arc));
                    }

                    int jump_count = static_cast<int>(jump_arcs.size());
                    oss.write(reinterpret_cast<const char *>(&jump_count), sizeof(jump_count));
                    for (const auto &[fst, snd]: jump_arcs) {
                        oss.write(reinterpret_cast<const char *>(&fst), sizeof(fst));
                        oss.write(reinterpret_cast<const char *>(&snd), sizeof(snd));
                    }
                }
            }

            // Write backward bucket graph data if symmetry is not enabled.
            if constexpr (!if_symmetry) {
                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < num_buckets_per_vertex; ++j) {
                        int bwdLabelSize = static_cast<int>(if_exist_extra_labels_in_backward_sense[i][j].first.size());
                        oss.write(reinterpret_cast<const char *>(&bwdLabelSize), sizeof(bwdLabelSize));

                        const auto &[bucket_arcs_backward, jump_arcs_backward, s_backward] = all_backward_buckets[i][j];
                        int arc_count_backward = static_cast<int>(bucket_arcs_backward.size());
                        oss.write(reinterpret_cast<const char *>(&arc_count_backward), sizeof(arc_count_backward));
                        for (int arc: bucket_arcs_backward) {
                            oss.write(reinterpret_cast<const char *>(&arc), sizeof(arc));
                        }

                        int jump_count_backward = static_cast<int>(jump_arcs_backward.size());
                        oss.write(reinterpret_cast<const char *>(&jump_count_backward), sizeof(jump_count_backward));
                        for (const auto &[fst, snd]: jump_arcs_backward) {
                            oss.write(reinterpret_cast<const char *>(&fst), sizeof(fst));
                            oss.write(reinterpret_cast<const char *>(&snd), sizeof(snd));
                        }
                    }
                }
            }

            // Compress and write the data to outFile using the general function.
            writeCompressedData(outFile, oss.str());
        }

        inline void writeEnuMatrix(const std::vector<bool> &deleted_columns_in_enumeration_pool,
                                   const RowVectorXT &index_columns_in_enumeration_column_pool,
                                   const RowVectorXd &cost_for_columns_in_enumeration_column_pool,
                                   const int *col_pool4_pricing) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");
            std::ostringstream oss(std::ios::binary);


            int cnt = 0;
            for (int i = 0; i < deleted_columns_in_enumeration_pool.size(); ++i) {
                if (deleted_columns_in_enumeration_pool.at(i)) continue;
                ++cnt;
            }

            oss.write(reinterpret_cast<const char *>(&cnt), sizeof(cnt));

            for (int i = 0; i < deleted_columns_in_enumeration_pool.size(); ++i) {
                if (deleted_columns_in_enumeration_pool.at(i)) continue;
                auto idx = index_columns_in_enumeration_column_pool[i];
                int curr_node = 0;
                oss.write(reinterpret_cast<const char *>(&curr_node), sizeof(curr_node));
                for (auto j = idx + 1; ; ++j) {
                    curr_node = col_pool4_pricing[j];
                    oss.write(reinterpret_cast<const char *>(&curr_node), sizeof(curr_node));
                    if (!curr_node) break;
                }
                oss.write(reinterpret_cast<const char *>(&cost_for_columns_in_enumeration_column_pool[i]),
                          sizeof(cost_for_columns_in_enumeration_column_pool[i]));
            }

            writeCompressedData(outFile, oss.str());
        }

        inline void writeEnumerationTryGap(double max_enumeration_success_gap,
                                           const std::pair<double, int> &success_enumeration_gap,
                                           double min_enumeration_fail_gap,
                                           double max_gap2try_enumeration) {
            if (!outFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::ostringstream oss(std::ios::binary);
            oss.write(reinterpret_cast<const char *>(&max_enumeration_success_gap),
                      sizeof(max_enumeration_success_gap));
            oss.write(reinterpret_cast<const char *>(&success_enumeration_gap.first),
                      sizeof(success_enumeration_gap.first));
            oss.write(reinterpret_cast<const char *>(&success_enumeration_gap.second),
                      sizeof(success_enumeration_gap.second));
            oss.write(reinterpret_cast<const char *>(&min_enumeration_fail_gap), sizeof(min_enumeration_fail_gap));
            oss.write(reinterpret_cast<const char *>(&max_gap2try_enumeration), sizeof(max_gap2try_enumeration));

            writeCompressedData(outFile, oss.str());
        }

        template<typename Map>
        void writeMap(const Map &mapData) {
            std::ostringstream oss(std::ios::binary);
            const int mapSize = static_cast<int>(mapData.size());
            oss.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
            for (const auto &pair: mapData) {
                oss.write(reinterpret_cast<const char *>(&pair.first.first), sizeof(pair.first.first));
                oss.write(reinterpret_cast<const char *>(&pair.first.second), sizeof(pair.first.second));
                oss.write(reinterpret_cast<const char *>(&pair.second.first), sizeof(pair.second.first));
                oss.write(reinterpret_cast<const char *>(&pair.second.second), sizeof(pair.second.second));
            }
            writeCompressedData(outFile, oss.str());
        }

        inline void writeBranchInfo(const Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                                    const Branching::BKF::BKFDataShared &bkf_data_shared) {
            writeMap(bkf_data_shared.getRStarDepth());
            writeMap(history.exact_improvement_down);
            writeMap(history.exact_improvement_up);
            writeMap(history.heuristic_improvement_down);
            writeMap(history.heuristic_improvement_up);
            writeMap(history.lp_testing_improvement_down);
            writeMap(history.lp_testing_improvement_up);
            writeMap(history.increase_depth);
        }
    }


    template<bool if_symmetry>
    void TwoStageController::writeNodeOut(const std::string &ins_name,
                                          const Solver &node_solver,
                                          const std::vector<SequenceInfo> &cols,
                                          const std::vector<Rcc> &rccs,
                                          const std::vector<R1c> &r1cs,
                                          const std::vector<Brc> &brcs,
                                          double val, int idx, bool if_enu,
                                          int dim, int num_buckets_per_vertex, res_int step_size,
                                          VecLabel *const *if_exist_extra_labels_in_forward_sense,
                                          Bucket *const *all_forward_buckets,
                                          VecLabel *const *if_exist_extra_labels_in_backward_sense,
                                          Bucket *const *all_backward_buckets,
                                          const std::vector<bool> &deleted_columns_in_enumeration_pool,
                                          const RowVectorXT &index_columns_in_enumeration_column_pool,
                                          const RowVectorXd &cost_for_columns_in_enumeration_column_pool,
                                          const int *col_pool4_pricing,
                                          double max_enumeration_success_gap,
                                          const std::pair<double, int> &success_enumeration_gap,
                                          double min_enumeration_fail_gap,
                                          double max_gap2try_enumeration,
                                          const Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                                          const Branching::BKF::BKFDataShared &bkf_data_shared) {
        int recommended_memory_usage;
        int col_pool_size = static_cast<int>(index_columns_in_enumeration_column_pool.size());
        if (if_enu) {
            if (col_pool_size > static_cast<int>(COL_POOL_TYPE::LARGE)) {
                recommended_memory_usage = static_cast<int>(OUT_NODE_MEMORY_USE::LARGE);
            } else if (col_pool_size > static_cast<int>(COL_POOL_TYPE::MID)) {
                recommended_memory_usage = static_cast<int>(OUT_NODE_MEMORY_USE::MID);
            } else {
                recommended_memory_usage = static_cast<int>(OUT_NODE_MEMORY_USE::SMALL);
            }
        } else {
            recommended_memory_usage = static_cast<int>(OUT_NODE_MEMORY_USE::MID);
        }
        OutNodeNameSpace::initFileStream(ins_name, recommended_memory_usage);

        OutNodeNameSpace::writeColSequence(node_solver, cols);
        OutNodeNameSpace::writeCutOut(rccs, r1cs, brcs);

        OutNodeNameSpace::writeEnuState(if_enu);
        OutNodeNameSpace::writeNodeInfo(val, idx);

        if (if_enu) {
            OutNodeNameSpace::writeEnuMatrix(
                deleted_columns_in_enumeration_pool,
                index_columns_in_enumeration_column_pool,
                cost_for_columns_in_enumeration_column_pool,
                col_pool4_pricing);
        } else {
            OutNodeNameSpace::writeBucketGraphOut<if_symmetry>(
                dim, num_buckets_per_vertex, step_size,
                if_exist_extra_labels_in_forward_sense,
                all_forward_buckets,
                if_exist_extra_labels_in_backward_sense,
                all_backward_buckets);
        }

        OutNodeNameSpace::writeEnumerationTryGap(
            max_enumeration_success_gap,
            success_enumeration_gap,
            min_enumeration_fail_gap,
            max_gap2try_enumeration);

        OutNodeNameSpace::writeBranchInfo(history, bkf_data_shared);

        //close the file stream
        OutNodeNameSpace::closeFileStream();
    }


    inline void TwoStageController::updateUB(const std::string &ub_name, double &ub) {
        std::filesystem::create_directories(std::string(UB_FOLDER));

        std::string filePath = std::string(UB_FOLDER) + "/" + ub_name + std::string(UB_FILE_SUFFIX);

        double fileUB = ub;
        bool fileExists = false;


        std::ifstream inFile(filePath);
        if (inFile.is_open()) {
            fileExists = true;
            inFile >> fileUB;
            inFile.close();
        }

        // Always use the smaller value
        double bestUB = std::min(ub, fileUB);
        ub = bestUB;

        // Update the file if it doesn't exist or if its value is larger than the new bestUB.
        if (!fileExists || fileUB != bestUB) {
            std::ofstream outFile(filePath);
            if (outFile.is_open()) {
                outFile << bestUB;
                outFile.close();
            }
        }
    }
}

#endif // ROUTE_OPT_T_NODE_OUT_HPP
