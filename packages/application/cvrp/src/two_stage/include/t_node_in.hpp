/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_T_NODE_IN_HPP
#define ROUTE_OPT_T_NODE_IN_HPP
#include <fstream>
#include <iostream>
#include <zlib.h>


#include "cvrp_pricing_controller.hpp"
#include "label.hpp"
#include "two_stage_macro.hpp"
#include "two_stage_controller.hpp"

namespace RouteOpt::Application::CVRP {
    namespace InNodeNameSpace {
        inline std::ifstream inFile{};

        inline void initInputFileStream(const std::string &node_file_name) {
            inFile.open(node_file_name, std::ios::binary | std::ios::in);
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("cannot open file stream.");
        }

        inline void closeInputFileStream() {
            inFile.clear(); // clear any error flags
            inFile.close();
        }

        inline void readCompressedData(std::istream &in, std::string &rawData) {
            // Read header information: uncompressed size and compressed size.
            size_t srcLen = 0;
            size_t compLen = 0;
            in.read(reinterpret_cast<char *>(&srcLen), sizeof(srcLen));
            in.read(reinterpret_cast<char *>(&compLen), sizeof(compLen));
            if (!in)
                THROW_RUNTIME_ERROR("error reading compressed header.");

            // Read the compressed block.
            std::vector<Bytef> compressedData(compLen);
            in.read(reinterpret_cast<char *>(compressedData.data()), static_cast<std::streamsize>(compLen));
            if (!in)
                THROW_RUNTIME_ERROR("error reading compressed data.");

            // Prepare the buffer for the uncompressed data.
            rawData.resize(srcLen, '\0');
            uLong rawDataSize = srcLen;
            int ret = uncompress(reinterpret_cast<Bytef *>(&rawData[0]), &rawDataSize,
                                 compressedData.data(), compLen);
            if (ret != Z_OK)
                THROW_RUNTIME_ERROR("decompression failed with error code: " + std::to_string(ret));
        }

        inline void readColSequence(std::vector<double> &cost, std::vector<SequenceInfo> &cols) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);

            std::istringstream iss(rawData, std::ios::binary);

            int num_cols = 0;
            iss.read(reinterpret_cast<char *>(&num_cols), sizeof(num_cols));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading number of columns.");

            cols.resize(num_cols);
            cost.resize(num_cols);

            for (int i = 0; i < num_cols; ++i) {
                int num_seq = 0;
                iss.read(reinterpret_cast<char *>(&num_seq), sizeof(num_seq));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading size of col_seq for column " + std::to_string(i));

                auto &seqInfo = cols[i];
                seqInfo.col_seq.resize(num_seq);
                for (int j = 0; j < num_seq; ++j) {
                    iss.read(reinterpret_cast<char *>(&seqInfo.col_seq[j]), sizeof(seqInfo.col_seq[j]));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading col_seq element " + std::to_string(j) +
                        " for column " + std::to_string(i));
                }

                iss.read(reinterpret_cast<char *>(&seqInfo.forward_concatenate_pos),
                         sizeof(seqInfo.forward_concatenate_pos));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading forward_concatenate_pos for column " + std::to_string(i));

                double obj = 0.0;
                iss.read(reinterpret_cast<char *>(&obj), sizeof(obj));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading objective value for column " + std::to_string(i));

                cost[i] = obj;
            }
        }

        inline void readCutOut(std::vector<Rcc> &rccs, std::vector<R1c> &r1cs, std::vector<Brc> &brcs) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);

            std::istringstream iss(rawData, std::ios::binary);

            // Read rccs data.
            int num_rccs = 0;
            iss.read(reinterpret_cast<char *>(&num_rccs), sizeof(num_rccs));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading number of rccs.");
            rccs.resize(num_rccs);
            for (int i = 0; i < num_rccs; ++i) {
                Rcc &rcc = rccs[i];
                iss.read(reinterpret_cast<char *>(&rcc.idx_rcc), sizeof(rcc.idx_rcc));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading rcc idx_rcc for rcc " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&rcc.rhs), sizeof(rcc.rhs));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading rcc rhs for rcc " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&rcc.form_rcc), sizeof(rcc.form_rcc));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading rcc form_rcc for rcc " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&rcc.if_keep), sizeof(rcc.if_keep));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading rcc if_keep for rcc " + std::to_string(i));

                int customer_count = 0;
                iss.read(reinterpret_cast<char *>(&customer_count), sizeof(customer_count));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading customer count for rcc " + std::to_string(i));
                rcc.info_rcc_customer.resize(customer_count);
                for (int j = 0; j < customer_count; ++j) {
                    iss.read(reinterpret_cast<char *>(&rcc.info_rcc_customer[j]), sizeof(rcc.info_rcc_customer[j]));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading rcc info_rcc_customer element " + std::to_string(j) +
                        " for rcc " + std::to_string(i));
                }

                int outside_count = 0;
                iss.read(reinterpret_cast<char *>(&outside_count), sizeof(outside_count));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading outside customer count for rcc " + std::to_string(i));
                rcc.info_rcc_outside_customer.resize(outside_count);
                for (int j = 0; j < outside_count; ++j) {
                    iss.read(reinterpret_cast<char *>(&rcc.info_rcc_outside_customer[j]),
                             sizeof(rcc.info_rcc_outside_customer[j]));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading rcc info_rcc_outside_customer element " + std::to_string(j) +
                        " for rcc " + std::to_string(i));
                }
            }

            // Read r1cs data.
            int num_r1cs = 0;
            iss.read(reinterpret_cast<char *>(&num_r1cs), sizeof(num_r1cs));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading number of r1cs.");
            r1cs.resize(num_r1cs);
            for (int i = 0; i < num_r1cs; ++i) {
                R1c &r1c = r1cs[i];
                iss.read(reinterpret_cast<char *>(&r1c.idx_r1c), sizeof(r1c.idx_r1c));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading r1c idx_r1c for r1c " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&r1c.rhs), sizeof(r1c.rhs));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading r1c rhs for r1c " + std::to_string(i));

                int info_count = 0;
                iss.read(reinterpret_cast<char *>(&info_count), sizeof(info_count));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading info count for r1c " + std::to_string(i));
                r1c.info_r1c.first.resize(info_count);
                for (int j = 0; j < info_count; ++j) {
                    iss.read(reinterpret_cast<char *>(&r1c.info_r1c.first[j]), sizeof(r1c.info_r1c.first[j]));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading r1c info element " + std::to_string(j) +
                        " for r1c " + std::to_string(i));
                }
                iss.read(reinterpret_cast<char *>(&r1c.info_r1c.second), sizeof(r1c.info_r1c.second));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading r1c info second for r1c " + std::to_string(i));

                int arc_mem_count = 0;
                iss.read(reinterpret_cast<char *>(&arc_mem_count), sizeof(arc_mem_count));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading arc_mem count for r1c " + std::to_string(i));
                r1c.arc_mem.resize(arc_mem_count);
                for (int j = 0; j < arc_mem_count; ++j) {
                    auto &arc_mem_pair = r1c.arc_mem[j];
                    int arc_size = 0;
                    iss.read(reinterpret_cast<char *>(&arc_size), sizeof(arc_size));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading arc size for arc_mem " + std::to_string(j) +
                        " for r1c " + std::to_string(i));
                    arc_mem_pair.first.resize(arc_size);
                    for (int k = 0; k < arc_size; ++k) {
                        iss.read(reinterpret_cast<char *>(&arc_mem_pair.first[k]), sizeof(arc_mem_pair.first[k]));
                        if (!iss)
                            THROW_RUNTIME_ERROR("error reading arc_mem element " + std::to_string(k) +
                            " for arc_mem " + std::to_string(j) +
                            " for r1c " + std::to_string(i));
                    }
                    iss.read(reinterpret_cast<char *>(&arc_mem_pair.second), sizeof(arc_mem_pair.second));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading arc_mem mem value for arc_mem " + std::to_string(j) +
                        " for r1c " + std::to_string(i));
                }
            }

            // Read brcs data.
            int num_brcs = 0;
            iss.read(reinterpret_cast<char *>(&num_brcs), sizeof(num_brcs));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading number of brcs.");
            brcs.resize(num_brcs);
            for (int i = 0; i < num_brcs; ++i) {
                Brc &brc = brcs[i];
                iss.read(reinterpret_cast<char *>(&brc.idx_brc), sizeof(brc.idx_brc));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading brc idx_brc for brc " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&brc.edge.first), sizeof(brc.edge.first));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading brc edge.first for brc " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&brc.edge.second), sizeof(brc.edge.second));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading brc edge.second for brc " + std::to_string(i));
                iss.read(reinterpret_cast<char *>(&brc.br_dir), sizeof(brc.br_dir));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading brc br_dir for brc " + std::to_string(i));
            }
        }

        inline void readEnuState(bool &if_enu) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);

            unsigned char u = 0;
            std::istringstream iss(rawData, std::ios::binary);
            iss.read(reinterpret_cast<char *>(&u), sizeof(u));
            if_enu = (u != 0);
        }

        inline void readNodeInfo(double &val, int &idx) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);
            std::istringstream iss(rawData, std::ios::binary);


            // read header values
            iss.read(reinterpret_cast<char *>(&val), sizeof(val));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading val");
            iss.read(reinterpret_cast<char *>(&idx), sizeof(idx));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading idx");
        }


        template<bool if_symmetry>
        void readBucketGraphOut(
            int dim, int &num_buckets_per_vertex, res_int &step_size,
            VecLabel **&if_exist_extra_labels_in_forward_sense,
            Bucket **&all_forward_buckets,
            VecLabel **&if_exist_extra_labels_in_backward_sense,
            Bucket **&all_backward_buckets) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);
            std::istringstream iss(rawData, std::ios::binary);


            iss.read(reinterpret_cast<char *>(&step_size), sizeof(step_size));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading step_size");

            iss.read(reinterpret_cast<char *>(&num_buckets_per_vertex), sizeof(num_buckets_per_vertex));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading num_buckets_per_vertex");

            // allocate forward arrays
            if_exist_extra_labels_in_forward_sense = new VecLabel *[dim];
            all_forward_buckets = new Bucket *[dim];
            for (int i = 0; i < dim; ++i) {
                if_exist_extra_labels_in_forward_sense[i] = new VecLabel[num_buckets_per_vertex];
                all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
            }

            // read forward bucket graph data
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < num_buckets_per_vertex; ++j) {
                    // read the size of the extra forward label array
                    int fwdLabelSize = 0;
                    iss.read(reinterpret_cast<char *>(&fwdLabelSize), sizeof(fwdLabelSize));
                    if (!iss)
                        THROW_RUNTIME_ERROR(
                        "error reading fwd label size for bucket (" + std::to_string(i) + "," + std::to_string(j) +
                        ")");
                    if_exist_extra_labels_in_forward_sense[i][j].first.resize(fwdLabelSize);

                    // read forward bucket arcs
                    int arc_count = 0;
                    iss.read(reinterpret_cast<char *>(&arc_count), sizeof(arc_count));
                    if (!iss)
                        THROW_RUNTIME_ERROR(
                        "error reading arc count for forward bucket (" + std::to_string(i) + "," + std::to_string(j) +
                        ")");
                    all_forward_buckets[i][j].bucket_arcs.resize(arc_count);
                    for (int a = 0; a < arc_count; ++a) {
                        int arc = 0;
                        iss.read(reinterpret_cast<char *>(&arc), sizeof(arc));
                        if (!iss)
                            THROW_RUNTIME_ERROR("error reading forward bucket arc " + std::to_string(a) +
                            " for bucket (" + std::to_string(i) + "," + std::to_string(j) + ")");
                        all_forward_buckets[i][j].bucket_arcs[a] = arc;
                    }

                    // read forward jump arcs
                    int jump_count = 0;
                    iss.read(reinterpret_cast<char *>(&jump_count), sizeof(jump_count));
                    if (!iss)
                        THROW_RUNTIME_ERROR(
                        "error reading jump count for forward bucket (" + std::to_string(i) + "," + std::to_string(j) +
                        ")");
                    all_forward_buckets[i][j].jump_arcs.resize(jump_count);
                    for (int k = 0; k < jump_count; ++k) {
                        res_int fst = 0;
                        int snd = 0;
                        iss.read(reinterpret_cast<char *>(&fst), sizeof(fst));
                        if (!iss)
                            THROW_RUNTIME_ERROR(
                            "error reading forward jump arc fst for bucket (" + std::to_string(i) + "," + std::to_string
                            (j) + ")");
                        iss.read(reinterpret_cast<char *>(&snd), sizeof(snd));
                        if (!iss)
                            THROW_RUNTIME_ERROR(
                            "error reading forward jump arc snd for bucket (" + std::to_string(i) + "," + std::to_string
                            (j) + ")");
                        all_forward_buckets[i][j].jump_arcs[k] = std::make_pair(fst, snd);
                    }
                }
            }

            // read backward bucket graph data if symmetry is not enabled
            if constexpr (!if_symmetry) {
                if_exist_extra_labels_in_backward_sense = new VecLabel *[dim];
                all_backward_buckets = new Bucket *[dim];
                for (int i = 0; i < dim; ++i) {
                    if_exist_extra_labels_in_backward_sense[i] = new VecLabel[num_buckets_per_vertex];
                    all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
                }
                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < num_buckets_per_vertex; ++j) {
                        int bwdLabelSize = 0;
                        iss.read(reinterpret_cast<char *>(&bwdLabelSize), sizeof(bwdLabelSize));
                        if (!iss)
                            THROW_RUNTIME_ERROR(
                            "error reading bwd label size for bucket (" + std::to_string(i) + "," + std::to_string(j) +
                            ")");

                        if_exist_extra_labels_in_backward_sense[i][j].first.resize(bwdLabelSize);

                        int arc_count_backward = 0;
                        iss.read(reinterpret_cast<char *>(&arc_count_backward), sizeof(arc_count_backward));
                        if (!iss)
                            THROW_RUNTIME_ERROR(
                            "error reading arc count for backward bucket (" + std::to_string(i) + "," + std::to_string(j
                            ) + ")");
                        all_backward_buckets[i][j].bucket_arcs.resize(arc_count_backward);
                        for (int a = 0; a < arc_count_backward; ++a) {
                            int arc = 0;
                            iss.read(reinterpret_cast<char *>(&arc), sizeof(arc));
                            if (!iss)
                                THROW_RUNTIME_ERROR("error reading backward bucket arc " + std::to_string(a) +
                                " for bucket (" + std::to_string(i) + "," + std::to_string(j) + ")");
                            all_backward_buckets[i][j].bucket_arcs[a] = arc;
                        }

                        int jump_count_backward = 0;
                        iss.read(reinterpret_cast<char *>(&jump_count_backward), sizeof(jump_count_backward));
                        if (!iss)
                            THROW_RUNTIME_ERROR(
                            "error reading jump count for backward bucket (" + std::to_string(i) + "," + std::to_string(
                                j) + ")");
                        all_backward_buckets[i][j].jump_arcs.resize(jump_count_backward);
                        for (int k = 0; k < jump_count_backward; ++k) {
                            res_int fst = 0;
                            int snd = 0;
                            iss.read(reinterpret_cast<char *>(&fst), sizeof(fst));
                            if (!iss)
                                THROW_RUNTIME_ERROR(
                                "error reading backward jump arc fst for bucket (" + std::to_string(i) + "," + std::
                                to_string(j) + ")");
                            iss.read(reinterpret_cast<char *>(&snd), sizeof(snd));
                            if (!iss)
                                THROW_RUNTIME_ERROR(
                                "error reading backward jump arc snd for bucket (" + std::to_string(i) + "," + std::
                                to_string(j) + ")");
                            all_backward_buckets[i][j].jump_arcs[k] = std::make_pair(fst, snd);
                        }
                    }
                }
            } else {
                // if symmetry is true, set backward pointers to nullptr.
                if_exist_extra_labels_in_backward_sense = nullptr;
                all_backward_buckets = nullptr;
            }
        }


        inline void readEnuMatrix(
            std::vector<bool> &deleted_columns_in_enumeration_pool,
            RowVectorXT &index_columns_in_enumeration_column_pool,
            RowVectorXd &cost_for_columns_in_enumeration_column_pool,
            int *&col_pool4_pricing,
            size_t &mem4_pricing,
            size_t &pool_beg4_pricing) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);
            std::istringstream iss(rawData, std::ios::binary);


            // read the count of non-deleted columns
            int cnt = 0;
            iss.read(reinterpret_cast<char *>(&cnt), sizeof(cnt));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading non-deleted columns count");

            // resize the deleted_columns_in_enumeration_pool to have cnt entries, all set to false
            deleted_columns_in_enumeration_pool.resize(cnt, false);
            // resize the index and cost vectors accordingly
            index_columns_in_enumeration_column_pool.resize(cnt);
            cost_for_columns_in_enumeration_column_pool.resize(cnt);

            std::vector<int> col_pool_vector;
            col_pool_vector.resize(static_cast<size_t>(EXPECTED_AVER_ROUTE_LENGTH * cnt * READ_MEMORY_RATIO));
            int pricing_index = 0;
            // for each non-deleted column, read its pricing chain and cost
            for (int i = 0; i < cnt; ++i) {
                int curr_node = 0;
                // the pricing chain always starts with a 0 (written by write function)
                iss.read(reinterpret_cast<char *>(&curr_node), sizeof(curr_node));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading starting curr_node for column " + std::to_string(i));
                if (curr_node != 0)
                    THROW_RUNTIME_ERROR("expected starting curr_node to be 0 for column " + std::to_string(i));


                index_columns_in_enumeration_column_pool[i] = pricing_index;
                col_pool_vector[pricing_index++] = curr_node;

                while (true) {
                    iss.read(reinterpret_cast<char *>(&curr_node), sizeof(curr_node));
                    if (!iss)
                        THROW_RUNTIME_ERROR("error reading pricing chain for column " + std::to_string(i));
                    col_pool_vector[pricing_index++] = curr_node;
                    if (pricing_index == col_pool_vector.size())
                        col_pool_vector.resize(
                            static_cast<size_t>(static_cast<double>(col_pool_vector.size()) * READ_MEMORY_RATIO));
                    if (curr_node == 0)
                        break;
                }

                double cost = 0.0;
                iss.read(reinterpret_cast<char *>(&cost), sizeof(cost));
                if (!iss)
                    THROW_RUNTIME_ERROR("error reading cost for column " + std::to_string(i));
                cost_for_columns_in_enumeration_column_pool[i] = cost;
            }
            col_pool4_pricing = new int[pricing_index];
            std::copy_n(col_pool_vector.begin(), pricing_index, col_pool4_pricing);
            mem4_pricing = pricing_index;
            pool_beg4_pricing = pricing_index;
        }

        inline void readEnumerationTryGap(double &max_enumeration_success_gap,
                                          std::pair<double, int> &success_enumeration_gap,
                                          double &min_enumeration_fail_gap,
                                          double &max_gap2try_enumeration) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);
            std::istringstream iss(rawData, std::ios::binary);

            iss.read(reinterpret_cast<char *>(&max_enumeration_success_gap), sizeof(max_enumeration_success_gap));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading max_enumeration_success_gap");

            iss.read(reinterpret_cast<char *>(&success_enumeration_gap.first), sizeof(success_enumeration_gap.first));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading success_enumeration_gap.first");

            iss.read(reinterpret_cast<char *>(&success_enumeration_gap.second), sizeof(success_enumeration_gap.second));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading success_enumeration_gap.second");

            iss.read(reinterpret_cast<char *>(&min_enumeration_fail_gap), sizeof(min_enumeration_fail_gap));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading min_enumeration_fail_gap");

            iss.read(reinterpret_cast<char *>(&max_gap2try_enumeration), sizeof(max_gap2try_enumeration));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading max_gap2try_enumeration");
        }


        template<typename T1, typename T2>
        void readPairFromStream(std::istream &is, std::pair<T1, T2> &p) {
            is.read(reinterpret_cast<char *>(&p.first), sizeof(p.first));
            if (!is)
                THROW_RUNTIME_ERROR("error reading first element of pair");
            is.read(reinterpret_cast<char *>(&p.second), sizeof(p.second));
            if (!is)
                THROW_RUNTIME_ERROR("error reading second element of pair");
        }

        template<typename Map, bool if_map = true>
        void readMap(Map &mapData) {
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);
            std::istringstream iss(rawData, std::ios::binary);

            int mapSize = 0;
            iss.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading map size");

            if constexpr (!if_map) {
                mapData.resize(mapSize);
                for (int i = 0; i < mapSize; ++i) {
                    auto &key = mapData[i].first;
                    auto &value = mapData[i].second;
                    readPairFromStream(iss, key);
                    readPairFromStream(iss, value);
                }
            } else {
                mapData.clear();
                for (int i = 0; i < mapSize; ++i) {
                    typename Map::key_type key;
                    typename Map::mapped_type value;
                    readPairFromStream(iss, key);
                    readPairFromStream(iss, value);
                    auto ret = mapData.insert(std::make_pair(key, value));
                    if (!ret.second)
                        THROW_RUNTIME_ERROR("duplicate key encountered during map read");
                }
            }
        }


        // Read branch info by reading each map in the same order as written.
        inline void readBranchInfo(Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                                   Branching::BKF::BKFDataShared &bkf_data_shared) {
            // read f;
            if (!inFile.is_open())
                THROW_RUNTIME_ERROR("file stream is not open.");

            std::string rawData;
            readCompressedData(inFile, rawData);
            std::istringstream iss(rawData, std::ios::binary);

            auto f = bkf_data_shared.getF();
            iss.read(reinterpret_cast<char *>(&f), sizeof(f));
            if (!iss)
                THROW_RUNTIME_ERROR("error reading max_enumeration_success_gap");
            bkf_data_shared.updateF(f);
            readMap<decltype(bkf_data_shared.refRStarDepth()), false>(bkf_data_shared.refRStarDepth());
            readMap(history.exact_improvement_down);
            readMap(history.exact_improvement_up);
            readMap(history.heuristic_improvement_down);
            readMap(history.heuristic_improvement_up);
            readMap(history.lp_testing_improvement_down);
            readMap(history.lp_testing_improvement_up);
            readMap<decltype(history.increase_depth), false>(history.increase_depth);
        }
    }

    template<bool if_symmetry>
    void TwoStageController::readNodeIn(const std::string &node_file_name,
                                        std::vector<double> &cost,
                                        std::vector<SequenceInfo> &cols,
                                        std::vector<Rcc> &rccs,
                                        std::vector<R1c> &r1cs,
                                        std::vector<Brc> &brcs,
                                        double &val, int &idx, bool &if_enu,
                                        int dim, int &num_buckets_per_vertex, res_int &step_size,
                                        VecLabel **&if_exist_extra_labels_in_forward_sense,
                                        Bucket **&all_forward_buckets,
                                        VecLabel **&if_exist_extra_labels_in_backward_sense,
                                        Bucket **&all_backward_buckets,
                                        std::vector<bool> &deleted_columns_in_enumeration_pool,
                                        RowVectorXT &index_columns_in_enumeration_column_pool,
                                        RowVectorXd &cost_for_columns_in_enumeration_column_pool,
                                        int *&col_pool4_pricing,
                                        size_t &mem4_pricing,
                                        size_t &pool_beg4_pricing,
                                        double &max_enumeration_success_gap,
                                        std::pair<double, int> &success_enumeration_gap,
                                        double &min_enumeration_fail_gap,
                                        double &max_gap2try_enumeration,
                                        Branching::BranchingHistory<std::pair<int, int>, PairHasher> &history,
                                        Branching::BKF::BKFDataShared &bkf_data_shared) {
        // Open the file in read mode.
        InNodeNameSpace::initInputFileStream(node_file_name);

        InNodeNameSpace::readColSequence(cost, cols);
        InNodeNameSpace::readCutOut(rccs, r1cs, brcs);

        InNodeNameSpace::readEnuState(if_enu);

        InNodeNameSpace::readNodeInfo(val, idx);

        // read enumeration data or bucket-graph data depending on if_enu.
        if (if_enu) {
            InNodeNameSpace::readEnuMatrix(deleted_columns_in_enumeration_pool,
                                           index_columns_in_enumeration_column_pool,
                                           cost_for_columns_in_enumeration_column_pool,
                                           col_pool4_pricing, mem4_pricing, pool_beg4_pricing);
        } else {
            InNodeNameSpace::readBucketGraphOut<if_symmetry>(dim, num_buckets_per_vertex, step_size,
                                                             if_exist_extra_labels_in_forward_sense,
                                                             all_forward_buckets,
                                                             if_exist_extra_labels_in_backward_sense,
                                                             all_backward_buckets);
        }

        // read enumeration try gap values
        InNodeNameSpace::readEnumerationTryGap(max_enumeration_success_gap,
                                               success_enumeration_gap,
                                               min_enumeration_fail_gap,
                                               max_gap2try_enumeration);

        // read branch info
        InNodeNameSpace::readBranchInfo(history, bkf_data_shared);

        // close the file stream
        InNodeNameSpace::closeInputFileStream();
    }

    inline void TwoStageController::deleteInFile(const std::string &node_file_name) {
        try {
            if (std::filesystem::exists(node_file_name)) {
                if (std::filesystem::remove(node_file_name)) {
                    std::cout << "File " << node_file_name << " deleted successfully." << std::endl;
                } else {
                    PRINT_WARNING("Failed to delete file " + node_file_name + ".");
                }
            } else {
                PRINT_WARNING("File " + node_file_name + " does not exist.");
            }
        } catch (const std::filesystem::filesystem_error &e) {
            std::cerr << "Filesystem error: " << e.what() << std::endl;
        }
    }
}

#endif // ROUTE_OPT_T_NODE_IN_HPP
