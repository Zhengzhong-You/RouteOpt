/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include <iostream>
#include "helper_rank1_arc_based_memory.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    void constructMemoryArcBased(const Rank1CutsDataShared &rank1CutsDataShared, DataShared &sharedData) {
        std::vector<Eigen::Triplet<int> > tripletList;
        const auto &sol = sharedData.getSol();
        const auto dim = rank1CutsDataShared.getDim();
        tripletList.reserve(static_cast<size_t>(static_cast<double>(dim * sol.size()) * RANK1_ROW_DENSITY));
        std::unordered_map<int, int> map_vertex_index;
        for (int i = 0; i < sol.size(); ++i) {
            map_vertex_index.clear();
            auto &s = sol[i];
            for (int j: s.col_seq) {
                ++map_vertex_index[j];
            }
            auto old_size = tripletList.size();
            tripletList.resize(old_size + map_vertex_index.size());
            for (auto &it: map_vertex_index) {
                tripletList[old_size++] = {it.first, i, it.second};
            }
        }
        sparseRowMatrixXI sol_matrix(dim, static_cast<int>(sol.size()));
        SAFE_EIGEN(sol_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
        )
        auto &cuts = sharedData.refCuts();
        for (auto it = cuts.begin(); it != cuts.end();) {
            auto &cut = *it;
            bool if_suc;
            findLeastMemoryArcBased(rank1CutsDataShared, sharedData, sol_matrix, cut, if_suc);
            if (!if_suc) {
                it = cuts.erase(it);
            } else ++it;
        }
    }


    void findLeastMemoryArcBased(const Rank1CutsDataShared &rank1CutsDataShared,
                                 const DataShared &sharedData,
                                 const sparseRowMatrixXI &sol_matrix, R1c &cut, bool &if_suc) {
        if_suc = true;
        auto &arc_mem = cut.arc_mem;
        std::unordered_set<std::pair<int, int>, PairHasher> existing_arcs;
        if (!arc_mem.empty()) {
            for (auto &info: arc_mem) {
                int end = info.second;
                for (int j: info.first) {
                    existing_arcs.emplace(j, end);
                }
            }
        }


        const auto &sol = sharedData.getSol();
        std::unordered_map<int, int> map_vertex_state;
        auto &multiplier = rank1CutsDataShared.getMultiplier(static_cast<int>(cut.info_r1c.first.size()),
                                                         cut.info_r1c.second);
        for (int i = 0; i < multiplier.size(); ++i) {
            map_vertex_state[cut.info_r1c.first[i]] = multiplier[i];
        }
        int denominator = rank1CutsDataShared.getDenominator(static_cast<int>(cut.info_r1c.first.size()),
                                                             cut.info_r1c.second);
        sparseRowVectorXI vec(sol_matrix.cols());
        for (auto i: cut.info_r1c.first) {
            vec += sol_matrix.row(i);
        }
        vec /= denominator;
        vec.prune(0);

        std::vector<Arcs> all_arcs;
        for (sparseRowVectorXI::InnerIterator it(vec, 0); it; ++it) {
            int col = static_cast<int>(it.col());
            std::vector<std::vector<int> > ps;
            std::vector<int> vertex_states;
            std::vector<std::unordered_set<std::pair<int, int>, PairHasher> > arcs;
            getVertexStates(sol[col].col_seq,
                            sol[col].forward_concatenate_pos,
                            map_vertex_state,
                            vertex_states,
                            arcs);
            findLeastPlans2MakeCoeffRight(vertex_states, denominator, ps, if_suc);
            if (!if_suc) {
                PRINT_WARNING("std::find least plans failed");
            } else {
                auto &arc = all_arcs.emplace_back();
                for (auto &p: ps) {
                    auto &arc2 = arc.arc_plan.emplace_back();
                    for (auto &it2: p) {
                        arc2.insert(arcs[it2].begin(), arcs[it2].end());
                    }
                }
            }
        }

        getLeastMemory(all_arcs, existing_arcs, if_suc);


        if (if_suc) {
            arc_mem.clear();
            std::unordered_map<int, std::vector<int> > map_vertex_arcs;
            cutLong tmp = 0;
            for (auto i: cut.info_r1c.first)tmp.set(i);
            for (auto &it: existing_arcs) {
                if (tmp.test(it.second)) continue;
                map_vertex_arcs[it.second].emplace_back(it.first);
            }
            for (auto &it: map_vertex_arcs) {
                arc_mem.emplace_back(it.second, it.first);
            }
        }
    }

    void getVertexStates(const std::vector<int> &sequence,
                         const int forward_pos,
                         const std::unordered_map<int, int> &mp_,
                         std::vector<int> &vertex_states,
                         std::vector<std::unordered_set<std::pair<int, int>, PairHasher> > &arcs) {
        vertex_states.clear();
        arcs.clear();
        std::vector<std::vector<int> > tmp_arcs;
        std::vector<int> tmp;
        std::pair<int, int> change_direction{-1, -1};
        for (int i = 0; i < sequence.size(); ++i) {
            int current = sequence[i];
            if (mp_.find(current) != mp_.end()) {
                vertex_states.emplace_back(mp_.at(current));
                tmp.emplace_back(current);
                tmp_arcs.emplace_back(tmp);
                tmp = {current};
            } else tmp.emplace_back(current);
            if (i == forward_pos) change_direction = {static_cast<int>(tmp_arcs.size()), static_cast<int>(tmp.size())};
        }
        tmp_arcs.emplace_back(tmp);

        int size = static_cast<int>(tmp_arcs.size()) - 1;
        if (change_direction.first != -1) {
            int idx = change_direction.first;
            for (int i = 1; i < idx; ++i) {
                auto &arc = tmp_arcs[i];
                auto &arc_set = arcs.emplace_back();
                for (int j = 0; j < static_cast<int>(arc.size()) - 1; ++j) {
                    arc_set.emplace(arc[j], arc[j + 1]);
                }
            } {
                int i = idx;
                if (i >= 1 && i <= size - 1) {
                    int idx2 = change_direction.second - 1;
                    auto &arc = tmp_arcs[i];
                    auto &arc_set = arcs.emplace_back();
                    for (int j = 0; j < idx2; ++j) {
                        arc_set.emplace(arc[j], arc[j + 1]);
                    }
                    for (int j = static_cast<int>(arc.size()) - 1; j > idx2 + 1; --j) {
                        arc_set.emplace(arc[j], arc[j - 1]);
                    }
                }
            }
            for (int i = idx + 1; i < size; ++i) {
                auto &arc = tmp_arcs[i];
                auto &arc_set = arcs.emplace_back();
                for (int j = static_cast<int>(arc.size()) - 1; j > 0; --j) {
                    arc_set.emplace(arc[j], arc[j - 1]);
                }
            }
        } else {
            for (int i = 1; i < size; ++i) {
                auto &arc = tmp_arcs[i];
                auto &arc_set = arcs.emplace_back();
                for (int j = 0; j < static_cast<int>(arc.size()) - 1; ++j) {
                    arc_set.emplace(arc[j], arc[j + 1]);
                }
            }
        }
    }

    void findLeastPlans2MakeCoeffRight(const std::vector<int> &vertex_states,
                                       int denominator,
                                       std::vector<std::vector<int> > &plans,
                                       bool &if_suc) {
        if_suc = true;
        if (vertex_states.empty() || vertex_states.size() >= MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN) {
            if_suc = false;
            return;
        }
        /**
         * dp: cut one column to several segments
         * e.g.: --------|------------|--------------|----------------
         * 	                 0            1
         * we have now 2 segments.
         * and we std::set initial state () to be 0 and we do extend the state where it must be larger than sparse_rep.back()
         * the extension ends when the coeff reaches the std::max coeff,
         * dominance rule 1(extension): used segments are subset; coeff is same; state is the same
         * dominance rule 2(finish): used segments are subset;
         */
        std::vector<int> segments(vertex_states.size() - 1);
        for (int i = 0; i < segments.size(); ++i) segments[i] = vertex_states[i + 1] + vertex_states[i];
        int max_coeff = std::accumulate(vertex_states.begin(), vertex_states.end(), 0) / denominator;
        std::vector<std::vector<std::pair<std::vector<State>, int> > > states(max_coeff + 1,
                                                                              std::vector<std::pair<std::vector<State>,
                                                                                  int> >(denominator));
        for (int i = 0; i < states.size(); ++i) {
            for (int j = 0; j < states[i].size(); ++j) {
                states[i][j].first.resize(static_cast<int>(std::pow(2, i)) < 100
                                              ? static_cast<int>(std::pow(2, i))
                                              : 100);
                states[i][j].second = 0;
            }
        }
        for (int i = 0; i < segments.size(); ++i) {
            auto &seg = segments[i];
            int coeff = seg / denominator;
            int state = seg % denominator;
            auto &state_vec = states[coeff][state].first;
            auto &num = states[coeff][state].second;
            if (num == state_vec.size()) state_vec.resize(state_vec.size() * 2);
            state_vec[num].state = state;
            state_vec[num].coeff = coeff;
            state_vec[num].end_segment = i;
            state_vec[num].bit.set(i);
            ++num;
        }
        for (int i = 0; i < static_cast<int>(states.size()) - 1; ++i) {
            for (int s = 0; s < states[i].size(); ++s) {
                auto &state = states[i][s];
                auto &state_vec = state.first;
                auto &num = state.second;
                for (int j = 0; j < num; ++j) {
                    auto &label = state_vec[j];
                    for (int k = label.end_segment + 1; k < segments.size(); ++k) {
                        State new_label{};
                        new_label.state =
                                k == label.end_segment + 1 ? label.state + segments[k] - vertex_states[k] : segments[k];
                        new_label.coeff = label.coeff + new_label.state / denominator;
                        if (new_label.coeff >= max_coeff) new_label.state = 0;
                        else new_label.state %= denominator;
                        if (new_label.state == label.state && new_label.coeff == label.coeff) continue;
                        new_label.bit = label.bit;
                        new_label.bit.set(k);
                        new_label.end_segment = k;
                        bool is_keep{true};
                        auto &dominance_vec = states[new_label.coeff][new_label.state].first;
                        auto &dominance_num = states[new_label.coeff][new_label.state].second;
                        for (int l = 0; l < dominance_num;) {
                            auto &dominance_label = dominance_vec[l];
                            if ((dominance_label.bit & new_label.bit) == dominance_label.bit) {
                                is_keep = false;
                                break;
                            } else if ((dominance_label.bit & new_label.bit) == new_label.bit) {
                                dominance_label = dominance_vec[--dominance_num];
                            } else ++l;
                        }
                        if (is_keep) {
                            if (dominance_num == dominance_vec.size()) {
                                dominance_vec.resize(dominance_vec.size() * 2);
                            }
                            dominance_vec[dominance_num++] = new_label;
                        }
                    }
                }
            }
        }

        auto &last_state = states[max_coeff][0];
        auto &last_state_vec = last_state.first;
        auto &last_num = last_state.second;
        plans.clear();
        plans.resize(last_num);
        for (int i = 0; i < last_num; ++i) {
            auto &label = last_state_vec[i];
            for (int j = 0; j < MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN; ++j) {
                if (label.bit.test(j)) {
                    plans[i].emplace_back(j);
                }
            }
        }
        if (plans.empty())
            THROW_RUNTIME_ERROR("plans.empty()");
    }

    void getLeastMemory(std::vector<Arcs> &all_arcs, std::unordered_set<std::pair<int, int>, PairHasher> &existing_arcs,
                        bool &if_suc) {
        bool loop = true;
        while (loop) {
            int old_size = static_cast<int>(existing_arcs.size());
            reduceArcs(all_arcs, existing_arcs, if_suc);
            if (old_size == existing_arcs.size()) loop = false;
            if (!if_suc) return;
        }

        if (all_arcs.empty()) return;
        std::vector<std::vector<arcBit> > bins(all_arcs.size());
        bins[0].resize(all_arcs[0].arc_plan.size());
        for (int i = 1; i < all_arcs.size(); ++i) {
            int tmp = 2 * bins[i - 1].size() < 100 ? static_cast<int>(2 * bins[i - 1].size()) : 100;
            bins[i].resize(tmp);
        }
        std::vector<int> bin_num(all_arcs.size(), 0);
        std::unordered_map<std::pair<int, int>, int, PairHasher> arc_bit_map;
        std::unordered_map<int, std::pair<int, int> > bit_arc_map;
        int bit_pos = 0;
        for (auto &it: all_arcs) {
            for (auto &arcs: it.arc_plan) {
                for (auto &arc: arcs) {
                    if (arc_bit_map.find(arc) == arc_bit_map.end()) {
                        bit_arc_map[bit_pos] = arc;
                        arc_bit_map[arc] = bit_pos++;
                        if (bit_pos > MAX_UNDETERMINED_ARC_NUMBER) {
                            if_suc = false;
                            return;
                        }
                    }
                }
            }
        }
        std::vector<std::vector<arcBit> > all_arcs_bit(all_arcs.size());
        for (int i = 0; i < all_arcs.size(); ++i) {
            auto &arcs = all_arcs[i].arc_plan;
            auto &arcs_bit = all_arcs_bit[i];
            arcs_bit.resize(arcs.size());
            for (int j = 0; j < arcs.size(); ++j) {
                auto &arc = arcs[j];
                for (auto &it: arc) {
                    arcs_bit[j].set(arc_bit_map[it]);
                }
            }
        }
        bins[0] = all_arcs_bit[0];
        bin_num[0] = static_cast<int>(all_arcs_bit[0].size());
        for (int i = 0; i < static_cast<int>(bins.size()) - 1; ++i) {
            auto &bin = bins[i];
            int num = bin_num[i];
            auto &bin_next = bins[i + 1];
            auto &bin_num_next = bin_num[i + 1];
            auto &add_arcs = all_arcs[i + 1].arc_plan;
            for (int j = 0; j < num; ++j) {
                auto &arc_bit_1 = bin[j];
                for (int k = 0; k < add_arcs.size(); ++k) {
                    arcBit new_arc_bit = arc_bit_1 | all_arcs_bit[i + 1][k];
                    bool if_keep = true;
                    for (int l = 0; l < bin_num_next;) {
                        if ((new_arc_bit & bin_next[l]) == bin_next[l]) {
                            if_keep = false;
                            break;
                        }
                        if ((new_arc_bit & bin_next[l]) == new_arc_bit) {
                            bin_next[l] = bin_next[--bin_num_next];
                        } else ++l;
                    }
                    if (if_keep) {
                        if (bin_num_next == bin_next.size()) {
                            bin_next.resize(bin_next.size() * 2);
                        }
                        bin_next[bin_num_next++] = new_arc_bit;
                        if (bin_num_next > MAX_LABELS) {
                            if_suc = false;
                            return;
                        }
                    }
                }
            }
        }

        auto &least_memory = bins.back();
        int least_memory_num = bin_num.back();
        arcBit best = least_memory[0];
        for (int i = 1; i < least_memory_num; ++i) {
            auto &arc_bit = least_memory[i];
            if (best.count() > arc_bit.count()) best = arc_bit;
        }
        for (int i = 0; i < bit_pos; ++i) {
            if (best.test(i)) {
                if (existing_arcs.find(bit_arc_map[i]) != existing_arcs.end())
                    THROW_RUNTIME_ERROR("existing_arcs.find(bit_arc_map[i]) != existing_arcs.end()");
                existing_arcs.emplace(bit_arc_map[i]);
            }
        }
    }

    void reduceArcs(std::vector<Arcs> &all_arcs, std::unordered_set<std::pair<int, int>, PairHasher> &existing_arcs,
                    bool &if_suc) {
        if_suc = true;
        std::unordered_map<std::pair<int, int>, std::unordered_set<int>, PairHasher> arc_map;
        for (auto &it: all_arcs) {
            auto &arc_plan = it.arc_plan;
            arc_map.clear();
            for (int i = 0; i < arc_plan.size(); ++i) {
                auto &arcs = arc_plan[i];
                for (auto &arc: arcs) {
                    arc_map[arc].emplace(i);
                }
            }
            for (auto &it2: arc_map) {
                if (it2.second.size() == arc_plan.size()) {
                    existing_arcs.emplace(it2.first);
                }
            }
        }
        for (auto it = all_arcs.begin(); it != all_arcs.end();) {
            auto &arc_plan = it->arc_plan;
            bool if_delete_all = false;
            for (auto it2 = arc_plan.begin(); it2 != arc_plan.end();) {
                auto &arcs = *it2;
                for (auto it3 = arcs.begin(); it3 != arcs.end();) {
                    if (existing_arcs.find(*it3) != existing_arcs.end()) {
                        it3 = arcs.erase(it3);
                    } else ++it3;
                }
                if (arcs.empty()) {
                    if_delete_all = true;
                    break;
                }
                ++it2;
            }
            if (if_delete_all) {
                it = all_arcs.erase(it);
            } else ++it;
        }
        std::unordered_map<std::pair<int, int>, int, PairHasher> arc_bit_map;
        int bit_pos = 0;
        for (auto &it: all_arcs) {
            auto &arc_plan = it.arc_plan;
            for (auto &arcs: arc_plan) {
                for (auto &arc: arcs) {
                    if (arc_bit_map.find(arc) == arc_bit_map.end()) {
                        arc_bit_map[arc] = bit_pos++;
                        if (bit_pos > MAX_UNDETERMINED_ARC_NUMBER) {
                            if_suc = false;
                            return;
                        }
                    }
                }
            }
        }

        for (auto it = all_arcs.begin(); it != all_arcs.end();) {
            auto &arc_plan = it->arc_plan;
            std::vector<arcBit> arc_bit(arc_plan.size(), 0);
            for (int i = 0; i < arc_plan.size(); ++i) {
                auto &arcs = arc_plan[i];
                for (auto &arc: arcs) {
                    arc_bit[i].set(arc_bit_map[arc]);
                }
            }
            for (int i = 0; i < arc_plan.size(); ++i) {
                for (int j = i + 1; j < arc_plan.size(); ++j) {
                    if ((arc_bit[i] & arc_bit[j]) == arc_bit[i]) {
                        arc_plan.erase(arc_plan.begin() + j);
                        arc_bit.erase(arc_bit.begin() + j);
                        --j;
                    } else if ((arc_bit[i] & arc_bit[j]) == arc_bit[j]) {
                        arc_plan.erase(arc_plan.begin() + i);
                        arc_bit.erase(arc_bit.begin() + i);
                        --i;
                        break;
                    }
                }
            }
            if (arc_plan.size() == 1) {
                auto &arcs = arc_plan[0];
                for (auto &arc: arcs) {
                    existing_arcs.emplace(arc);
                }
                it = all_arcs.erase(it);
            } else ++it;
        }
    }
}
