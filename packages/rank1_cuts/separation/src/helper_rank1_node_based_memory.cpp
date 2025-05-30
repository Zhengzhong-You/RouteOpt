/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include <iostream>
#include "rank1_macro.hpp"
#include "helper_rank1_node_based_memory.hpp"
#include "rank1_memory_finder.hpp"


namespace RouteOpt::Rank1Cuts::Separation {
    void findMem4R1CMulti(
        const Rank1CutsDataShared &rank1CutsDataShared,
        const DataShared &sharedData,
        const PRICING_HARD_LEVEL pricing_hard_level,
        Solver &solver,
        std::unordered_map<std::vector<int>, std::vector<std::vector<int> >, VectorHashInRank1> &
        rank1_multi_mem_plan_map,
        const std::pair<std::vector<int>, int> &cut_pair,
        std::unordered_set<int> &mem,
        bool &if_suc) {
        if_suc = true;
        auto &cut = cut_pair.first;
        auto plan_idx = cut_pair.second;
        int size = static_cast<int>(cut.size());
        const auto &sol = sharedData.getSol();
        const auto &v_r_map = sharedData.getVRMap();
        const auto dim = rank1CutsDataShared.getDim();
        const auto &multi = rank1CutsDataShared.getMultiplier(size, plan_idx);
        auto denominator = rank1CutsDataShared.getDenominator(size, plan_idx);
        std::vector<std::vector<std::vector<int> > > vec_data;
        std::vector<std::vector<std::unordered_set<int> > > vec_segment_route;
        std::unordered_map<int, int> map_cut_mul;
        std::vector<int> num_vis_times(sol.size(), 0);
        for (int i = 0; i < cut.size(); ++i) {
            map_cut_mul[cut[i]] = multi[i];
            for (auto &pr: v_r_map[cut[i]]) {
                num_vis_times[pr.first] += multi[i] * pr.second;
            }
        }
        std::transform(num_vis_times.begin(), num_vis_times.end(), num_vis_times.begin(),
                       [denominator](const int x) { return static_cast<int>(x / denominator); });
        cutLong mem_long = 0;
        int num = 0;
        for (auto &route: sol) {
            auto &i = route.col_seq;
            if (num_vis_times[num++] == 0) continue;
            std::vector<std::vector<int> > data;
            std::vector<int> vis;
            std::vector<std::unordered_set<int> > segment_route;
            std::unordered_set<int> tmp_seg;
            for (auto &j: i) {
                if (map_cut_mul.find(j) != map_cut_mul.end()) {
                    vis.emplace_back(map_cut_mul[j]);
                    segment_route.emplace_back(tmp_seg);
                    tmp_seg.clear();
                } else {
                    tmp_seg.insert(j);
                }
            }
            if (!segment_route.empty())
                segment_route.erase(segment_route.begin()); //std::remove the first one
            findPlan4R1CMulti(rank1CutsDataShared, rank1_multi_mem_plan_map, vis, denominator, mem_long, segment_route,
                              data);
            if (!data.empty()) {
                vec_data.emplace_back(data);
                vec_segment_route.emplace_back(segment_route);
            }
        }

        size_t cnt = 1;
        for (int i = 0; i < vec_data.size();) {
            bool if_clear = false;
            for (auto &j: vec_data[i]) {
                bool if_all_satis = true;
                for (auto k: j) {
                    for (auto l: vec_segment_route[i][k]) {
                        if (!mem_long[l]) {
                            if_all_satis = false;
                            goto outside;
                        }
                    }
                }
            outside:
                if (if_all_satis) {
                    vec_data.erase(vec_data.begin() + i);
                    vec_segment_route.erase(vec_segment_route.begin() + i);
                    if_clear = true;
                    break;
                }
            }
            if (!if_clear) {
                cnt *= static_cast<int>(vec_data[i].size());
                ++i;
            }
        }

        if (cnt != 1) {
            for (auto &r: vec_segment_route) {
                for (auto &s: r) {
                    for (auto i = s.begin(); i != s.end();) {
                        if (mem_long[*i]) {
                            i = s.erase(i);
                        } else ++i;
                    }
                }
            }
        }

        if (mem_long.count() > static_cast<int>(dim * MAX_CUT_MEM_FACTOR)) {
            if_suc = false;
            return;
        }

        for (int i = 1; i < dim; ++i) {
            if (mem_long[i]) {
                mem.emplace(i);
            }
        }

        if (pricing_hard_level == PRICING_HARD_LEVEL::EASY) {
            if (cnt == 1) return;
            findMem(vec_data, vec_segment_route, mem);
            return;
        }

        cnt = 1;
        for (auto &i: vec_data) {
            cnt *= i.size();
            if (cnt >= FIND_MEM_USE_ENUMERATION_OR_MIP) break;
        }

        if (cnt == 1) return;

        if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
            std::vector<int> tmp;
            std::unordered_set<int> new_mem;
            int record_min = std::numeric_limits<int>::max();
            combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
            mem = new_mem;
        } else {
            getMemoryByMIP(rank1CutsDataShared, solver, vec_data, vec_segment_route, mem, if_suc);
        }
    }

    void findMem4R1C3(const Rank1CutsDataShared &rank1CutsDataShared,
                      const DataShared &sharedData,
                      const PRICING_HARD_LEVEL pricing_hard_level,
                      Solver &solver,
                      const cutLong &v_comb,
                      std::unordered_set<int> &mem,
                      bool &if_suc) {
        if_suc = true;
        int times;
        std::unordered_set<int> aux;
        std::vector<std::vector<std::vector<int> > > vec_data;
        std::vector<std::vector<std::unordered_set<int> > > vec_segment_route;
        cutLong mem_long = 0;
        const auto &sol = sharedData.getSol();
        const auto dim = rank1CutsDataShared.getDim();
        for (const auto &route: sol) {
            const auto &i = route.col_seq;
            std::vector<std::unordered_set<int> >
                    tmp_segment_route;
            times = 0;
            aux.clear();
            for (auto v: i) {
                if (v_comb[v]) {
                    ++times;
                    if (times == 2) {
                        tmp_segment_route.emplace_back(aux);
                        times = 1; //change to 1! Notice!
                        aux.clear();
                    }
                } else if (times) {
                    aux.emplace(v);
                }
            }
            if (!tmp_segment_route.empty()) {
                if (tmp_segment_route.size() % 2) {
                    for (int j = 0; j < tmp_segment_route.size(); j += 2) {
                        for (auto v: tmp_segment_route[j]) {
                            mem_long.set(v);
                        }
                    }
                } else {
                    std::vector<int> arr(tmp_segment_route.size());
                    int n = static_cast<int>(arr.size());
                    int r = static_cast<int>(std::floor((n + 1) / 2) + RANK1_TOLERANCE);
                    std::vector<std::vector<int> > plans;
                    std::vector<int> tmp(r);
                    std::iota(arr.begin(), arr.end(), 0);
                    combinationUtil(arr, tmp, plans, 0, n - 1, 0, r);
                    vec_data.emplace_back(plans);
                    vec_segment_route.emplace_back(tmp_segment_route);
                    auto bs = cutLong{}.set();
                    for (auto &p: plans) {
                        cutLong tmp_bs = 0;
                        for (auto j: p) {
                            for (auto v: tmp_segment_route[j]) {
                                tmp_bs.set(v);
                            }
                        }
                        bs &= tmp_bs;
                    }
                    if (bs.any()) {
                        for (int v = 0; v < dim; ++v) {
                            if (bs[v]) {
                                mem_long.set(v);
                            }
                        }
                    }
                }
            }
        }

        int cnt = 1;
        for (int i = 0; i < vec_data.size();) {
            bool if_clear = false;
            for (auto &j: vec_data[i]) {
                bool if_all_satis = true;
                for (auto k: j) {
                    for (auto l: vec_segment_route[i][k]) {
                        if (!mem_long[l]) {
                            if_all_satis = false;
                            goto outside;
                        }
                    }
                }
            outside:
                if (if_all_satis) {
                    vec_data.erase(vec_data.begin() + i);
                    vec_segment_route.erase(vec_segment_route.begin() + i);
                    if_clear = true;
                    break;
                }
            }
            if (!if_clear) {
                if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) cnt *= static_cast<int>(vec_data[i].size());
                ++i;
            }
        }

        if (cnt != 1) {
            for (auto &r: vec_segment_route) {
                for (auto &s: r) {
                    for (auto i = s.begin(); i != s.end();) {
                        if (mem_long[*i]) {
                            i = s.erase(i);
                        } else ++i;
                    }
                }
            }
        }

        if (mem_long.count() > static_cast<int>(dim * MAX_CUT_MEM_FACTOR)) {
            if_suc = false;
            return;
        }

        for (int i = 1; i < dim; ++i) {
            if (mem_long[i]) mem.emplace(i);
        }

        if (pricing_hard_level == PRICING_HARD_LEVEL::EASY) {
            if (cnt == 1) return;
            findMem(vec_data, vec_segment_route, mem);
            return;
        }

        if (cnt == 1) return;
        if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
            std::vector<int> tmp;
            std::unordered_set<int> new_mem;
            int record_min = std::numeric_limits<int>::max();
            combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
            mem = new_mem;
        } else {
            getMemoryByMIP(rank1CutsDataShared, solver, vec_data, vec_segment_route, mem, if_suc);
        }
    }

    void translateMemInt2MemPair(const std::unordered_set<int> &mem_set,
                                 std::vector<std::pair<int, int> > &mem_pair) {
        std::vector<int> mem(mem_set.begin(), mem_set.end());
        std::sort(mem.begin(), mem.end());
        //write mem_pair
        mem_pair.clear();
        for (int i = 0; i < static_cast<int>(mem.size()) - 1; ++i) {
            for (int j = i + 1; j < static_cast<int>(mem.size()); ++j) {
                mem_pair.emplace_back(mem[i], mem[j]);
            }
        }
    }

    void constructMemoryVertexBased(const Rank1CutsDataShared &rank1CutsDataShared, DataShared &sharedData,
                                    const PRICING_HARD_LEVEL pricing_hard_level, Solver &solver,
                                    std::unordered_map<std::vector<int>, std::vector<std::vector<int> >,
                                        VectorHashInRank1> &
                                    rank1_multi_mem_plan_map) {
        const auto dim = rank1CutsDataShared.getDim();
        auto &cuts = sharedData.refCuts();
        std::vector<int> tmp_fill(dim);
        std::iota(tmp_fill.begin(), tmp_fill.end(), 0);
        for (auto it = cuts.begin(); it != cuts.end();) {
            auto &c = *it;
            std::unordered_set<int> mem = {};
            if (c.idx_r1c != INITIAL_IDX_R1C) {
                auto &arc_mem = c.arc_mem;
                for (auto &m: arc_mem) {
                    mem.emplace(m.first);
                    mem.emplace(m.second);
                }
            }
            bool if_suc;
            auto &cut = c.info_r1c;
            if (cut.second == 0) {
                cutLong v_comb = 0;
                for (const auto &i: cut.first) v_comb.set(i);
                findMem4R1C3(rank1CutsDataShared, sharedData, pricing_hard_level, solver, v_comb, mem, if_suc);
            } else {
                findMem4R1CMulti(rank1CutsDataShared, sharedData, pricing_hard_level, solver, rank1_multi_mem_plan_map,
                                 cut, mem,
                                 if_suc);
            }
            if (if_suc) {
                //add c.info_r1c.first to mem;
                std::transform(cut.first.begin(), cut.first.end(), std::inserter(mem, mem.end()),
                               [](int i) { return i; });
                translateMemInt2MemPair(mem, c.arc_mem);
                ++it;
            } else {
                it = cuts.erase(it);
            }
        }
    }

    void findMem(const std::vector<std::vector<std::vector<int> > > &arr,
                 const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                 std::unordered_set<int> &mem) {
        for (int i = 0; i < arr.size(); ++i) {
            auto &r = arr[i];
            std::vector<int> ele_size(r.size(), 0);
            for (int j = 0; j < r.size(); ++j) {
                for (auto k: r[j]) {
                    ele_size[j] += static_cast<int>(vec_segment[i][k].size());
                }
            }
            auto min_idx = distance(ele_size.begin(), std::min_element(ele_size.begin(), ele_size.end()));
            for (auto k: r[min_idx]) {
                for (auto l: vec_segment[i][k]) {
                    mem.emplace(l);
                }
            }
        }
    }

    void getMemoryByMIP(const Rank1CutsDataShared &rank1CutsDataShared,
                        Solver &solver,
                        const std::vector<std::vector<std::vector<int> > > &arr,
                        const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                        std::unordered_set<int> &mem, bool &if_suc) {
        std::unordered_map<int, int> new_mem_map;
        std::unordered_map<int, int> re_new_mem_map;
        const auto dim = rank1CutsDataShared.getDim();
        new_mem_map.reserve(dim);

        int idx = 0;
        for (auto &r: vec_segment) {
            for (auto &s: r) {
                for (auto i: s) {
                    if (new_mem_map.find(i) == new_mem_map.end()) {
                        new_mem_map.emplace(i, idx);
                        re_new_mem_map.emplace(idx, i);
                        ++idx;
                    }
                }
            }
        }

        int num_r_p = 0;
        for (auto &r: arr) {
            num_r_p += static_cast<int>(r.size());
        }

        int half_num_row = static_cast<int>(new_mem_map.size());
        int last_half = static_cast<int>(arr.size());
        int local_num_row = half_num_row + last_half;
        std::vector<char> sense(local_num_row, SOLVER_LESS_EQUAL);
        std::fill(sense.begin() + half_num_row, sense.end(), SOLVER_EQUAL);
        std::vector<double> rhs(local_num_row, 0);
        std::fill(rhs.begin() + half_num_row, rhs.end(), 1);

        Solver local_solver{};
        local_solver.getEnv(&solver); //need load environment
        RANK1_VERBOSE_EXEC(std::cout << "we build model= getMemByMIP_.lp to std::get mem!" << std::endl;
        )
        SAFE_SOLVER(local_solver.newModel("getMemByMIP_.lp", 0, nullptr, nullptr, nullptr, nullptr, nullptr))
        SAFE_SOLVER(local_solver.addConstraints(local_num_row,
            0,
            nullptr,
            nullptr,
            nullptr,
            sense.data(),
            rhs.data(),
            nullptr))

        int last_num_idx = half_num_row;
        std::vector<size_t> solver_beg;
        std::vector<int> solver_ind;
        std::vector<double> solver_val, solver_obj;
        int numRe = 10000;
        solver_beg.reserve(numRe);
        solver_ind.reserve(numRe);
        solver_val.reserve(numRe);
        solver_obj.reserve(numRe);
        for (int i = 0; i < arr.size(); ++i, ++last_num_idx) {
            for (auto &p: arr[i]) {
                solver_beg.emplace_back(static_cast<int>(solver_ind.size()));
                cutLong tmp = 0;
                for (auto n: p) {
                    for (auto j: vec_segment[i][n]) {
                        if (tmp[j]) continue;
                        tmp.set(j);
                        solver_ind.emplace_back(new_mem_map[j]);
                    }
                }
                solver_ind.emplace_back(last_num_idx);
            }
        }
        int old_ccnt = static_cast<int>(solver_beg.size());
        int old_nzcnt = static_cast<int>(solver_ind.size());
        solver_obj.assign(old_ccnt, 0);
        solver_val.assign(old_nzcnt, 1);
        solver_obj.resize(old_ccnt + half_num_row, 1);
        solver_beg.resize(old_ccnt + half_num_row + 1);
        std::iota(solver_beg.begin() + old_ccnt, solver_beg.end(), old_nzcnt);
        old_ccnt += half_num_row;
        solver_ind.resize(old_nzcnt + half_num_row);
        std::iota(solver_ind.begin() + old_nzcnt, solver_ind.end(), 0);
        old_nzcnt += half_num_row;
        solver_val.resize(old_nzcnt, -num_r_p);
        std::vector<char> xtype(old_ccnt, SOLVER_BINARY);
        SAFE_SOLVER(local_solver.XaddVars(old_ccnt,
            old_nzcnt,
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            solver_obj.data(),
            nullptr,
            nullptr,
            xtype.data(),
            nullptr))
        SAFE_SOLVER(local_solver.setEnvTimeLimit(TIME_LIMIT_FOR_MIP_FIND_MEM))
        SAFE_SOLVER(local_solver.optimize())
        SAFE_SOLVER(local_solver.setEnvTimeLimit(std::numeric_limits<int>::max()))
        int status, left_ = old_ccnt - num_r_p;
        SAFE_SOLVER(local_solver.getStatus(&status))
        std::vector<double> X(left_);
        if (status == SOLVER_TIME_LIMIT) {
            RANK1_VERBOSE_EXEC(std::cout << "time limit for getMemoryByMIP" << std::endl;
            )
            if_suc = false;
            goto HERE;
        }
        SAFE_SOLVER(local_solver.getX(num_r_p, left_, X.data()))
        for (int i = 0; i < left_; ++i) {
            if (X[i] > 0.5) {
                mem.emplace(re_new_mem_map[i]);
            }
        }
    HERE:
        SAFE_SOLVER(local_solver.freeModel())
    }

    void findPlan4R1CMulti(
        const Rank1CutsDataShared &rank1CutsDataShared,
        std::unordered_map<std::vector<int>, std::vector<std::vector<int> >, VectorHashInRank1> &
        rank1_multi_mem_plan_map,
        const std::vector<int> &vis, const int denominator, cutLong &mem,
        std::vector<std::unordered_set<int> > &segment,
        std::vector<std::vector<int> > &plan) {
        int sum = std::accumulate(vis.begin(), vis.end(), 0);
        int mod = sum % denominator;
        std::vector<int> key = vis;
        key.emplace_back(mod);
        const auto dim = rank1CutsDataShared.getDim();
        auto &other2 = rank1_multi_mem_plan_map[key];
        if (other2.empty()) {
            std::list<other_> other; //we can not use
            other.emplace_back(0, mod, vis, std::vector<int>{});
            for (auto i = other.begin(); i != other.end(); ++i) {
                auto &o = *i;
                int cnt = 0;
                int tor = o.tor;
                int beg = o.beg;
                for (int j = 0; j < o.left_c.size(); ++j) {
                    cnt += o.left_c[j];
                    if (cnt >= denominator) {
                        cnt -= denominator;
                    }
                    if (cnt) {
                        if (cnt <= tor && o.left_c.begin() + j + 1 < o.left_c.end()) {
                            other.emplace_back(beg + j + 1,
                                               tor - cnt,
                                               std::vector<int>(o.left_c.begin() + j + 1, o.left_c.end()),
                                               o.mem_c);
                        }
                        int rem = beg + j;
                        if (rem != static_cast<int>(vis.size()) - 1)
                            o.mem_c.emplace_back(beg + j); //can not change the sequence!
                    }
                }
            }
            other2.resize(other.size());
            std::transform(other.begin(), other.end(), other2.begin(), [](const other_ &o) {
                return o.mem_c;
            });
        }

        for (int i = 1; i < dim; ++i) {
            if (mem[i]) {
                for (auto &j: segment) {
                    j.erase(i);
                }
            }
        }
        std::vector<std::unordered_set<int> >
                num_vis(dim);
        for (int i = 0; i < other2.size(); ++i) {
            for (auto j: other2[i]) {
                for (auto k: segment[j]) {
                    num_vis[k].emplace(i);
                }
            }
        }
        for (int i = 1; i < dim; ++i) {
            if (num_vis[i].size() == other2.size()) {
                mem.set(i);
                for (auto &j: segment) {
                    j.erase(i);
                }
            }
        }

        std::vector<std::pair<cutLong, std::vector<int> > > mem_other;
        for (const auto &i: other2) {
            cutLong p_mem = 0;
            for (auto j: i) {
                for (auto k: segment[j]) {
                    p_mem.set(k);
                }
            }
            if (p_mem.none()) {
                plan.clear();
                goto QUIT;
            }
            for (auto &j: mem_other) {
                if (((p_mem & j.first) ^ p_mem).none()) {
                    j = {p_mem, i};
                    goto NEXT;
                }
            }
            mem_other.emplace_back(p_mem, i);
        NEXT:;
        }
        plan.resize(mem_other.size());
        std::transform(mem_other.begin(), mem_other.end(), plan.begin(),
                       [](const auto &i) { return i.second; });
    QUIT:;
    }

    void combinationUtil(const std::vector<int> &arr,
                         std::vector<int> &tmp,
                         std::vector<std::vector<int> > &data,
                         int start,
                         int end,
                         int index,
                         int r) {
        if (index == r) {
            data.emplace_back(tmp);
            return;
        }

        for (int i = start; end - i >= 2 * (r - index - 1); i++) {
            tmp[index] = arr[i];
            combinationUtil(arr, tmp, data, i + 2,
                            end, index + 1, r);
        }
    }


    void combinations(const std::vector<std::vector<std::vector<int> > > &arr,
                      const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                      const int i,
                      const std::vector<int> &accum,
                      const std::unordered_set<int> &mem,
                      int &record_min,
                      std::unordered_set<int> &new_mem) {
        if (i == arr.size()) {
            int num = 0;
            auto tmp_mem = mem;
            for (int j = 0; j < arr.size(); ++j) {
                for (auto k: arr[j][accum[j]]) {
                    for (auto l: vec_segment[j][k]) {
                        if (tmp_mem.find(l) == tmp_mem.end()) {
                            tmp_mem.emplace(l);
                            num++;
                        }
                    }
                }
            }
            if (num < record_min) {
                record_min = num;
                new_mem = tmp_mem;
            }
        } else {
            for (int j = 0; j < arr[i].size(); ++j) {
                std::vector<int> tmp(accum);
                tmp.emplace_back(j);
                combinations(arr, vec_segment, i + 1, tmp, mem, record_min, new_mem);
            }
        }
    }
}
