//
// Created by Ricky You on 10/12/24.
//


#include "rank1_cuts_separator.hpp"

#include <bitset>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <list>


using namespace std;
using namespace Eigen;

struct State {
    int coeff{};
    int state{};
    int end_segment{};
    bitset<MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN> bit{};
};

struct Arcs {
    vector<unordered_set<pair<int, int>, PairHasher> > arc_plan;
};

struct other_ {
    int beg{};
    int tor{};
    vector<int> left_c{};
    vector<int> mem_c{};

    other_(int beg, int tor, vector<int> left_c, vector<int> mem_c) : beg(beg), tor(tor), left_c(std::move(left_c)),
                                                                      mem_c(std::move(mem_c)) {
    }

    other_() = default;
};

void reduceArcs(std::vector<Arcs> &all_arcs, std::unordered_set<std::pair<int, int>, PairHasher> &existing_arcs,
                bool &if_suc);

void getLeastMemory(std::vector<Arcs> &all_arcs, std::unordered_set<std::pair<int, int>, PairHasher> &existing_arcs,
                    bool &if_suc);

void findLeastPlans2MakeCoeffRight(const std::vector<int> &vertex_states,
                                   int denominator,
                                   std::vector<std::vector<int> > &plans,
                                   bool &if_suc);

void getVertexStates(const vector<int> &sequence,
                     int forward_pos,
                     const std::unordered_map<int, int> &map,
                     std::vector<int> &vertex_states,
                     std::vector<std::unordered_set<std::pair<int, int>, PairHasher> > &arcs);


void Rank1CutsSeparator::constructMemoryVertexBased() {
    vector<int> tmp_fill(dim);
    iota(tmp_fill.begin(), tmp_fill.end(), 0);
    for (auto &c: cuts) {
        unordered_set<int> mem = {};
        if (c.idx_r1c != INITIAL_IDX_R1C) {
            auto &arc_mem = c.arc_mem;
            for (auto &m: arc_mem) {
                mem.emplace(m.second);
            }
        }
        bool if_suc;
        auto &cut = c.info_r1c;
        if (cut.second == 0) {
            cutLong v_comb = 0;
            for (const auto &i: cut.first) v_comb.set(i);
            findMemoryForR1CsMode3ByEnumerationAndMIP(v_comb, mem, if_suc);
        } else {
            findMemoryForRank1Multi(cut, mem, if_suc);
        }
        if (if_suc) {
            auto &arc_mem = c.arc_mem;
            arc_mem.assign(mem.size(), make_pair(tmp_fill, 0));
            int cnt = 0;
            for (auto &m: mem) {
                arc_mem[cnt++].second = m;
            }
        }
    }
}


void Rank1CutsSeparator::findMemoryForR1CsMode3ByEnumerationAndMIP(const cutLong &v_comb,
                                                                   std::unordered_set<int> &mem,
                                                                   bool &if_suc) {
    if_suc = true;
    int times;
    unordered_set<int> aux;
    vector<vector<vector<int> > > vec_data;
    vector<vector<unordered_set<int> > > vec_segment_route;
    cutLong mem_long = 0;
    for (const auto &route: sol) {
        const auto &i = route.col_seq;
        vector<unordered_set<int> >
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
                vector<int> arr(tmp_segment_route.size());
                int n = static_cast<int>(arr.size());
                int r = static_cast<int>(floor((n + 1) / 2) + TOLERANCE);
                vector<vector<int> > plans;
                vector<int> tmp(r);
                iota(arr.begin(), arr.end(), 0);
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

    if (mem_long.count() > static_cast<int>(dim * max_cut_mem_factor)) {
        if_suc = false;
        return;
    }

    for (int i = 1; i < dim; ++i) {
        if (mem_long[i]) mem.emplace(i);
    }

    if (pricing_hard_level == 0) {
        if (cnt == 1) return;
        findMemAggressively(vec_data, vec_segment_route, mem);
        return;
    }

    if (cnt == 1) {
        return;
    } else if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
        vector<int> tmp;
        unordered_set<int> new_mem;
        int record_min = numeric_limits<int>::max();
        combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
        mem = new_mem;
    } else {
        getMemoryByMIP(vec_data, vec_segment_route, mem, if_suc);
    }
}

void Rank1CutsSeparator::findMemAggressively(const std::vector<std::vector<std::vector<int> > > &array,
                                             const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                                             std::unordered_set<int> &mem) {
    for (int i = 0; i < array.size(); ++i) {
        auto &r = array[i];
        vector<int> ele_size(r.size(), 0);
        for (int j = 0; j < r.size(); ++j) {
            for (auto k: r[j]) {
                ele_size[j] += static_cast<int>(vec_segment[i][k].size());
            }
        }
        auto min_idx = distance(ele_size.begin(), min_element(ele_size.begin(), ele_size.end()));
        for (auto k: r[min_idx]) {
            for (auto l: vec_segment[i][k]) {
                mem.emplace(l);
            }
        }
    }
}

void Rank1CutsSeparator::combinationUtil(const std::vector<int> &arr,
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


void Rank1CutsSeparator::combinations(const std::vector<std::vector<std::vector<int> > > &array,
                                      const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                                      const int i,
                                      const std::vector<int> &accum,
                                      const std::unordered_set<int> &mem,
                                      int &record_min,
                                      std::unordered_set<int> &new_mem) {
    if (i == array.size()) {
        int num = 0;
        auto tmp_mem = mem;
        for (int j = 0; j < array.size(); ++j) {
            for (auto k: array[j][accum[j]]) {
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
        for (int j = 0; j < array[i].size(); ++j) {
            vector<int> tmp(accum);
            tmp.emplace_back(j);
            combinations(array, vec_segment, i + 1, tmp, mem, record_min, new_mem);
        }
    }
}

void Rank1CutsSeparator::getMemoryByMIP(const std::vector<std::vector<std::vector<int> > > &array,
                                        const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                                        std::unordered_set<int> &mem, bool &if_suc) {
    unordered_map<int, int> new_mem_map;
    unordered_map<int, int> re_new_mem_map;
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
    for (auto &r: array) {
        num_r_p += static_cast<int>(r.size());
    }

    int half_num_row = static_cast<int>(new_mem_map.size());
    int last_half = static_cast<int>(array.size());
    int local_num_row = half_num_row + last_half;
    vector<char> sense(local_num_row, SOLVER_LESS_EQUAL);
    fill(sense.begin() + half_num_row, sense.end(), SOLVER_EQUAL);
    vector<double> rhs(local_num_row, 0);
    fill(rhs.begin() + half_num_row, rhs.end(), 1);

    Solver local_solver{};
    local_solver.getEnv(&solver); //need load environment
    cout << "we build model= getMemByMIP_.lp to get mem!" << endl;
    safe_solver(local_solver.newModel("getMemByMIP_.lp", 0, nullptr, nullptr, nullptr, nullptr, nullptr))
    safe_solver(local_solver.addConstraints(local_num_row,
        0,
        nullptr,
        nullptr,
        nullptr,
        sense.data(),
        rhs.data(),
        nullptr))

    int last_num_idx = half_num_row;
    vector<size_t> solver_beg;
    vector<int> solver_ind;
    vector<double> solver_val, solver_obj;
    int numRe = 10000;
    solver_beg.reserve(numRe);
    solver_ind.reserve(numRe);
    solver_val.reserve(numRe);
    solver_obj.reserve(numRe);
    for (int i = 0; i < array.size(); ++i, ++last_num_idx) {
        for (auto &p: array[i]) {
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
    iota(solver_beg.begin() + old_ccnt, solver_beg.end(), old_nzcnt);
    old_ccnt += half_num_row;
    solver_ind.resize(old_nzcnt + half_num_row);
    iota(solver_ind.begin() + old_nzcnt, solver_ind.end(), 0);
    old_nzcnt += half_num_row;
    solver_val.resize(old_nzcnt, -num_r_p);
    vector<char> xtype(old_ccnt, SOLVER_BINARY);
    safe_solver(local_solver.XaddVars(old_ccnt,
        old_nzcnt,
        solver_beg.data(),
        solver_ind.data(),
        solver_val.data(),
        solver_obj.data(),
        nullptr,
        nullptr,
        xtype.data(),
        nullptr))
    safe_solver(local_solver.setEnvTimeLimit(TIME_LIMIT_FOR_MIP_FIND_MEM))
    safe_solver(local_solver.optimize())
    safe_solver(local_solver.setEnvTimeLimit(MAX_TIME_LIMIT_FOR_MIP))
    int status, left = old_ccnt - num_r_p;
    safe_solver(local_solver.getStatus(&status))
    vector<double> X(left);
    if (status == SOLVER_TIME_LIMIT) {
        cout << "time limit for getMemoryByMIP" << endl;
        if_suc = false;
        goto HERE;
    }
    safe_solver(local_solver.getX(num_r_p, left, X.data()))
    for (int i = 0; i < left; ++i) {
        if (X[i] > 0.5) {
            mem.emplace(re_new_mem_map[i]);
        }
    }
HERE:
    safe_solver(local_solver.freeModel())
}


void Rank1CutsSeparator::findMemoryForRank1Multi(
    const std::pair<std::vector<int>, int> &cut_pair,
    std::unordered_set<int> &mem,
    bool &if_suc) {
    if_suc = true;
    auto &cut = cut_pair.first;
    auto plan_idx = cut_pair.second;
    int size = static_cast<int>(cut.size());
    const auto &multi = get<0>(map_rank1_multiplier[size][plan_idx]);
    auto denominator = get<1>(map_rank1_multiplier[size][plan_idx]);
    vector<vector<vector<int> > > vec_data;
    vector<vector<unordered_set<int> > > vec_segment_route;
    unordered_map<int, int> map_cut_mul;
    vector<int> num_vis_times(sol.size(), 0);
    for (int i = 0; i < cut.size(); ++i) {
        map_cut_mul[cut[i]] = multi[i];
        for (auto &pr: v_r_map[cut[i]]) {
            num_vis_times[pr.first] += multi[i] * pr.second;
        }
    }
    transform(num_vis_times.begin(), num_vis_times.end(), num_vis_times.begin(),
              [denominator](const int x) { return static_cast<int>(x / denominator); });
    cutLong mem_long = 0;
    int num = 0;
    for (auto &route: sol) {
        auto &i = route.col_seq;
        if (num_vis_times[num++] == 0) continue;
        vector<vector<int> > data;
        vector<int> vis;
        vector<unordered_set<int> > segment_route;
        unordered_set<int> tmp_seg;
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
            segment_route.erase(segment_route.begin()); //remove the first one
        findPlanForRank1Multi(vis, denominator, mem_long, segment_route, data);
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

    if (mem_long.count() > static_cast<int>(dim * max_cut_mem_factor)) {
        if_suc = false;
        return;
    }

    for (int i = 1; i < dim; ++i) {
        if (mem_long[i]) {
            mem.emplace(i);
        }
    }

    if (pricing_hard_level == 0) {
        if (cnt == 1) return;
        findMemAggressively(vec_data, vec_segment_route, mem);
        return;
    }

    cnt = 1;
    for (auto &i: vec_data) {
        cnt *= i.size();
        if (cnt >= FIND_MEM_USE_ENUMERATION_OR_MIP) break;
    }

    if (cnt == 1) {
        return;
    } else if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
        vector<int> tmp;
        unordered_set<int> new_mem;
        int record_min = numeric_limits<int>::max();
        combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
        mem = new_mem;
    } else {
        getMemoryByMIP(vec_data, vec_segment_route, mem, if_suc);
    }
}

void Rank1CutsSeparator::findPlanForRank1Multi(const std::vector<int> &vis, const int denominator, cutLong &mem,
                                               std::vector<std::unordered_set<int> > &segment,
                                               std::vector<std::vector<int> > &plan) {
    int sum = accumulate(vis.begin(), vis.end(), 0);
    int mod = sum % denominator;
    vector<int> key = vis;
    key.emplace_back(mod);
    auto &other2 = rank1_multi_mem_plan_map[key];
    if (other2.empty()) {
        list<other_> other; //we can not use
        other.emplace_back(0, mod, vis, vector<int>{});
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
                                           vector<int>(o.left_c.begin() + j + 1, o.left_c.end()),
                                           o.mem_c);
                    }
                    int rem = beg + j;
                    if (rem != static_cast<int>(vis.size()) - 1)
                        o.mem_c.emplace_back(beg + j); //can not change the sequence!
                }
            }
        }
        other2.resize(other.size());
        transform(other.begin(), other.end(), other2.begin(), [](const other_ &o) {
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
    vector<unordered_set<int> >
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

    vector<pair<cutLong, vector<int> > > mem_other;
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
    transform(mem_other.begin(), mem_other.end(), plan.begin(),
              [](const auto &i) { return i.second; });
QUIT:;
}


void Rank1CutsSeparator::constructMemoryArcBased() {
    vector<Eigen::Triplet<int> > tripletList;
    tripletList.reserve(static_cast<size_t>(static_cast<double>(dim * sol.size()) * 0.05));
    unordered_map<int, int> map_vertex_index;
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
    safe_eigen(sol_matrix.setFromTriplets(tripletList.begin(), tripletList.end());)

    for (auto it = cuts.begin(); it != cuts.end();) {
        auto &cut = *it;
        bool if_suc;
        findLeastMemoryArcBased(sol_matrix, cut, if_suc);
        if (!if_suc) {
            it = cuts.erase(it);
        } else ++it;
    }
}


void Rank1CutsSeparator::findLeastMemoryArcBased(const sparseRowMatrixXI &sol_matrix, R1c &cut, bool &if_suc) {
    if_suc = true;
    auto &arc_mem = cut.arc_mem;
    unordered_set<pair<int, int>, PairHasher> existing_arcs;
    if (!arc_mem.empty()) {
        for (auto &info: arc_mem) {
            int end = info.second;
            for (int j: info.first) {
                existing_arcs.emplace(j, end);
            }
        }
    }

    unordered_map<int, int> map_vertex_state;
    auto &multiplier = get<0>(map_rank1_multiplier[static_cast<int>(cut.info_r1c.first.size())][cut.info_r1c.second]);
    for (int i = 0; i < multiplier.size(); ++i) {
        map_vertex_state[cut.info_r1c.first[i]] = multiplier[i];
    }
    int denominator = get<1>(map_rank1_multiplier[static_cast<int>(cut.info_r1c.first.size())][cut.info_r1c.second]);
    sparseRowVectorXI vec(sol_matrix.cols());
    for (auto i: cut.info_r1c.first) {
        vec += sol_matrix.row(i);
    }
    vec /= denominator;
    vec.prune(0);

    vector<Arcs> all_arcs;
    for (sparseRowVectorXI::InnerIterator it(vec, 0); it; ++it) {
        int col = static_cast<int>(it.col());
        vector<vector<int> > ps;
        vector<int> vertex_states;
        vector<unordered_set<pair<int, int>, PairHasher> > arcs;
        getVertexStates(sol[col].col_seq,
                        sol[col].forward_concatenate_pos,
                        map_vertex_state,
                        vertex_states,
                        arcs);
        findLeastPlans2MakeCoeffRight(vertex_states, denominator, ps, if_suc);
        if (!if_suc) cout << "WARNING: findLeastPlans2MakeCoeffRight failed!" << endl;
        else {
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
        unordered_map<int, vector<int> > map_vertex_arcs;
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

void reduceArcs(vector<Arcs> &all_arcs, unordered_set<pair<int, int>, PairHasher> &existing_arcs, bool &if_suc) {
    if_suc = true;
    unordered_map<pair<int, int>, unordered_set<int>, PairHasher> arc_map;
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
            } else ++it2;
        }
        if (if_delete_all) {
            it = all_arcs.erase(it);
        } else ++it;
    }
    unordered_map<pair<int, int>, int, PairHasher> arc_bit_map;
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
        vector<ARCBIT> arc_bit(arc_plan.size(), 0);
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


void getLeastMemory(vector<Arcs> &all_arcs, unordered_set<pair<int, int>, PairHasher> &existing_arcs, bool &if_suc) {
    bool loop = true;
    while (loop) {
        int old_size = static_cast<int>(existing_arcs.size());
        reduceArcs(all_arcs, existing_arcs, if_suc);
        if (old_size == existing_arcs.size()) loop = false;
        if (!if_suc) return;
    }

    if (all_arcs.empty()) return;
    vector<vector<ARCBIT> > bins(all_arcs.size());
    bins[0].resize(all_arcs[0].arc_plan.size());
    for (int i = 1; i < all_arcs.size(); ++i) {
        int tmp = 2 * bins[i - 1].size() < 100 ? static_cast<int>(2 * bins[i - 1].size()) : 100;
        bins[i].resize(tmp);
    }
    vector<int> bin_num(all_arcs.size(), 0);
    unordered_map<pair<int, int>, int, PairHasher> arc_bit_map;
    unordered_map<int, pair<int, int> > bit_arc_map;
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
    vector<vector<ARCBIT> > all_arcs_bit(all_arcs.size());
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
            auto &arc_bit = bin[j];
            for (int k = 0; k < add_arcs.size(); ++k) {
                ARCBIT new_arc_bit = arc_bit | all_arcs_bit[i + 1][k];
                bool if_keep = true;
                for (int l = 0; l < bin_num_next;) {
                    if ((new_arc_bit & bin_next[l]) == bin_next[l]) {
                        if_keep = false;
                        break;
                    } else if ((new_arc_bit & bin_next[l]) == new_arc_bit) {
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
    ARCBIT best = least_memory[0];
    for (int i = 1; i < least_memory_num; ++i) {
        auto &arc_bit = least_memory[i];
        if (best.count() > arc_bit.count()) best = arc_bit;
    }
    for (int i = 0; i < bit_pos; ++i) {
        if (best.test(i)) {
            if (existing_arcs.find(bit_arc_map[i]) != existing_arcs.end())
                throw runtime_error("existing_arcs.find(bit_arc_map[i]) != existing_arcs.end()");
            existing_arcs.emplace(bit_arc_map[i]);
        }
    }
}

void getVertexStates(const vector<int> &sequence,
                     const int forward_pos,
                     const unordered_map<int, int> &map,
                     vector<int> &vertex_states,
                     vector<unordered_set<pair<int, int>, PairHasher> > &arcs) {
    vertex_states.clear();
    arcs.clear();
    vector<vector<int> > tmp_arcs;
    vector<int> tmp;
    pair<int, int> change_direction{-1, -1};
    for (int i = 0; i < sequence.size(); ++i) {
        int current = sequence[i];
        if (map.find(current) != map.end()) {
            vertex_states.emplace_back(map.at(current));
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
     * and we set initial state () to be 0 and we do extend the state where it must be larger than sparse_rep.back()
     * the extension ends when the coeff reaches the max coeff,
     * dominance rule 1(extension): used segments are subset; coeff is same; state is the same
     * dominance rule 2(finish): used segments are subset;
     */
    vector<int> segments(vertex_states.size() - 1);
    for (int i = 0; i < segments.size(); ++i) segments[i] = vertex_states[i + 1] + vertex_states[i];
    int max_coeff = accumulate(vertex_states.begin(), vertex_states.end(), 0) / denominator;
    vector<vector<pair<vector<State>, int> > > states(max_coeff + 1, vector<pair<vector<State>, int> >(denominator));
    for (int i = 0; i < states.size(); ++i) {
        for (int j = 0; j < states[i].size(); ++j) {
            states[i][j].first.resize(static_cast<int>(pow(2, i)) < 100 ? static_cast<int>(pow(2, i)) : 100);
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
    if (plans.empty()) throw runtime_error("plans.empty()");
}


void Rank1CutsSeparator::selectR1CsByVioNMemory() {
    unordered_map<int, vector<pair<int, int> > > map_v_cut; // key: v; val: idx & multiplier
    map_v_cut.reserve(cuts.size());
    vector<int> rhs(cuts.size());
    for (int i = 0; i < cuts.size(); ++i) {
        auto &c = cuts[i].info_r1c;
        auto &multipliers = get<0>(map_rank1_multiplier.at(static_cast<int>(c.first.size()))[c.second]);
        for (int j = 0; j < c.first.size(); ++j) {
            map_v_cut[c.first[j]].emplace_back(i, multipliers[j]);
        }
        rhs[i] = get<2>(map_rank1_multiplier.at(static_cast<int>(c.first.size()))[c.second]);
    }
    vector<Eigen::RowVectorXd> cuts_coeffs(cuts.size(), Eigen::RowVectorXd::Zero(static_cast<int>(sol.size())));
    for (int i = 0; i < sol.size(); ++i) {
        for (const auto &j: sol[i].col_seq) {
            for (auto &k: map_v_cut[j]) {
                cuts_coeffs[k.first][i] += k.second;
            }
        }
    }

    for (int i = 0; i < cuts.size(); ++i) {
        auto &c = cuts[i].info_r1c;
        const auto denominator = get<1>(map_rank1_multiplier[static_cast<int>(c.first.size())][c.second]);
        for (auto &j: cuts_coeffs[i]) j = static_cast<int>(j / denominator + TOLERANCE);
    }

    vector<pair<double, int> > cuts_vio(cuts.size());
    for (int i = 0; i < cuts.size(); ++i)cuts_vio[i] = {0, i};
    Eigen::RowVectorXd frac_routes_vio(sol.size());
    for (int i = 0; i < sol.size(); ++i)frac_routes_vio[i] = sol[i].frac_x;
    for (int i = 0; i < cuts.size(); ++i)cuts_vio[i].first = cuts_coeffs[i].dot(frac_routes_vio) - rhs[i];

    sort(cuts_vio.begin(), cuts_vio.end(), [](const pair<double, int> &a, const pair<double, int> &b) {
        return a.first > b.first;
    });

    int left_cuts = MAX_NUM_R1CS_IN_PRICING - static_cast<int>(old_cuts.size()) - 10;
    if (left_cuts < 0) {
        cuts.clear();
        return;
    }

    cuts_vio.resize(min(left_cuts, static_cast<int>(cuts_vio.size())));

    vector<int> vertex_related_r1c(dim, MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX - 1);

    for (auto &r1c: old_cuts) {
        for (auto &it: r1c.info_r1c.first) {
            --vertex_related_r1c[it];
        }
        for (auto &it: r1c.arc_mem) {
            --vertex_related_r1c[it.second];
        }
    }

    vector<R1c> select_cut(cuts_vio.size());
    int num = 0;
    for (auto &i: cuts_vio) {
        int idx = i.second;
        auto &c = cuts[idx];
        bool if_keep = true;
        auto tmp = vertex_related_r1c;
        for (auto j: c.info_r1c.first) {
            if (tmp[j] <= 0) {
                if_keep = false;
                break;
            }
            --tmp[j];
        }
        for (auto &j: c.arc_mem) {
            if (tmp[j.second] <= 0) {
                if_keep = false;
                break;
            }
            --tmp[j.second];
        }
        if (!if_keep) continue;
        vertex_related_r1c = tmp;
        select_cut[num++] = c;
    }
    select_cut.resize(num);
    cuts = select_cut;
}

void Rank1CutsSeparator::getSeparatedCuts(std::vector<R1c> &cuts) {
    cuts = Rank1CutsSeparator::cuts;
}
