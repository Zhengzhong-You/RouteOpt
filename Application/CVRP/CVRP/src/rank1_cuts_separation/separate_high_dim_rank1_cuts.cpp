//
// Created by Ricky You on 10/12/24.
//

#include "rank1_cuts_separator.hpp"

#include <bitset>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <stdexcept>


using namespace std;

void Rank1CutsSeparator::initialSupportVector() {
    cut_record.clear();
    v_r_map.clear();
    c_N_noC.clear();
    map_cut_plan_vio.clear();
    map_cut_plan_vio.reserve(4096);
    generated_rank1_multi_pool.clear();
    rank1_multi_label_pool.clear();
    num_label = 0;
    rank1_multi_mem_plan_map.clear();
    cuts.clear();
}


void Rank1CutsSeparator::getHighDimCuts() {
    constructVRMapAndSeedCrazy();

    startSeedCrazy();

    for (int i = 0; i < num_label;) {
        operationsCrazy(rank1_multi_label_pool[i], i);
    }
    constructCutsCrazy();
}

void Rank1CutsSeparator::startSeedCrazy() {
    rank1_multi_label_pool.resize(INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE);
    num_label = 0;
    for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
        for (auto &p: c_N_noC) {
            auto &c = p.first;
            auto &wc = p.second;
            double vio, best_vio;
            exactFindBestPermutationForOnePlan(c, plan_idx, vio);

            if (vio < TOLERANCE) continue;
            best_vio = vio;
            char best_oper = 'o';

            int add_j, remove_j;
            pair<int, int> swap_i_j;
            addSearchCrazy(plan_idx, c, wc, vio, add_j);

            if (vio > best_vio) {
                best_vio = vio;
                best_oper = 'a';
            }

            removeSearchCrazy(plan_idx, c, vio, remove_j);

            if (vio > best_vio) {
                best_vio = vio;
                best_oper = 'r';
            }

            swapSearchCrazy(plan_idx, c, wc, vio, swap_i_j);

            if (vio > best_vio) {
                best_vio = vio;
                best_oper = 's';
            }
            cutLong tmp;
            vector<int> new_c, new_w_no_c;
            switch (best_oper) {
                case 'o': if (c.size() < 2 || c.size() > max_row_rank1) break;
                    tmp = 0;
                    for (auto &i: c) {
                        tmp.set(i);
                    }
                    generated_rank1_multi_pool[static_cast<int>(c.size())].emplace_back(tmp, plan_idx, best_vio);
                    break;
                case 'a': new_c = c;
                    new_c.emplace_back(add_j);
                    new_w_no_c.resize(static_cast<int>(wc.size()) - 1);
                    for (int i = 0, j = 0; i < wc.size(); ++i) {
                        if (wc[i] != add_j) {
                            new_w_no_c[j++] = wc[i];
                        }
                    }
                    break;
                case 'r': new_c.resize(static_cast<int>(c.size()) - 1);
                    for (int i = 0, j = 0; i < c.size(); ++i) {
                        if (c[i] != remove_j) {
                            new_c[j++] = c[i];
                        }
                    }
                    new_w_no_c = wc;
                    break;
                case 's': new_c.resize(c.size());
                    new_c.resize(c.size());
                    for (int i = 0; i < c.size(); ++i) {
                        if (c[i] != swap_i_j.first) {
                            new_c[i] = c[i];
                        } else {
                            new_c[i] = swap_i_j.second;
                        }
                    }
                    new_w_no_c.resize(wc.size());
                    for (int i = 0; i < wc.size(); ++i) {
                        if (wc[i] != swap_i_j.second) {
                            new_w_no_c[i] = wc[i];
                        } else {
                            new_w_no_c[i] = swap_i_j.first;
                        }
                    }
                    break;
                default: ;
            }
            if (best_oper != 'o') {
                rank1_multi_label_pool[num_label++] = Rank1MultiLabel{new_c, new_w_no_c, plan_idx, best_vio, best_oper};
                if (num_label == rank1_multi_label_pool.size()) {
                    rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
                }
            }
        }
    }
}

void Rank1CutsSeparator::exactFindBestPermutationForOnePlan(std::vector<int> &cut, const int plan_idx, double &vio) {
    const int cut_size = static_cast<int>(cut.size());
    const auto &plan = map_rank1_multiplier[cut_size][plan_idx];
    if (!get<1>(plan)) {
        vio = -numeric_limits<double>::max();
        return;
    }
    cutLong tmp = 0;
    for (const auto &it: cut) {
        tmp.set(it);
    }
    if (map_cut_plan_vio.find(tmp) != map_cut_plan_vio.end()) {
        auto &vec = map_cut_plan_vio[tmp][plan_idx];
        if (!vec.first.empty()) {
            vio = vec.second;
            return;
        }
    } else {
        map_cut_plan_vio[tmp].resize(7);
    }

    const int denominator = get<1>(plan);
    const double rhs = get<2>(plan); //must be double
    const auto &coeffs = record_map_rank1_combinations[cut_size][plan_idx];

    unordered_map<int, int> map_r_numbers;
    for (const auto &i: cut) {
        for (const auto &[fst, snd]: v_r_map[i]) {
            map_r_numbers[fst] += snd;
        }
    }

    vector<double> valid_routes;
    valid_routes.reserve(sol.size());

    unordered_map<int, int> map_old_new_routes;
    map_old_new_routes.reserve(sol.size());

    for (auto &pr: map_r_numbers) {
        if (pr.second > 1) {
            valid_routes.emplace_back(sol[pr.first].frac_x);
            map_old_new_routes[pr.first] = static_cast<int>(valid_routes.size()) - 1;
        }
    }

    vector<vector<int> > cut_num_times_vis_routes(cut_size);

    for (int i = 0; i < cut_size; ++i) {
        int c = cut[i];
        for (auto &pr: v_r_map[c]) {
            if (map_old_new_routes.find(pr.first) != map_old_new_routes.end()) {
                for (int j = 0; j < pr.second; ++j) {
                    cut_num_times_vis_routes[i].emplace_back(map_old_new_routes[pr.first]);
                }
            }
        }
    }

    vector<double> num_times_vis_routes(map_old_new_routes.size());

    double best_vio = -numeric_limits<double>::max();
    int best_idx = -1;
    int cnt = 0;
    for (auto &cof: coeffs) {
        memset(num_times_vis_routes.data(), 0, sizeof(double) * num_times_vis_routes.size());
        for (int i = 0; i < cut_size; ++i) {
            for (auto &j: cut_num_times_vis_routes[i]) {
                num_times_vis_routes[j] += cof[i];
            }
        }
        transform(num_times_vis_routes.begin(),
                  num_times_vis_routes.end(),
                  valid_routes.begin(),
                  num_times_vis_routes.begin(),
                  [denominator](const double a, const double b) {
                      return static_cast<int>(a / denominator + tolerance) * b;
                  });
        const double vio_tmp = accumulate(num_times_vis_routes.begin(),
                                          num_times_vis_routes.end(),
                                          -rhs);
        if (vio_tmp > best_vio) {
            best_vio = vio_tmp;
            best_idx = cnt;
        }
        ++cnt;
    }


    vio = best_vio;
    vector<pair<int, int> > cut_coeff(cut_size);
    for (int i = 0; i < cut_size; ++i) {
        cut_coeff[i] = make_pair(cut[i], coeffs[best_idx][i]);
    }
    sort(cut_coeff.begin(), cut_coeff.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
        return a.second > b.second;
    });
    vector<int> new_cut(cut_size);
    transform(cut_coeff.begin(), cut_coeff.end(), new_cut.begin(), [](const pair<int, int> &a) {
        return a.first;
    });
    map_cut_plan_vio[tmp][plan_idx] = make_pair(std::move(new_cut), vio);
}


void Rank1CutsSeparator::constructVRMapAndSeedCrazy() {
    v_r_map.assign(dim, {});
    for (auto &i: v_r_map) i.reserve(sol.size());
    for (int r = 0; r < sol.size(); ++r) {
        const auto &route = sol[r].col_seq;
        for (const auto &i: route) {
            ++v_r_map[i][r];
        }
    }
    unordered_map<cutLong, cutLong> seed_map; // c ij, and w_no_c
    seed_map.reserve(4096);
    for (int i = 1; i < dim; ++i) {
        if (v_r_map[i].empty()) continue;
        cutLong wc = 0; // c within i
        for (auto &[fst, snd]: v_r_map[i]) {
            for (const auto &v: sol[fst].col_seq) {
                if (rank1_sep_heur_mem4_vertex[i][v]) {
                    wc.set(v);
                }
            }
        }

        for (int r = 0; r < sol.size(); ++r) {
            if (v_r_map[i].find(r) != v_r_map[i].end()) continue;
            cutLong tmp_c = 0;
            for (const auto &v: sol[r].col_seq) {
                if (wc[v]) {
                    tmp_c.set(v); //r no i but c within i
                }
            }
            tmp_c.set(i);
            if (int c_size = static_cast<int>(tmp_c.count());
                c_size < 3 || c_size > max_heuristic_initial_seed_set_size_row_rank1c)
                continue;
            if (seed_map.find(tmp_c) == seed_map.end()) {
                seed_map[tmp_c] = wc ^ tmp_c;
            } else {
                seed_map[tmp_c] |= wc ^ tmp_c; //union
            }
        }
    }
    c_N_noC.assign(seed_map.size(), {});
    int cnt = 0;
    for (auto &[fst, snd]: seed_map) {
        int cnt2 = 0, cnt3 = 0;
        auto &[tmp_fst, tmp_snd] = c_N_noC[cnt];
        tmp_fst.resize(fst.count());
        tmp_snd.resize(snd.count());
        for (int i = 1; i < dim; ++i) {
            if (fst.test(i)) {
                tmp_fst[cnt2++] = i;
            } else if (snd.test(i)) {
                tmp_snd[cnt3++] = i;
            }
        }
        ++cnt;
    }
}

void Rank1CutsSeparator::addSearchCrazy(int plan_idx,
                                        const std::vector<int> &c,
                                        const std::vector<int> &w_no_c,
                                        double &new_vio,
                                        int &add_j) {
    const auto new_c_size = static_cast<int>(c.size()) + 1;
    if (const auto &plan = map_rank1_multiplier[new_c_size][plan_idx]; new_c_size > max_row_rank1 || !get<1>(plan)) {
        new_vio = -numeric_limits<double>::max();
        return;
    }

    vector<int> tmp_c = c;
    tmp_c.emplace_back();
    double vio, best_vio = -numeric_limits<double>::max();
    for (const auto &cplus: w_no_c) {
        tmp_c.back() = cplus;
        exactFindBestPermutationForOnePlan(tmp_c, plan_idx, vio);
        if (vio > best_vio) {
            best_vio = vio;
            add_j = cplus;
        }
    }
    new_vio = best_vio - tolerance; // penalty
}


void Rank1CutsSeparator::removeSearchCrazy(int plan_idx,
                                           const std::vector<int> &c,
                                           double &new_vio,
                                           int &remove_j

) {
    const auto new_c_size = static_cast<int>(c.size()) - 1;
    const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
    if (new_c_size < 3 || !get<1>(plan)) {
        new_vio = -numeric_limits<double>::max();
        return;
    }

    vector<int> tmp_c(static_cast<int>(c.size()) - 1);
    double vio, best_vio = -numeric_limits<double>::max();
    for (int i = 0; i < c.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            tmp_c[j] = c[j];
        }
        for (int j = i + 1; j < c.size(); ++j) {
            tmp_c[j - 1] = c[j];
        }
        exactFindBestPermutationForOnePlan(tmp_c, plan_idx, vio);
        if (vio > best_vio) {
            best_vio = vio;
            remove_j = c[i];
        }
    }
    new_vio = best_vio + tolerance; // good penalty
}

void Rank1CutsSeparator::swapSearchCrazy(int plan_idx,
                                         const std::vector<int> &c,
                                         const std::vector<int> &w_no_c,
                                         double &new_vio,
                                         std::pair<int, int> &swap_i_j) {
    const auto &plan = map_rank1_multiplier[static_cast<int>(c.size())][plan_idx];
    int new_c_size = static_cast<int>(c.size());
    if ((new_c_size < 3 || new_c_size > max_row_rank1) || !get<1>(plan)) {
        new_vio = -numeric_limits<double>::max();
        return;
    }

    vector<int> tmp_c = c;
    double vio, best_vio = -numeric_limits<double>::max();
    for (int i = 0; i < c.size(); ++i) {
        for (int j: w_no_c) {
            tmp_c[i] = j;
            exactFindBestPermutationForOnePlan(tmp_c, plan_idx, vio);
            if (vio > best_vio) {
                best_vio = vio;
                swap_i_j = {c[i], j};
            }
        }
        tmp_c[i] = c[i];
    }
    new_vio = best_vio;
}


void Rank1CutsSeparator::operationsCrazy(
    Rank1MultiLabel &label,
    int &i) {
    auto &vio = label.vio;
    auto &new_cij = label.c;
    auto &w_no_cij = label.w_no_c;
    auto &plan_idx = label.plan_idx;

    auto dir = label.search_dir;
    int add_j, remove_j;
    pair<int, int> swap_i_j;
    vector<pair<int, double> > move_vio(4);
    for (int j = 1; j < 4; ++j) {
        move_vio[j] = {j, -numeric_limits<double>::max()};
    }
    move_vio[0] = {0, vio};
    double new_vio;
    if (dir == 'a' || dir == 's') {
        addSearchCrazy(plan_idx, new_cij, w_no_cij, new_vio, add_j);
        move_vio[1].second = new_vio;
    }
    if (dir == 'r' || dir == 's') {
        removeSearchCrazy(plan_idx, new_cij, new_vio, remove_j);
        move_vio[2].second = new_vio;
    }
    if (dir == 'a' || dir == 'r') {
        swapSearchCrazy(plan_idx, new_cij, w_no_cij, new_vio, swap_i_j);
        move_vio[3].second = new_vio;
    }
    sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                              const pair<int, double> &b) {
        return a.second > b.second;
    });
    int best_move = move_vio[0].first;
    double best_move_vio = move_vio[0].second;

    cutLong tmp;
    switch (best_move) {
        case 0: tmp = 0;
            for (auto j: new_cij) {
                tmp.set(j);
            }
            generated_rank1_multi_pool[static_cast<int>(new_cij.size())].emplace_back(tmp, plan_idx, best_move_vio);
            ++i;
            break;
        case 1: new_cij.emplace_back(add_j);
            w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), add_j));
            break;
        case 2: new_cij.erase(find(new_cij.begin(), new_cij.end(), remove_j));
            break;
        case 3: *find(new_cij.begin(), new_cij.end(), swap_i_j.first) = swap_i_j.second;
            w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), swap_i_j.second));
            break;
        default: throw std::runtime_error("best move error");
    }
    vio = best_move_vio;
}


void Rank1CutsSeparator::constructCutsCrazy() {
    for (auto &i: generated_rank1_multi_pool) {
        sort(i.second.begin(), i.second.end(),
             [](auto &a, auto &b) {
                 return get<2>(a) > get<2>(b);
             });
    }

    for (auto &i: generated_rank1_multi_pool) {
        unordered_set<cutLong> cut_set;
        unordered_set<int> p_set;
        const double vio_std = get<2>(i.second[0]) * cut_vio_factor;
        vector<R1c> tmp_cuts;
        for (auto &j: i.second) {
            if (get<2>(j) < vio_std) break;
            auto &key = get<0>(j);
            if (cut_set.find(key) != cut_set.end() && p_set.find(get<1>(j)) != p_set.end()) continue;
            tmp_cuts.emplace_back(R1c{make_pair(map_cut_plan_vio.at(key)[get<1>(j)].first, get<1>(j))});
            cut_set.insert(key);
            p_set.insert(get<1>(j));
        }
        chooseCuts(tmp_cuts, cuts, max_num_r1c_per_round);
    }
}
