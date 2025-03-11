/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <numeric>
#include <iostream>
#include "rank1_macro.hpp"
#include "rank1_cuts_generator.hpp"
#include "rank1_data_shared.hpp"


namespace RouteOpt::Rank1Cuts::Separation {
    void chooseCuts(std::vector<R1c> &tmp_cuts,
                    int numCuts,
                    const Rank1CutsDataShared &rank1CutsDataShared,
                    DataShared &sharedData,
                    std::vector<R1c> &chosen_cuts);

    void generatePermutations(std::unordered_map<int, int> &count_map,
                              std::vector<int> &result,
                              std::vector<std::vector<int> > &results,
                              int remaining);

    void CutGenerator::generateSepHeurMem4Vertex() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        const auto &sharedData = sharedData_ref.get();
        int dim = rank1CutsDataShared.getDim();
        const auto &cost_mat4_vertex = sharedData.getCostMat4Vertex();
        rank1_sep_heur_mem4_vertex.resize(dim);
        std::vector<std::pair<int, double> > cost(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 1; j < dim; ++j) {
                cost[j].first = j;
                cost[j].second = cost_mat4_vertex[i][j];
            }
            cost[0].first = 0;
            cost[0].second = std::numeric_limits<int>::max();
            std::stable_sort(cost.begin(), cost.end(),
                             [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
                                 return a.second < b.second;
                             });
            cutLong &vst2 = rank1_sep_heur_mem4_vertex[i];
            int max_heuristic_sep_mem4_row_rank1 = std::min(dim - 1, MAX_HEURISTIC_SEP_ROW_RANK1);
            for (int k = 0; k < max_heuristic_sep_mem4_row_rank1; ++k) {
                vst2.set(cost[k].first);
            }
        }
    }

    void CutGenerator::generateRecordMapRank1Combinations() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        std::unordered_map<int, std::bitset<8> > support;
        for (int i = 1; i <= max_heuristic_initial_seed_set_size_row_rank1c; ++i) {
            auto &sup = support[i];
            for (int j = 0; j < 7; ++j) {
                if (rank1CutsDataShared.getDenominator(i, j) == 0) continue;
                sup.set(j);
            }
        }

        for (int i = 1; i <= max_heuristic_initial_seed_set_size_row_rank1c; ++i) {
            for (int j = 1; j <= max_heuristic_initial_seed_set_size_row_rank1c; ++j) {
                if (i == j)continue;
                if (support[i].count() < support[j].count()) continue;
                if (i > j && map_rank1_multiplier_dominance[j].find(i) != map_rank1_multiplier_dominance[j].end())
                    continue;
                if (((support[i] & support[j]) ^ support[j]).none()) {
                    map_rank1_multiplier_dominance[i].emplace(j);
                }
            }
        }
    }

    void CutGenerator::generateMapRank1MultiplierDominance() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        record_map_rank1_combinations.resize(max_heuristic_initial_seed_set_size_row_rank1c + 1,
                                             std::vector<std::vector<std::vector<int> > >(
                                                 7, std::vector<std::vector<int> >()));
        for (int i = 1; i <= max_heuristic_initial_seed_set_size_row_rank1c; ++i) {
            for (int j = 0; j < 7; ++j) {
                if (rank1CutsDataShared.getDenominator(i, j) == 0) continue;
                std::unordered_map<int, int> count_map;
                for (const auto &it: rank1CutsDataShared.getMultiplier(i, j)) {
                    count_map[it]++;
                }
                std::vector<int> result;
                std::vector<std::vector<int> > results;
                generatePermutations(count_map, result, results, i);
                record_map_rank1_combinations[i][j] = std::move(results);
            }
        }
    }


    void CutGenerator::generateR1C3() {
        /**
         * we only search for cuts that could potentially make the cuts have a small memory
         * for the routes, we can construct a residual graph
         * the cuts std::set {i,j,k}. for i, (i,j) or (i,k) must exist in the graph
         */
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        const auto &sol = sharedData.getSol();
        const auto dim = rank1CutsDataShared.getDim();
        const auto all_num_routes = static_cast<int>(sol.size());
        std::vector<std::unordered_map<int, int> >
                visited_times_by_routes(dim);
        std::vector<std::unordered_set<int> >
                i_connections(dim); //j,k that connects i
        int routeNumber = 0;
        for (const auto &[frac_x, col_seq, forward_concatenate_pos]: sol) {
            const auto &route = col_seq;
            for (int j = 0; j < static_cast<int>(route.size()) - 1; ++j) {
                i_connections[route[j]].emplace(route[j + 1]);
                i_connections[route[j + 1]].emplace(route[j]);
                ++visited_times_by_routes[route[j]][routeNumber];
            }
            ++visited_times_by_routes[route.back()][routeNumber];
            ++routeNumber;
        }
        std::vector<std::vector<int> > i_connections_vec(dim);
        std::transform(i_connections.begin(), i_connections.end(), i_connections_vec.begin(),
                       [](const std::unordered_set<int> &ele) -> std::vector<int> {
                           return {ele.begin(), ele.end()};
                       });
        for (auto &i: i_connections_vec) {
            std::sort(i.begin(), i.end());
        }
        std::vector<int> v_tmp;
        std::vector<int> v_tmp2(all_num_routes); //std::find the number of times visited the cut
        std::vector<std::pair<std::vector<int>, double> >
                cuts_pool;
        cuts_pool.reserve(2048);
        for (int i = 1; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                v_tmp.clear();
                if (i_connections[i].find(j) != i_connections[i].end()) {
                    set_union(i_connections_vec[i].begin(), i_connections_vec[i].end(),
                              i_connections_vec[j].begin(), i_connections_vec[j].end(),
                              std::back_inserter(v_tmp));
                } else {
                    set_intersection(i_connections_vec[i].begin(), i_connections_vec[i].end(),
                                     i_connections_vec[j].begin(), i_connections_vec[j].end(),
                                     std::back_inserter(v_tmp));
                }
                for (int k: v_tmp) {
                    if (k <= i || k <= j) continue;
                    double vio = -1;
                    memset(v_tmp2.data(), 0, sizeof(int) * all_num_routes);
                    for (auto &[fst, snd]: visited_times_by_routes[i]) {
                        v_tmp2[fst] += snd;
                    }
                    for (auto &[fst, snd]: visited_times_by_routes[j]) {
                        v_tmp2[fst] += snd;
                    }
                    for (auto &[fst, snd]: visited_times_by_routes[k]) {
                        v_tmp2[fst] += snd;
                    }
                    for (int l = 0; l < all_num_routes; ++l) {
                        vio += static_cast<int>(v_tmp2[l] / 2.) * sol[l].frac_x;
                    }
                    if (vio > RANK1_TOLERANCE) {
                        cuts_pool.emplace_back(std::vector<int>{i, j, k}, vio);
                    }
                }
            }
        }

        std::sort(cuts_pool.begin(), cuts_pool.end(),
                  [](const std::pair<std::vector<int>, double> &a,
                     const std::pair<std::vector<int>, double> &b) -> bool {
                      return a.second > b.second;
                  });
        std::vector<R1c> tmp_cut(cuts_pool.size());
        std::transform(cuts_pool.begin(), cuts_pool.end(), tmp_cut.begin(),
                       [](const std::pair<std::vector<int>, double> &ele) -> R1c {
                           return R1c{std::make_pair(ele.first, 0)};
                       });
        auto &cuts = sharedData.refCuts();
        chooseCuts(tmp_cut, max_num_r1c3_per_round, rank1CutsDataShared, sharedData, cuts);
    }

    void CutGenerator::generateR1C1() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        const auto &sol = sharedData.getSol();
        std::vector<std::pair<double, R1c> > tmp_cuts;
        std::unordered_map<int, int> vis_map;
        for (const auto &r: sol) {
            vis_map.clear();
            for (const auto i: r.col_seq) {
                ++vis_map[i];
            }
            for (auto &[v, times]: vis_map) {
                if (times > 1) {
                    tmp_cuts.emplace_back();
                    tmp_cuts.back().first = std::floor(times / 2. + RANK1_TOLERANCE) * r.frac_x;
                    tmp_cuts.back().second.info_r1c = std::make_pair(std::vector{v}, 0);
                }
            }
        }
        if (tmp_cuts.empty()) return;

        std::sort(tmp_cuts.begin(), tmp_cuts.end(),
                  [](const std::pair<double, R1c> &a, const std::pair<double, R1c> &b) {
                      return a.first > b.first;
                  });
        std::vector<R1c> pure_cuts(tmp_cuts.size());
        std::transform(tmp_cuts.begin(), tmp_cuts.end(), pure_cuts.begin(), [](const std::pair<double, R1c> &a) {
            return a.second;
        });
        auto &cuts = sharedData.refCuts();
        chooseCuts(pure_cuts, max_num_r1c_per_round, rank1CutsDataShared, sharedData, cuts);
    }

    void CutGenerator::getHighDimCuts() {
        constructVRMapNSeed();

        startSeed();

        for (int i = 0; i < num_label;) {
            operationsControl(rank1_multi_label_pool[i], i);
        }
        constructCuts();
    }

    void CutGenerator::startSeed() {
        rank1_multi_label_pool.resize(INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE);
        num_label = 0;
        for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
            for (auto &p: c_N_noC) {
                auto &c = p.first;
                auto &wc = p.second;
                double vio, best_vio;
                exactFindBestPermutationForOnePlan(c, plan_idx, vio);

                if (vio < RANK1_TOLERANCE) continue;
                best_vio = vio;
                char best_oper = 'o';

                int add_j, remove_j;
                std::pair<int, int> swap_i_j;
                addSearch(plan_idx, c, wc, vio, add_j);

                if (vio > best_vio) {
                    best_vio = vio;
                    best_oper = 'a';
                }

                removeSearch(plan_idx, c, vio, remove_j);

                if (vio > best_vio) {
                    best_vio = vio;
                    best_oper = 'r';
                }

                swapSearch(plan_idx, c, wc, vio, swap_i_j);

                if (vio > best_vio) {
                    best_vio = vio;
                    best_oper = 's';
                }
                cutLong tmp;
                std::vector<int> new_c, new_w_no_c;
                switch (best_oper) {
                    case 'o': if (c.size() < 2 || c.size() > max_row_rank1) break;
                        tmp = 0;
                        for (auto &i: c) {
                            tmp.set(i);
                        }
                        generated_rank1_multi_pool[static_cast<int>(c.size())].
                                emplace_back(tmp, plan_idx, best_vio);
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
                    rank1_multi_label_pool[num_label++] = Rank1MultiLabel{
                        new_c, new_w_no_c, plan_idx, best_vio, best_oper
                    };
                    if (num_label == rank1_multi_label_pool.size()) {
                        rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
                    }
                }
            }
        }
    }

    void CutGenerator::exactFindBestPermutationForOnePlan(
        std::vector<int> &cut,
        int plan_idx,
        double &vio) {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        const auto &sharedData = sharedData_ref.get();
        const int cut_size = static_cast<int>(cut.size());
        if (rank1CutsDataShared.getDenominator(cut_size, plan_idx) == 0) {
            vio = -std::numeric_limits<double>::max();
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

        const int denominator = rank1CutsDataShared.getDenominator(cut_size, plan_idx);
        const double rhs = rank1CutsDataShared.getRhs(cut_size, plan_idx);
        const auto &coeffs = record_map_rank1_combinations[cut_size][plan_idx];
        const auto &v_r_map = sharedData.getVRMap();

        std::unordered_map<int, int> map_r_numbers;
        for (const auto &i: cut) {
            for (const auto &[fst, snd]: v_r_map[i]) {
                map_r_numbers[fst] += snd;
            }
        }
        const auto &sol = sharedData.getSol();
        std::vector<double> valid_routes;
        valid_routes.reserve(sol.size());

        std::unordered_map<int, int> map_old_new_routes;
        map_old_new_routes.reserve(sol.size());

        for (auto &pr: map_r_numbers) {
            if (pr.second > 1) {
                valid_routes.emplace_back(sol[pr.first].frac_x);
                map_old_new_routes[pr.first] = static_cast<int>(valid_routes.size()) - 1;
            }
        }

        std::vector<std::vector<int> > cut_num_times_vis_routes(cut_size);

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

        std::vector<double> num_times_vis_routes(map_old_new_routes.size());

        double best_vio = -std::numeric_limits<double>::max();
        int best_idx = -1;
        int cnt = 0;
        for (auto &cof: coeffs) {
            memset(num_times_vis_routes.data(), 0, sizeof(double) * num_times_vis_routes.size());
            for (int i = 0; i < cut_size; ++i) {
                for (auto &j: cut_num_times_vis_routes[i]) {
                    num_times_vis_routes[j] += cof[i];
                }
            }
            std::transform(num_times_vis_routes.begin(),
                           num_times_vis_routes.end(),
                           valid_routes.begin(),
                           num_times_vis_routes.begin(),
                           [denominator](const double a, const double b) {
                               return static_cast<int>(a / denominator + RANK1_TOLERANCE) * b;
                           });
            const double vio_tmp = std::accumulate(num_times_vis_routes.begin(),
                                                   num_times_vis_routes.end(),
                                                   -rhs);
            if (vio_tmp > best_vio) {
                best_vio = vio_tmp;
                best_idx = cnt;
            }
            ++cnt;
        }


        vio = best_vio;
        std::vector<std::pair<int, int> > cut_coeff(cut_size);
        for (int i = 0; i < cut_size; ++i) {
            cut_coeff[i] = std::make_pair(cut[i], coeffs[best_idx][i]);
        }
        std::sort(cut_coeff.begin(), cut_coeff.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
            return a.second > b.second;
        });
        std::vector<int> new_cut(cut_size);
        std::transform(cut_coeff.begin(), cut_coeff.end(), new_cut.begin(), [](const std::pair<int, int> &a) {
            return a.first;
        });
        map_cut_plan_vio[tmp][plan_idx] = std::make_pair(std::move(new_cut), vio);
    }


    void CutGenerator::constructVRMapNSeed() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        const auto &sol = sharedData.getSol();
        const auto dim = rank1CutsDataShared.getDim();
        std::vector<std::unordered_map<int, int> > v_r_map(dim, std::unordered_map<int, int>());
        for (auto &i: v_r_map) i.reserve(sol.size());
        for (int r = 0; r < sol.size(); ++r) {
            const auto &route = sol[r].col_seq;
            for (const auto &i: route) {
                ++v_r_map[i][r];
            }
        }
        std::unordered_map<cutLong, cutLong> seed_map; // c ij, and w_no_c
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
                    c_size < 3 || c_size > max_heuristic_initial_seed_set_size_row_rank1c
                )
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
        sharedData.setVRMap(v_r_map);
    }

    void CutGenerator::addSearch(int plan_idx,
                                 const std::vector<int> &c,
                                 const std::vector<int> &w_no_c,
                                 double &new_vio,
                                 int &add_j) {
        const auto new_c_size = static_cast<int>(c.size()) + 1;
        if (
            new_c_size > max_row_rank1 || (rank1CutsDataShared_ref.get().getDenominator(new_c_size, plan_idx) == 0)
        ) {
            new_vio = -std::numeric_limits<double>::max();
            return;
        }

        std::vector<int> tmp_c = c;
        tmp_c.emplace_back();
        double vio, best_vio = -std::numeric_limits<double>::max();
        for (const auto &cplus: w_no_c) {
            tmp_c.back() = cplus;
            exactFindBestPermutationForOnePlan(tmp_c, plan_idx, vio);
            if (vio > best_vio) {
                best_vio = vio;
                add_j = cplus;
            }
        }
        new_vio = best_vio - RANK1_TOLERANCE; // penalty
    }


    void CutGenerator::removeSearch(int plan_idx,
                                    const std::vector<int> &c,
                                    double &new_vio,
                                    int &remove_j

    ) {
        const auto new_c_size = static_cast<int>(c.size()) - 1;
        if (new_c_size < 3 || (rank1CutsDataShared_ref.get().getDenominator(new_c_size, plan_idx) == 0)) {
            new_vio = -std::numeric_limits<double>::max();
            return;
        }

        std::vector<int> tmp_c(static_cast<int>(c.size()) - 1);
        double vio, best_vio = -std::numeric_limits<double>::max();
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
        new_vio = best_vio + RANK1_TOLERANCE; // good penalty
    }

    void CutGenerator::swapSearch(int plan_idx,
                                  const std::vector<int> &c,
                                  const std::vector<int> &w_no_c,
                                  double &new_vio,
                                  std::pair<int, int> &swap_i_j) {
        int new_c_size = static_cast<int>(c.size());
        if ((new_c_size < 3 || new_c_size > max_row_rank1) || (
                rank1CutsDataShared_ref.get().getDenominator(new_c_size, plan_idx) == 0)) {
            new_vio = -std::numeric_limits<double>::max();
            return;
        }

        std::vector<int> tmp_c = c;
        double vio, best_vio = -std::numeric_limits<double>::max();
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


    void CutGenerator::operationsControl(
        Rank1MultiLabel &label,
        int &i) {
        auto &vio = label.vio;
        auto &new_cij = label.c;
        auto &w_no_cij = label.w_no_c;
        auto &plan_idx = label.plan_idx;

        auto dir = label.search_dir;
        int add_j, remove_j;
        std::pair<int, int> swap_i_j;
        std::vector<std::pair<int, double> > move_vio(4);
        for (int j = 1; j < 4; ++j) {
            move_vio[j] = {j, -std::numeric_limits<double>::max()};
        }
        move_vio[0] = {0, vio};
        double new_vio;
        if (dir == 'a' || dir == 's') {
            addSearch(plan_idx, new_cij, w_no_cij, new_vio, add_j);
            move_vio[1].second = new_vio;
        }
        if (dir == 'r' || dir == 's') {
            removeSearch(plan_idx, new_cij, new_vio, remove_j);
            move_vio[2].second = new_vio;
        }
        if (dir == 'a' || dir == 'r') {
            swapSearch(plan_idx, new_cij, w_no_cij, new_vio, swap_i_j);
            move_vio[3].second = new_vio;
        }
        std::sort(move_vio.begin(), move_vio.end(), [](const std::pair<int, double> &a,
                                                       const std::pair<int, double> &b) {
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
                generated_rank1_multi_pool[static_cast<int>(new_cij.size())].emplace_back(
                    tmp, plan_idx, best_move_vio);
                ++i;
                break;
            case 1: new_cij.emplace_back(add_j);
                w_no_cij.erase(std::find(w_no_cij.begin(), w_no_cij.end(), add_j));
                break;
            case 2: new_cij.erase(std::find(new_cij.begin(), new_cij.end(), remove_j));
                break;
            case 3: *std::find(new_cij.begin(), new_cij.end(), swap_i_j.first) = swap_i_j.second;
                w_no_cij.erase(std::find(w_no_cij.begin(), w_no_cij.end(), swap_i_j.second));
                break;
            default: THROW_RUNTIME_ERROR("best move error");
        }
        vio = best_move_vio;
    }


    void CutGenerator::constructCuts() {
        const auto &rank1CutsDataShared = rank1CutsDataShared_ref.get();
        auto &sharedData = sharedData_ref.get();
        for (auto &i: generated_rank1_multi_pool) {
            std::sort(i.second.begin(), i.second.end(),
                      [](auto &a, auto &b) {
                          return std::get<2>(a) > std::get<2>(b);
                      });
        }

        for (auto &i: generated_rank1_multi_pool) {
            std::unordered_set<cutLong> cut_set;
            std::unordered_set<int> p_set;
            const double vio_std = std::get<2>(i.second[0]) * CUT_VIO_FACTOR;
            std::vector<R1c> tmp_cuts;
            for (auto &j: i.second) {
                if (std::get<2>(j) < vio_std) break;
                auto &key = std::get<0>(j);
                if (cut_set.find(key) != cut_set.end() && p_set.find(std::get<1>(j)) != p_set.end()) continue;
                tmp_cuts.emplace_back(R1c{
                    std::make_pair(map_cut_plan_vio.at(key)[std::get<1>(j)].first, std::get<1>(j))
                });
                cut_set.insert(key);
                p_set.insert(std::get<1>(j));
            }
            auto &cuts = sharedData.refCuts();
            chooseCuts(tmp_cuts, max_num_r1c_per_round, rank1CutsDataShared, sharedData, cuts);
        }
    }

    inline void sortCuts(const Rank1CutsDataShared &rank1CutsDataShared,
                         std::vector<int> &fst,
                         int snd) {
        auto &multi = rank1CutsDataShared.getMultiplier(static_cast<int>(fst.size()), snd);
        std::vector<std::pair<int, int> > pr(fst.size());
        std::transform(fst.begin(), fst.end(), multi.begin(), pr.begin(),
                       [](int a, int b) {
                           return std::make_pair(a, b);
                       });
        std::sort(pr.begin(), pr.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
            return (a.second > b.second) || (a.second == b.second && a.first > b.first);
        });
        std::transform(pr.begin(), pr.end(), fst.begin(),
                       [](const std::pair<int, int> &p) {
                           return p.first;
                       });
    }

    void chooseCuts(std::vector<R1c> &tmp_cuts,
                    int numCuts,
                    const Rank1CutsDataShared &rank1CutsDataShared,
                    DataShared &sharedData,
                    std::vector<R1c> &chosen_cuts) {
        /**
         * this function select valid cuts from all cuts
         */
        numCuts = std::min(numCuts, static_cast<int>(tmp_cuts.size()));
        if (numCuts == 0) return;

        auto &cut_record = sharedData.refCutRecord();
        for (auto &cut: tmp_cuts) {
            auto &fst = cut.info_r1c.first;
            auto &snd = cut.info_r1c.second;
            sortCuts(rank1CutsDataShared, fst, snd);
            if (cut_record[fst].find(snd) != cut_record[fst].end()) continue;
            cut_record[fst].insert(snd);
            chosen_cuts.emplace_back();
            chosen_cuts.back().info_r1c = std::make_pair(fst, snd);
            numCuts--;
            if (numCuts == 0) break;
        }
    }

    void generatePermutations(std::unordered_map<int, int> &count_map,
                              std::vector<int> &result,
                              std::vector<std::vector<int> > &results,
                              const int remaining) {
        if (remaining == 0) {
            results.emplace_back(result);
            return;
        }

        for (auto &pr_: count_map) {
            if (pr_.second > 0) {
                result.emplace_back(pr_.first);
                pr_.second--;
                generatePermutations(count_map, result, results, remaining - 1);
                pr_.second++; // go back
                result.pop_back(); // erase and go next
            }
        }
    }
} // namespace RouteOpt::Rank1Cuts::Separation
