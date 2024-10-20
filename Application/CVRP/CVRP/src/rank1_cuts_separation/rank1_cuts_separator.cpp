//
// Created by Ricky You on 10/11/24.
//

#include "rank1_cuts_separator.hpp"

#include <bitset>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <iostream>
#include <stdexcept>

int Rank1CutsSeparator::max_row_rank1{};
double Rank1CutsSeparator::cut_vio_factor{};
int Rank1CutsSeparator::max_heuristic_sep_mem4_row_rank1{};
int Rank1CutsSeparator::max_heuristic_initial_seed_set_size_row_rank1c{};
int Rank1CutsSeparator::max_num_r1c3_per_round{};
int Rank1CutsSeparator::max_num_r1c_per_round{};
double Rank1CutsSeparator::max_cut_mem_factor{};
const double Rank1CutsSeparator::tolerance{1e-6};
std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, int> > > Rank1CutsSeparator::map_rank1_multiplier
        {};
std::unordered_map<int, std::unordered_set<int> > Rank1CutsSeparator::map_rank1_multiplier_dominance{};
int Rank1CutsSeparator::dim{};
int Rank1CutsSeparator::pricing_hard_level{};
Solver Rank1CutsSeparator::solver{};
std::vector<cutLong> Rank1CutsSeparator::rank1_sep_heur_mem4_vertex{};
std::vector<std::vector<double> > Rank1CutsSeparator::cost_mat4_vertex{};
std::vector<std::vector<std::vector<std::vector<int> > > > Rank1CutsSeparator::record_map_rank1_combinations{};


int Rank1CutsSeparator::limited_memory_type{NODE_MEMORY};
std::vector<RouteInfo> Rank1CutsSeparator::sol{};
std::vector<R1c> Rank1CutsSeparator::cuts{};
std::vector<R1c> Rank1CutsSeparator::old_cuts{};
std::unordered_map<std::vector<int>, std::unordered_set<int>, VectorHashInRank1> Rank1CutsSeparator::cut_record{};
std::vector<std::unordered_map<int, int> > Rank1CutsSeparator::v_r_map{};
std::vector<std::pair<std::vector<int>, std::vector<int> > > Rank1CutsSeparator::c_N_noC{};
std::unordered_map<cutLong, std::vector<std::pair<std::vector<int>, double> > > Rank1CutsSeparator::map_cut_plan_vio{};
std::unordered_map<int, std::vector<std::tuple<cutLong, int, double> > > Rank1CutsSeparator::generated_rank1_multi_pool
        {};
std::vector<Rank1MultiLabel> Rank1CutsSeparator::rank1_multi_label_pool{};
int Rank1CutsSeparator::num_label{};
std::unordered_map<std::vector<int>, std::vector<std::vector<int> >, VectorHashInRank1>
Rank1CutsSeparator::rank1_multi_mem_plan_map{};


using namespace std;

void Rank1CutsSeparator::setInitialInfo(
    const int max_row_rank1,
    const int max_num_r1c3_per_round,
    const int max_num_r1c_per_round,
    const int dim,
    Solver &solver,
    const std::vector<std::vector<double> > &cost_mat4_vertex) {
    Rank1CutsSeparator::max_row_rank1 = max_row_rank1;
    Rank1CutsSeparator::dim = dim;
    Rank1CutsSeparator::cost_mat4_vertex = cost_mat4_vertex;

    Rank1CutsSeparator::max_num_r1c3_per_round = max_num_r1c3_per_round;
    Rank1CutsSeparator::max_num_r1c_per_round = max_num_r1c_per_round;

    Rank1CutsSeparator::cut_vio_factor = 0.1;
    Rank1CutsSeparator::max_heuristic_sep_mem4_row_rank1 = min(dim - 1, 16);
    Rank1CutsSeparator::max_heuristic_initial_seed_set_size_row_rank1c = Rank1CutsSeparator::max_row_rank1 + 1;
    Rank1CutsSeparator::max_cut_mem_factor = 0.15;
    Rank1CutsSeparator::solver.getEnv(&solver);

    generateOptimalMultiplier();
    generateSepHeurMem4Vertex();
}

void Rank1CutsSeparator::updateInfo(const int limited_memory_type,
                                    const int pricing_hard_level,
                                    const std::vector<RouteInfo> &sol,
                                    const std::vector<R1c> &old_cuts) {
    initialSupportVector();
    Rank1CutsSeparator::pricing_hard_level = pricing_hard_level;
    Rank1CutsSeparator::limited_memory_type = limited_memory_type;
    Rank1CutsSeparator::sol = sol;
    Rank1CutsSeparator::old_cuts = old_cuts;
    std::sort(Rank1CutsSeparator::sol.begin(), Rank1CutsSeparator::sol.end(),
              [](const RouteInfo &a, const RouteInfo &b) {
                  return a.frac_x > b.frac_x;
              });
}



void Rank1CutsSeparator::separateRank1Cuts() {
    fillMemory();
    generateR1C1();
    generateR1C3();
    getHighDimCuts();
    if (limited_memory_type == NODE_MEMORY) {
        constructMemoryVertexBased();
    } else if (limited_memory_type == ARC_MEMORY) {
        constructMemoryArcBased();
    }
}

void Rank1CutsSeparator::findMemory4Cuts(std::vector<R1c> &cuts, bool is_node_memory) {
    initialSupportVector();
    Rank1CutsSeparator::cuts = cuts;
    is_node_memory ? constructMemoryVertexBased() : constructMemoryArcBased();
    cuts = Rank1CutsSeparator::cuts;
}


void Rank1CutsSeparator::generateR1C3() {
    /**
     * we only search for cuts that could potentially make the cuts have a small memory
     * for the routes, we can construct a residual graph
     * the cuts set {i,j,k}. for i, (i,j) or (i,k) must exist in the graph
     */
    const auto all_num_routes = static_cast<int>(sol.size());
    vector<unordered_map<int, int> >
            visited_times_by_routes(dim);
    vector<unordered_set<int> >
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
    vector<vector<int> > i_connections_vec(dim);
    transform(i_connections.begin(), i_connections.end(), i_connections_vec.begin(),
              [](const unordered_set<int> &ele) -> vector<int> {
                  return {ele.begin(), ele.end()};
              });
    for (auto &i: i_connections_vec) {
        sort(i.begin(), i.end());
    }
    vector<int> v_tmp;
    vector<int> v_tmp2(all_num_routes); //find the number of times visited the cut
    vector<pair<vector<int>, double> >
            cuts_pool;
    cuts_pool.reserve(2048);
    for (int i = 1; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            v_tmp.clear();
            if (i_connections[i].find(j) != i_connections[i].end()) {
                set_union(i_connections_vec[i].begin(), i_connections_vec[i].end(),
                          i_connections_vec[j].begin(), i_connections_vec[j].end(),
                          back_inserter(v_tmp));
            } else {
                set_intersection(i_connections_vec[i].begin(), i_connections_vec[i].end(),
                                 i_connections_vec[j].begin(), i_connections_vec[j].end(),
                                 back_inserter(v_tmp));
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
                if (vio > tolerance) {
                    cuts_pool.emplace_back(vector<int>{i, j, k}, vio);
                }
            }
        }
    }

    sort(cuts_pool.begin(), cuts_pool.end(),
         [](const pair<vector<int>, double> &a, const pair<vector<int>, double> &b) -> bool {
             return a.second > b.second;
         });
    vector<R1c> tmp_cut(cuts_pool.size());
    transform(cuts_pool.begin(), cuts_pool.end(), tmp_cut.begin(),
              [](const pair<vector<int>, double> &ele) -> R1c {
                  return R1c{make_pair(ele.first, 0)};
              });
    chooseCuts(tmp_cut, cuts, max_num_r1c3_per_round);
}

void Rank1CutsSeparator::generateR1C1() {
    vector<std::pair<double, R1c> > tmp_cuts;

    unordered_map<int, int> vis_map;
    for (const auto &r: sol) {
        vis_map.clear();
        for (const auto i: r.col_seq) {
            ++vis_map[i];
        }
        for (auto &[v, times]: vis_map) {
            if (times > 1) {
                tmp_cuts.emplace_back();
                tmp_cuts.back().first = floor(times / 2. + tolerance) * r.frac_x;
                tmp_cuts.back().second.info_r1c = make_pair(vector{v}, 0);
            }
        }
    }
    if (tmp_cuts.empty()) return;

    sort(tmp_cuts.begin(), tmp_cuts.end(), [](const pair<double, R1c> &a, const pair<double, R1c> &b) {
        return a.first > b.first;
    });
    vector<R1c> pure_cuts(tmp_cuts.size());
    transform(tmp_cuts.begin(), tmp_cuts.end(), pure_cuts.begin(), [](const pair<double, R1c> &a) {
        return a.second;
    });
    chooseCuts(pure_cuts, cuts, max_num_r1c_per_round);
}

void Rank1CutsSeparator::chooseCuts(const std::vector<R1c> &tmp_cuts,
                                    std::vector<R1c> &chosen_cuts,
                                    int numCuts) {
    /**
     * this function select valid cuts from all cuts
     */
    numCuts = min(numCuts, static_cast<int>(tmp_cuts.size()));
    if (numCuts == 0) return;

    for (auto &cut: tmp_cuts) {
        auto &fst = cut.info_r1c.first;
        auto &snd = cut.info_r1c.second;
        if (cut_record[fst].find(snd) != cut_record[fst].end()) continue;
        cut_record[fst].insert(snd);
        /**
         * make the cuts exist one and only one (very important for avoiding repeated cuts)
         */
        int size = static_cast<int>(fst.size());
        const auto &coeff = get<0>(map_rank1_multiplier[size][snd]);
        vector<vector<int> > tmp_cut(coeff[0] + 1);
        for (int i = 0; i < size; ++i) {
            tmp_cut[coeff[i]].emplace_back(fst[i]);
        }
        for (int i = 0; i <= coeff[0]; ++i) {
            if (tmp_cut[i].size() > 1)
                sort(tmp_cut[i].begin(), tmp_cut[i].end());
        }
        vector<int> new_cut(size);
        int index = 0;
        for (int i = coeff[0]; i >= 0; --i) {
            for (auto j: tmp_cut[i]) {
                new_cut[index++] = j;
            }
        }
        chosen_cuts.emplace_back();
        chosen_cuts.back().info_r1c = make_pair(new_cut, snd);
        numCuts--;
        if (numCuts == 0) break;
    }
}

void Rank1CutsSeparator::fillMemory() {
    if (limited_memory_type == NO_MEMORY) return;
    if (!cuts.empty()) throw runtime_error("cuts should be empty at first!");
    std::vector<std::unordered_map<int, int> > v_r_map(dim);
    int r_idx = 0;
    for (const auto &i: sol) {
        for (const auto &j: i.col_seq) {
            ++v_r_map[j][r_idx];
        }
        ++r_idx;
    }
    int num_add_mem = 0;
    vector<double> vis_times(sol.size());
    int idx = 0;


    for (auto &r1c: old_cuts) {
        memset(&vis_times[0], 0, sizeof(double) * sol.size());
        const auto &plan = map_rank1_multiplier[static_cast<int>(r1c.info_r1c.first.size())][r1c.info_r1c.second];
        const auto &coeff = get<0>(plan);
        const auto deno = get<1>(plan);
        const auto rhs = get<2>(plan);
        int cnt = 0;
        for (const auto &v: r1c.info_r1c.first) {
            for (auto &[fst, snd]: v_r_map[v]) {
                vis_times[fst] += snd * coeff[cnt];
            }
            ++cnt;
        }

        transform(vis_times.begin(), vis_times.end(), sol.begin(), vis_times.begin(), [deno](auto &a, auto &b) {
            return static_cast<int>(a / deno + tolerance) * b.frac_x;
        });

        if (auto vio = accumulate(vis_times.begin(), vis_times.end(), -static_cast<double>(rhs)); vio > tolerance) {
            cuts.emplace_back(r1c);
            cut_record[r1c.info_r1c.first].insert(r1c.info_r1c.second);
            ++num_add_mem;
        }
        ++idx;
    }
    print_cuts(cout<< "num_add_mem= " << num_add_mem << endl;)
}

void Rank1CutsSeparator::generateOptimalMultiplier() {
    map_rank1_multiplier[max_heuristic_initial_seed_set_size_row_rank1c + 1].resize(7, {});
    unordered_map<int, bitset<8> > support;
    for (int i = 1; i <= max_heuristic_initial_seed_set_size_row_rank1c; ++i) {
        auto &it = map_rank1_multiplier[i];
        auto &sup = support[i];
        it.resize(7, {});
        if (i % 2) {
            get<0>(it[0]) = vector<int>(i, 1);
            get<1>(it[0]) = 2;
            get<2>(it[0]) = static_cast<int>(i / 2);
            sup.set(0);
        }
        if ((i - 2) % 3 == 0 && i >= 5) {
            get<0>(it[1]) = vector<int>(i, 1);
            get<1>(it[1]) = 3;
            get<2>(it[1]) = static_cast<int>(i / 3);
            sup.set(1);
        }
        if (i >= 5) {
            auto &tmp = get<0>(it[2]);
            tmp.resize(i, 1);
            tmp[0] = i - 3;
            tmp[1] = i - 3;
            tmp[2] = 2;
            get<1>(it[2]) = i - 2;
            get<2>(it[2]) = 2;
            sup.set(2);
            auto &tmp2 = get<0>(it[3]);
            tmp2.resize(i, 1);
            tmp2[0] = i - 2;
            tmp2[1] = i - 2;
            tmp2[2] = 2;
            tmp2[3] = 2;
            get<1>(it[3]) = i - 1;
            get<2>(it[3]) = 2;
            sup.set(3);
            auto &tmp3 = get<0>(it[4]);
            tmp3.resize(i, 1);
            tmp3[0] = i - 3;
            tmp3[1] = 2;
            get<1>(it[4]) = i - 1;
            get<2>(it[4]) = 1;
            sup.set(4);
            auto &tmp4 = get<0>(it[6]);
            tmp4.resize(i, 1);
            tmp4[0] = i - 2;
            tmp4[1] = 2;
            tmp4[2] = 2;
            get<1>(it[6]) = i;
            get<2>(it[6]) = 1;
            sup.set(6);
        }
        if (i >= 4) {
            auto &tmp5 = get<0>(it[5]);
            tmp5.resize(i, 1);
            tmp5[0] = i - 2;
            get<1>(it[5]) = i - 1;
            get<2>(it[5]) = 1;
            sup.set(5);
        }
    }


    for (int i = 1; i <= max_heuristic_initial_seed_set_size_row_rank1c; ++i) {
        for (int j = 1; j <= max_heuristic_initial_seed_set_size_row_rank1c; ++j) {
            if (i == j)continue;
            if (support[i].count() < support[j].count()) continue;
            if (i > j && map_rank1_multiplier_dominance[j].find(i) != map_rank1_multiplier_dominance[j].end())continue;
            if (((support[i] & support[j]) ^ support[j]).none()) {
                map_rank1_multiplier_dominance[i].emplace(j);
            }
        }
    }

    record_map_rank1_combinations.resize(max_heuristic_initial_seed_set_size_row_rank1c + 1,
                                         vector<vector<vector<int> > >(7, vector<vector<int> >()));
    for (int i = 1; i <= max_heuristic_initial_seed_set_size_row_rank1c; ++i) {
        for (int j = 0; j < 7; ++j) {
            if (get<1>(map_rank1_multiplier[i][j]) == 0) continue;
            unordered_map<int, int> count_map;
            for (auto &it: get<0>(map_rank1_multiplier[i][j])) {
                count_map[it]++;
            }
            vector<int> result;
            vector<vector<int> > results;
            generatePermutations(count_map, result, results, i);
            record_map_rank1_combinations[i][j] = std::move(results);
        }
    }
}

void Rank1CutsSeparator::generateSepHeurMem4Vertex() {
    rank1_sep_heur_mem4_vertex.resize(dim);
    vector<pair<int, double> > cost(dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 1; j < dim; ++j) {
            cost[j].first = j;
            cost[j].second = cost_mat4_vertex[i][j];
        }
        cost[0].first = 0;
        cost[0].second = INFINITY;
        std::stable_sort(cost.begin(), cost.end(),
                         [](const pair<int, double> &a, const pair<int, double> &b) {
                             return a.second < b.second;
                         });
        cutLong &vst2 = rank1_sep_heur_mem4_vertex[i];
        for (int k = 0; k < max_heuristic_sep_mem4_row_rank1; ++k) {
            vst2.set(cost[k].first);
        }
    }
}


void Rank1CutsSeparator::generatePermutations(std::unordered_map<int, int> &count_map,
                                              std::vector<int> &result,
                                              std::vector<std::vector<int> > &results,
                                              const int remaining) {
    if (remaining == 0) {
        results.emplace_back(result);
        return;
    }

    for (auto &pair: count_map) {
        if (pair.second > 0) {
            result.emplace_back(pair.first);
            pair.second--;
            generatePermutations(count_map, result, results, remaining - 1);
            pair.second++; // go back
            result.pop_back(); // erase and go next
        }
    }
}


void Rank1CutsSeparator::getMapPlanInfo(std::vector<int> &states, int &denominator, int &rhs, int cut_dim,
                                        int plan_idx) {
    const auto &plan = map_rank1_multiplier.at(cut_dim).at(plan_idx);
    states = get<0>(plan);
    denominator = get<1>(plan);
    rhs = get<2>(plan);
}

void Rank1CutsSeparator::getMapRhs(int &rhs, int cut_dim, int plan_idx) {
    rhs = get<2>(map_rank1_multiplier.at(cut_dim).at(plan_idx));
}

void Rank1CutsSeparator::getMapDenominator(int &denominator, int cut_dim, int plan_idx) {
    denominator = get<1>(map_rank1_multiplier.at(cut_dim).at(plan_idx));
}


bool Rank1CutsSeparator::tellIfCutsAreChangedByMem(const std::vector<R1c> &current_cuts) {
    return current_cuts.size() != old_cuts.size() ||
           !equal(current_cuts.begin(), current_cuts.end(), old_cuts.begin(),
                  [](const R1c &a, const R1c &b) {
                      return a.arc_mem == b.arc_mem;
                  });
}
