#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace std::chrono;

/**
 * compare the yzzLong seed and find union of the w_no_c set! do not need the dominance check!
 * and also record the information like, map[yzzLong][plan]= {best_vec, best_vio},
 * if match, then avoid calculating the same information!
 */
void CVRP::construct_v_r_map_n_seed_crazy(const vector<vector<int>> &routes,
                                          std::vector<std::unordered_map<int, int>> &v_r_map,
                                          vector<pair<vector<int>, vector<int>>> &c_N_noC) {
  v_r_map.resize(Dim);
  for (auto &i : v_r_map) i.reserve(routes.size());
  for (int r = 0; r < routes.size(); ++r) {
    for (auto &i : routes[r]) {
      ++v_r_map[i][r];
    }
  }
  int route_size = (int) routes.size();
  unordered_map<yzzLong, yzzLong> seed_map;// c ij, and w_no_c
  seed_map.reserve(4096);
  for (int i = 1; i < Dim; ++i) {
    if (v_r_map[i].empty()) continue;
    yzzLong wc = 0;// c within i
    for (auto &j : v_r_map[i]) {
      for (auto v : routes[j.first]) {
        if (Rank1SepHeurMem4Vertex[i][v]) {
          wc.set(v);
        }
      }
    }

    for (int r = 0; r < routes.size(); ++r) {
      if (v_r_map[i].find(r) != v_r_map[i].end()) continue;
      yzzLong tmp_c = 0;
      for (auto v : routes[r]) {
        if (wc[v]) {
          tmp_c.set(v);//r no i but c within i
        }
      }
      tmp_c.set(i);
      int c_size = (int) tmp_c.count();// basic c set
      if (c_size < 3 || c_size > CONFIG::MaxHeurInitialCsetSize4RowRank1C) continue;
      if (seed_map.find(tmp_c) == seed_map.end()) {
        seed_map[tmp_c] = wc ^ tmp_c;
      } else {
        seed_map[tmp_c] |= wc ^ tmp_c;//union
      }
    }
  }
  c_N_noC.resize(seed_map.size());
  int cnt = 0;
  for (auto &pr : seed_map) {
    int cnt2 = 0, cnt3 = 0;
    auto &tmp = c_N_noC[cnt];
    tmp.first.resize(pr.first.count());
    tmp.second.resize(pr.second.count());
    for (int i = 1; i < Dim; ++i) {
//      if (pr.first.test(i) && pr.second.test(i)) {
//        throw runtime_error("construct_v_r_map_n_seed_crazy: error!");
//      }
      if (pr.first.test(i)) {
        tmp.first[cnt2++] = i;
      } else if (pr.second.test(i)) {
        tmp.second[cnt3++] = i;
      }
    }
    ++cnt;
  }
//  cout << "c_N_noC.size()=" << c_N_noC.size() << endl;
//  for (auto &i : c_N_noC) {
//    cout << "----------------------" << endl;
//    cout << i.first.size() << " " << i.second.size() << endl;
//    for (auto j : i.first) cout << j << " ";
//    cout << endl;
//    for (auto j : i.second) cout << j << " ";
//    cout << endl;
//    cout << "----------------------" << endl;
//  }
}

void CVRP::generateOptimalMultiplier4R1C_multi() {
//  map_rank1_multiplier[3].resize(7, {});
  map_rank1_multiplier[CONFIG::MaxHeurInitialCsetSize4RowRank1C + 1].resize(7, {});
  unordered_map<int, yzzLong> support;
  for (int i = 1; i <= CONFIG::MaxHeurInitialCsetSize4RowRank1C; ++i) {
    auto &it = map_rank1_multiplier[i];
    auto &sup = support[i];
    it.resize(7, {});
    if (i % 2) {
      //plan 0
      get<0>(it[0]) = vector<int>(i, 1);
      get<1>(it[0]) = 2;
      get<2>(it[0]) = int(i / 2);
      sup.set(0);
    }
    if ((i - 2) % 3 == 0 && i >= 5) {
      //plan 1
      get<0>(it[1]) = vector<int>(i, 1);
      get<1>(it[1]) = 3;
      get<2>(it[1]) = int(i / 3);
      sup.set(1);
    }
    if (i >= 5) {
      //plan 2
      auto &tmp = get<0>(it[2]);
      tmp.resize(i, 1);
      tmp[0] = i - 3;
      tmp[1] = i - 3;
      tmp[2] = 2;
      get<1>(it[2]) = i - 2;
      get<2>(it[2]) = 2;
      sup.set(2);
      //plan 3
      auto &tmp2 = get<0>(it[3]);
      tmp2.resize(i, 1);
      tmp2[0] = i - 2;
      tmp2[1] = i - 2;
      tmp2[2] = 2;
      tmp2[3] = 2;
      get<1>(it[3]) = i - 1;
      get<2>(it[3]) = 2;
      sup.set(3);
      //plan 4
      auto &tmp3 = get<0>(it[4]);
      tmp3.resize(i, 1);
      tmp3[0] = i - 3;
      tmp3[1] = 2;
      get<1>(it[4]) = i - 1;
      get<2>(it[4]) = 1;
      sup.set(4);
      //plan 6
      auto &tmp4 = get<0>(it[6]);
      //note here we got 6
      tmp4.resize(i, 1);
      tmp4[0] = i - 2;
      tmp4[1] = 2;
      tmp4[2] = 2;
      get<1>(it[6]) = i;
      get<2>(it[6]) = 1;
      sup.set(6);
    }
    if (i >= 4) {
      //plan 5
      auto &tmp5 = get<0>(it[5]);
      tmp5.resize(i, 1);
      tmp5[0] = i - 2;
      get<1>(it[5]) = i - 1;
      get<2>(it[5]) = 1;
      sup.set(5);
    }
  }
  for (int i = 1; i <= CONFIG::MaxHeurInitialCsetSize4RowRank1C; ++i) {
    for (int j = 1; j <= CONFIG::MaxHeurInitialCsetSize4RowRank1C; ++j) {
      if (i == j)continue;
      if (support[i].count() < support[j].count()) continue;
      if (i > j && map_rank1_multiplier_dominance[j].find(i) != map_rank1_multiplier_dominance[j].end())continue;
      if (((support[i] & support[j]) ^ support[j]).none()) {
        map_rank1_multiplier_dominance[i].emplace(j);
      }
    }
//    cout << "i= " << i <<"| ";
//    for (auto &it : map_rank1_multiplier_dominance[i]) {
//      cout << it << " ";
//    }
//    cout << endl;
  }

  record_map_rank1_combinations.resize(CONFIG::MaxHeurInitialCsetSize4RowRank1C + 1,
                                       vector<vector<vector<int>>>(7, vector<vector<int>>()));
  for (int i = 1; i <= CONFIG::MaxHeurInitialCsetSize4RowRank1C; ++i) {
    for (int j = 0; j < 7; ++j) {
      if (get<1>(map_rank1_multiplier[i][j]) == 0) continue;
      unordered_map<int, int> count_map;
      for (auto &it : get<0>(map_rank1_multiplier[i][j])) {
        count_map[it]++;
      }
      vector<int> result;
      vector<vector<int>> results;
      generatePermutations(count_map, result, results, i);
      record_map_rank1_combinations[i][j] = std::move(results);
    }
  }

//  for (auto &it : record_map_rank1_combinations) {
//    for (auto &it2 : it) {
//      for (auto &it3 : it2) {
//        for (auto &it4 : it3) {
//          cout << it4 << " ";
//        }
//        cout << endl;
//      }
//      cout << endl;
//    }
//    cout << endl;
//  }
}

// Function to generate all permutations of the map
void CVRP::generatePermutations(std::unordered_map<int, int> &count_map,
                                std::vector<int> &result,
                                std::vector<std::vector<int>> &results,
                                int remaining) {
  if (remaining == 0) { // if all used, add them
    results.emplace_back(result);
    return;
  }

  for (auto &pair : count_map) {
    if (pair.second > 0) { // if left, add and --
      result.emplace_back(pair.first);
      pair.second--;
      generatePermutations(count_map, result, results, remaining - 1);
      pair.second++; // go back
      result.pop_back(); // erase and go next
    }
  }
}

void CVRP::exactFindBestPermuatation4Oneplan(const std::vector<std::unordered_map<int, int>> &v_r_map,
                                             const std::vector<double> &frac_routes,
                                             const std::vector<int> &cut,
                                             int plan_idx,
                                             double &vio,
                                             std::unordered_map<
                                                 yzzLong, std::vector<
                                                     std::pair<std::vector<int>, double>>> &map_cut_plan_vio) {
  int cut_size = (int) cut.size();
  auto &plan = map_rank1_multiplier[cut_size][plan_idx];
  if (!get<1>(plan)) {
    vio = -numeric_limits<double>::max();
    return;
  }
  yzzLong tmp = 0;
  for (auto &it : cut) {
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

  int denominator = get<1>(plan);
  double rhs = get<2>(plan);//must be double
  const auto &coeffs = record_map_rank1_combinations[cut_size][plan_idx];

  unordered_map<int, int> map_r_numbers;
  for (auto &i : cut) {
    for (auto &pr : v_r_map[i]) {
      map_r_numbers[pr.first] += pr.second;
    }
  }

  vector<double> valid_routes;
  valid_routes.reserve(frac_routes.size());

  unordered_map<int, int> map_old_new_routes;
  map_old_new_routes.reserve(frac_routes.size());

  for (auto &pr : map_r_numbers) {
    if (pr.second > 1) {
      valid_routes.emplace_back(frac_routes[pr.first]);
      map_old_new_routes[pr.first] = (int) valid_routes.size() - 1;
    }
  }

  vector<vector<int>> cut_num_times_vis_routes(cut_size);

  for (int i = 0; i < cut_size; ++i) {
    int c = cut[i];
    for (auto &pr : v_r_map[c]) {
      if (map_old_new_routes.find(pr.first) != map_old_new_routes.end()) {
        for (int j = 0; j < pr.second; ++j) {
          cut_num_times_vis_routes[i].emplace_back(map_old_new_routes[pr.first]);
        }
      }
    }
  }

  /// now we can calculate all the violation of all the combinations
  vector<double> num_times_vis_routes(map_old_new_routes.size());

  double best_vio = -numeric_limits<double>::max();
  int best_idx = -1;
  int cnt = 0;
  for (auto &cof : coeffs) {
    memset(num_times_vis_routes.data(), 0, sizeof(double) * num_times_vis_routes.size());
    for (int i = 0; i < cut_size; ++i) {
      for (auto &j : cut_num_times_vis_routes[i]) {
        num_times_vis_routes[j] += cof[i];
      }
    }
    transform(num_times_vis_routes.begin(),
              num_times_vis_routes.end(),
              valid_routes.begin(),
              num_times_vis_routes.begin(),
              [denominator](double a, double b) {
                return int(a / denominator + TOLERANCE) * b;
              });
    double vio_tmp = accumulate(num_times_vis_routes.begin(),
                                num_times_vis_routes.end(),
                                -rhs);
    if (vio_tmp > best_vio) {
      best_vio = vio_tmp;
      best_idx = cnt;
    }
    ++cnt;
  }

  vio = best_vio;
  vector<pair<int, int>> cut_coeff(cut_size);
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
//  yzzLong tmp1 = 0;
//  tmp1.set(133);
//  tmp1.set(42);
//  tmp1.set(62);
//  tmp1.set(284);
//
//  if (map_cut_plan_vio.find(tmp1) != map_cut_plan_vio.end() && map_cut_plan_vio[tmp1][5].second != 0) {
//    if (abs(map_cut_plan_vio[tmp1][5].second - 0.166667) > 0.01) {
//      throw runtime_error("wrong");
//    }
//  }
//  if (tmp == tmp1) {
//    cout << "best_vio= " << best_vio << endl;
//    cout << "best_idx= " << best_idx << endl;
//    cout << "vio= " << vio << endl;
//    unordered_set<int> set_tmp;
//    for (auto &i : new_cut) {
//      for (auto &pr : v_r_map[i]) {
//        set_tmp.insert(pr.first);
//      }
//    }
//    for (auto &i : set_tmp) {
//      for (auto &pr : v_r_map[i]) {
//        if (pr.second > 1) {
//          cout << i << " " << pr.second << endl;
//        }
//      }
//    }
//  }
//  if (new_cut == vector<int>{133, 42, 62, 284}) {
//    cout << "best_vio= " << best_vio << endl;
//    cout << "best_idx= " << best_idx << endl;
//    cout << "vio= " << vio << endl;
//    unordered_set<int> set_tmp;
//    for (auto &i : new_cut) {
//      for (auto &pr : v_r_map[i]) {
//        set_tmp.insert(pr.first);
//      }
//    }
//    for (auto &i : set_tmp) {
//      for (auto &pr : v_r_map[i]) {
//        if (pr.second > 1) {
//          cout << i << " " << pr.second << endl;
//        }
//      }
//    }
//  }
  map_cut_plan_vio[tmp][plan_idx] = make_pair(std::move(new_cut), vio);
}

void CVRP::start_seed_crazy(
    const std::vector<std::unordered_map<int, int>> &v_r_map,
    const vector<double> &frac_routes,
    const vector<pair<vector<int>, vector<int>>> &c_N_noC,
    std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
    vector<rank1_multi_label> &rank1_multi_label_pool,
    std::unordered_map<
        yzzLong, std::vector<
            std::pair<std::vector<int>, double>>> &map_cut_plan_vio,
    int &num_label
) {
  rank1_multi_label_pool.resize(Initial_rank1_multi_label_pool_size);
  num_label = 0;
  for (int plan_idx = 0; plan_idx < 7; ++plan_idx) {
    for (auto &p : c_N_noC) {
      auto &c = p.first;
      auto &wc = p.second;
      double vio, best_vio;
      char best_oper;
      exactFindBestPermuatation4Oneplan(v_r_map, frac_routes, c, plan_idx, vio, map_cut_plan_vio);
      if (vio < TOLERANCE) continue;
      best_vio = vio;
      best_oper = 'o';

      int add_j, remove_j;
      pair<int, int> swap_i_j;
      addSearch_crazy(plan_idx, c, wc, vio, add_j, v_r_map, frac_routes, map_cut_plan_vio);

      if (vio > best_vio) {
        best_vio = vio;
        best_oper = 'a';
      }

      removeSearch_crazy(plan_idx, c, vio, remove_j, v_r_map, frac_routes, map_cut_plan_vio);

      if (vio > best_vio) {
        best_vio = vio;
        best_oper = 'r';
      }

      swapSearch_crazy(plan_idx, c, wc, vio, swap_i_j, v_r_map, frac_routes, map_cut_plan_vio);

      if (vio > best_vio) {
        best_vio = vio;
        best_oper = 's';
      }
      yzzLong tmp;
      vector<int> new_c, new_w_no_c;
      switch (best_oper) {
        case 'o': if (c.size() < 2 || c.size() > CONFIG::MaxRowRank1) break;
          tmp = 0;
          for (auto &i : c) {
            tmp.set(i);
          }
          generated_rank1_multi_pool[(int) c.size()].emplace_back(tmp, plan_idx, best_vio);
          break;
        case 'a':new_c = c;
          new_c.emplace_back(add_j);
          new_w_no_c.resize(wc.size() - 1);
          for (int i = 0, j = 0; i < wc.size(); ++i) {
            if (wc[i] != add_j) {
              new_w_no_c[j++] = wc[i];
            }
          }
          break;
        case 'r':new_c.resize(c.size() - 1);
          for (int i = 0, j = 0; i < c.size(); ++i) {
            if (c[i] != remove_j) {
              new_c[j++] = c[i];
            }
          }
          new_w_no_c = wc;
          break;
        case 's':new_c.resize(c.size());
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
        default:throw runtime_error("error in start_seed_crazy");
      }
      if (best_oper != 'o') {
        rank1_multi_label_pool[num_label++] = rank1_multi_label{new_c, new_w_no_c, plan_idx, best_vio, best_oper};
        if (num_label == rank1_multi_label_pool.size()) {
          rank1_multi_label_pool.resize(rank1_multi_label_pool.size() * 2);
        }
      }
    }
  }
}

void CVRP::addSearch_crazy(int plan_idx,
                           const vector<int> &c,
                           const vector<int> &w_no_c,
                           double &new_vio,
                           int &add_j,
                           const vector<unordered_map<int, int>> &v_r_map,
                           const vector<double> &frac_routes,
                           unordered_map<yzzLong, vector<pair<vector<int>, double>>> &map_cut_plan_vio
) {
  int new_c_size = (int) c.size() + 1;
  const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
  if (new_c_size > CONFIG::MaxRowRank1 || !get<1>(plan)) {
    new_vio = -numeric_limits<double>::max();
    return;
  }

  vector<int> tmp_c = c;
  tmp_c.emplace_back();
  double vio, best_vio = -numeric_limits<double>::max();
  for (auto cplus : w_no_c) {
    tmp_c.back() = cplus;
    exactFindBestPermuatation4Oneplan(v_r_map, frac_routes, tmp_c, plan_idx, vio, map_cut_plan_vio);
    if (vio > best_vio) {
      best_vio = vio;
      add_j = cplus;
    }
  }
  new_vio = best_vio - TOLERANCE;// penalty
}

void CVRP::removeSearch_crazy(int plan_idx,
                              const vector<int> &c,
                              double &new_vio,
                              int &remove_j,
                              const vector<unordered_map<int, int>> &v_r_map,
                              const vector<double> &frac_routes,
                              unordered_map<yzzLong, vector<pair<vector<int>, double>>> &map_cut_plan_vio
) {
  int new_c_size = (int) c.size() - 1;
  const auto &plan = map_rank1_multiplier[new_c_size][plan_idx];
  if (new_c_size < 3 || !get<1>(plan)) {
    new_vio = -numeric_limits<double>::max();
    return;
  }

  vector<int> tmp_c(c.size() - 1);
  double vio, best_vio = -numeric_limits<double>::max();
  for (int i = 0; i < c.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      tmp_c[j] = c[j];
    }
    for (int j = i + 1; j < c.size(); ++j) {
      tmp_c[j - 1] = c[j];
    }
    exactFindBestPermuatation4Oneplan(v_r_map, frac_routes, tmp_c, plan_idx, vio, map_cut_plan_vio);
    if (vio > best_vio) {
      best_vio = vio;
      remove_j = c[i];
    }
  }
  new_vio = best_vio + TOLERANCE;// good penalty
}

void CVRP::swapSearch_crazy(int plan_idx,
                            const vector<int> &c,
                            const vector<int> &w_no_c,
                            double &new_vio,
                            pair<int, int> &swap_i_j,
                            const vector<unordered_map<int, int>> &v_r_map,
                            const vector<double> &frac_routes,
                            unordered_map<yzzLong, vector<pair<vector<int>, double>>> &map_cut_plan_vio
) {
  const auto &plan = map_rank1_multiplier[(int) c.size()][plan_idx];
  int new_c_size = (int) c.size();
  if ((new_c_size < 3 || new_c_size > CONFIG::MaxRowRank1) || !get<1>(plan)) {
    new_vio = -numeric_limits<double>::max();
    return;
  }

  vector<int> tmp_c = c;
  double vio, best_vio = -numeric_limits<double>::max();
  for (int i = 0; i < c.size(); ++i) {
    for (int j : w_no_c) {
      tmp_c[i] = j;
      exactFindBestPermuatation4Oneplan(v_r_map, frac_routes, tmp_c, plan_idx, vio, map_cut_plan_vio);
      if (vio > best_vio) {
        best_vio = vio;
        swap_i_j = {c[i], j};
      }
    }
    tmp_c[i] = c[i];
  }
  new_vio = best_vio;
}

void CVRP::crazySearch(const std::vector<std::vector<int>> &routes,
                       const std::vector<double> &frac_routes,
                       std::vector<std::unordered_map<int, int>> &v_r_map,
                       std::vector<std::pair<std::vector<int>, int>> &cuts) {
  cout << "crazySearch!" << endl;
  vector<pair<vector<int>, vector<int>>> c_N_noC;
  construct_v_r_map_n_seed_crazy(routes, v_r_map, c_N_noC);

  unordered_map<yzzLong, vector<pair<vector<int>, double>>> map_cut_plan_vio;// cut and vio
  map_cut_plan_vio.reserve(4096);

  std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>>
      generated_rank1_multi_pool;

  std::vector<rank1_multi_label> rank1_multi_label_pool;
  int num_label;
  start_seed_crazy(v_r_map,
                   frac_routes,
                   c_N_noC,
                   generated_rank1_multi_pool,
                   rank1_multi_label_pool,
                   map_cut_plan_vio,
                   num_label);

  for (int i = 0; i < num_label;) {
    operations_crazy(rank1_multi_label_pool[i],
                     v_r_map,
                     frac_routes,
                     i,
                     generated_rank1_multi_pool,
                     map_cut_plan_vio);
  }
  construct_cuts_crazy(map_cut_plan_vio, generated_rank1_multi_pool, cuts);
}

void CVRP::construct_cuts_crazy(
    const std::unordered_map<yzzLong, std::vector<std::pair<std::vector<int>, double>>> &map_cut_plan_vio,
    std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
    vector<pair<vector<int>, int>> &cuts) {

  for (auto &i : generated_rank1_multi_pool) {
    sort(i.second.begin(), i.second.end(),
         [](auto &a, auto &b) {
           return get<2>(a) > get<2>(b);
         });
  }

  for (auto &i : generated_rank1_multi_pool) {
//    if (i.first == 4) continue;
    unordered_set<yzzLong> cut_set;
    unordered_set<int> p_set;
    double vio_std = get<2>(i.second[0]) * CONFIG::CutVioFactor;
    vector<pair<vector<int>, int>> tmp_cuts;
    for (auto &j : i.second) {
      if (get<2>(j) < vio_std) break;
      auto &key = get<0>(j);
      if (cut_set.find(key) != cut_set.end() && p_set.find(get<1>(j)) != p_set.end()) continue;
#ifdef debugVio
      if (abs(get<2>(j) - map_cut_plan_vio.at(key)[get<1>(j)].second) > 0.01) {
        cout << "vio= " << get<2>(j) << endl;
        cout << "map_vio= " << map_cut_plan_vio.at(key)[get<1>(j)].second << endl;
        auto &cut = map_cut_plan_vio.at(key)[get<1>(j)].first;
        for (auto &k : cut) {
          cout << k << " ";
        }
        cout << endl;
        throw runtime_error("vio error!");
      }
#endif
      tmp_cuts.emplace_back(map_cut_plan_vio.at(key)[get<1>(j)].first, get<1>(j));
      cut_set.insert(key);
      p_set.insert(get<1>(j));
    }
#ifdef debugVio
    for (auto &j : i.second) {
      if (get<2>(j) < vio_std) break;
      auto &key = get<0>(j);
      if (vio_map.find(key) == vio_map.end()) vio_map[key].resize(7, -numeric_limits<double>::max());
//      if (vio_map[key][get<1>(j)] != -numeric_limits<double>::max()) continue;
      vio_map[key][get<1>(j)] = get<2>(j);
    }
#endif
    selectCuts(tmp_cuts, cuts, CONFIG::MaxNumR1CPerRound);
  }
//  cout << "-------------------" << endl;
//
//  vector<pair<int, double>> map_vec(map_plan_vio.begin(), map_plan_vio.end());
//  sort(map_vec.begin(), map_vec.end(), [](auto &a, auto &b) {
//    return a.first < b.first;
//  });
//  for (auto &i : map_vec) {
//    cout << i.first << " : " << i.second << endl;
//  }
}

void CVRP::operations_crazy(
    rank1_multi_label &label,
    const vector<unordered_map<int, int>> &v_r_map,
    const vector<double> &frac_routes,
    int &i,
    std::unordered_map<int, std::vector<std::tuple<yzzLong, int, double>>> &generated_rank1_multi_pool,
    unordered_map<yzzLong, vector<pair<vector<int>, double>>> &map_cut_plan_vio
) {
  auto &vio = label.vio;
  auto &new_cij = label.c;
  auto &w_no_cij = label.w_no_c;
  auto &plan_idx = label.plan_idx;

  auto dir = label.search_dir;
  int best_move;
  double best_move_vio;
  int add_j, remove_j;
  pair<int, int> swap_i_j;
  vector<pair<int, double>> move_vio(4);
  for (int j = 1; j < 4; ++j) {
    move_vio[j] = {j, -numeric_limits<double>::max()};
  }
  move_vio[0] = {0, vio};
  double new_vio;
  if (dir == 'a' || dir == 's') {
    addSearch_crazy(plan_idx, new_cij, w_no_cij, new_vio, add_j, v_r_map, frac_routes, map_cut_plan_vio);
    move_vio[1].second = new_vio;
  }
  if (dir == 'r' || dir == 's') {
    removeSearch_crazy(plan_idx, new_cij, new_vio, remove_j, v_r_map, frac_routes, map_cut_plan_vio);
    move_vio[2].second = new_vio;
  }
  if (dir == 'a' || dir == 'r') {
    swapSearch_crazy(plan_idx, new_cij, w_no_cij, new_vio, swap_i_j, v_r_map, frac_routes, map_cut_plan_vio);
    move_vio[3].second = new_vio;
  }
  sort(move_vio.begin(), move_vio.end(), [](const pair<int, double> &a,
                                            const pair<int, double> &b) {
    return a.second > b.second;
  });
  best_move = move_vio[0].first;
  best_move_vio = move_vio[0].second;

  yzzLong tmp;
  switch (best_move) {
    case 0:tmp = 0;
      for (auto j : new_cij) {
        tmp.set(j);
      }
      generated_rank1_multi_pool[(int) new_cij.size()].emplace_back(tmp, plan_idx, best_move_vio);
      ++i;
      break;
    case 1:new_cij.emplace_back(add_j);
      w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), add_j));
      break;
    case 2:new_cij.erase(find(new_cij.begin(), new_cij.end(), remove_j));
      break;
    case 3:*find(new_cij.begin(), new_cij.end(), swap_i_j.first) = swap_i_j.second;
      w_no_cij.erase(find(w_no_cij.begin(), w_no_cij.end(), swap_i_j.second));
      break;
    default:throw std::runtime_error("best_move error");
  }
  vio = best_move_vio;
}

void CVRP::sepHybridCutsInEnu(BBNODE *const node) {
  //For CVRP, we consider RCC as well as other robust cuts and SR3
  throw std::runtime_error("sepHybridCutsInEnu has been banned! use sepHybridCuts instead!");
  int cnt_tail_off = 0;
  int count = 0;
  double standard;
  int RCCnum, R1Cnum;
  double prior_val;
  //strategy for adding cuts
  //add RCC in first round and add R1Cs cuts

#ifdef VERBOSE
  cout << "Begin Hybrid separation...\n";
#endif

  while (true) {

    ++count;

    prior_val = node->Val;

    int oldNum = NumRow;

    auto beg = high_resolution_clock::now();

    generateRCCsInEnu(node);

    generateAllR1CsInEnu(node, false);

//    NewgenerateAllR1Cs(node, false);

    safe_solver(node->solver.SOLVERreoptimize())

//    deleteNonActiveCutsInEnu(node, oldNum);

    deleteNonActiveCutsSafely(node, oldNum);

    auto end = high_resolution_clock::now();
    auto eps = (double) duration_cast<milliseconds>(end - beg).count() * 1e-3;
    cout << "time for sep is " << eps << endl;

    int cuts_sum = NumRow - oldNum;
    if (!cuts_sum) {
      cout << "sep fails due to zero cuts" << endl;
      goto QUIT;
    }

    printCutsInfo(node);

    beg = high_resolution_clock::now();

    solveLPByInspection(node, false, false, true);

    findNonactiveCuts(node);

    end = high_resolution_clock::now();
    eps = (double) duration_cast<milliseconds>(end - beg).count() * 1e-3;
    cout << "time for solveLPByInspection= " << eps << endl;

    if (node->if_terminated)goto QUIT;

    if (node->SizeEnuColPool + NumCol <= MaxNumRoute4MIP) goto QUIT;

    standard = calculateGapImprovement(node->Val, prior_val);

    cout << "we delete cuts, attention! the lp formation is changed! some information maybe updated\n";

    cout << "local gap= "
         << (double(UB - node->Val) / UB > TOLERANCE ? double(UB - node->Val) / UB * 100 : 0)
         << "%\n";
    cout << "gap_improved= " << standard << endl;
    cout << SMALL_PHASE_SEPARATION;

    if (ceil_transformed_number_related(node->Val - TOLERANCE) + TOLERANCE >= UB) {
      node->if_terminated = true;
      goto QUIT;
    }

    if (standard < TOLERANCE) goto QUIT;
    else if (standard < CONFIG::CutsTailOff[2])if (++cnt_tail_off >= CONFIG::CutsToleranceInEnu)goto QUIT;
  }

  QUIT:
#ifdef VERBOSE
  cout << "Hybrid tali-off! Stop cut-separation! Last column reduction... ncol= " << NumCol
       << endl;
  cout << BIG_PHASE_SEPARATION;
#endif
}

void CVRP::addR1C_atOnceInEnu(BBNODE *node,
                              const std::vector<std::pair<std::vector<int>, int>> &cuts) {
  size_t num_nz, numnzP;
  vector<int> ai_col(NumCol);

  for (auto &cut : cuts) {
    num_nz = 0;
    memset(ai_col.data(), 0, sizeof(int) * NumCol);
    const auto &plan = map_rank1_multiplier[(int) cut.first.size()][cut.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto i : cut.first) {
      safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, i - 1, 1))
      for (size_t j = 0; j < numnzP; ++j) {
        ai_col[solver_ind[j]] += multi[count];
      }
      ++count;
    }
    int rhs = get<2>(plan);
    if (rhs > 0.1) {// is not 0
      solver_ind[num_nz] = 0;
      solver_val[num_nz++] = rhs;
    }
    for (int i = 1; i < NumCol; ++i) {
      int tmp = int(ai_col[i] / denominator);
      if (tmp) {
        solver_ind[num_nz] = i;
        solver_val[num_nz++] = tmp;
      }
    }

    safe_solver(node->solver.SOLVERaddconstr((int) num_nz,
                                             solver_ind,
                                             solver_val,
                                             SOLVER_LESS_EQUAL,
                                             rhs,
                                             nullptr))
    if (cut.second == 0) {
      R1C r1c3;
      r1c3.InfoR1C = cut.first;
      r1c3.IdxR1C = NumRow++;
      r1c3.RHS = rhs;
      node->R1Cs.emplace_back(r1c3);
    } else {
      R1C_multi r1c;
      r1c.InfoR1C = cut;
      r1c.IdxR1C = NumRow++;
      r1c.RHS = rhs;
      node->R1Cs_multi.emplace_back(r1c);
    }
  }
}