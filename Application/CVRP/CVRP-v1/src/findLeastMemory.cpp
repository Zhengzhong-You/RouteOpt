

#include "CVRP.hpp"

using namespace std;

void CVRP::constructMemoryVertexBased(BbNode *node, const vector<vector<int>> &routes,
									  const std::vector<std::unordered_map<int, int>> &v_r_map,
									  const vector<pair<vector<int>, int>> &cuts,
									  vector<R1c> &full_cuts) {//cut self, plan idx, cut idx, mem
  unordered_map<pair<vector<int>, int>, std::pair<std::set<int>, int>, Rank1MultiPairHasher>
	  R1C_multi_Pool;//mem & cut_index
  for (int i = 0; i < node->r1cs.size(); ++i) {
	R1C_multi_Pool[node->r1cs[i].info_r1c] =
		{set < int > (node->r1cs[i].mem.begin(), node->r1cs[i].mem.end()), i};
  }
  full_cuts.reserve(cuts.size());
  int num_r1c =  MAX_NUM_R1CS_IN_PRICING - node->r1cs.size() - 1;
  for (auto &cut : cuts) {
	set<int> mem = {};
	int index = numeric_limits<int>::max();
	auto it_find = R1C_multi_Pool.find(cut);
	bool if_suc;
	if (it_find != R1C_multi_Pool.end()) {
	  mem = it_find->second.first;
	  index = it_find->second.second;
	}
	if (!cut.second) {
	  if (--num_r1c < 0) break;
	  yzzLong v_comb = 0;
	  for (auto i : cut.first) v_comb.set(i);
	  findMemoryForR1CsMode3ByEnumerationAndMIP(v_comb, mem, routes, if_suc);
	} else {
	  if (--num_r1c < 0) break;
	  findMemoryForRank1Multi(routes, v_r_map, cut, mem, if_suc);
	}
	if (if_suc) {
	  full_cuts.emplace_back();
	  auto &tmp = full_cuts.back();
	  tmp.info_r1c = cut;
	  tmp.mem = vector<int>(mem.begin(), mem.end());
	  tmp.idx_r1c = index;
	}
  }
}



void CVRP::findMemoryForRank1Multi(const vector<vector<int>> &routes,
                               const std::vector<std::unordered_map<int, int>> &v_r_map,
                               const pair<vector<int>, int> &cut_pair,
                               set<int> &mem,
                               bool &if_suc) {
  if_suc = true;
  auto &cut = cut_pair.first;
  auto plan_idx = cut_pair.second;
  int size = (int) cut.size();
  const auto &multi = get<0>(map_rank1_multiplier[size][plan_idx]);
  auto denominator = get<1>(map_rank1_multiplier[size][plan_idx]);
  vector<vector<vector<int>>> vec_data;
  vector<vector<set<int>>> vec_segment_route;
  unordered_map<int, int> map_cut_mul;
  vector<int> num_vis_times(routes.size(), 0);
  for (int i = 0; i < cut.size(); ++i) {
    map_cut_mul[cut[i]] = multi[i];
    for (auto &pr : v_r_map[cut[i]]) {
      num_vis_times[pr.first] += multi[i] * pr.second;
    }
  }
  transform(num_vis_times.begin(), num_vis_times.end(), num_vis_times.begin(),
            [denominator](int x) { return int(x / denominator); });
  yzzLong mem_long = 0;
  int num = 0;
  for (auto &i : routes) {
    if (num_vis_times[num++] == 0) continue;
    vector<vector<int>> data;
    vector<int> vis;
    vector<set<int>> segment_route;
    set<int> tmp_seg;
    for (auto &j : i) {
      if (map_cut_mul.find(j) != map_cut_mul.end()) {
        vis.emplace_back(map_cut_mul[j]);
        segment_route.emplace_back(tmp_seg);
        tmp_seg.clear();
      } else {
        tmp_seg.insert(j);
      }
    }
    if (!segment_route.empty())
      segment_route.erase(segment_route.begin());//remove the first one
    findPlanForRank1Multi(vis, denominator, mem_long, segment_route, data);
    if (!data.empty()) {
      vec_data.emplace_back(data);
      vec_segment_route.emplace_back(segment_route);
    }
  }

  size_t cnt = 1;
  for (int i = 0; i < vec_data.size();) {
    bool if_clear = false;
    for (auto &j : vec_data[i]) {
      bool if_all_satis = true;
      for (auto k : j) {
        for (auto l : vec_segment_route[i][k]) {
          if (!mem_long[l]) {
            if_all_satis = false;
            goto outside;
          }
        }
      }
      outside:
      if (if_all_satis) {//clear this vec
        vec_data.erase(vec_data.begin() + i);
        vec_segment_route.erase(vec_segment_route.begin() + i);
        if_clear = true;
        break;
      }
    }
    if (!if_clear) {
      cnt *= (int) vec_data[i].size();
      ++i;
    }
  }

  if (cnt != 1) {
    for (auto &r : vec_segment_route) {
      for (auto &s : r) {
        for (auto i = s.begin(); i != s.end();) {
          if (mem_long[*i]) {
            i = s.erase(i);
          } else ++i;
        }
      }
    }
  }

  if (mem_long.count() > rank1_mem_size_limit) {
    if_suc = false;
    return;
  }

  for (int i = 1; i < dim; ++i) {
    if (mem_long[i]) {
      mem.emplace(i);
    }
  }

  cnt = 1;
  for (auto &i : vec_data) {// in this way we don't have to worry about overflow!
    cnt *= i.size();
    if (cnt >= FIND_MEM_USE_ENUMERATION_OR_MIP) break;
  }

  if (cnt == 1) {
    return;
  } else if (cnt < FIND_MEM_USE_ENUMERATION_OR_MIP) {
    vector<int> tmp;
    set<int> new_mem;
    int record_min = MAX_INT;
    combinations(vec_data, vec_segment_route, 0, tmp, mem, record_min, new_mem);
    mem = new_mem;
  } else {
    getMemoryByMIP(vec_data, vec_segment_route, mem, if_suc);
  }
}

void CVRP::combinationUtilAddOne(const std::vector<int> &arr,
								 std::vector<int> &tmp,
								 std::vector<std::vector<int>> &data,
								 int start,
								 int end,
								 int index,
								 int r) {
  if (index == r) {
	data.emplace_back(tmp);
	return;
  }

  for (int i = start; end - i >= (r - index - 1); i++) {
	tmp[index] = arr[i];
	combinationUtilAddOne(arr, tmp, data, i + 1,
						  end, index + 1, r);
  }
}

void CVRP::combineAll(int n, int r, vector<vector<int>> &data) {
  vector<int> arr(n);
  iota(arr.begin(), arr.end(), 0);
  vector<int> tmp(r);
  combinationUtilAddOne(arr, tmp, data, 0, (int)arr.size() - 1, 0, r);
}

void CVRP::combinationUtil(const std::vector<int> &arr,
						   std::vector<int> &tmp,
						   std::vector<std::vector<int>> &data,
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

void CVRP::combinations(const vector<vector<vector<int>>> &array,
						const vector<vector<set<int>>> &vec_segment,
						int i,
						const vector<int> &accum,
						const set<int> &mem,
						int &record_min,
						set<int> &new_mem) {
  if (i == array.size()) {
	int num = 0;
	auto tmp_mem = mem;
	for (int j = 0; j < array.size(); ++j) {
	  for (auto k : array[j][accum[j]]) {
		for (auto l : vec_segment[j][k]) {
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


void CVRP::getMemoryByMIP(const vector<vector<vector<int>>> &array,
                       const vector<vector<set<int>>> &vec_segment,
                       set<int> &mem, bool &if_suc) {
  unordered_map<int, int> new_mem_map;
  unordered_map<int, int> re_new_mem_map;
  new_mem_map.reserve(dim);

  int idx = 0;
  for (auto &r : vec_segment) {
    for (auto &s : r) {
      for (auto i : s) {
        if (new_mem_map.find(i) == new_mem_map.end()) {
          new_mem_map.emplace(i, idx);
          re_new_mem_map.emplace(idx, i);
          ++idx;
        }
      }
    }
  }

  int num_r_p = 0;
  for (auto &r : array) {
    num_r_p += (int) r.size();
  }

  int half_num_row = (int) new_mem_map.size();
  int last_half = (int) array.size();
  int local_num_row = half_num_row + last_half;
  vector<char> sense(local_num_row, SOLVER_LESS_EQUAL);
  fill(sense.begin() + half_num_row, sense.end(), SOLVER_EQUAL);
  vector<double> rhs(local_num_row, 0);
  fill(rhs.begin() + half_num_row, rhs.end(), 1);

  Solver local_solver{};
  local_solver.getEnv(&solver);//need load environment
  cout << "we build model= getMemByMIP_.lp to get mem!" << endl;
  safe_solver(local_solver.newModel("getMemByMIP_.lp", 0, nullptr, nullptr, nullptr, nullptr, nullptr))
  safe_solver(local_solver.addConstraints(local_num_row, 0, nullptr, nullptr, nullptr, sense.data(), rhs.data(), nullptr))

  int last_num_idx = half_num_row;
  vector<size_t> solver_beg;
  vector<int> solver_ind;
  vector<double> solver_val, solver_obj;
  int numRe=10000;
  solver_beg.reserve(numRe);
  solver_ind.reserve(numRe);
  solver_val.reserve(numRe);
  solver_obj.reserve(numRe);
  for (int i = 0; i < array.size(); ++i, ++last_num_idx) {
    for (auto &p : array[i]) {
	  solver_beg.emplace_back((int) solver_ind.size());
      yzzLong tmp = 0;
      for (auto n : p) {
        for (auto j : vec_segment[i][n]) {
          if (tmp[j]) continue;
          tmp.set(j);
		  solver_ind.emplace_back(new_mem_map[j]);
        }
      }
	  solver_ind.emplace_back(last_num_idx);
    }
  }
  int old_ccnt=(int) solver_beg.size();
  int old_nzcnt= (int)solver_ind.size();
  solver_obj.assign(old_ccnt, 0);
  solver_val.assign(old_nzcnt, 1);
  solver_obj.resize(old_ccnt+ half_num_row, 1);
  solver_beg.resize(old_ccnt+ half_num_row+1);
  iota(solver_beg.begin()+old_ccnt, solver_beg.end(), old_nzcnt);
  old_ccnt += half_num_row;
  solver_ind.resize(old_nzcnt+ half_num_row);
  iota(solver_ind.begin()+ old_nzcnt, solver_ind.end(), 0);
  old_nzcnt+= half_num_row;
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
  int status, left= old_ccnt-num_r_p;
  safe_solver(local_solver.getStatus(&status))
  vector<double> X(left);
  if (status == SOLVER_TIME_LIMIT) {
    cout << "time limit for getMemoryByMIP" << endl;
    if_suc = false;
    goto here;
  }
  safe_solver(local_solver.getX(num_r_p, left, X.data()))
  for (int i = 0; i < left; ++i) {
    if (X[i] > 0.5) {
      mem.emplace(re_new_mem_map[i]);
    }
  }
  here:
  safe_solver(local_solver.freeModel())
}



