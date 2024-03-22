


#include "CVRP.hpp"

#define MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN 16
#define MAX_UNDETERMINED_ARC_NUMBER 128
#define MAX_LABELS 1024
using ARCBIT = std::bitset<MAX_UNDETERMINED_ARC_NUMBER>;

using namespace std;
using namespace Eigen;
using sparseRowVectorXI = Eigen::SparseVector<int, Eigen::RowMajor>;

struct State {
  int coeff{};
  int state{};
  int end_segment{};
  bitset<MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN> bit{};
};

struct Arcs {
  vector<unordered_set<pair<int, int>, PairHasher>> arc_plan;
};

void getLeastMemory(vector<Arcs> &all_arcs, unordered_set<pair<int, int>, PairHasher> &existing_arcs, bool &if_suc);

void CVRP::addLimitedMemoryR1Cs(BbNode *node,
								std::vector<R1c> &full_cuts) {
  vector<int> idx;
  writeIntoNode(node, full_cuts, idx);
  addR1CAtOnce(node, idx);
}

void CVRP::writeIntoNode(BbNode *node, vector<R1c> &full_cuts, vector<int> &idx) {
  idx.clear();
  int numRow = num_row;
  for (auto &cut : full_cuts) {
	int cut_index = cut.idx_r1c;
	const auto &big_plan = map_rank1_multiplier[(int)get<0>(cut.info_r1c).size()][get<1>(cut.info_r1c)];
	const auto &plan = get<0>(big_plan);
	int rhs = get<2>(big_plan);
	if (cut_index == numeric_limits<int>::max()) {
	  cut.idx_r1c = numRow++;
	  cut.rhs = rhs;
	  node->r1cs.emplace_back(cut);
	  idx.emplace_back((int)node->r1cs.size() - 1);
	} else {
	  auto &r1c = node->r1cs[cut_index];
	  r1c.arc_mem = cut.arc_mem;
	  idx.emplace_back(cut_index);
	}
  }
}

void findLeastPlans2MakeCoeffRight(const vector<int> &vertex_states,
								   int denominator,
								   vector<vector<int>> &plans,
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
  vector<vector<pair<vector<State>, int>>> states(max_coeff + 1, vector<pair<vector<State>, int>>(denominator));
  for (int i = 0; i < states.size(); ++i) {
	for (int j = 0; j < states[i].size(); ++j) {
	  states[i][j].first.resize((int)pow(2, i) < 100 ? (int)pow(2, i) : 100);
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
  for (int i = 0; i < states.size() - 1; ++i) {
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

void getVertexStates(const vector<int> &sequence,
					 int forward_pos,
					 const unordered_map<int, int> &map,
					 vector<int> &vertex_states,
					 vector<unordered_set<pair<int, int>, PairHasher>> &arcs) {
  vertex_states.clear();
  arcs.clear();
  vector<vector<int>> tmp_arcs;
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
	if (i == forward_pos) change_direction = {(int)tmp_arcs.size(), (int)tmp.size()};
  }
  tmp_arcs.emplace_back(tmp);

  int size = (int)tmp_arcs.size() - 1;
  if (change_direction.first != -1) {
	int idx = change_direction.first;
	for (int i = 1; i < idx; ++i) {
	  auto &arc = tmp_arcs[i];
	  auto &arc_set = arcs.emplace_back();
	  for (int j = 0; j < arc.size() - 1; ++j) {
		arc_set.emplace(arc[j], arc[j + 1]);
	  }
	}
	{
	  int i = idx;
	  if (i >= 1 && i <= size - 1) {
		int idx2 = change_direction.second - 1;
		auto &arc = tmp_arcs[i];
		auto &arc_set = arcs.emplace_back();
		for (int j = 0; j < idx2; ++j) {
		  arc_set.emplace(arc[j], arc[j + 1]);
		}
		for (int j = (int)arc.size() - 1; j > idx2 + 1; --j) {
		  arc_set.emplace(arc[j], arc[j - 1]);
		}
	  }
	}
	for (int i = idx + 1; i < size; ++i) {
	  auto &arc = tmp_arcs[i];
	  auto &arc_set = arcs.emplace_back();
	  for (int j = (int)arc.size() - 1; j > 0; --j) {
		arc_set.emplace(arc[j], arc[j - 1]);
	  }
	}
  } else {
	for (int i = 1; i < size; ++i) {
	  auto &arc = tmp_arcs[i];
	  auto &arc_set = arcs.emplace_back();
	  for (int j = 0; j < arc.size() - 1; ++j) {
		arc_set.emplace(arc[j], arc[j + 1]);
	  }
	}
  }

}

void CVRP::findLeastMemoryArcBased(BbNode *node, const sparseRowMatrixXI &sol_matrix, R1c &cut, bool &if_suc) {
  if_suc = true;
  auto &arc_mem = cut.arc_mem;
  unordered_set<pair<int, int>, PairHasher> existing_arcs;
  if (!arc_mem.empty()) {
	for (auto &info : arc_mem) {
	  int end = info.second;
	  for (int j : info.first) {
		existing_arcs.emplace(j, end);
	  }
	}
  }
#ifdef EXTRA_ARC_MEMORY
  auto old_existing = existing_arcs;
#endif
  unordered_map<int, int> map_vertex_state;
  auto &multiplier = get<0>(map_rank1_multiplier[(int)cut.info_r1c.first.size()][cut.info_r1c.second]);
  for (int i = 0; i < multiplier.size(); ++i) {
	map_vertex_state[cut.info_r1c.first[i]] = multiplier[i];
  }
  int denominator = get<1>(map_rank1_multiplier[(int)cut.info_r1c.first.size()][cut.info_r1c.second]);
  sparseRowVectorXI vec(sol_matrix.cols());
  for (auto i : cut.info_r1c.first) {
	vec += sol_matrix.row(i);
  }
  vec /= denominator;
  vec.prune(0);
  int cnt = 0;
  vector<Arcs> all_arcs;
  for (sparseRowVectorXI::InnerIterator it(vec, 0); it; ++it) {
	int col = (int)it.col();
	vector<vector<int>> ps;
	vector<int> vertex_states;
	vector<unordered_set<pair<int, int>, PairHasher>> arcs;
	getVertexStates(node->only_frac_sol.first[col].col_seq,
					node->only_frac_sol.first[col].forward_concatenate_pos,
					map_vertex_state,
					vertex_states,
					arcs);
	findLeastPlans2MakeCoeffRight(vertex_states, denominator, ps, if_suc);
	if (!if_suc) cout << "WARNING: findLeastPlans2MakeCoeffRight failed!" << endl;
	else {
	  auto &arc = all_arcs.emplace_back();
	  for (auto &p : ps) {
		auto &arc2 = arc.arc_plan.emplace_back();
		for (auto &it2 : p) {
		  arc2.insert(arcs[it2].begin(), arcs[it2].end());
		}
	  }
	}
  }

  getLeastMemory(all_arcs, existing_arcs, if_suc);


  if (if_suc) {
	arc_mem.clear();
	unordered_map<int, vector<int>> map_vertex_arcs;
	yzzLong tmp = 0;
	for (auto i : cut.info_r1c.first)tmp.set(i);
	for (auto &it : existing_arcs) {
	  if (tmp.test(it.second)) continue;
	  map_vertex_arcs[it.second].emplace_back(it.first);
	}
	for (auto &it : map_vertex_arcs) {
	  arc_mem.emplace_back(it.second, it.first);
	}
#ifdef EXTRA_ARC_MEMORY
	crossOverExtraArcMem(arc_mem);
	auto &arc_memory = node->arc_memory_vector[cut.info_r1c];
	auto &arc_memory_vector = arc_memory.first;
	auto &arc_memory_bar = arc_memory.second;
	for (auto &it : arc_mem) {
	  int in_v = it.second;
	  for (auto &arc : it.first) {
		if (old_existing.find({arc, in_v}) == old_existing.end()) {
		  arc_memory_vector[{arc, in_v}] = arc_memory_bar + 1;
		}
	  }
	}
	if (!arc_memory_vector.empty()) {
	  vector<int> tmp_vec(arc_memory_vector.size());
	  transform(arc_memory_vector.begin(),
				arc_memory_vector.end(),
				tmp_vec.begin(),
				[](const pair<pair<int, int>, int> &it) {
				  return it.second;
				});
	  auto size = int(Config::extra_arc_memory_percentage * (int)arc_memory_vector.size());
	  nth_element(tmp_vec.begin(), tmp_vec.begin() + size, tmp_vec.end(), greater<>());
	  arc_memory_bar = tmp_vec[size];
	} else arc_memory_bar = 0;
#endif
  }
}

void reduceArcs(vector<Arcs> &all_arcs, unordered_set<pair<int, int>, PairHasher> &existing_arcs, bool &if_suc) {
  if_suc = true;
  unordered_map<pair<int, int>, unordered_set<int>, PairHasher> arc_map;
  for (auto &it : all_arcs) {
	auto &arc_plan = it.arc_plan;
	arc_map.clear();
	for (int i = 0; i < arc_plan.size(); ++i) {
	  auto &arcs = arc_plan[i];
	  for (auto &arc : arcs) {
		arc_map[arc].emplace(i);
	  }
	}
	for (auto &it2 : arc_map) {
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
  for (auto &it : all_arcs) {
	auto &arc_plan = it.arc_plan;
	for (auto &arcs : arc_plan) {
	  for (auto &arc : arcs) {
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
	  for (auto &arc : arcs) {
		arc_bit[i].set(arc_bit_map[arc]);
	  }
	}
	for (int i = 0; i < arc_plan.size(); ++i) {
	  auto &arcs = arc_plan[i];
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
	  for (auto &arc : arcs) {
		existing_arcs.emplace(arc);
	  }
	  it = all_arcs.erase(it);
	} else ++it;
  }
}

void getLeastMemory(vector<Arcs> &all_arcs, unordered_set<pair<int, int>, PairHasher> &existing_arcs, bool &if_suc) {
  bool loop = true;
  while (loop) {
	int old_size = (int)existing_arcs.size();
	reduceArcs(all_arcs, existing_arcs, if_suc);
	if (old_size == existing_arcs.size()) loop = false;
	if (!if_suc) return;
  }

  if (all_arcs.empty()) return;
  vector<vector<ARCBIT >> bins(all_arcs.size());
  bins[0].resize(all_arcs[0].arc_plan.size());
  for (int i = 1; i < all_arcs.size(); ++i) {
	int tmp = 2 * bins[i - 1].size() < 100 ? int(2 * bins[i - 1].size()) : 100;
	bins[i].resize(tmp);
  }
  vector<int> bin_num(all_arcs.size(), 0);
  unordered_map<pair<int, int>, int, PairHasher> arc_bit_map;
  unordered_map<int, pair<int, int>> bit_arc_map;
  int bit_pos = 0;
  for (auto &it : all_arcs) {
	for (auto &arcs : it.arc_plan) {
	  for (auto &arc : arcs) {
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
  vector<vector<ARCBIT>> all_arcs_bit(all_arcs.size());
  for (int i = 0; i < all_arcs.size(); ++i) {
	auto &arcs = all_arcs[i].arc_plan;
	auto &arcs_bit = all_arcs_bit[i];
	arcs_bit.resize(arcs.size());
	for (int j = 0; j < arcs.size(); ++j) {
	  auto &arc = arcs[j];
	  for (auto &it : arc) {
		arcs_bit[j].set(arc_bit_map[it]);
	  }
	}
  }
  bins[0] = all_arcs_bit[0];
  bin_num[0] = (int)all_arcs_bit[0].size();
  for (int i = 0; i < bins.size() - 1; ++i) {
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

void CVRP::constructMemoryArcBased(BbNode *node, std::vector<R1c> &new_cuts) {//cut self, plan idx, cut idx, mem
  vector<Triplet<int>> tripletList;
  tripletList.reserve(size_t(double(dim * node->only_frac_sol.first.size()) * 0.05));
  unordered_map<int, int> map_vertex_index;
  for (int i = 0; i < node->only_frac_sol.first.size(); ++i) {
	map_vertex_index.clear();
	auto &sol = node->only_frac_sol.first[i];
	for (int j : sol.col_seq) {
	  ++map_vertex_index[j];
	}
	auto old_size = tripletList.size();
	tripletList.resize(old_size + map_vertex_index.size());
	for (auto &it : map_vertex_index) {
	  tripletList[old_size++] = {it.first, i, it.second};
	}
  }
  sparseRowMatrixXI sol_matrix(dim, (int)node->only_frac_sol.first.size());
  sol_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
  unordered_map<pair<vector<int>, int>, int, Rank1MultiPairHasher> R1C_Pool;//cut_index
  for (int i = 0; i < node->r1cs.size(); ++i) {
	R1C_Pool[node->r1cs[i].info_r1c] = i;
  }
  int num_r1c = MAX_NUM_R1CS_IN_PRICING - node->r1cs.size() - 1;
  if (new_cuts.size() > num_r1c) new_cuts.resize(num_r1c);

  for (auto it = new_cuts.begin(); it != new_cuts.end();) {
	auto &cut = *it;
	int index = numeric_limits<int>::max();
	auto it_find = R1C_Pool.find(cut.info_r1c);
	bool if_suc;
	if (it_find != R1C_Pool.end()) {
	  cut.arc_mem = node->r1cs[it_find->second].arc_mem;
	  index = it_find->second;
	}
#ifdef EXTRA_ARC_MEMORY
	node->populateExtraArcMem(cut.info_r1c, cut.arc_mem);
#endif
	cut.idx_r1c = index;//index, here, is the index in node->r1cs, not row index
	findLeastMemoryArcBased(node, sol_matrix, cut, if_suc);
	if (!if_suc) {
	  it = new_cuts.erase(it);
	} else ++it;
  }

}

void CVRP::selectR1CsByVioNMemory(BbNode *node,
								  const vector<vector<int>> &routes,
								  const vector<double> &frac_routes,
								  vector<R1c> &cuts) {
  unordered_map<int, vector<pair<int, int>>> map_v_cut;// key: v; val: idx & multiplier
  map_v_cut.reserve(cuts.size());
  vector<int> rhs(cuts.size());
  for (int i = 0; i < cuts.size(); ++i) {
	auto &c = cuts[i].info_r1c;
	auto &multipliers = get<0>(map_rank1_multiplier[(int)c.first.size()][c.second]);
	for (int j = 0; j < c.first.size(); ++j) {
	  map_v_cut[c.first[j]].emplace_back(i, multipliers[j]);
	}
	rhs[i] = get<2>(map_rank1_multiplier[(int)c.first.size()][c.second]);
  }
  vector<Eigen::RowVectorXd> cuts_coeffs(cuts.size(), Eigen::RowVectorXd::Zero((int)routes.size()));// use xd
  for (int i = 0; i < routes.size(); ++i) {
	for (auto j : routes[i]) {
	  for (auto &k : map_v_cut[j]) {
		cuts_coeffs[k.first][i] += k.second;
	  }
	}
  }

  for (int i = 0; i < cuts.size(); ++i) {
	auto &c = cuts[i].info_r1c;
	auto denominator = get<1>(map_rank1_multiplier[(int)c.first.size()][c.second]);
	for (auto &j : cuts_coeffs[i]) j = int(j / denominator + TOLERANCE);
  }

  vector<pair<double, int>> cuts_vio(cuts.size());
  for (int i = 0; i < cuts.size(); ++i)cuts_vio[i] = {0, i};
  Eigen::RowVectorXd frac_routes_vio(routes.size());
  for (int i = 0; i < routes.size(); ++i)frac_routes_vio[i] = frac_routes[i];
  for (int i = 0; i < cuts.size(); ++i)cuts_vio[i].first = cuts_coeffs[i].dot(frac_routes_vio) - rhs[i];

#ifdef MEMORY_SELECTION
  vector<double> mem_size(cuts.size());
#if LIMITED_MEMORY_TYPE == 1
  for (int i = 0; i < cuts.size(); ++i) {
	mem_size[i] = max((double)cuts[i].mem.size(), 1.);
  }
#elif LIMITED_MEMORY_TYPE == 2
  for (int i = 0; i < cuts.size(); ++i) {
	mem_size[i] = 0;
	for (auto &it : cuts[i].arc_mem) {
	  mem_size[i] += (double)it.first.size();
	}
	mem_size[i] = max(mem_size[i], 1.);
  }
#endif
  for (int i = 0; i < cuts.size(); ++i) cuts_vio[i].first /= log(mem_size[i] + 1);
#endif
  sort(cuts_vio.begin(), cuts_vio.end(), [](const pair<double, int> &a, const pair<double, int> &b) {
	return a.first > b.first;
  });

  vector<int> vertex_related_r1c(dim, MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX);
#if LIMITED_MEMORY_TYPE == 1
  for (auto &r1c : node->r1cs) {
	for (auto &it : r1c.info_r1c.first) {
	  --vertex_related_r1c[it];
	}
	for (auto &it : r1c.mem) {
	  --vertex_related_r1c[it];
	}
  }
#elif LIMITED_MEMORY_TYPE == 2
  for (auto &r1c : node->r1cs) {
	for (auto &it : r1c.info_r1c.first) {
	  --vertex_related_r1c[it];
	}
	for (auto &it : r1c.arc_mem) {
	  --vertex_related_r1c[it.second];
	}
  }
#endif

  vector<R1c> select_cut(cuts.size());
  int num = 0;
  for (int i = 0; i < cuts.size(); ++i) {
	int idx = cuts_vio[i].second;
	auto &c = cuts[idx];
	bool if_keep = true;
	auto tmp = vertex_related_r1c;
	for (auto j : c.info_r1c.first) {
	  if (tmp[j] <= 0) {
		if_keep = false;
		break;
	  }
	  --tmp[j];
	}
#if LIMITED_MEMORY_TYPE == 1
	for (auto j : c.mem) {
	  if (tmp[j] <= 0) {
		if_keep = false;
		break;
	  }
	  --tmp[j];
	}
#elif LIMITED_MEMORY_TYPE == 2
	for (auto &j : c.arc_mem) {
	  if (tmp[j.second] <= 0) {
		if_keep = false;
		break;
	  }
	  --tmp[j.second];
	}
#endif
	if (!if_keep) {
	  continue;
	}
	vertex_related_r1c = tmp;
	select_cut[num++] = c;
  }
  select_cut.resize(num);
  cuts = select_cut;
}

#ifdef EXTRA_ARC_MEMORY
void CVRP::crossOverExtraArcMem(std::vector<std::pair<std::vector<int>, int>> &arc_mem) {
  unordered_map<int, unordered_set<int>> in_flow_vertex;
  unordered_set<int> all_source_vertex;
  for (auto &it : arc_mem) {
	in_flow_vertex[it.second] = unordered_set<int>();
	for (auto &it2 : it.first) {
	  in_flow_vertex[it.second].emplace(it2);
	  all_source_vertex.emplace(it2);
	}
  }
  for (auto &it : in_flow_vertex) {
	int in_v = it.first;
	for (auto &it2 : all_source_vertex) {
	  if (neighborhood_indicator_4_extra_arc_mem[in_v].test(it2)) {
		it.second.emplace(it2);
	  }
	}
  }
  arc_mem.clear();
  for (auto &it : in_flow_vertex) {
	arc_mem.emplace_back(vector<int>(it.second.begin(), it.second.end()), it.first);
  }
}

void BbNode::populateExtraArcMem(const pair<vector<int>, int> &cut_info,
								 std::vector<std::pair<std::vector<int>, int>> &arc_mem) {
  if (arc_memory_vector.find(cut_info) == arc_memory_vector.end()) return;
  unordered_map<int, int> in_flow_vertex;
  unordered_set<pair<int, int>, PairHasher> in_flow_edge;
  in_flow_vertex.reserve(arc_mem.size());
  for (int i = 0; i < arc_mem.size(); ++i) {
	auto &it = arc_mem[i];
	in_flow_vertex[it.second] = i;
	for (auto &it2 : it.first) {
	  in_flow_edge.emplace(it2, it.second);
	}
  }
  auto &arc_mem_tuple = arc_memory_vector[cut_info];
  auto &arc_mem_vector = arc_mem_tuple.first;
  auto bar = arc_mem_tuple.second;
  for (auto &it : arc_mem_vector) {
	if (it.second < bar) continue;
	if (in_flow_vertex.find(it.first.second) != in_flow_vertex.end()) {
	  if (in_flow_edge.find(it.first) != in_flow_edge.end()) continue;
	  in_flow_edge.emplace(it.first);
	  arc_mem[in_flow_vertex[it.first.second]].first.emplace_back(it.first.first);
	} else {
	  in_flow_vertex[it.first.second] = (int)arc_mem.size();
	  in_flow_edge.emplace(it.first);
	  arc_mem.emplace_back(vector<int>{it.first.first}, it.first.second);
	}
  }
}
#endif