
#include "CVRP.hpp"
using namespace std;
void CVRP::getRank1DualsInCG(BbNode *node, const std::vector<double> &pi_vector) {
  std::fill(cg_v_cut_map.begin(), cg_v_cut_map.end(), R1CUseStates());
  vector<std::pair<int, double>> cut_dual;
  for (int i = 0; i < node->r1cs.size(); ++i) {
	auto &r1c = node->r1cs[i];
	double dual = abs(pi_vector[r1c.idx_r1c]);
	if (dual < DUAL_TOLERANCE)continue;
	cut_dual.emplace_back(i, dual);
  }

  sort(cut_dual.begin(), cut_dual.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
	return a.second > b.second;
  });

  for (int i = 0; i < dim; ++i) {// 0 is necessary
	for (int j = 1; j < dim; ++j) {
	  cg_v_v_use_states[i][j].assign(cut_dual.size(), RANK1_INVALID);
	}
  }

  rank1_dual.resize(cut_dual.size());
  cg_r1c_denominator.resize(cut_dual.size());
  revised_rank1_dual.resize(cut_dual.size());

  int num = 0;
  for (auto &pr : cut_dual) {
	int i = pr.first;
	auto &r1c = node->r1cs[i];
	const auto &plan = map_rank1_multiplier[(int)r1c.info_r1c.first.size()][r1c.info_r1c.second];
	const auto &multi = get<0>(plan);
	int denominator = get<1>(plan);
	cg_r1c_denominator.emplace_back(denominator);
	rank1_dual[num] = pi_vector[r1c.idx_r1c];
	if (if_exact_labeling_cg) revised_rank1_dual[num] = rank1_dual[num];
	else revised_rank1_dual[num] = rank1_dual[num] * r1c.convert_dual_coeff;
	cg_r1c_denominator[num] = denominator;
	for (int j = 0; j < r1c.info_r1c.first.size(); ++j) {
	  int n = r1c.info_r1c.first[j];
	  auto &tmp_n = cg_v_cut_map[n];
	  int add = multi[j];
	  tmp_n.sparse.emplace_back(num, add);
	  tmp_n.sparse_map.set(num);
	  tmp_n.v_union_mem.emplace_back(num);
	  tmp_n.union_map.set(num);
	  for (int k = 0; k < dim; ++k) {
		cg_v_v_use_states[k][n][num] = add;
	  }
	}
	for (auto &m : r1c.arc_mem) {
	  cg_v_cut_map[m.second].v_union_mem.emplace_back(num);
	  cg_v_cut_map[m.second].union_map.set(num);
	  for (auto &k : m.first) {
		cg_v_v_use_states[k][m.second][num] = 0;
	  }
	}
	++num;
  }

}

void CVRP::getVCutMapLP(BbNode *node) {
  lp_r1c_map.resize(node->r1cs.size());
  vector<int> idx(node->r1cs.size());
  for (int i = 0; i < node->r1cs.size(); ++i) {
	idx[i] = i;
	lp_r1c_map[i] = node->r1cs[i].idx_r1c;
  }
  getLimitedR1CPre(node, idx);
}