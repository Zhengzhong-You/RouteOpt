
#include "CVRP.hpp"

using namespace std;

void CVRP::updateR1CStates(double &rc, R1CPricingStat &out_states,
						   const R1CPricingStat &in_states,
						   int from,
						   int to) {
  auto &in_idx = in_states.valid_cut_idx;
  auto &in_map = in_states.cut_map;
  auto in_num = in_states.num;

  auto &out_idx = out_states.valid_cut_idx;
  auto &out_map = out_states.cut_map;
  auto &out_num = out_states.num;

  auto &cut_map = cg_v_cut_map[to];
  for (auto &i : cut_map.v_union_mem) out_map[i] = 0;
  auto record = cut_map.sparse_map;
  auto &v_cut_map = cg_v_v_use_states[from][to];
  out_num = 0;
  for (int i = 0; i < in_num; ++i) {
	int idx = in_idx[i];
	int add = v_cut_map[idx];
	if (add == 0) {
	  out_idx[out_num++] = idx;
	  out_map[idx] = in_map[idx];
	} else if (add > 0) {
	  int state = in_map[idx] + add;
	  record.reset(idx);
	  if (state == cg_r1c_denominator[idx]) {
		rc -= rank1_dual[idx];
	  } else if (state < cg_r1c_denominator[idx]) {
		out_idx[out_num++] = idx;
		out_map[idx] = state;
	  } else {
		rc -= rank1_dual[idx];
		state -= cg_r1c_denominator[idx];
		out_idx[out_num++] = idx;
		out_map[idx] = state;
	  }
	}
  }

  for (auto &pr : cut_map.sparse) {
	if (record.test(pr.first)) {
	  out_idx[out_num++] = pr.first;
	  out_map[pr.first] = pr.second;
	}
  }

  /** could have been useful?
   * sort(out_idx, out_idx + out_num, [&](int a, int b) {
   * return out_map[a] < out_map[b];  });
   */
}


bool CVRP::doR1CDominance(double out_rc, double in_rc,
						  const R1CPricingStat &out_states,
						  const R1CPricingStat &in_states) {
  auto &out_num = out_states.num;
  auto &out_idx = out_states.valid_cut_idx;
  auto &out_map = out_states.cut_map;
  auto &in_map = in_states.cut_map;
  double gap = in_rc - out_rc;//out_rc is smaller
  for (int l = 0; l < out_num; ++l) {
	int cut = out_idx[l];
	if (out_map[cut] > in_map[cut]) {
	  gap += revised_rank1_dual[cut];
	  if (gap < RC_TOLERANCE) return false;
	}
  }
  return true;
}

bool CVRP::concatenateR1CStates(double &rc, double req, R1CPricingStat &out_states,
								const R1CPricingStat &in_states, int out, int in) {
  auto &out_num = out_states.num;
  auto &in_num = in_states.num;

  if (out_num < in_num) {
	auto &out_idx = out_states.valid_cut_idx;
	auto &out_map = out_states.cut_map;
	auto &in_map = in_states.cut_map;
	for (int i = 0; i < out_num; ++i) {
	  int idx = out_idx[i];
	  if (!cg_v_cut_map[in].union_map.test(idx)) continue;
	  int add = in_map[idx];
	  if (add == 0) continue;
	  int state = out_map[idx] + add;
	  if (state >= cg_r1c_denominator[idx]) {
		rc -= rank1_dual[idx];
		if (rc > req) return false;
	  }
	}
  } else {
	auto &in_idx = in_states.valid_cut_idx;
	auto &in_map = in_states.cut_map;
	auto &out_map = out_states.cut_map;
	for (int i = 0; i < in_num; ++i) {
	  int idx = in_idx[i];
	  if (!cg_v_cut_map[out].union_map.test(idx)) continue;
	  int add = out_map[idx];
	  if (add == 0) continue;
	  int state = in_map[idx] + add;
	  if (state >= cg_r1c_denominator[idx]) {
		rc -= rank1_dual[idx];
		if (rc > req) return false;
	  }
	}
  }
  return true;
}
