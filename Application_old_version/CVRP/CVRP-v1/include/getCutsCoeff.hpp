
#ifndef GET_CUTS_COEFF_HPP
#define GET_CUTS_COEFF_HPP

#include "CVRP.hpp"

template<typename MatrixType>
void CVRP::getLimitedR1CCoeffs(const std::vector<SequenceInfo> &seq_info,
							   MatrixType &mat) {
  /**
   * lp_r1c_denominator, limited_v_idx_map_lp, forward_out_mem, backward_out_mem
   */
  mat.resize(lp_r1c_denominator.size(), seq_info.size());
  mat.setZero();

  std::vector<int> state(lp_r1c_denominator.size());
  std::vector<int> state_back(lp_r1c_denominator.size());
  std::vector<int> sparse_rep(lp_r1c_denominator.size());
  int valid_sparse_num;
  std::vector<int> sparse_rep_back(lp_r1c_denominator.size());
  int valid_sparse_num_back;

  std::unordered_map<int, int> cnt;
  cnt.reserve(lp_r1c_denominator.size());

  std::vector<Eigen::Triplet<int>> triplets;
  triplets.reserve(lp_r1c_denominator.size() * seq_info.size());

  for (int col_idx = 0; col_idx < seq_info.size(); ++col_idx) {
	const auto &seq = seq_info[col_idx];
	fill(state.begin(), state.end(), 0);
	valid_sparse_num = 0;
	cnt.clear();
	int old = 0;
	for (int i = 0; i <= seq.forward_concatenate_pos; ++i) {
	  int n = seq.col_seq[i];
	  getCoefficientExtendR1C(state, sparse_rep, cnt, valid_sparse_num, old, n);
	  old = n;
	}
	fill(state_back.begin(), state_back.end(), 0);
	valid_sparse_num_back = 0;
	old = 0;
	for (int i = (int)seq.col_seq.size() - 1; i > seq.forward_concatenate_pos; --i) {
	  int n = seq.col_seq[i];
	  getCoefficientExtendR1C(state_back, sparse_rep_back, cnt, valid_sparse_num_back, old, n);
	  old = n;
	}
	for (int i = 0; i < valid_sparse_num; ++i) {
	  int idx = sparse_rep[i];
	  if (state[idx] + state_back[idx] >= lp_r1c_denominator[idx]) {
		++cnt[idx];
	  }
	}
	int old_size = (int)triplets.size();
	triplets.resize(old_size + cnt.size());
	for (auto &pr : cnt) {
	  triplets[old_size++] = {pr.first, col_idx, pr.second};
	}
  }
  mat.setFromTriplets(triplets.begin(), triplets.end());
}

#endif //GET_CUTS_COEFF_HPP
