
#include "CVRP.hpp"
#include "getCutsCoeff.hpp"
using namespace std;

void CVRP::addLimitedMemoryR1CsNodeBased(BbNode *node,
										 const std::vector<R1c> &full_cuts) {
  vector<int> idx;
  giveMemInNode(node, full_cuts, idx);
  addR1CAtOnce(node, idx);
}

void CVRP::giveMemInNode(BbNode *node,
						 const std::vector<R1c> &full_cuts,
						 std::vector<int> &idx) {
  idx.clear();
  int numRow = num_row;
  vector<int> tmp(dim);
  iota(tmp.begin(), tmp.end(), 0);
  for (const auto &full_cut : full_cuts) {
	int cut_index = full_cut.idx_r1c;
	const auto &big_plan = map_rank1_multiplier[(int)full_cut.info_r1c.first.size()][full_cut.info_r1c.second];
	const auto &plan = get<0>(big_plan);
	int rhs = get<2>(big_plan);
	if (cut_index == numeric_limits<int>::max()) {
	  R1c r1c;
	  r1c.info_r1c = full_cut.info_r1c;
	  r1c.mem.assign(full_cut.mem.begin(), full_cut.mem.end());
	  r1c.arc_mem.resize(r1c.mem.size(), make_pair(tmp, 0));
	  for (int m = 0; m < r1c.mem.size(); ++m) {
		r1c.arc_mem[m].second = r1c.mem[m];
	  }
	  r1c.idx_r1c = numRow++;
	  r1c.rhs = rhs;
	  node->r1cs.emplace_back(r1c);
	  idx.emplace_back((int)node->r1cs.size() - 1);
	} else {
	  auto &r1c = node->r1cs[cut_index];
	  r1c.mem.assign(full_cut.mem.begin(), full_cut.mem.end());
	  r1c.arc_mem.assign(r1c.mem.size(), make_pair(tmp, 0));
	  for (int m = 0; m < r1c.mem.size(); ++m) {
		r1c.arc_mem[m].second = r1c.mem[m];
	  }
	  idx.emplace_back(cut_index);
	}
  }
}

void CVRP::getLimitedR1CPre(BbNode *node,
							const std::vector<int> &idx) {// idx: idx in r1c
  /**
   * get std::vector<std::vector<int>> lp_v_union_mem{};
  		std::vector<std::vector<R1CUseStates>> lp_v_v_use_states{};
  		std::vector<int> lp_r1c_denominator{};
   * to satisfy the requirement of getLimitedR1CCoeffs
   */
  lp_r1c_denominator.resize(idx.size());
  fill(lp_v_cut_map.begin(), lp_v_cut_map.end(), make_pair(R1CINDEX(), vector<int>()));

  for (int i = 0; i < dim; ++i) {// 0 is necessary
	for (int j = 1; j < dim; ++j) {
	  lp_v_v_use_states[i][j].assign(idx.size(), RANK1_INVALID);
	}
  }

  int num = 0;
  for (auto i : idx) {
	auto &r1c = node->r1cs[i];
	const auto &plan = map_rank1_multiplier[(int)r1c.info_r1c.first.size()][r1c.info_r1c.second];
	const auto &multi = get<0>(plan);
	int denominator = get<1>(plan);
	lp_r1c_denominator[num] = denominator;
	for (int j = 0; j < r1c.info_r1c.first.size(); ++j) {
	  int n = r1c.info_r1c.first[j];
	  auto &tmp_n = lp_v_cut_map[n];
	  int add = multi[j];
	  tmp_n.first.set(num);
	  tmp_n.second.emplace_back(num);
	  for (int k = 0; k < dim; ++k) {
		lp_v_v_use_states[k][n][num] = add;
	  }
	}
	for (auto &m : r1c.arc_mem) {
	  for (auto &k : m.first) {
		lp_v_v_use_states[k][m.second][num] = 0;
	  }
	}
	++num;
  }
}

void CVRP::addR1CAtOnce(BbNode *node,
						const std::vector<int> &idx) {
  getLimitedR1CPre(node, idx);

  sparseRowMatrixXd mat;
  getLimitedR1CCoeffs(node->cols, mat);

  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
  int old_num_row = num_row;
  for (int i = 0; i < idx.size(); ++i) {
	int numnz = 0;
	int row_index = node->r1cs[idx[i]].idx_r1c;
	double rhs = node->r1cs[idx[i]].rhs;
	double coeff = mat.coeff(i, 0);
	if (abs(coeff - rhs) > TOLERANCE) {
	  if (abs(coeff) > TOLERANCE) {
		solver_ind[numnz] = 0;
		solver_val[numnz] = rhs;
		++numnz;
		sparseRowMatrixXd::InnerIterator it(mat, i);
		++it;
		for (; it; ++it) {
		  solver_ind[numnz] = (int)it.col();
		  solver_val[numnz] = it.value();
		  ++numnz;
		}
	  } else {
		solver_ind[numnz] = 0;
		solver_val[numnz] = rhs;
		++numnz;
		goto HERE;
	  }
	} else {
	  HERE:
	  for (sparseRowMatrixXd::InnerIterator it(mat, i); it; ++it) {
		solver_ind[numnz] = (int)it.col();
		solver_val[numnz] = it.value();
		++numnz;
	  }
	}
	if (row_index >= old_num_row) {
	  safe_solver(node->solver.addConstraint(numnz,
											 solver_ind.data(),
											 solver_val.data(),
											 SOLVER_LESS_EQUAL,
											 rhs,
											 nullptr))
	  int vind = 0;
	  safe_solver(node->solver.changeCoeffs(1, &row_index, &vind, &rhs))
	  ++num_row;
	} else {
	  vector<int> solver_ind2(numnz);
	  fill(solver_ind2.begin(), solver_ind2.end(), row_index);
	  safe_solver(node->solver.XchangeCoeffs(numnz, solver_ind2.data(), solver_ind.data(), solver_val.data()))
	}
  }
}

void CVRP::getCoefficientExtendR1C(std::vector<int> &states,
								   std::vector<int> &sparse_rep,
								   std::unordered_map<int, int> &cnt,
								   int &valid_sparse_num,
								   int from,
								   int to
) {
  auto &v_cut_map = lp_v_v_use_states[from][to];
  auto record = lp_v_cut_map[to].first;
  for (int i = 0; i < valid_sparse_num;) {
	int idx = sparse_rep[i];
	int add = v_cut_map[idx];
	bool is_add = true;
	if (add > 0) {
	  int state = states[idx] + add;
	  record.reset(idx);
	  if (state == lp_r1c_denominator[idx]) {
		++cnt[idx];
		states[idx] = 0;
		is_add = false;
	  } else if (state < lp_r1c_denominator[idx]) {
		is_add = true;
		states[idx] = state;
	  } else {
		++cnt[idx];
		is_add = true;
		states[idx] = state - lp_r1c_denominator[idx];
	  }
	} else if (add == RANK1_INVALID) {
	  is_add = false;
	  states[idx] = 0;
	}
	if (is_add) ++i;
	else {
	  sparse_rep[i] = sparse_rep[--valid_sparse_num];
	}
  }
  for (auto &pr : lp_v_cut_map[to].second) {
	if (record.test(pr)) {
	  sparse_rep[valid_sparse_num++] = pr;
	  states[pr] = v_cut_map[pr];
	}
  }
}
