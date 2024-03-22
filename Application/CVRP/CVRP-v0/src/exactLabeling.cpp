
#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace Eigen;
using namespace chrono;

void CVRP::addColumns(BbNode *node, int &ccnt) {

  int ccnt_cnt;
  int curr_node, past_node;
  double cost_sum;

  vector<int> r1c_eff(node->r1cs.size(), 0);
  vector<int> r1c_multi_eff(node->r1cs_multi.size(), 0);
  vector<int> r1c_multi_state(node->r1cs_multi.size(), 0);
  unordered_map<int, vector<int>> map_node_lp;
  map_node_lp.reserve(dim * dim);
  R1CMem rank1_cut_mem = 0;
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(num_row, ccnt);
  Eigen::RowVectorXd cost(ccnt);

  ccnt_cnt = 0;
  past_node = 0;
  cost_sum = 0;
  for (size_t start = prior_pool_beg4_pricing + 1; start < pool_beg4_pricing; ++start) {
	curr_node = col_pool4_pricing[start];
	cost_sum += cost_mat4_vertex[past_node][curr_node];
	map_node_lp[past_node * dim + curr_node].emplace_back(ccnt_cnt);

	for (auto l : vertex2_all_in_one_lp_r1cs[curr_node].first) {
	  if (rank1_cut_mem[l]) {
		rank1_cut_mem[l] = false;
		++r1c_eff[l];
	  } else rank1_cut_mem[l] = true;
	}
	rank1_cut_mem &= vertex2_all_in_one_lp_r1cs[curr_node].second;

	for (auto &l : Vertex2AllInOneLPR1C_multi[curr_node].first) {
	  int tmp_cut = l.first;
	  r1c_multi_state[tmp_cut] += l.second;
	  if (r1c_multi_state[tmp_cut] >= r1c_multi_denominator_in_lp[tmp_cut]) {
		r1c_multi_state[tmp_cut] -= r1c_multi_denominator_in_lp[tmp_cut];
		++r1c_multi_eff[tmp_cut];
	  }
	}

	for (auto l : Vertex2AllInOneLPR1C_multi[curr_node].second) r1c_multi_state[l] = 0;

	if (!curr_node) {
	  cost(ccnt_cnt) = cost_sum;
	  mat(real_dim, ccnt_cnt) = 1;
	  for (int j = 0; j < r1c_eff.size(); ++j) {
		if (r1c_eff[j]) {
		  mat(node->r1cs[j].idx_r1c, ccnt_cnt) = r1c_eff[j];
		}
	  }
	  for (int j = 0; j < r1c_multi_eff.size(); ++j) {
		if (r1c_multi_eff[j]) {
		  mat(node->r1cs_multi[j].idx_r1c, ccnt_cnt) = r1c_multi_eff[j];
		}
	  }
	  ++start;
	  cost_sum = 0;
	  memset(r1c_eff.data(), 0, sizeof(int) * r1c_eff.size());
	  memset(r1c_multi_eff.data(), 0, sizeof(int) * r1c_multi_eff.size());
	  memset(r1c_multi_state.data(), 0, sizeof(int) * r1c_multi_state.size());
	  rank1_cut_mem = 0;
	  ++ccnt_cnt;
	} else {
	  ++mat(curr_node - 1, ccnt_cnt);
	}
	past_node = curr_node;
  }

  for (auto &rcc : node->rccs) {
	if (rcc.form_rcc) {
	  auto &info = rcc.info_rcc_customer;
	  int idx = rcc.idx_rcc;
	  for (auto it = info.begin(); it != info.end(); ++it) {
		int ai = *it;
		auto it_inner = it;
		++it_inner;
		for (; it_inner != info.end(); ++it_inner) {
		  int aj = *it_inner;
		  for (auto it_map : map_node_lp[ai * dim + aj]) ++mat(idx, it_map);
		  for (auto it_map : map_node_lp[aj * dim + ai]) ++mat(idx, it_map);
		}
	  }
	} else {
	  auto &customer_info = rcc.info_rcc_customer;
	  auto &outside_customer_info = rcc.info_rcc_outside_customer;
	  int idx = rcc.idx_rcc;
	  for (auto it = outside_customer_info.begin(); it != outside_customer_info.end(); ++it) {
		int ai = *it;
		auto it_inner = it;
		++it_inner;
		for (; it_inner != outside_customer_info.end(); ++it_inner) {
		  int aj = *it_inner;
		  for (auto it_map : map_node_lp[ai * dim + aj]) ++mat(idx, it_map);
		  for (auto it_map : map_node_lp[aj * dim + ai]) ++mat(idx, it_map);
		}
	  }
	  for (int aj : outside_customer_info) {
		for (auto it_map : map_node_lp[aj]) mat(idx, it_map) += 0.5;
		for (auto it_map : map_node_lp[aj * dim]) mat(idx, it_map) += 0.5;
	  }
	  for (int aj : customer_info) {
		for (auto it_map : map_node_lp[aj]) mat(idx, it_map) -= 0.5;
		for (auto it_map : map_node_lp[aj * dim]) mat(idx, it_map) -= 0.5;
	  }
	}
  }

  for (auto &br : node->brcs) {
	int ai = br.edge.first;
	int aj = br.edge.second;
	int idx = br.idx_br_c;
	for (auto it_map : map_node_lp[ai * dim + aj]) ++mat(idx, it_map);
	for (auto it_map : map_node_lp[aj * dim + ai]) ++mat(idx, it_map);
  }

  if (ccnt != ccnt_cnt) {
	cerr << "Wrong in exactLabeling 1122" << endl;
	exit(0);
  }

  Map<RowVectorXd> local_pi(pi4_labeling.data(), num_row);
  RowVectorXd rc = cost - local_pi * mat;

  ccnt_cnt = 0;
  size_t nzcnt = 0;
  vector<size_t> solver_beg(ccnt + 1);
  vector<int> solver_ind(num_row * ccnt);
  vector<double> solver_val(num_row * ccnt);
  vector<double> solver_obj(ccnt);
  for (int i = 0; i < ccnt; ++i) {
	if (rc(i) < RC_TOLERANCE) {
	  solver_beg[ccnt_cnt] = nzcnt;
	  solver_obj[ccnt_cnt] = cost(i);
	  ++ccnt_cnt;
	  for (int row = 0; row < num_row; ++row) {
		if (mat(row, i) != 0) {
		  solver_ind[nzcnt] = row;
		  solver_val[nzcnt++] = mat(row, i);
		}
	  }
	} else {
	  cout << "col " << i << " is not allowed! The rc= " << rc(i) << endl;
	}
  }

  solver_beg[ccnt_cnt] = nzcnt;

  ccnt = ccnt_cnt;
  if (!ccnt) return;

  safe_solver(node->solver.XaddVars(ccnt_cnt,
									nzcnt,
									solver_beg.data(),
									solver_ind.data(),
									solver_val.data(),
									solver_obj.data(),
									nullptr,
									nullptr,
									nullptr,
									nullptr))
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getNumCol(&num_col))
}

void CVRP::writeColumnsInPricingPool(BbNode *const node, int &index) {
  Label *p;
  if (checkPricingPool()) reallocatePricingPool();
  auto seq = new int[int(max_main_resource) + 3];// because it is not ng!
  int cnt;
  vector<size_t> tmp;
  tmp.reserve(Config::MaxNumRoutesInExact);
  for (auto &i : negative_rc_label_tuple) {
	auto &ki = get<0>(i);
	auto &kj = get<1>(i);
	tmp.emplace_back(pool_beg4_pricing);
	cnt = 0;
	p = ki;
	while (p) {
	  seq[cnt++] = p->end_vertex;
	  p = p->p_label;
	}
	for (int k = 0; k < cnt; ++k) {
	  col_pool4_pricing[pool_beg4_pricing++] = seq[cnt - 1 - k];
	}
	if (kj) {
	  p = kj;
	  while (p) {
		col_pool4_pricing[pool_beg4_pricing++] = p->end_vertex;
		p = p->p_label;
	  }
	} else {
	  col_pool4_pricing[pool_beg4_pricing++] = 0;
	}
  }
  int mem = index + (int)tmp.size() + 1;
  if (node->index_columns.size() < mem) node->index_columns.resize(mem);
  copy(tmp.begin(), tmp.end(), node->index_columns.begin() + index + 1);
  index += (int)tmp.size();
  delete[]seq;
}

void CVRP::initializeLabels(BbNode *const node,
							int mode,
							bool if_resetLabelPoint,
							tuple<bool, int, bool> control_cleanAllPtr) {

  if (if_resetLabelPoint) {
#ifdef SYMMETRY_PROHIBIT
	idx_glo = 2 * dim - 1;
#else
	idx_glo = dim;
#endif
	rc_std = RC_TOLERANCE;
	num_dominance_checks = 0;
	negative_rc_label_tuple.clear();
	map4_each_negative_rc_route.clear();
  }

  if (get<0>(control_cleanAllPtr)) {
	cleanAllPointers(node, get<1>(control_cleanAllPtr), get<2>(control_cleanAllPtr));
  }

  if (mode == 1) {
	unordered_set<int> depot_set;
	for (auto j : node->all_forward_buckets[0][0].bucket_arcs) depot_set.emplace(j);
	for (int i = 1; i < dim; ++i) {
	  if (depot_set.find(i) == depot_set.end()) continue;
	  all_label[i].rc = chg_cost_mat4_vertex[0][i];
	  all_label[i].is_extended = false;
	  all_label[i].rank1_cut_mem = 0;
	  all_label[i].num_valid_rank1_cut = 0;
	  for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[i])) {
		all_label[i].rank1_cut_mem.set(l);
		all_label[i].valid_rank1_cut[all_label[i].num_valid_rank1_cut++] = l;
	  }
	  memset(all_label[i].rank1_cut_mem_multi, 0, sizeof(int) * num_valid_r1c_multi_in_cg);
	  all_label[i].num_valid_rank1_cut_multi = 0;
	  for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[i])) {
		int tmp_cut = get<0>(l);
		all_label[i].rank1_cut_mem_multi[tmp_cut] = get<1>(l);
		all_label[i].valid_rank1_cut_multi[all_label[i].num_valid_rank1_cut_multi++] = tmp_cut;
	  }
	  int bin = int(all_label[i].sum_main_resource / step_size);
	  auto &bucket = label_array_in_forward_sense[i][bin];
	  bucket.first[bucket.second++] = all_label + i;
	  auto &bucket2 = if_exist_extra_labels_in_forward_sense[i][bin];
	  bucket2.first[bucket2.second++] = all_label + i;
	}
  }
#ifdef SYMMETRY_PROHIBIT
  else if (mode == 2) {
	unordered_set<int> depot_set;
	for (auto j : node->all_backward_buckets[0][0].bucket_arcs) depot_set.emplace(j);
	int max_num = 2 * dim - 1;
	for (int i = dim; i < max_num; ++i) {
	  int point = i - dim + 1;
	  if (depot_set.find(point) == depot_set.end()) continue;
	  all_label[i].rc = chg_cost_mat4_vertex[0][point];
	  all_label[i].is_extended = false;
	  all_label[i].rank1_cut_mem = 0;
	  all_label[i].num_valid_rank1_cut = 0;
	  for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[point])) {
		all_label[i].rank1_cut_mem.set(l);
		all_label[i].valid_rank1_cut[all_label[i].num_valid_rank1_cut++] = l;
	  }
	  memset(all_label[i].rank1_cut_mem_multi, 0, sizeof(int) * num_valid_r1c_multi_in_cg);
	  all_label[i].num_valid_rank1_cut_multi = 0;
	  for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[point])) {
		int tmp_cut = get<0>(l);
		all_label[i].rank1_cut_mem_multi[tmp_cut] = get<1>(l);
		all_label[i].valid_rank1_cut_multi[all_label[i].num_valid_rank1_cut_multi++] = tmp_cut;
	  }
	  int bin = int(all_label[i].sum_main_resource / step_size);
	  auto &bucket = label_array_in_backward_sense[point][bin];
	  bucket.first[bucket.second++] = all_label + i;
	  auto &bucket2 = if_exist_extra_labels_in_backward_sense[point][bin];
	  bucket2.first[bucket2.second++] = all_label + i;
	}
  }
#endif
}

void CVRP::addPathByRC(double path_rc, Label *ki, Label *kj, int num) {
  if (path_rc < rc_std) {
	auto p = ki;
	yzzLong tmp_seq = 0;
	while (p) {
	  tmp_seq.set(p->end_vertex);
	  p = p->p_label;
	}
	p = kj;
	while (p) {
	  tmp_seq.set(p->end_vertex);
	  p = p->p_label;
	}
	if (map4_each_negative_rc_route.find(tmp_seq) == map4_each_negative_rc_route.end()) {
	  map4_each_negative_rc_route[tmp_seq] = {ki, kj, path_rc};
	  auto it = std::lower_bound(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(), path_rc,
								 [](const std::tuple<Label *, Label *, double> &a, double b) {
								   return get<2>(a) < b;
								 });
	  negative_rc_label_tuple.insert(it, {ki, kj, path_rc});
	  if (negative_rc_label_tuple.size() > num) negative_rc_label_tuple.resize(num);
	  rc_std = get<2>(negative_rc_label_tuple.back());
	} else {
	  auto &tmp_route = map4_each_negative_rc_route[tmp_seq];
	  if (path_rc < get<2>(tmp_route)) {
		auto tmp_old = std::find(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(), tmp_route);
		if (tmp_old != negative_rc_label_tuple.end()) {
		  auto &_old = *tmp_old;
		  get<0>(_old) = ki;
		  get<1>(_old) = kj;
		  get<2>(_old) = path_rc;
		  std::sort(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(),
					[](const std::tuple<Label *, Label *, double> &a, const std::tuple<Label *, Label *, double> &b) {
					  return get<2>(a) < get<2>(b);
					});
		} else {
		  auto it = std::lower_bound(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(), path_rc,
									 [](const std::tuple<Label *, Label *, double> &a, double b) {
									   return get<2>(a) < b;
									 });
		  negative_rc_label_tuple.insert(it, {ki, kj, path_rc});
		  negative_rc_label_tuple.pop_back();
		}
		get<0>(tmp_route) = ki;
		get<1>(tmp_route) = kj;
		get<2>(tmp_route) = path_rc;
		rc_std = get<2>(negative_rc_label_tuple.back());
	  }
	}
  }
}

