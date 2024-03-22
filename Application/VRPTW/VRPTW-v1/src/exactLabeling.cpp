
#include "CVRP.hpp"
#include "getCutsCoeff.hpp"

using namespace std;
using namespace Eigen;
using namespace chrono;

void CVRP::addColumns(BbNode *node, int &ccnt) {
  if (ccnt == 0) return;

  int ccnt_cnt;
  int past_node;
  double cost_sum;

  unordered_map<pair<int, int>, vector<int>, PairHasher> edge_map;
  edge_map.reserve(dim * dim);
  sparseColMatrixXd mat(num_row, ccnt);
  mat.setZero();
  Eigen::RowVectorXd cost(ccnt);

  unordered_map<pair<int, int>, double, PairHasher> coeff_map;
  coeff_map.reserve(dim * dim);

  ccnt_cnt = 0;
  for (int i = num_col; i < node->cols.size(); ++i) {
	auto &col = node->cols[i].col_seq;
	past_node = 0;
	cost_sum = 0;
	for (int curr_node : col) {
	  cost_sum += cost_mat4_vertex[past_node][curr_node];
	  ++coeff_map[{curr_node - 1, ccnt_cnt}];
	  auto pr = past_node < curr_node ? make_pair(past_node, curr_node) : make_pair(curr_node, past_node);
	  edge_map[pr].emplace_back(ccnt_cnt);
	  past_node = curr_node;
	}
	edge_map[make_pair(0, past_node)].emplace_back(ccnt_cnt);
	cost_sum += cost_mat4_vertex[past_node][0];
	cost(ccnt_cnt) = cost_sum;
	++ccnt_cnt;
  }

  safe_Hyperparameter(ccnt != ccnt_cnt)

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
		  auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
		  for (auto col : edge_map[pr]) ++coeff_map[{idx, col}];
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
		  auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
		  for (auto it_map : edge_map[pr]) ++coeff_map[{idx, it_map}];
		}
	  }
	  for (int aj : outside_customer_info) {
		for (auto it_map : edge_map[make_pair(0, aj)]) coeff_map[{idx, it_map}] += 0.5;
	  }
	  for (int aj : customer_info) {
		for (auto it_map : edge_map[make_pair(0, aj)]) coeff_map[{idx, it_map}] -= 0.5;
	  }
	}
  }

  for (auto &br : node->brcs) {
	int ai = br.edge.first;
	int aj = br.edge.second;
	int idx = br.idx_br_c;
	auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
	for (auto it_map : edge_map[pr]) ++coeff_map[{idx, it_map}];
  }

  sparseColMatrixXd mat_r1c;

  vector<SequenceInfo> seq_info(ccnt_cnt);
  transform(node->cols.begin() + num_col, node->cols.end(), seq_info.begin(),
			[](const SequenceInfo &a) { return a; });
  getLimitedR1CCoeffs(seq_info, mat_r1c);

  size_t n = coeff_map.size() + mat_r1c.nonZeros() + ccnt;
  vector<Eigen::Triplet<double>> triplet(n);
  n = 0;

  for (auto &it : coeff_map) {
	triplet[n++] = {it.first.first, it.first.second, it.second};
  }

  for (int i = 0; i < mat_r1c.outerSize(); ++i) {
	for (sparseColMatrixXd::InnerIterator it(mat_r1c, i); it; ++it) {
	  triplet[n++] = {lp_r1c_map[it.row()], (int)it.col(), it.value()};
	}
  }
  for (int i = 0; i < ccnt; ++i) {
	triplet[n++] = {real_dim, i, 1};
  }

  mat.setFromTriplets(triplet.begin(), triplet.end());

  Map<RowVectorXd> local_pi(pi4_labeling.data(), num_row);
  RowVectorXd rc = cost - local_pi * mat;

  ccnt_cnt = 0;
  size_t nzcnt = 0;
  vector<size_t> solver_beg(ccnt + 1);
  auto numNz = mat.nonZeros();
  vector<int> solver_ind(numNz);
  vector<double> solver_val(numNz);
  vector<double> solver_obj(ccnt);
  for (int i = 0; i < ccnt; ++i) {
	if (rc(i) < -RC_TOLERANCE) {
	  solver_beg[ccnt_cnt] = nzcnt;
	  solver_obj[ccnt_cnt] = cost(i);
	  ++ccnt_cnt;
	  for (sparseColMatrixXd::InnerIterator it(mat, i); it; ++it) {
		solver_ind[nzcnt] = (int)it.row();
		solver_val[nzcnt++] = it.value();
	  }
	} else if (rc(i) > -RC_TOLERANCE * MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX) {
	  auto &col = node->cols[num_col + i];
	  for (auto j : col.col_seq) cout << j << " ";
	  cout << " | " << col.forward_concatenate_pos << endl;
	  cout << "col " << i << " is not allowed! The rc= " << rc(i) << endl;
#ifdef CHECK_PRICING_LABELS
	  if (abs(seq_rc[node->cols[num_col + i].col_seq] - rc(i)) > TOLERANCE) {
		cout << "seq_rc= " << seq_rc[node->cols[num_col + i].col_seq] << " rc= " << rc(i) << endl;
		throw runtime_error("seq_rc is not equal to rc");
	  }
#endif
	  throw runtime_error("infeasible solution");
	}
  }

#ifdef CHECK_PRICING_LABELS
  seq_rc.clear();
#endif

  if (ccnt != ccnt_cnt) {
	solver_beg[ccnt_cnt] = nzcnt;
	solver_beg.resize(ccnt_cnt + 1);
	solver_ind.resize(nzcnt);
	solver_val.resize(nzcnt);
	solver_obj.resize(ccnt_cnt);
	auto &cols = node->cols;
	int keep = num_col;
	for (int i = 0; i < ccnt; ++i) {
	  if (rc(i) < -RC_TOLERANCE) {
		cols[keep++] = cols[num_col + i];
	  }
	}
	cols.resize(keep);
  }

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

#ifdef CHANGE_DUAL
  if (sub_pricing_solver.model) {
	safe_solver(sub_pricing_solver.XaddVars(ccnt_cnt,
											nzcnt,
											solver_beg.data(),
											solver_ind.data(),
											solver_val.data(),
											solver_obj.data(),
											nullptr,
											nullptr,
											nullptr,
											nullptr))
	safe_solver(sub_pricing_solver.updateModel())
  }
#endif
}

void cleanNegativeTuple(vector<tuple<Label *, Label *, double>> &negative_rc_label_tuple) {
  if (negative_rc_label_tuple.empty()) return;
  std::sort(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(),
			[](const std::tuple<Label *, Label *, double> &a, const std::tuple<Label *, Label *, double> &b) {
			  return get<2>(a) < get<2>(b);
			});
  int cnt = 0;
  for (auto it = negative_rc_label_tuple.begin() + 1; it != negative_rc_label_tuple.end();) {
	if (abs(get<2>(*it) - get<2>(*(it - 1))) < TOLERANCE) {
	  it = negative_rc_label_tuple.erase(it);
	} else {
	  ++it;
	  ++cnt;
	  if (cnt >= Config::MaxNumRoutesInExact) break;
	}
  }
  if (negative_rc_label_tuple.size() > Config::MaxNumRoutesInExact)
	negative_rc_label_tuple.resize(Config::MaxNumRoutesInExact);
}

void CVRP::writeColumnsInPricingPool(BbNode *node, int &index) {
  Label *p;
  if (checkPricingPool()) reallocatePricingPool();

  cleanNegativeTuple(negative_rc_label_tuple);

  auto col_idx = node->cols.size();
  node->cols.resize(node->cols.size() + negative_rc_label_tuple.size());

  for (auto &i : negative_rc_label_tuple) {
	auto &ki = get<0>(i);
	auto &kj = get<1>(i);
	auto &col = node->cols[col_idx].col_seq;
	auto &res = node->cols[col_idx].main_res;
	col.reserve(dim);
	res.reserve(dim);
	p = ki;
	while (p && p->end_vertex) {
	  col.emplace_back(p->end_vertex);
	  res.emplace_back(p->res.first_res);
	  p = p->p_label;
	}
	reverse(col.begin(), col.end());
	reverse(res.begin(), res.end());

	bool if_reverse = false;
	if (!col.empty()) {
	  node->cols[col_idx].forward_concatenate_pos = (int)col.size() - 1;
	} else {
	  if_reverse = true;
	}

	if (kj) {
	  p = kj;
	  while (p && p->end_vertex) {
		col.emplace_back(p->end_vertex);
		res.emplace_back(p->res.first_res);
		p = p->p_label;
	  }
	}
	if (if_reverse) node->cols[col_idx].forward_concatenate_pos = (int)col.size() - 1;
#ifdef CHECK_PRICING_LABELS
	seq_rc[col] = get<2>(i);
#endif
	++col_idx;
  }
  index += (int)negative_rc_label_tuple.size();

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
	for (int i = 1; i < dim; ++i) {
	  if (!node->canLeaveDepot_forward.test(i)) continue;
	  auto &new_label = all_label[i];
	  new_label.rc = chg_cost_mat4_vertex[0][i];
	  new_label.is_extended = false;
	  updateR1CStates(new_label.rc, new_label.r1c, all_label->r1c, 0, i);
	  int bin = int(new_label.res.first_res / step_size);
	  auto &bucket = label_array_in_forward_sense[i][bin];
	  bucket.push_front(all_label + i);
	  auto &bucket2 = if_exist_extra_labels_in_forward_sense[i][bin];
	  bucket2.first[bucket2.second++] = all_label + i;
	}
  }
#ifdef SYMMETRY_PROHIBIT
  else if (mode == 2) {
	int max_num = 2 * dim - 1;
	for (int i = dim; i < max_num; ++i) {
	  int point = i - dim + 1;
	  if (!node->canLeaveDepot_backward.test(point)) continue;
	  auto &new_label = all_label[i];
	  new_label.rc = chg_cost_mat4_vertex[0][point];
	  new_label.is_extended = false;
	  updateR1CStates(new_label.rc, new_label.r1c, all_label->r1c, 0, point);
	  int bin = int(new_label.res.first_res / step_size);
	  auto &bucket = label_array_in_backward_sense[point][bin];
	  bucket.push_front(all_label + i);
	  auto &bucket2 = if_exist_extra_labels_in_backward_sense[point][bin];
	  bucket2.first[bucket2.second++] = all_label + i;
	}
  }
#endif
}



void CVRP::addPathByRC(double path_rc, Label *ki, Label *kj, int num) {
  if (path_rc < rc_std) {
	negative_rc_label_tuple.emplace_back(ki, kj, path_rc);
	if (negative_rc_label_tuple.size() >= num)
	  rc_std = get<2>(negative_rc_label_tuple[negative_rc_label_tuple.size() - num]);
  }
}