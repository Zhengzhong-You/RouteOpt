
#include "BbNode.hpp"
#include "CVRP.hpp"

using namespace std;

BbNode::BbNode(int num, int p_col, CVRP *cvrp) {
  allocateMem(num);
  num_parent_cols = p_col;
  num_rows_in_bucket_graph = cvrp->dim;
  all_forward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_forward_buckets[i] = new Bucket[cvrp->num_buckets_per_vertex];
  }
#ifdef SYMMETRY_PROHIBIT
  all_backward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_backward_buckets[i] = new Bucket[cvrp->num_buckets_per_vertex];
  }
#endif
}

#ifdef READ_ENUMERATION_TREES
BbNode::BbNode(int num) {
  allocateMem(num);
}
#endif

BbNode::BbNode(BbNode *node,
			   BidirectionalLinkedList *ptr,
			   int p_col,
			   int idx,
			   const Brc &bf,
			   int num_buckets_per_vertex) {
  ptr = ptr;
  num_parent_cols = p_col;
  index = idx;

  index_columns.assign(node->index_columns.begin(), node->index_columns.begin() + num_parent_cols);

  solver.getSolver(&node->solver);
  tree_level = node->tree_level;
  rccs = node->rccs;
  r1cs = node->r1cs;
  r1cs_multi = node->r1cs_multi;
  brcs = node->brcs;
  brcs.emplace_back(bf);
  value = node->value;
  num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
  num_forward_bucket_arcs = node->num_forward_bucket_arcs;
  num_forward_jump_arcs = node->num_forward_jump_arcs;
#ifdef SYMMETRY_PROHIBIT
  num_backward_bucket_arcs = node->num_backward_bucket_arcs;
  num_backward_jump_arcs = node->num_backward_jump_arcs;
#endif
  last_gap = node->last_gap;

  all_forward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  all_forward_buckets[i][j] = node->all_forward_buckets[i][j];
	}
  }
#ifdef USE_M_DYNAMICS
  t_for_one_lp = node->t_for_one_lp;
  geo_r_star = node->geo_r_star;
  c = node->c;
  objective_change = node->objective_change;
  l_r_ratio = node->l_r_ratio;
#endif
#ifdef SYMMETRY_PROHIBIT
  all_backward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  all_backward_buckets[i][j] = node->all_backward_buckets[i][j];
	}
  }
#endif
}

BbNode::BbNode(BbNode *node, int p_col, int idx, const Brc &bf) {
  index = idx;
  num_parent_cols = p_col;
  allocateMem((int)node->edge_head.size());
  tree_level = node->tree_level;
  rccs = node->rccs;
  r1cs = node->r1cs;
  r1cs_multi = node->r1cs_multi;
  brcs = node->brcs;
  last_gap = node->last_gap;
  brcs.emplace_back(bf);

  solver.getSolver(&node->solver);
  size_enumeration_col_pool = node->size_enumeration_col_pool;
  valid_size = node->valid_size;

  index_columns_in_enumeration_column_pool = node->index_columns_in_enumeration_column_pool;
  cost_for_columns_in_enumeration_column_pool = node->cost_for_columns_in_enumeration_column_pool;
  deleted_columns_in_enumeration_pool.assign(node->deleted_columns_in_enumeration_pool.begin(),
											 node->deleted_columns_in_enumeration_pool.begin()
												 + size_enumeration_col_pool);
  index_columns.assign(node->index_columns.begin(), node->index_columns.begin() + num_parent_cols);
  value = node->value;
  num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
#ifdef USE_M_DYNAMICS
  t_for_one_lp = node->t_for_one_lp;
  geo_r_star = node->geo_r_star;
  c = node->c;
  objective_change = node->objective_change;
  l_r_ratio = node->l_r_ratio;
#endif
}

BbNode::~BbNode() {
  if (all_forward_buckets) {
	for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	  delete[]all_forward_buckets[i];
	}
	delete[]all_forward_buckets;
#ifdef SYMMETRY_PROHIBIT
	for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	  delete[]all_backward_buckets[i];
	}
	delete[]all_backward_buckets;
#endif
	if (ptr)
	  this->ptr->deleteSelf();
  }

  solver.freeModel();
}

void BbNode::allocateMem(int num) {
  edge_head.resize(num);
  edge_tail.resize(num);
  edge_value.resize(num);
}

void BidirectionalLinkedList::becomeParent(BbNode *node, CVRP *cvrp) {
  auto &tmp_p = node->index_columns;
  int *pool = cvrp->col_pool4_mem;
  int fx_Dimension = cvrp->dim;
  int num_col = cvrp->num_col;

  if (!edge_to_cols.empty()) {
	cout << "edge_to_cols is not empty" << endl;
	exit(0);
  }
  edge_to_cols.clear();

  for (int i = 1; i < fx_Dimension; ++i) {
	edge_to_cols[i].emplace_back(i);
	edge_to_cols[i].emplace_back(i);
  }

  for (int i = fx_Dimension; i < num_col; ++i) {
	for (size_t j = tmp_p[i];;) {
	  int ai = pool[j], aj = pool[++j];
	  if (ai > aj) {
		edge_to_cols[aj * fx_Dimension + ai].emplace_back(i);
	  } else {
		edge_to_cols[ai * fx_Dimension + aj].emplace_back(i);
	  }
	  if (!aj) break;
	}
  }

  if (this == p_node->l_node) {
	p_node->l_node = nullptr;
  } else {
	p_node->r_node = nullptr;
  }
  p_node->deleteSelf();
  p_node = nullptr;
  node->num_parent_cols = num_col;
}

void BidirectionalLinkedList::giveBirth(BidirectionalLinkedList *&lnode, BidirectionalLinkedList *&rnode) {
  lnode = new BidirectionalLinkedList(this);
  rnode = new BidirectionalLinkedList(this);
  l_node = lnode;
  r_node = rnode;
}

void BidirectionalLinkedList::deleteSelf() {
  if ((!l_node) && (!r_node)) {
	if (p_node) {
	  if (this == p_node->l_node) {
		p_node->l_node = nullptr;
	  } else {
		p_node->r_node = nullptr;
	  }
	  p_node->deleteSelf();
	}
	delete this;
  }
}

PtrAllR1CS::PtrAllR1CS(BbNode *node, CVRP *cvrp) {
  int size_r1c = int(node->r1cs.size());
  r1c_to_pi = new double[size_r1c];
  num_r1c_nzero = 0;

  auto &pi = cvrp->pi4_labeling;

  for (int i = 0; i < size_r1c; ++i) {
	if (abs(pi[node->r1cs[i].idx_r1c]) < DUAL_TOLERANCE) continue;
	r1c_to_pi[num_r1c_nzero++] = pi[node->r1cs[i].idx_r1c];
  }

  int size_r1c_multi = int(node->r1cs_multi.size());
  r1c_multi_to_pi = new double[size_r1c_multi];
  num_r1c_multi_nzero = 0;

  for (int i = 0; i < size_r1c_multi; ++i) {
	if (abs(pi[node->r1cs_multi[i].idx_r1c]) < DUAL_TOLERANCE) continue;
	r1c_multi_to_pi[num_r1c_multi_nzero++] = pi[node->r1cs_multi[i].idx_r1c];
  }
  cvrp->convertVertexToR1CsInPricing(node);
}

PtrAllR1CS::~PtrAllR1CS() {
  delete[]r1c_to_pi;
  delete[]r1c_multi_to_pi;
}

BbNode::BbNode(BbNode *node, int num_buckets_per_vertex, int num_col, const bool *if_use_arc) {
  num_parent_cols = node->num_parent_cols;
  index = node->index;

  index_columns.assign(node->index_columns.begin(), node->index_columns.begin() + num_col);
  allocateMem((int)node->edge_head.size());

  solver.getSolver(&node->solver);
  safe_solver(solver.reoptimize())
  tree_level = node->tree_level;
  rccs = node->rccs;
  r1cs = node->r1cs;
  r1cs_multi = node->r1cs_multi;
  brcs = node->brcs;
  value = node->value;
  num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
  num_forward_bucket_arcs = node->num_forward_bucket_arcs;
  num_forward_jump_arcs = node->num_forward_jump_arcs;
#ifdef SYMMETRY_PROHIBIT
  num_backward_bucket_arcs = node->num_backward_bucket_arcs;
  num_backward_jump_arcs = node->num_backward_jump_arcs;
#endif
  last_gap = node->last_gap;

  all_forward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
	for (int b = 0; b < num_buckets_per_vertex; ++b) {
	  auto &bucket = all_forward_buckets[i][b];
	  for (auto &j : node->all_forward_buckets[i][b].bucket_arcs) {
		if (if_use_arc[i * num_rows_in_bucket_graph + j]) {
		  bucket.bucket_arcs.emplace_back(j);
		}
	  }
	  for (auto &arc : node->all_forward_buckets[i][b].jump_arcs) {
		if (if_use_arc[i * num_rows_in_bucket_graph + arc.second]) {
		  bucket.jump_arcs.emplace_back(arc);
		}
	  }
	}
  }
#ifdef USE_M_DYNAMICS
  t_for_one_lp = node->t_for_one_lp;
  geo_r_star = node->geo_r_star;
  c = node->c;
  objective_change = node->objective_change;
  l_r_ratio = node->l_r_ratio;
#endif
#ifdef SYMMETRY_PROHIBIT
  all_backward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
	for (int b = 0; b < num_buckets_per_vertex; ++b) {
	  auto &bucket = all_backward_buckets[i][b];
	  for (auto &j : node->all_backward_buckets[i][b].bucket_arcs) {
		if (if_use_arc[i * num_rows_in_bucket_graph + j]) {
		  bucket.bucket_arcs.emplace_back(j);
		}
	  }
	  for (auto &arc : node->all_backward_buckets[i][b].jump_arcs) {
		if (if_use_arc[i * num_rows_in_bucket_graph + arc.second]) {
		  bucket.jump_arcs.emplace_back(arc);
		}
	  }
	}
  }
#endif
}

#ifdef USE_M_DYNAMICS
void BbNode::updateState(double new_value, double &old_value, int n) {
  if (n == 0) {
	old_value = new_value;
  } else {
	old_value = old_value * pow(new_value, 1.0 / (n + 1)) * pow(1.0 / old_value, 1.0 / (n + 1));
  }
}

void BbNode::calculateRStar(double lift, double &new_r_star, CVRP *cvrp) {
  if (index == 0) throw runtime_error("calculateRStar: index==0");
  auto &edge = brcs.back().edge;
  auto dir = brcs.back().br_dir;
  auto &info = objective_change[edge];

  if (dir) {
	if (get<0>(info) == 0) {
	  get<0>(info) = max(lift * l_r_ratio, TOLERANCE);
	} else {
	  if (get<1>(info) < TOLERANCE) {
		get<0>(info) = TOLERANCE;
	  } else {
		get<0>(info) *= max(lift / get<1>(info), TOLERANCE);
	  }
	}
	get<1>(info) = max(lift, TOLERANCE);
  } else {
	get<0>(info) = lift;
	if (get<1>(info) == 0) {
	  get<1>(info) = max(lift / l_r_ratio, TOLERANCE);
	} else {
	  if (get<0>(info) < TOLERANCE) {
		get<1>(info) = TOLERANCE;
	  } else {
		get<1>(info) *= max(lift / get<0>(info), TOLERANCE);
	  }
	}
  }

  cvrp->is_use_full_k = false;
  if (get<0>(info) == TOLERANCE || get<1>(info) == TOLERANCE) {
	cvrp->is_use_full_k = true;
	new_r_star = geo_r_star;
	return;
  }

  auto l = get<0>(info);
  auto r = get<1>(info);

  auto bar_r = sqrt(l * r);
  double m = cvrp->esti_m;
  int n = cvrp->ml.give_initial_screening_num_candidates();
  auto k = min(get<2>(info) + 0., m);// generalize the root node branching decision
  new_r_star = bar_r / ((m + 1.0) / n * k / (k + 1.0) + 1.0 - m / n);
  cout << "new_r_star = " << new_r_star << endl;
  cout << "l_r_ratio1 = " << l_r_ratio << endl;
  auto new_l_r = get<0>(info) / get<1>(info);
  updateState(new_l_r, l_r_ratio, (int)brcs.size() - 1);
  cout << "l_r_ratio2 = " << l_r_ratio << endl;
}
#endif
