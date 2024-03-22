
#include "BbNode.hpp"
#include "CVRP.hpp"

using namespace std;

BbNode::BbNode(int num, CVRP *cvrp) {
  allocateMem(num);
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
			   int idx,
			   const Brc &bf,
			   int num_buckets_per_vertex) {
  index = idx;
  cols = node->cols;
  solver.getSolver(&node->solver);
  tree_level = node->tree_level;
  rccs = node->rccs;
  r1cs = node->r1cs;
  brcs = node->brcs;
  brcs.emplace_back(bf);
  value = node->value;
  num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
  num_forward_bucket_arcs = node->num_forward_bucket_arcs;
  num_forward_jump_arcs = node->num_forward_jump_arcs;
  topological_order_forward = node->topological_order_forward;
  canLeaveDepot_forward = node->canLeaveDepot_forward;
#ifdef SYMMETRY_PROHIBIT
  topological_order_backward = node->topological_order_backward;
  canLeaveDepot_backward=node->canLeaveDepot_backward;
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

#ifdef EXTRA_ARC_MEMORY
  arc_memory_vector = node->arc_memory_vector;
#endif
}

BbNode::BbNode(BbNode *node, int idx, const Brc &bf) {
  index = idx;
  allocateMem((int)node->edge_head.size());
  tree_level = node->tree_level;
  rccs = node->rccs;
  r1cs = node->r1cs;
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
  cols = node->cols;
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
  }

  solver.freeModel();
}

void BbNode::allocateMem(int num) {
  edge_head.resize(num);
  edge_tail.resize(num);
  edge_value.resize(num);
}

BbNode::BbNode(BbNode *node, int num_buckets_per_vertex, const bool *if_use_arc) {
  index = node->index;

  cols = node->cols;
  allocateMem((int)node->edge_head.size());

  solver.getSolver(&node->solver);
  safe_solver(solver.reoptimize())
  tree_level = node->tree_level;
  rccs = node->rccs;
  r1cs = node->r1cs;
  brcs = node->brcs;
  value = node->value;
  num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
  num_forward_bucket_arcs = node->num_forward_bucket_arcs;
  num_forward_jump_arcs = node->num_forward_jump_arcs;
  topological_order_forward = node->topological_order_forward;
  canLeaveDepot_forward = node->canLeaveDepot_forward;
#ifdef SYMMETRY_PROHIBIT
  topological_order_backward = node->topological_order_backward;
  canLeaveDepot_backward=node->canLeaveDepot_backward;
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

#ifdef EXTRA_ARC_MEMORY
  arc_memory_vector = node->arc_memory_vector;
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
  double alpha = cvrp->alpha;
  auto k = min(get<2>(info) + 0., cvrp->est_m);// generalize the root node branching decision
  new_r_star = bar_r / (1 - alpha / (k + 1));
#if VERBOSE_MODE == 1
  cout << "new_r_star = " << new_r_star << endl;
  cout << "l_r_ratio1 = " << l_r_ratio << endl;
#endif
  auto new_l_r = get<0>(info) / get<1>(info);
  updateState(new_l_r, l_r_ratio, (int)brcs.size() - 1);
#if VERBOSE_MODE == 1
  cout << "l_r_ratio2 = " << l_r_ratio << endl;
#endif
}
#endif
