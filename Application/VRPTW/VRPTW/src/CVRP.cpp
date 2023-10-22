

#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace chrono;


CVRP::CVRP(const InstanceData &instanceData) {
  dim = instanceData.dim;
  real_dim = dim - 1;
  num_vehicle = instanceData.k;
  cap = instanceData.cap;
  info_vertex = instanceData.info_vertex;
  file_name = instanceData.name;
  safe_Hyperparameter(checkMaximumNumberCustomers())
}

void CVRP::lateProcessing() {
  count4_tolerance4_try_enumeration_when_arc_elimination_fails =
      Config::InitialTolerance4tryEnumerationWhenArcEliminationFails;
  max_num_route4_mip = Config::max_num_route4_mip;
  rank1_mem_size_limit = int(dim * Config::MaxCutMemFactor);
  cut_gen_time_threshold_in_pricing = Config::CutGenTimeThresholdInPricingInitial;
  size_ng_mem = Config::InitialNGSize;
  gap_tolerance4_arc_elimination_n_enumeration = Config::InitGapTolerance4ArcEliminationNEnumeration;
  num_buckets_per_vertex = Config::InitialNumBuckets;
  cost_mat4_vertex.resize(dim, vector<double>(dim, 0));
  vertex2_all_in_one_lp_r1cs.resize(dim);
  Vertex2ActiveInOnePricingR1Cs.resize(dim);
  ng_mem4_vertex.resize(dim, 0);
  rank1_sep_heur_mem4_vertex.resize(dim, 0);
  size_ng_mem4_vertex.resize(dim, Config::InitialNGSize);
  prior_mem_sets.resize(dim, Config::InitialNGSize);
  chg_cost_mat4_vertex.resize(dim, vector<double>(dim));
  round_up_tolerance = -1. / pow_self(10, transformed_number) + TOLERANCE;
  Config::MaxHeurSepMem4RowRank1=min(Config::MaxHeurSepMem4RowRank1, dim-1);

#ifdef READ_ENUMERATION_TREES
    tree_path = TREE_FOLDER + "/" + Config::tree_path;
    col_pool_path = COL_POOL_FOLDER + "/" + Config::col_pool_path;
	auto pos = tree_path.find_last_of('/')+1;
	auto end=tree_path.find_last_of('.');
	string new_path = tree_path.substr(pos, end);
	file_name= new_path;
#endif

  for (int i = 0; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      cost_mat4_vertex[i][j] = transformCost(
          sqrt_self(float((info_vertex[i][1] - info_vertex[j][1]) * (info_vertex[i][1] - info_vertex[j][1]) +
              (info_vertex[i][2] - info_vertex[j][2]) * (info_vertex[i][2] - info_vertex[j][2]))));
    }
  }
  for (int i = 0; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      cost_mat4_vertex[j][i] = cost_mat4_vertex[i][j];
    }
  }

  vector<pair<int, double>> cost(dim);
  vector<vector<int>> tex(dim, vector<int>(size_ng_mem, 0));
  for (int i = 1; i < dim; ++i) {
    for (int j = 1; j < dim; ++j) {
      cost[j].first = j;
      cost[j].second = cost_mat4_vertex[i][j];
    }
    cost[0].second = INFINITY;
    std::stable_sort(cost.begin(), cost.end(),
                     [](const pair<int, double> &a, const pair<int, double> &b) {
                       return a.second < b.second;
                     });
    yzzLong &vst = ng_mem4_vertex[i];
    for (int k = 0; k < size_ng_mem; ++k) {
      vst.set(cost[k].first);
    }
    yzzLong &vst2 = rank1_sep_heur_mem4_vertex[i];
    for (int k = 0; k < Config::MaxHeurSepMem4RowRank1; ++k) {
      vst2.set(cost[k].first);
    }
  }
  demand = new double[dim];
  for (int i = 1; i < dim; ++i) {
    demand[i] = info_vertex[i][3];
  }

  generateOptimalMultiplierForR1CMulti();
  setResourceInBucketGraph();

  step_size = 2 * int(ceil(max_main_resource / num_buckets_per_vertex / 2));
  num_buckets_per_vertex = (int) floor(max_main_resource / step_size) + 1;

  assignMemory();

#ifndef READ_ENUMERATION_TREES
  initializeBucketGraph();
  initializeLabels();
#endif

  int candi_size = dim * dim / 2;
  real_improvement_up.reserve(candi_size);
  real_improvement_down.reserve(candi_size);
  heuristic_improvement_up.reserve(candi_size);
  heuristic_improvement_down.reserve(candi_size);
  lp_testing_improvement_up.reserve(candi_size);
  lp_testing_improvement_down.reserve(candi_size);

  getLowerBoundofMinimumNumberCars();

#ifdef MASTER_VALVE_ML
  ml.calculate_prerequisites();
#endif
}

void CVRP::initializeLabels() {
  all_label[0].end_vertex = 0;
  all_label[0].p_label = nullptr;
  for (int i = 1; i < dim; ++i) {
    all_label[i].pi.set(i);
    all_label[i].cost = cost_mat4_vertex[0][i];
    all_label[i].end_vertex = i;
    all_label[i].p_label = all_label;
    increaseMainResourceConsumption(0, all_label[i].sum_main_resource, 0, i);
  }
#ifdef SYMMETRY_PROHIBIT
  int max_num = 2 * dim - 1;
  for (int i = dim; i < max_num; ++i) {
    int point = i - dim + 1;
    all_label[i].pi.set(point);
    all_label[i].cost = cost_mat4_vertex[0][point];
    all_label[i].end_vertex = point;
    all_label[i].p_label = all_label;
    decreaseMainResourceConsumption(max_main_resource, all_label[i].sum_main_resource, 0, point);
  }
#endif
}


void CVRP::buildModel() {
  auto *node = new BbNode(500, dim, this);
  initializeBucketGraphForNode(node);
  auto p0 = new BidirectionalLinkedList(nullptr);
  node->ptr = new BidirectionalLinkedList(p0);
  p0->l_node = node->ptr;

  for (int i = 1; i < dim; ++i) {
    p0->edge_to_cols[i].emplace_back(i);
    p0->edge_to_cols[i].emplace_back(i);
  }

  int tmp_cnt = 0;

  pool_beg4_mem = 0;
  if(node->index_columns.size()< dim) node->index_columns.resize(dim);
  node->index_columns[tmp_cnt++] = pool_beg4_mem;

  for (int i = 0; i < dim; ++i) {
    col_pool4_mem[i] = i;
  }
  col_pool4_mem[dim] = 0;
  pool_beg4_mem = dim + 1;
  for (int i = 1; i < dim; ++i) {
    node->index_columns[tmp_cnt++] = pool_beg4_mem;
    col_pool4_mem[pool_beg4_mem++] = 0;
    col_pool4_mem[pool_beg4_mem++] = i;
    col_pool4_mem[pool_beg4_mem++] = 0;
  }

  int nz=2 * real_dim + dim;
  vector<int> solver_beg(dim+1), solver_ind(nz);
  vector<double> solver_val(nz),solver_obj(dim);

  double obj_sum = 0;
  for (int i = 1; i < dim; ++i) {
	solver_obj[i]=2 * cost_mat4_vertex[0][i];
    obj_sum += solver_obj[i];
  }

  ub = Config::ub < obj_sum ? Config::ub : obj_sum;
  cout << "ub= " << ub << endl;

  double increase_UB = int(ub * 1.5) + 1;
  obj4_first_col = min(increase_UB, obj_sum);
  solver_obj[0] = obj4_first_col;

  safe_solver(solver.loadEnv(nullptr))
  safe_solver(solver.setEnvThreads(NUM_THREADS_LP, true))
  safe_solver(solver.setEnvOutputFlag(0, true))
  safe_solver(solver.setEnvInfUnbdInfo(1, true))
  safe_solver(solver.setEnvMIPGap(MIP_GAP_TOLERANCE, true))
  node->solver.getEnv(&solver);
  const char *model_name = "CVRP.lp";

  cout << SMALL_PHASE_SEPARATION;
  cout << "<Instance  " << file_name << "  Capacity  " << cap << ">" << endl;

  auto *rhs = new double[dim];
  char *sense = new char[dim];
  for (int i = 0; i < real_dim; ++i) {
    rhs[i] = 1;
    sense[i] = SOLVER_EQUAL;
    solver_beg[i] = 2 * i;
    solver_ind[2 * i] = 0;
    solver_ind[2 * i + 1] = i + 1;
    solver_val[2 * i] = 1;
    solver_val[2 * i + 1] = 1;
  }
  rhs[real_dim] = num_vehicle;
  sense[real_dim] = SOLVER_GREATER_EQUAL;
  solver_beg[real_dim] = 2 * real_dim;
  auto i_sup_it = solver_beg[real_dim];
  solver_ind[i_sup_it] = 0;
  solver_val[i_sup_it++] = rhs[real_dim];
  for (int i = 1; i < dim; ++i, ++i_sup_it) {
    solver_ind[i_sup_it] = i;
    solver_val[i_sup_it] = 1;
  }

  safe_solver(node->solver.newModel(model_name, dim, solver_obj.data(), nullptr, nullptr, nullptr, nullptr))
  safe_solver(node->solver.addConstraints(dim,
                                             2 * real_dim + dim,
                                             solver_beg.data(),
                                             solver_ind.data(),
                                             solver_val.data(),
                                             sense,
                                             rhs,
                                             nullptr))
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getNumRow(&num_row))
  safe_Hyperparameter(checkCSTLimit())
  safe_solver(node->solver.getNumCol(&num_col))
  safe_solver(node->solver.optimize())

  idx_node = 0;
  node->tree_level = 0;
  node->value = 0;
  node->index = idx_node;
  node->num_parent_cols = dim;

  bbt.push(node);
  lb = 0;
  lb_transformed = 0;

  num_strong_arti_vars = 1;
  num_arti_vars = dim;
#ifndef SYMMETRY_PROHIBIT
  seq_size_arti_vars = 2 * real_dim;
#else
  seq_size_arti_vars = 4 * real_dim;
#endif

  delete[] rhs;
  delete[] sense;
}


void CVRP::solveLPInLabeling(BbNode *node, bool if_open_heur, bool if_open_exact, bool if_record_sol) {
  node->is_integer = false;

  if (!node->index) {
    ratio_dominance_checks_non_dominant = {};
  }

  optimizeLPForOneIteration(node, node->value);

  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.getEnvMethod(&env_method))
  if (env_method != SOLVER_PRIMAL_SIMPLEX) {
    safe_solver(node->solver.setEnvMethod(SOLVER_PRIMAL_SIMPLEX))
    if_changed = true;
  }

  if (if_open_heur) {
    if (runColumnGenerationType(node, 1) || runColumnGenerationType(node, 2))
      goto QUIT;
  }

  if (if_open_exact) {
    if (Config::If_ExactIsMixedStyle) {
      vector<double> rank1_dual(node->r1cs.size() + node->r1cs_multi.size(), 0);
      int cnt = 0;
      for (auto &r1c : node->r1cs) {
        if (-pi4_labeling[r1c.idx_r1c] > TOLERANCE)
          rank1_dual[cnt++] = -pi4_labeling[r1c.idx_r1c];
      }
      for (auto &r1c : node->r1cs_multi) {
        if (-pi4_labeling[r1c.idx_r1c] > TOLERANCE)
          rank1_dual[cnt++] = -pi4_labeling[r1c.idx_r1c];
      }
      if (cnt) {
        gap_between_last_smallest_rc_and_rc_threshold = accumulate(rank1_dual.begin(), rank1_dual.end(), 0.0) / cnt;
      } else gap_between_last_smallest_rc_and_rc_threshold = 0;
    } else {
      gap_between_last_smallest_rc_and_rc_threshold = 0;
    }

    if (runColumnGenerationType(node, 3)) {
	  goto QUIT;
	}
    if (num_col > LP_COL_FINAL_LIMIT) cleanIndexColForNode(node, node->num_parent_cols);
    if_can_arc_elimination_by_exact_cg = true;
  }

  if (if_record_sol && rollback != 1) {
    recordOptimalColumn(node);
  }

  QUIT:
  if (if_changed) {
    safe_solver(node->solver.setEnvMethod(env_method))
  }
  if (node->is_terminated) return;
  if (!node->index && if_open_exact) {
    double aver_ratio = ratio_dominance_checks_non_dominant.first / ratio_dominance_checks_non_dominant.second;
    if (aver_ratio > Config::BucketResizeFactorRatioDominanceChecksNonDominant) {
      auto numArcs = node->num_forward_bucket_arcs + node->num_forward_jump_arcs;
#ifdef SYMMETRY_PROHIBIT
      numArcs += node->num_backward_bucket_arcs + node->num_backward_jump_arcs;
#endif
      if (double(numArcs) / dim < Config::BucketResizeFactorNumBucketArcsPerVertex
          && !if_force_not_regenerate_bucket_graph) {
        regenerateGraphBucket(node);
      }
    }
  }
}



void CVRP::changeModelForBetterDual(BbNode *node) {
  safe_solver(node->solver.changeObj(0, 1, &lp_val))

  vector<double> old_rhs(num_row);
  safe_solver(node->solver.getRhs(0, num_row, old_rhs.data()))
  vector<double> rhs(num_row, 0);
  for (auto &r1c : node->r1cs) {
    rhs[r1c.idx_r1c] = int(r1c.info_r1c.size() + r1c.mem.size());
  }
  for (auto &r1c : node->r1cs_multi) {
    rhs[r1c.idx_r1c] = int(r1c.info_r1c.first.size() + r1c.mem.size());
  }
  safe_solver(node->solver.setRhs(0, num_row, rhs.data()))
  safe_solver(node->solver.setColUpper(0, 0.))
  safe_solver(node->solver.removeColLower(0))
  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.getEnvMethod(&env_method))
  if (env_method != SOLVER_DUAL_SIMPLEX) {
    safe_solver(node->solver.setEnvMethod(SOLVER_DUAL_SIMPLEX))
    if_changed = true;
  }
  safe_solver(node->solver.reoptimize())
  bool if_new_dual = true;
  int status;
  safe_solver(node->solver.getStatus(&status))
  if (status == SOLVER_UNBOUNDED || status == SOLVER_INF_OR_UNBD) {
    if_new_dual = false;
  }
  if (if_new_dual) {
    safe_solver(node->solver.getDual(0, num_row, pi4_labeling.data()))
    safe_solver(node->solver.getSlack(0, num_row, slack4_labeling.data()))
  }
  safe_solver(node->solver.changeObj(0, 1, &obj4_first_col))
  safe_solver(node->solver.setRhs(0, num_row, old_rhs.data()))
  safe_solver(node->solver.setColLower(0, 0.))
  safe_solver(node->solver.removeColUpper(0))
  safe_solver(node->solver.reoptimize())
  if (!if_new_dual) {
    safe_solver(node->solver.getDual(0, num_row, pi4_labeling.data()))
    safe_solver(node->solver.getSlack(0, num_row, slack4_labeling.data()))
  }
  if (if_changed) {
    safe_solver(node->solver.setEnvMethod(env_method))
  }
}

void CVRP::priceLabeling(BbNode *node) {
  pi4_labeling.resize(num_row);
  slack4_labeling.resize(num_row);

  if (
      if_exact_labeling_cg && !if_ban_convert_dual && gap_between_last_smallest_rc_and_rc_threshold < 1
          && (!node->r1cs.empty() || !node->r1cs_multi.empty())
      ) {
    if (time_pricing4_iter > Config::ChangeDualByTimeRatio * time_resolve_lp4_iter) {
      changeModelForBetterDual(node);
    } else {
      if (time_resolve_lp4_iter > time_pricing4_iter) {
        if_ban_convert_dual = true;
      }
      goto here;
    }
  } else {
    here:
    safe_solver(node->solver.getDual(0, num_row, pi4_labeling.data()))
    safe_solver(node->solver.getSlack(0, num_row, slack4_labeling.data()))
  }

  for (int i = 1; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      chg_cost_mat4_vertex[i][j] = cost_mat4_vertex[i][j] - 0.5 * (pi4_labeling[i - 1] + pi4_labeling[j - 1]);
    }
  }
  for (int i = 1; i < dim; ++i) {
    chg_cost_mat4_vertex[0][i] = cost_mat4_vertex[0][i] - 0.5 * (pi4_labeling[i - 1] + pi4_labeling[real_dim]);
  }
  for (int i = 1; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      chg_cost_mat4_vertex[j][i] = chg_cost_mat4_vertex[i][j];
    }
  }
  for (int i = 1; i < dim; ++i) {
    chg_cost_mat4_vertex[i][0] = chg_cost_mat4_vertex[0][i];
  }

  double rc;
  for (auto &rcc : node->rccs) {
    if (rcc.form_rcc) {
      auto &info = rcc.info_rcc_customer;
      rc = pi4_labeling[rcc.idx_rcc];
      for (auto i = info.begin(); i != info.end(); ++i) {
        auto j = i;
        ++j;
        for (; j != info.end(); ++j) {
          chg_cost_mat4_vertex[*i][*j] -= rc;
          chg_cost_mat4_vertex[*j][*i] -= rc;
        }
      }
    } else {
      auto &outside_customer_info = rcc.info_rcc_outside_customer;
      auto &customer_info = rcc.info_rcc_customer;
      rc = pi4_labeling[rcc.idx_rcc];
      for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
        auto j = i;
        ++j;
        for (; j != outside_customer_info.end(); ++j) {
          chg_cost_mat4_vertex[*i][*j] -= rc;
          chg_cost_mat4_vertex[*j][*i] -= rc;
        }
      }
      double half_rc = 0.5 * rc;
      for (auto it : outside_customer_info) {
        chg_cost_mat4_vertex[0][it] -= half_rc;
        chg_cost_mat4_vertex[it][0] -= half_rc;
      }
      for (auto it : customer_info) {
        chg_cost_mat4_vertex[0][it] += half_rc;
        chg_cost_mat4_vertex[it][0] += half_rc;
      }
    }
  }
  for (auto &brc : node->brcs) {
    if (!brc.br_dir) {
      chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] = numeric_limits<float>::max();
      chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] =
          numeric_limits<float>::max();//do not use double since the number will overflow
    } else {
      chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] -= pi4_labeling[brc.idx_br_c];
      chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] -= pi4_labeling[brc.idx_br_c];
    }
  }
}






void CVRP::assignMemory() {
  route_in_pricing_assign = MAX_ROUTE_PRICING;
  aver_route_length = int(ceil(dim / num_vehicle)) + 5;
  mem4_pricing = aver_route_length * (route_in_pricing_assign);
  col_pool4_pricing = new int[mem4_pricing];

#ifndef READ_ENUMERATION_TREES
  label_assign = LABEL_ASSIGN;
    route_in_mem_assign = MAX_ROUTE_MEMORY;
    all_label = new Label[label_assign];
    mem4_mem = aver_route_length * (route_in_mem_assign);
    col_pool4_mem = new int[mem4_mem];
    all_valid_r1cs = new int[(label_assign) * MAX_NUM_R1CS];
    all_valid_r1c_multi = new int[(label_assign) * MAX_NUM_R1C_MULTI];
    all_r1c_multi_mem = new int[(label_assign) * MAX_NUM_R1C_MULTI];

    for (size_t i = 0; i < label_assign; ++i) {
      all_label[i].valid_rank1_cut = all_valid_r1cs + i * MAX_NUM_R1CS;
      all_label[i].valid_rank1_cut_multi = all_valid_r1c_multi + i * MAX_NUM_R1C_MULTI;
      all_label[i].rank1_cut_mem_multi = all_r1c_multi_mem + i * MAX_NUM_R1C_MULTI;
    }
#endif

  pi = new double[CST_LIMIT];
  slack = new double[CST_LIMIT];
  arc_graph = new double *[dim];
  for (int i = 0; i < dim; ++i) {
    arc_graph[i] = new double[dim];
  }
  arc_graph_revised = new double *[dim];//
  for (int i = 0; i < dim; ++i) {
    arc_graph_revised[i] = new double[dim];
  }
}


CVRP::~CVRP() {
  delete[]all_label;
  delete[]pi;
  delete[]slack;
  delete[]col_pool4_pricing;
  delete[]copy_col_pool4_pricing;
  delete[]col_pool4_mem;
  delete[]demand;
  delete[]all_valid_r1cs;
  delete[]all_valid_r1c_multi;
  delete[]all_r1c_multi_mem;

#ifndef READ_ENUMERATION_TREES
    for (int i = 0; i < dim; ++i) {
      delete[]rc2_till_this_bin_in_forward_sense[i];
    }
    delete[]rc2_till_this_bin_in_forward_sense;
    for (int i = 0; i < dim; ++i) {
      delete[]rc2_bin_in_forward_sense[i];
    }
    delete[]rc2_bin_in_forward_sense;
    for (int i = 0; i < dim; ++i) {
      delete[]label_array_in_forward_sense[i];
    }
    delete[]label_array_in_forward_sense;

    for (int i = 0; i < dim; ++i) {
      delete[]if_exist_extra_labels_in_forward_sense[i];
    }
    delete[]if_exist_extra_labels_in_forward_sense;
#endif

#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < dim; ++i) {
    delete[]rc2_till_this_bin_in_backward_sense[i];
  }
  delete[]rc2_till_this_bin_in_backward_sense;
  for (int i = 0; i < dim; ++i) {
    delete[]rc2_bin_in_backward_sense[i];
  }
  delete[]rc2_bin_in_backward_sense;
  for (int i = 0; i < dim; ++i) {
    delete[]label_array_in_backward_sense[i];
  }
  delete[]label_array_in_backward_sense;
  for (int i = 0; i < dim; ++i) {
    delete[]if_exist_extra_labels_in_backward_sense[i];
  }
  delete[]if_exist_extra_labels_in_backward_sense;
#endif

  for (int i = 0; i < dim; ++i) {
    delete[]arc_graph[i];
  }
  delete[]arc_graph;
  for (int i = 0; i < dim; ++i) {
    delete[]arc_graph_revised[i];
  }
  delete[]arc_graph_revised;
}

bool operator==(const Rcc &lhs, const Rcc &rhs) {
  if (lhs.rhs != rhs.rhs) return false;
  if (lhs.form_rcc != rhs.form_rcc) {
    return false;
  } else if (lhs.form_rcc) {
    if (lhs.info_rcc_customer != rhs.info_rcc_customer) {
      return false;
    }
  } else {
    if (lhs.info_rcc_outside_customer != rhs.info_rcc_outside_customer) {
      return false;
    }
  }
  return true;
}

bool CmpLabelRCLess(const Label *l1, const Label *l2) {
  return (l1->rc < l2->rc);
}

void CVRP::recordOptimalColumn(BbNode *node, bool if_force_rewrite) {
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getObjVal(&lp_val))
  node->value = lp_val;
  if (ceilTransformedNumberRelated(node->value - TOLERANCE) + TOLERANCE >= ub) {
    node->is_terminated = true;
    return;
  }
  vector<double> X(num_col);
  safe_solver(node->solver.getX(0, num_col, X.data()))

  int parent_column = if_in_enu_state ? dim : node->num_parent_cols;
  vector<pair<size_t, double>> tmp_solindex;

  for (int i = 0; i < parent_column; ++i) {
    if (X[i] > TOLERANCE) {
      tmp_solindex.emplace_back(node->index_columns[i], X[i]);
    }
  }

  int break_record = (int) tmp_solindex.size();
  node->num_parent_cols_in_lp_solutions = break_record;

  for (int i = parent_column; i < num_col; ++i) {
    if (X[i] > TOLERANCE) {
      tmp_solindex.emplace_back(node->index_columns[i], X[i]);
    }
  }

  if (node->index_for_lp_solutions_in_column_pool == tmp_solindex) {
    cout << "sol remains the same" << endl;
    if (if_force_rewrite) {
      goto REWRITE;
    } else {
      return;
    }
  }

  node->index_for_lp_solutions_in_column_pool = tmp_solindex;

  REWRITE:
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      arc_graph[i][j] = 0;
      arc_graph_revised[i][j] = 0;
    }
  }

  int beg = 0;
  int end = (int) tmp_solindex.size();
  int *colpool = col_pool4_pricing;
  if (!if_in_enu_state) {
    end = break_record;
    colpool = col_pool4_mem;
  }

  AGAIN:
  for (int i = beg; i < end; ++i) {
    for (size_t j = tmp_solindex[i].first;;) {
      arc_graph[colpool[j]][colpool[j + 1]] += tmp_solindex[i].second;
      if (!colpool[++j]) break;
    }
  }

  if (end != tmp_solindex.size()) {
    beg = end;
    end = (int) tmp_solindex.size();
    colpool = col_pool4_pricing;
    goto AGAIN;
  }

  for (int i = 0; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      arc_graph[i][j] += arc_graph[j][i];
      arc_graph_revised[i][j] = arc_graph[i][j];
    }
  }


  vector<pair<int, double>> tmp_first_dim;
  for (int i = 0; i < dim; ++i) if (X[i] > TOLERANCE) tmp_first_dim.emplace_back(i, X[i]);
  if (if_in_enu_state) {
    for (auto &i : tmp_first_dim) {
      if (i.first) {
        int valid_length = 0, curr_node;
        for (size_t j = node->index_columns[i.first] + 1;; ++j) {
          curr_node = col_pool4_pricing[j];
          if (!curr_node) break;
          ++valid_length;
        }
        curr_node = col_pool4_pricing[node->index_columns[i.first] + 1];
        if (valid_length == 1)arc_graph_revised[0][curr_node] -= i.second;
        else break;
      }
    }
  } else {
    for (auto &i : tmp_first_dim)
      if (i.first) arc_graph_revised[0][i.first] -= i.second;
  }

  bool if_reformulate = true;
  for (int i = 0; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      if (arc_graph_revised[i][j] > TOLERANCE && abs(arc_graph_revised[i][j] - 1) > TOLERANCE) {
        if_reformulate = false;
        goto outside;
      }
    }
  }

  outside:
  if (if_reformulate) {
    reformulateIntegerProgramSolution(node);
  }

  num_edge = 0;

  if (if_force_rewrite) {
    for (int i = 0; i < dim; ++i) {
      for (int j = i + 1; j < dim; ++j) {
        if (arc_graph_revised[i][j] > TOLERANCE) {
          ++num_edge;
        }
      }
    }
  } else {
    for (int i = 0; i < dim; ++i) {
      for (int j = i + 1; j < dim; ++j) {
        if (arc_graph_revised[i][j] > TOLERANCE) {
          ++num_edge;
        }
#ifdef MASTER_VALVE_ML
        auto &tmp = ml.edge_long_info[{i, j}].aver_edge_lp;
        tmp.first += arc_graph_revised[i][j];
        ++tmp.second;
#endif
      }
    }
  }
}

int CVRP::optimizeLPForOneIteration(BbNode *node, double prior_value) {
  START_OVER:
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumCol(&num_col))
  safe_solver(node->solver.getObjVal(&lp_val))
  if (num_col - node->num_parent_cols > LIMIT_NEW_ADDED_COLUMNS && abs(prior_value - lp_val) > TOLERANCE)
    cleanIndexColForNode(node, node->num_parent_cols);
  safe_solver(node->solver.getObjVal(&lp_val))
  vector<double> X(num_col);
  safe_solver(node->solver.getX(0, num_col, X.data()))

  priceLabeling(node);

  node->is_integer = true;
  for (int i = 0; i < num_col; ++i) {
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      node->is_integer = false;
      break;
    }
  }
  if (node->is_integer) {
    if (ceilTransformedNumberRelated(lp_val - TOLERANCE) + TOLERANCE < ub) {
      ip_opt_sol.clear();
      int tmp_p_col = node->num_parent_cols;
      for (int i = 0; i < tmp_p_col; ++i) {
        if (X[i] > TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(dim);
          tmp.emplace_back(0);
          for (size_t j = node->index_columns[i] + 1;; ++j) {
            if (!col_pool4_mem[j])
              break;
            tmp.emplace_back(col_pool4_mem[j]);
          }
          tmp.emplace_back(0);
          ip_opt_sol.emplace_back(std::move(tmp));
        }
      }

      for (int i = tmp_p_col; i < num_col; ++i) {
        if (X[i] > TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(dim);
          tmp.emplace_back(0);
          for (size_t j = node->index_columns[i] + 1;; ++j) {
            if (!col_pool4_pricing[j])
              break;
            tmp.emplace_back(col_pool4_pricing[j]);
          }
          tmp.emplace_back(0);
          ip_opt_sol.emplace_back(std::move(tmp));
        }
      }

#ifdef SOLVER_VRPTW
	  bool if_feasible;
	  checkSolutionFeasibleByCapacity(if_feasible);
	  if (if_feasible){
		ub = ceilTransformedNumberRelated(lp_val - TOLERANCE);
	  }else{
		node->is_integer=false;
		ip_opt_sol.clear();
		/**
		 * add rcc cuts and they are non-removable before enumeration!
		 */
		if_force_keep_rcc=true;
		generateRCCs(node);
		if_force_keep_rcc=false;
		goto START_OVER;
	  }
#else
	  ub = ceilTransformedNumberRelated(lp_val - TOLERANCE);
#endif


#ifdef MASTER_VALVE_ML
      updateUpperBoundEdgeSolution();
#endif

      if (ceilTransformedNumberRelated(lb_transformed - TOLERANCE) + TOLERANCE >= ub) {
        node->value = lp_val;
        node->is_terminated = true;
        cout << TERMINATED_MESSAGE_PROMISING_UPDATE_UB;
        return 1;
      }
    }
  }

  return 0;
}

void CVRP::cleanIndexColForNode(BbNode *node, int beg, bool if_only_rcfixing) {
  if (beg >= num_col) return;
  safe_solver(node->solver.reoptimize())
  vector<double> rc(num_col);
  safe_solver(node->solver.getRC(0, num_col, rc.data()))
  vector<double> rc_copy(rc.data() + beg, rc.data()+ num_col);
  opt_gap = calculateOptimalGap(node);
  double threshold;
  if (if_only_rcfixing) {
    threshold = opt_gap;
  } else {
    int n = int((num_col - beg) * COL_KEEP_FRAC);
    nth_element(rc_copy.begin(), rc_copy.begin() + n, rc_copy.end());
    threshold = max(min(rc_copy[n], opt_gap), TOLERANCE);
  }
  vector<int> solver_ind(num_col - beg);
  int len = 0, keep = beg;
  for (int i = beg; i < num_col; ++i) {
    if (rc[i] > threshold) {
      solver_ind[len++] = i;
    } else {
      node->index_columns[keep++] = node->index_columns[i];
    }
  }
  safe_solver(node->solver.delVars(len, solver_ind.data()))
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumCol(&num_col))
}

void CVRP::getNewConstraintCoefficientByEdge(BbNode *node, std::pair<int, int> edge, int *c_ind, double *c_val, int &num_nz) {
  if (if_in_enu_state) {
    vector<int> cols(num_col, 0);
    for (auto it : map_edge_col_idx_in_enu[edge]) ++cols[it];//will be automatically populated!
    int ccnt = 0;
    int numnz = 0;
    for (auto it : cols) {
      if (it) {
        c_ind[numnz] = ccnt;
        c_val[numnz++] = it;
      }
      ++ccnt;
    }
    num_nz = numnz;
  } else {
    auto p = node->ptr;
    int code = edge.first * dim + edge.second;
    int cnt = 0;
    int times;
    int record;
    while (p) {
      if (!p->edge_to_cols[code].empty()) {
        times = 0;
        record = p->edge_to_cols[code][0];
        for (int i : p->edge_to_cols[code]) {
          if (record == i) {
            ++times;
          } else {
            c_ind[cnt] = record;
            if (record < num_arti_vars) {
              c_val[cnt++] = times - 1;
              safe_Hyperparameter(times - 2)
            } else
              c_val[cnt++] = times;
            record = i;
            times = 1;
          }
        }
        c_ind[cnt] = record;
        if (record < num_arti_vars) {
          c_val[cnt++] = 1;
          safe_Hyperparameter(times - 2)
        } else
          c_val[cnt++] = times;
      }
      p = p->p_node;
    }
    num_nz = cnt;
  }
}

void CVRP::getCoefficientRCC(BbNode *node, Rcc &rcc, int *c_ind, double *c_val, int &num_nz) const {
  unordered_map<int, vector<int>> seq_map;
  seq_map.reserve(dim * dim);
  int curr_node, next_node;
  vector<double> supp(num_col, 0);
  for (int n = node->num_parent_cols; n < num_col; ++n) {
    for (size_t m = node->index_columns[n];;) {
      curr_node = col_pool4_pricing[m];
      next_node = col_pool4_pricing[++m];
      seq_map[curr_node * dim + next_node].emplace_back(n);
      if (!next_node) break;
    }
  }

  if (rcc.form_rcc) {//normal case
    auto &customer_info = rcc.info_rcc_customer;
	vector<int> code( customer_info.size() * customer_info.size() / 2 + 1, 0);
    int cnt = 0;
    for (auto i = customer_info.begin(); i != customer_info.end(); ++i) {
      auto j = i;
      ++j;
      for (; j != customer_info.end(); ++j) {
        int tmp1 = (*i) * dim + (*j);
        int tmp2 = (*j) * dim + (*i);
        for (auto it : seq_map[tmp1])++supp[it];
        for (auto it : seq_map[tmp2])++supp[it];
        code[cnt++] = min(tmp1, tmp2);
      }
    }

    auto p = node->ptr->p_node;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->edge_to_cols[code[i]]) {
          ++supp[j];
        }
      }
      p = p->p_node;
    }
  } else {//another case
    auto &customer_info = rcc.info_rcc_customer;
    auto &outside_customer_info = rcc.info_rcc_outside_customer;
	vector<int> code( dim*dim, 0);
    int cnt = 0;
    for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
      auto j = i;
      ++j;
      for (; j != outside_customer_info.end(); ++j) {
        int tmp1 = (*i) * dim + (*j);
        int tmp2 = (*j) * dim + (*i);
        for (auto it : seq_map[tmp1])++supp[it];
        for (auto it : seq_map[tmp2])++supp[it];
        code[cnt++] = min(tmp1, tmp2);
      }
    }

    auto p = node->ptr->p_node;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->edge_to_cols[code[i]]) {
          ++supp[j];
        }
      }
      p = p->p_node;
    }

    cnt = 0;
    for (auto customer_it : outside_customer_info) {
      int tmp1 = customer_it;
      int tmp2 = customer_it * dim;
      for (auto it : seq_map[tmp1])supp[it] += 0.5;
      for (auto it : seq_map[tmp2])supp[it] += 0.5;
      code[cnt++] = min(tmp1, tmp2);
    }
    p = node->ptr->p_node;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->edge_to_cols[code[i]]) {
          supp[j] += 0.5;
        }
      }
      p = p->p_node;
    }

    cnt = 0;
    for (auto customer_it : customer_info) {
      int tmp1 = customer_it;
      int tmp2 = customer_it * dim;
      for (auto it : seq_map[tmp1])supp[it] -= 0.5;
      for (auto it : seq_map[tmp2])supp[it] -= 0.5;
      code[cnt++] = min(tmp1, tmp2);
    }
    p = node->ptr->p_node;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->edge_to_cols[code[i]]) {
          supp[j] -= 0.5;
        }
      }
      p = p->p_node;
    }
  }

	int cnt = 0;
	supp[0] = rcc.rhs;
	for (int i = 0; i < num_col; ++i) {
	  if (supp[i] != 0) {
		c_val[cnt] = supp[i];
		c_ind[cnt++] = i;
	  }
	}
	num_nz = cnt;
}

void CVRP::writeColumnToMemory(BbNode *node) {
  if (checkMemoryPool()) reallocateMemoryPool();
  for (int i = node->num_parent_cols; i < num_col; ++i) {
    size_t start = pool_beg4_mem;
    col_pool4_mem[pool_beg4_mem++] = 0;
    for (size_t j = node->index_columns[i] + 1;; ++j) {
      if (!(col_pool4_pricing[j]))break;
      col_pool4_mem[pool_beg4_mem++] = col_pool4_pricing[j];
    }
    col_pool4_mem[pool_beg4_mem++] = 0;
    node->index_columns[i] = start;
  }
  safe_solver(node->solver.updateModel())
}

void CVRP::constructMap(BbNode *node, int beg) const {
  auto edge_map = &(node->ptr->edge_to_cols);
  for (int i = beg; i < num_col; ++i) {
    for (size_t j = node->index_columns[i];; ++j) {
      int ai = col_pool4_mem[j], aj = col_pool4_mem[j + 1];
      if (ai > aj) {
        (*edge_map)[aj * dim + ai].emplace_back(i);
      } else {
        (*edge_map)[ai * dim + aj].emplace_back(i);
      }
      if (!aj) break;
    }
  }
  node->num_parent_cols = num_col;
}

void CVRP::reformulateIntegerProgramSolution(BbNode *node) {
  cout << "reformulateIntegerProgramSolution" << endl;
#ifdef SOLVER_VRPTW
  bool if_feasible;
  checkSolutionFeasibleByCapacity(if_feasible);
  if(!if_feasible) {
	cout<<"infeasible solution, done nothing, continue to solve..."<<endl;
	return;
  }
#endif
  node->is_integer = true;
  node->is_terminated = true;
  safe_solver(node->solver.getObjVal(&lp_val))
  if (ceilTransformedNumberRelated(lp_val - TOLERANCE) + TOLERANCE < ub) {
    ub = ceilTransformedNumberRelated(lp_val - TOLERANCE);
    ip_opt_sol.clear();
    for (int i = 0; i < node->num_parent_cols_in_lp_solutions; ++i) {
      if ((node->index_for_lp_solutions_in_column_pool[i].second - 1) < TOLERANCE) {
        vector<int> tmp;
        tmp.reserve(dim);
        tmp.emplace_back(0);
        for (auto j = node->index_for_lp_solutions_in_column_pool[i].first + 1;; ++j) {
          if (!col_pool4_mem[j]) break;
          tmp.emplace_back(col_pool4_mem[j]);
        }
        tmp.emplace_back(0);
        ip_opt_sol.emplace_back(std::move(tmp));
        for (auto j = node->index_for_lp_solutions_in_column_pool[i].first;;) {
          int p = col_pool4_mem[j], q = col_pool4_mem[j + 1];
          if (p > q)
            arc_graph_revised[q][p] = 0;
          else
            arc_graph_revised[p][q] = 0;
          if (!col_pool4_mem[++j]) break;
        }
      }
    }
    for (int i = node->num_parent_cols_in_lp_solutions; i < node->index_for_lp_solutions_in_column_pool.size(); ++i) {
      if ((node->index_for_lp_solutions_in_column_pool[i].second - 1) < TOLERANCE) {
        vector<int> tmp;
        for (auto j = node->index_for_lp_solutions_in_column_pool[i].first + 1;; ++j) {
          if (!col_pool4_pricing[j]) break;
          tmp.emplace_back(col_pool4_pricing[j]);
        }
        tmp.emplace_back(0);
        ip_opt_sol.emplace_back(std::move(tmp));
        for (auto j = node->index_for_lp_solutions_in_column_pool[i].first;;) {
          int p = col_pool4_pricing[j], q = col_pool4_pricing[j + 1];
          if (p > q)
            arc_graph_revised[q][p] = 0;
          else
            arc_graph_revised[p][q] = 0;
          if (!col_pool4_pricing[++j]) break;
        }
      }
    }
    for (int j = 1; j < dim; ++j) {
      if (arc_graph_revised[0][j] == 0) continue;
      vector<int> tmp;
      tmp.reserve(dim);
      tmp.emplace_back(0);
      tmp.emplace_back(j);
      arc_graph_revised[0][j] = 0;
      int tmp_j = j;
      while (true) {
        bool ifbreak = false;
        for (int k = 0; k < tmp_j; ++k) {
          if (arc_graph_revised[k][tmp_j] != 0) {
            arc_graph_revised[k][tmp_j] = 0;
            tmp.emplace_back(k);
            tmp_j = k;
            ifbreak = true;
            break;
          }
        }
        if (!ifbreak) {
          for (int k = tmp_j + 1; k < dim; ++k) {
            if (arc_graph_revised[tmp_j][k] != 0) {
              arc_graph_revised[tmp_j][k] = 0;
              tmp.emplace_back(k);
              tmp_j = k;
              break;
            }
          }
        }
        if (tmp_j == 0) {
          break;
        }
      }
      ip_opt_sol.emplace_back(std::move(tmp));
    }
#ifdef MASTER_VALVE_ML
    updateUpperBoundEdgeSolution();
#endif
  }
}


void CVRP::convertVertexToR1CsInPricing(BbNode *node) {
  Vertex2ActiveInOnePricingR1Cs.assign(dim, tuple<vector<int>, R1CMem, vector<int >>());
  int cnt = 0;
  unordered_map<int, double> cnt_dual;
  for (auto &r1c : node->r1cs) {
    if (abs(pi4_labeling[r1c.idx_r1c]) < DUAL_TOLERANCE) continue;
    cnt_dual.emplace(cnt, pi4_labeling[r1c.idx_r1c]);
    for (auto i : r1c.info_r1c) {
      get<0>(Vertex2ActiveInOnePricingR1Cs[i]).emplace_back(cnt);
      get<1>(Vertex2ActiveInOnePricingR1Cs[i]).set(cnt);
      get<2>(Vertex2ActiveInOnePricingR1Cs[i]).emplace_back(cnt);
    }
    for (auto &v : r1c.mem) {
      get<1>(Vertex2ActiveInOnePricingR1Cs[v]).set(cnt);
      get<2>(Vertex2ActiveInOnePricingR1Cs[v]).emplace_back(cnt);
    }
    ++cnt;
  }
  for (int i = 1; i < dim; ++i) {
    vector<pair<int, double>> tmp(get<2>(Vertex2ActiveInOnePricingR1Cs[i]).size());
    transform(get<2>(Vertex2ActiveInOnePricingR1Cs[i]).begin(),
              get<2>(Vertex2ActiveInOnePricingR1Cs[i]).end(), tmp.begin(),
              [&cnt_dual](const auto &a) { return make_pair(a, cnt_dual[a]); });
    sort(tmp.begin(), tmp.end(), [](const auto &a, const auto &b) {
      return a.second < b.second;
    });
    transform(tmp.begin(), tmp.end(), get<2>(Vertex2ActiveInOnePricingR1Cs[i]).begin(),
              [](const auto &a) { return a.first; });
  }

  cnt_dual.clear();
  cnt = 0;
  Vertex2ActiveInOnePricingR1C_multi.assign(dim, tuple<vector<tuple<int, int, int>>, vector<int>, vector<int>>());
  r1c_multi_denominator_in_cg.clear();
  for (auto &r1c : node->r1cs_multi) {
    if (abs(pi4_labeling[r1c.idx_r1c]) < DUAL_TOLERANCE) continue;
    const auto &plan = map_rank1_multiplier[(int) r1c.info_r1c.first.size()][r1c.info_r1c.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    r1c_multi_denominator_in_cg.emplace_back(denominator);
    int count = 0;
    for (auto &i : r1c.info_r1c.first) {
      get<0>(Vertex2ActiveInOnePricingR1C_multi[i]).emplace_back(cnt, multi[count], denominator);
      get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).emplace_back(cnt);
      ++count;
    }
    for (auto &v : r1c.mem) {
      get<2>(Vertex2ActiveInOnePricingR1C_multi[v]).emplace_back(cnt);
    }
    cnt_dual.emplace(cnt, pi4_labeling[r1c.idx_r1c]);
    ++cnt;
  }
  num_valid_r1c_multi_in_cg = cnt;
  vector<bool> if_use(num_valid_r1c_multi_in_cg);
  for (int i = 1; i < dim; ++i) {
	fill(if_use.begin(), if_use.end(), false);
    for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[i])) {
      if_use[l] = true;
    }
    for (int j = 0; j < num_valid_r1c_multi_in_cg; ++j) {
      if (!if_use[j]) get<1>(Vertex2ActiveInOnePricingR1C_multi[i]).emplace_back(j);
    }
  }

  for (int i = 1; i < dim; ++i) {
    vector<pair<int, double>> tmp(get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).size());
    transform(get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).begin(),
              get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).end(), tmp.begin(),
              [&cnt_dual](const auto &a) { return make_pair(a, cnt_dual[a]); });
    sort(tmp.begin(), tmp.end(), [](const auto &a, const auto &b) {
      return a.second < b.second;
    });
    transform(tmp.begin(), tmp.end(), get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).begin(),
              [](const auto &a) { return a.first; });
  }
}

void CVRP::convertVertexToR1CsInOneLP(BbNode *node) {
  vertex2_all_in_one_lp_r1cs.assign(dim, pair<vector<int>, R1CMem>());
  int cnt = 0;
  for (auto &r1c : node->r1cs) {
    for (auto i : r1c.info_r1c) {
      vertex2_all_in_one_lp_r1cs[i].first.emplace_back(cnt);
      vertex2_all_in_one_lp_r1cs[i].second.set(cnt);
    }
    for (auto &v : r1c.mem) {
      vertex2_all_in_one_lp_r1cs[v].second.set(cnt);
    }
    ++cnt;
  }
  cnt = 0;
  Vertex2AllInOneLPR1C_multi.assign(dim, pair<vector<pair<int, int>>, vector<int>>());
  r1c_multi_denominator_in_lp.clear();
  vector<vector<int>> tmp(dim, vector<int>());
  for (auto &r1c : node->r1cs_multi) {
    const auto &plan = map_rank1_multiplier[(int) r1c.info_r1c.first.size()][r1c.info_r1c.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    r1c_multi_denominator_in_lp.emplace_back(denominator);
    int count = 0;
    for (auto &i : r1c.info_r1c.first) {
      Vertex2AllInOneLPR1C_multi[i].first.emplace_back(cnt, multi[count]);
      tmp[i].emplace_back(cnt);
      ++count;
    }
    for (auto &v : r1c.mem) {
      tmp[v].emplace_back(cnt);
    }
    ++cnt;
  }
  auto if_use = new bool[cnt];
  for (int i = 1; i < dim; ++i) {
    memset(if_use, 0, sizeof(bool) * cnt);
    for (auto l : tmp[i]) {
      if_use[l] = true;
    }
    for (int j = 0; j < cnt; ++j) {
      if (!if_use[j]) Vertex2AllInOneLPR1C_multi[i].second.emplace_back(j);
    }
  }
  delete[] if_use;
}

void CVRP::cleanAllPointers(BbNode *const node, int mode, bool if_clear_concatenate) {
  if (mode == 1) {
    for (int i = 0; i < dim; ++i) {
      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        label_array_in_forward_sense[i][b].second = 0;
        if_exist_extra_labels_in_forward_sense[i][b].second = 0;
      }
    }
  }
#ifdef SYMMETRY_PROHIBIT
  else if (mode == 2) {
    for (int i = 0; i < dim; ++i) {
      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        label_array_in_backward_sense[i][b].second = 0;
        if_exist_extra_labels_in_backward_sense[i][b].second = 0;
      }
    }
  }
#endif
  if (if_clear_concatenate) {
    concatenate_labels_in_forward_cg.clear();
#ifdef SYMMETRY_PROHIBIT
    concatenate_labels_in_backward_cg.clear();
#endif
  }
}

double CVRP::calculateOptimalGap(BbNode *const node) const {
  return (ub - node->value + round_up_tolerance);
}

void CVRP::determineIfArcElimination(BbNode *node) {
  double now_gap, required_gap;
  if (!if_can_arc_elimination_by_exact_cg) {
    final_decision_4_arc_elimination = false;
    cerr << "Error: if_can_arc_elimination_by_exact_cg is false, but tries to call arcElimination!" << endl;
    exit(-1);
  } else if_can_arc_elimination_by_exact_cg = false;
  if (if_stop_arc_elimination) {
    final_decision_4_arc_elimination = false;
    if_stop_arc_elimination = false;
    cout << "ArcEliminationNEnumeration is banned!" << endl;
    goto QUIT;
  }
  if (final_decision_4_arc_elimination) {
    cout << "ArcEliminationNEnumeration must be performed!" << endl;
    goto QUIT;
  }
  if (abs(old_ub - ub) < TOLERANCE) {
    cout << "old_ub is updated from " << old_ub << " to " << ub << endl;
    cout << "ArcEliminationNEnumeration must be performed!" << endl;
    old_ub = ub;
    final_decision_4_arc_elimination = true;
    goto QUIT;
  }
  now_gap = (ub - node->value) / ub;
  if (node->last_gap / now_gap > Config::GapImproved4ArcEliminationNEnumeration
      || now_gap < gap_tolerance4_arc_elimination_n_enumeration) {
    node->last_gap = now_gap;
    cout << "gap updated as " << now_gap << endl;
    final_decision_4_arc_elimination = true;
    final_decision_4_enumeration = true;
    goto QUIT;
  }
  QUIT:
  return;
}

void CVRP::determineIfEnumeration(BbNode *node) {
  double gap;
  if (final_decision_4_enumeration) {
    cout << "Enumeration must be performed!" << endl;
    enumeration_mode = !if_arc_elimination_succeed;
    if (!if_arc_elimination_succeed) {
      cout << "But ArcElimination is failed! So skip enumeration!" << endl;
      final_decision_4_enumeration = false;
    }
    goto QUIT;
  }
  if (!if_arc_elimination_succeed && !if_arc_elimination_tried_but_failed) {
    goto QUIT;
  }

  if (abs(old_ub - ub) < TOLERANCE) {
    cout << "old_ub is updated from " << old_ub << " to " << ub << endl;
    cout << "ArcEliminationNEnumeration must be performed!" << endl;
    old_ub = ub;
    final_decision_4_enumeration = true;
    enumeration_mode = !if_arc_elimination_succeed;
    goto QUIT;
  }

  gap = (ub - node->value) / ub;
  if (gap > Config::MaxGap2TryEnumeration || gap > last_enumeration_fail_gap) {
    goto QUIT;
  }
  final_decision_4_enumeration = true;
  enumeration_mode = !if_arc_elimination_succeed;
  if (if_arc_elimination_tried_but_failed) {
    if (count4_tolerance4_try_enumeration_when_arc_elimination_fails
        <= Config::InitialTolerance4tryEnumerationWhenArcEliminationFails) {
      ++count4_tolerance4_try_enumeration_when_arc_elimination_fails;//if succeed in enumeration, this will be reset!
    } else {
      final_decision_4_enumeration = false;
    }
  }
  QUIT:
  return;
}

void CVRP::rollbackEasyWay(BbNode *const node, int old_num) {
  if (if_in_enu_state) throw runtime_error("Error: rollbackEasyWay is called in Enu_State!");
#if SETTING==2
  throw runtime_error("Error: CutRollBackNotAllowed is defined, but rollbackToPreviousState is called!");
#endif
  cout << "rollbackEasyWay" << endl;
  hard_rollback_factor = Config::HardTimeThresholdInPricing / last_max_time_labeling * 1.1;
  cut_gen_time_threshold_in_pricing *= Config::CutGenTimeThresholdInPricingReduceFactor;
  cout << "hard_rollback_factor: " << hard_rollback_factor << " cut_gen_time_threshold_in_pricing: "
       << cut_gen_time_threshold_in_pricing << endl;
  vector<int> delete_cuts(num_row - old_num);
  iota(delete_cuts.begin(), delete_cuts.end(), old_num);
  deleteNonactiveCuts(node, delete_cuts);

  for (auto &i : reset_cut_mem) {
    auto &d = get<2>(i);
    std::vector<int> tmp(d.count());
    int idx = 0;
    for (int j = 1; j < dim; ++j) {
      if (d.test(j)) tmp[idx++] = j;
    }
    if (get<0>(i)) {
      node->r1cs[get<1>(i)].mem = std::move(tmp);
    } else {
      node->r1cs_multi[get<1>(i)].mem = std::move(tmp);
    }
  }

  reset_cut_mem.clear();

  convertVertexToR1CsInOneLP(node);
  safe_solver(node->solver.optimize())
  rollback = 0;
  solveLPInLabeling(node);
  safe_Hyperparameter(rollback == 1)
}

void CVRP::initializeBucketGraphForNode(BbNode *node) {
  vector<int> bucketArc(real_dim - 1);
  for (int i = 1; i < dim; ++i) {
    int tmp = 0;
    for (int j = 1; j < i; ++j) bucketArc[tmp++] = j;
    for (int j = i + 1; j < dim; ++j) bucketArc[tmp++] = j;
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      auto &bucket = node->all_forward_buckets[i][b];
      bucket.bucket_arcs = bucketArc;
      bucket.i = i;
    }
  }
  auto &bucket = node->all_forward_buckets[0][0];
  bucket.bucket_arcs.resize(real_dim);
  iota(bucket.bucket_arcs.begin(), bucket.bucket_arcs.end(), 1);

  max_num_forward_graph_arc = num_buckets_per_vertex * (real_dim - 1) * real_dim;
  node->num_forward_bucket_arcs = max_num_forward_graph_arc;
#ifdef SYMMETRY_PROHIBIT
  for (int i = 1; i < dim; ++i) {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      node->all_backward_buckets[i][b] = node->all_forward_buckets[i][0];
    }
  }
  auto &bucket_back = node->all_backward_buckets[0][0];
    bucket_back=bucket;
  max_num_backward_graph_arc = num_buckets_per_vertex * (real_dim - 1) * real_dim;
  node->num_backward_bucket_arcs = max_num_backward_graph_arc;
#endif
}

double CVRP::ceilTransformedNumberRelated(double x) const {
  for (int i = 0; i < transformed_number; ++i) {
    x *= 10;
  }
  x = std::ceil(x);
  for (int i = 0; i < transformed_number; ++i) {
    x /= 10;
  }
  return x;
}

void splitString(const string &str, const string &splits, vector<string> &res) {
  if (str.empty()) return;
  auto strs = str + splits;
  size_t pos = strs.find(splits);
  int step = (int) splits.size();
  while (pos != std::string::npos) {
    string tmp = strs.substr(0, pos);
    res.push_back(tmp);
    strs = strs.substr(pos + step, strs.size());
    pos = strs.find(splits);
  }
}

string readInstanceFile(const std::string &file_name, int line) {
  safe_Hyperparameter(line == READ_NO_LINE)
  fstream file(file_name);
  if (!file.is_open()) {
    cout << "File not found!" << endl;
    exit(EXIT_FAILURE);
  }
  file.seekg(std::ios::beg);
  for (int i = 0; i < line; ++i) {
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  string line_str;
  getline(file, line_str);
  vector<string> strList;
  splitString(line_str, " ", strList);
  if (strList.size() == 2) {
    Config::ub = strtod(strList[1].c_str(), nullptr);
  }
#ifdef READ_ENUMERATION_TREES
  else if (strList.size() == 4) {
    Config::ub = strtod(strList[1].c_str(), nullptr);
    Config::col_pool_path = strList[2].c_str();
    Config::tree_path = strList[3].c_str();
  }
#endif
  file.close();
  return strList[0];
}

string generateInstancePath(int argc, char *argv[]) {
  int n = READ_NO_LINE;
  string d;
  bool if_type1 = false;
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-d" && i + 1 < argc) {
      d = argv[++i];
      if_type1 = true;
    } else if (arg == "-n" && i + 1 < argc) {
      stringstream convert(argv[++i]);
      if (!(convert >> n)) {
        printf("Invalid number: %s\n", argv[i]);
        continue;
      }
      if_type1 = true;
    } else if (arg == "-u" && i + 1 < argc) {
      stringstream convert(argv[++i]);
      if (!(convert >> Config::ub)) {
        printf("Invalid number: %s\n", argv[i]);
        continue;
      }
    }
#ifdef READ_ENUMERATION_TREES
      else if (arg == "-c" && i + 1 < argc) {
        stringstream convert(argv[++i]);
        if (!(convert >> Config::col_pool_path)) {
          printf("Invalid string: %s\n", argv[i]);
          continue;
        }
      } else if (arg == "-t" && i + 1 < argc) {
        stringstream convert(argv[++i]);
        if (!(convert >> Config::tree_path)) {
          printf("Invalid string: %s\n", argv[i]);
          continue;
        }
      }
#endif
    else {
      printf("Unknown option: %s\n", argv[i]);
    }
  }

  if (if_type1) {
    return readInstanceFile(d, n);
  } else {
    if (argc > 1) {
      return argv[1];
    } else {
      throw std::invalid_argument("Not enough arguments");
    }
  }
}

int CVRP::inverseLastBranchConstraint(char sense, double rhs, Solver &solver) {
  int error = solver.updateModel();
  if (error) return error;
  int idxconstrs;
  error += solver.getNumRow(&idxconstrs);
  if (error) return error;
  --idxconstrs;
  error += solver.setRhs(idxconstrs, 1, &sense, &rhs);
  if (error) return error;
  int vind = 0;
  error += solver.XchangeCoeffs(1, &idxconstrs, &vind, &rhs);
  return error;
}

int CVRP::addBranchConstraint(int numnz,
                          int *cind,
                          double *cval,
                          char sense,
                          double rhs,
                          const char *constrname,
                          Solver &solver) {
  int error = solver.addConstraint(numnz, cind, cval, sense, rhs, constrname);
  if (error) return error;
  error += solver.updateModel();
  if (error) return error;
  int idxconstrs;
  error += solver.getNumRow(&idxconstrs);
  if (error) return error;
  --idxconstrs;
  int vind = 0;
  error += solver.XchangeCoeffs(1, &idxconstrs, &vind, &rhs);
  return error;
}

void CVRP::changeBranchConstraint(double *val,
                           int *cind,
                           int *vind,
                           int numnz,
                           const int *old_ind,
                           const double *old_val,
                           char sense,
                           double rhs,
                           Solver &local_solver) {
  memset(val, 0, sizeof(double) * num_col);
  for (int i = 0; i < numnz; ++i) {
    val[old_ind[i]] = old_val[i];
  }
  safe_solver(local_solver.changeCoeffs(num_col, cind, vind, val))
  int BeforeNumRow = cind[0];
  safe_solver(local_solver.setRhs(BeforeNumRow, 1, &sense, &rhs))
  int vind2 = 0;
  safe_solver(local_solver.XchangeCoeffs(1, &BeforeNumRow, &vind2, &rhs))
}

float sqrt_self(float x) {
  float xhalf = 0.5f * x;
  int i = *(int *) &x; // get bits for floating VALUE
  i = 0x5f375a86 - (i >> 1); // gives initial guess y0
  x = *(float *) &i; // convert bits BACK to float
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  return 1 / x;
}

double pow_self(double x, int n) {
  if (n == 0) return 1;
  double y = x;
  for (int i = 1; i < n; ++i) {
    y *= x;
  }
  return y;
}


bool CVRP::increaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end) {
  newMainResource = nowMainResource + main_resource_across_arcs_in_forward_sense[start][end];
#ifdef CAPACITY_AS_MAIN_RESOURCE
#else
  if (newMainResource > ub4_vertex[end]) return false;
  newMainResource = newMainResource > lb4_vertex[end] ? newMainResource : lb4_vertex[end];
#endif
  if (newMainResource + main_resource_across_arcs_in_forward_sense[end][0] > max_main_resource) return false;
  return true;
}

bool CVRP::decreaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end) {
  newMainResource = nowMainResource - main_resource_across_arcs_in_backward_sense[start][end];
#ifdef CAPACITY_AS_MAIN_RESOURCE
#else
  if (newMainResource < lb4_vertex[end]) return false;
  newMainResource = newMainResource < ub4_vertex[end] ? newMainResource : ub4_vertex[end];
#endif
  if (newMainResource - main_resource_across_arcs_in_backward_sense[end][0] < 0) return false;
  return true;
}

bool CVRP::runColumnGenerationType(BbNode *node, int mode) {
  if (node->is_terminated) return true;
  int status = 0;
  bool switch_is_on = true;
  switch (mode) {
    case 1:cout << "LighterHeur phase begin...\n";
      break;
    case 2:cout << "HeavierHeur phase begin...\n";
      break;
    case 3:cout << "Exact phase begin...\n";
      break;
    default:cerr << "None of these modes are used!\n";
      exit(-1);
  }
  int iter = 0, ccnt = 0, old_ncol = num_col, tag;
  double eps_CG = 0, eps_LP = 0, b4_node_val = node->value;
  auto beg = chrono::high_resolution_clock::now();
  auto end = beg;
  bool if_tail_off = false;
  bool old_force;
  bool if_gap_0 = false;
  if (mode == 3) {
    if_exact_labeling_cg = true;
    old_force = force_not_rollback;
    if_ban_convert_dual = false;
  }

  while (true) {

    ++iter;

    beg = chrono::high_resolution_clock::now();

    switch (mode) {
      case 1:ccnt = generateColumnsByLighterHeuristic(node);
        break;
      case 2:ccnt = generateColumnsByHeavierHeuristic(node);
        break;
      case 3:
        if (abs(gap_between_last_smallest_rc_and_rc_threshold) < TOLERANCE) {
          if_gap_0 = true;
        }
        ccnt = generateColsByBidir<
#ifdef SYMMETRY_PROHIBIT
            false
#else
            true
#endif
        >(node);
        if (if_gap_0) {
          force_not_rollback = true;
        }
        if (rollback == 3) {
          if_tail_off = true;
        }

        break;
      default:cerr << "None of these modes are used!\n";
        exit(-1);
    }

    end = chrono::high_resolution_clock::now();
    time_pricing4_iter = chrono::duration<double>(end - beg).count();
    eps_CG += time_pricing4_iter;

    if (!ccnt) {
      if (mode != 3 || rollback == 1 || if_gap_0) {
        --iter;
        break;
      }
    }
    if (node->is_terminated) {
      --iter;
      break;
    }

    beg = chrono::high_resolution_clock::now();

    tag = optimizeLPForOneIteration(node, b4_node_val);
    b4_node_val = lp_val;

    end = chrono::high_resolution_clock::now();
    time_resolve_lp4_iter = chrono::duration<double>(end - beg).count();
    eps_LP += time_resolve_lp4_iter;

    if (tag) {
      status = 1;//QUIT
      goto QUIT;
    }

    if (!(iter % PRINT_LABELING_STEP_SIZE)) {
      REPORT:
      glo_end = chrono::high_resolution_clock::now();
      glo_eps = chrono::duration<double>(glo_end - glo_beg).count();
      printInfoLabeling(iter, num_col - old_ncol, num_col, num_row, eps_LP,
                        eps_CG, glo_eps,
                        lp_val, lb, ub);
      if (!switch_is_on) {
        goto QUIT;
      }
      eps_CG = 0;
      eps_LP = 0;
      old_ncol = num_col;
    }
  }
  if (!iter || (switch_is_on && iter % PRINT_LABELING_STEP_SIZE)) {
    switch_is_on = false;
    goto REPORT;
  }
  QUIT:
  if (mode == 3) {
    if_exact_labeling_cg = false;
    force_not_rollback = old_force;
    if (if_tail_off && !rollback)
      rollback = 3;
  }
  return status;
}

double CVRP::transformCost(double x) {
  return std::floor(x + 0.5);
}

void CVRP::setTailOffStandardAndRollBackStandard() const {
  for (int b = 0; b < num_buckets_per_vertex; ++b) {
    int min_num_labels = MAX_INT;
    for (int i = 0; i < dim; ++i) {
      if (label_array_in_forward_sense[i][b].second && min_num_labels > label_array_in_forward_sense[i][b].second) {
        min_num_labels = label_array_in_forward_sense[i][b].second;
      }
    }
    if (min_num_labels == MAX_INT) {
      min_num_labels = 1;
    }
    for (int i = 0; i < dim; ++i) {
      int hard_max = max(label_array_in_forward_sense[i][b].second, min_num_labels) * FACTOR_NUM_LABEL;
      hard_max = int(pow_self(2.0, int(ceil(log(hard_max) / log(2)) + TOLERANCE)));
      label_array_in_forward_sense[i][b].first.resize(hard_max);
      if_exist_extra_labels_in_forward_sense[i][b].first.resize(hard_max);
    }
#ifdef SYMMETRY_PROHIBIT
    min_num_labels = MAX_INT;
    for (int i = 0; i < dim; ++i) {
      if (label_array_in_backward_sense[i][b].second && min_num_labels > label_array_in_backward_sense[i][b].second) {
        min_num_labels = label_array_in_backward_sense[i][b].second;
      }
    }
    if (min_num_labels == MAX_INT) {
      min_num_labels = 1;
    }
    for (int i = 0; i < dim; ++i) {
      int hard_max = max(label_array_in_backward_sense[i][b].second, min_num_labels) * FACTOR_NUM_LABEL;
      hard_max = (int) (pow_self(2.0, (int) ceil(log(hard_max) / log(2))) + TOLERANCE);
      label_array_in_backward_sense[i][b].first.resize(hard_max);
      if_exist_extra_labels_in_backward_sense[i][b].first.resize(hard_max);
    }
#endif
  }
}

void CVRP::getEdgeInfo(BbNode *node, bool if_br) const {
  if (node->is_integer) return;
  double **local_arc_graph;
  if (if_br) local_arc_graph = arc_graph_revised;
  else local_arc_graph = arc_graph;
  int numEdge=num_edge+1;//as required
  if (node->edge_head.size() < numEdge) {
	node->allocateMem(numEdge);
  }
  node->num_edges = 0;
  for (int i = 0; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      if (local_arc_graph[i][j] > TOLERANCE) {
        ++node->num_edges;
        node->edge_tail[node->num_edges] = i;
        node->edge_head[node->num_edges] = j;
        node->edge_value[node->num_edges] = local_arc_graph[i][j];
      }
    }
  }
}


#ifdef DELUXING_APPLIED
void CVRP::applyRCF(BbNode *node, int round, bool if_verbose) {
  cout << "WARNING: RCF can only be applied by GRB model!" << endl;
  throw runtime_error("Error: RCF is not seen by far!");
}
#endif

void CVRP::regenerateGraphBucket(BbNode *node) {
  if (abs(step_size / 2) < 1 - TOLERANCE) return;
  if_stop_arc_elimination = true;
  step_size /= 2;
  num_buckets_per_vertex *= 2;
  max_num_forward_graph_arc *= 2;
  node->num_forward_bucket_arcs *= 2;
  node->num_forward_jump_arcs *= 2;
#ifdef SYMMETRY_PROHIBIT
  node->num_backward_bucket_arcs *= 2;
  node->num_backward_jump_arcs *= 2;
  max_num_backward_graph_arc *= 2;
#endif
  auto new_LabelArrayInForwardSense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
    new_LabelArrayInForwardSense[i] = new VecLabel[num_buckets_per_vertex];
    for (int j = 0; j < num_buckets_per_vertex; ++j) {
      new_LabelArrayInForwardSense[i][j].first.resize(label_array_in_forward_sense[i][j / 2].first.size() / 2);
    }
  }

  auto new_IfExistExtraLabelsInForwardSense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
    new_IfExistExtraLabelsInForwardSense[i] = new VecLabel[num_buckets_per_vertex];
    for (int j = 0; j < num_buckets_per_vertex; ++j) {
      new_IfExistExtraLabelsInForwardSense[i][j].first.resize(
          label_array_in_forward_sense[i][j / 2].first.size() / 2);
    }
  }

  for (int i = 0; i < dim; ++i) {
    delete[]rc2_till_this_bin_in_forward_sense[i];
    rc2_till_this_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
  }

  for (int i = 0; i < dim; ++i) {
    delete[]rc2_bin_in_forward_sense[i];
    rc2_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
  }
#ifdef SYMMETRY_PROHIBIT
  auto new_LabelArrayInBackwardSense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
    new_LabelArrayInBackwardSense[i] = new VecLabel[num_buckets_per_vertex];
    for (int j = 0; j < num_buckets_per_vertex; ++j) {
      new_LabelArrayInBackwardSense[i][j].first.resize(label_array_in_backward_sense[i][j / 2].first.size() / 2);
    }
  }

  auto new_IfExistExtraLabelsInBackwardSense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
    new_IfExistExtraLabelsInBackwardSense[i] = new VecLabel[num_buckets_per_vertex];
    for (int j = 0; j < num_buckets_per_vertex; ++j) {
      new_IfExistExtraLabelsInBackwardSense[i][j].first.resize(label_array_in_backward_sense[i][j / 2].first.size() / 2);
    }
  }
  for (int i = 0; i < dim; ++i) {
    delete[]rc2_till_this_bin_in_backward_sense[i];
    rc2_till_this_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
  }
  for (int i = 0; i < dim; ++i) {
    delete[]rc2_bin_in_backward_sense[i];
    rc2_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
  }
#endif
  for (int i = 0; i < dim; ++i) {
    delete[]label_array_in_forward_sense[i];
  }
  delete[]label_array_in_forward_sense;
  label_array_in_forward_sense = new_LabelArrayInForwardSense;

  for (int i = 0; i < dim; ++i) {
    delete[]if_exist_extra_labels_in_forward_sense[i];
  }
  delete[]if_exist_extra_labels_in_forward_sense;
  if_exist_extra_labels_in_forward_sense = new_IfExistExtraLabelsInForwardSense;

#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < dim; ++i) {
    delete[]label_array_in_backward_sense[i];
  }
  delete[]label_array_in_backward_sense;
  label_array_in_backward_sense = new_LabelArrayInBackwardSense;
  for (int i = 0; i < dim; ++i) {
    delete[]if_exist_extra_labels_in_backward_sense[i];
  }
  delete[]if_exist_extra_labels_in_backward_sense;
  if_exist_extra_labels_in_backward_sense = new_IfExistExtraLabelsInBackwardSense;
#endif

  auto new_AllForwardBuckets = new Bucket *[dim];
  for (int i = 0; i < dim; ++i) {
    new_AllForwardBuckets[i] = new Bucket[num_buckets_per_vertex];
    for (int j = 0; j < num_buckets_per_vertex; ++j) {
      new_AllForwardBuckets[i][j] = node->all_forward_buckets[i][j / 2];
    }
  }
#ifdef SYMMETRY_PROHIBIT
  auto new_AllBackwardBuckets = new Bucket *[dim];
  for (int i = 0; i < dim; ++i) {
    new_AllBackwardBuckets[i] = new Bucket[num_buckets_per_vertex];
    for (int j = 0; j < num_buckets_per_vertex; ++j) {
      new_AllBackwardBuckets[i][j] = node->all_backward_buckets[i][j / 2];
    }
  }
#endif
  for (int i = 0; i < dim; ++i) {
    delete[]node->all_forward_buckets[i];
  }
  delete[]node->all_forward_buckets;
  node->all_forward_buckets = new_AllForwardBuckets;
#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < dim; ++i) {
    delete[]node->all_backward_buckets[i];
  }
  delete[]node->all_backward_buckets;
  node->all_backward_buckets = new_AllBackwardBuckets;
#endif
  tell_which_bin4_arc_elimination_in_forward_sense.clear();
  populateTellWhichBin4ArcElimination<true>();

#ifdef SYMMETRY_PROHIBIT
  tell_which_bin4_arc_elimination_in_backward_sense.clear();
  populateTellWhichBin4ArcElimination<false>();
#endif

  cout << "new generated bucket graph: num_buckets_per_vertex= " << num_buckets_per_vertex << " step_size= " << step_size
       << endl;
  cout << "we cannot use arc elimination next!" << endl;
}

void CVRP::getLowerBoundofMinimumNumberCars() {
  double sum_demand = accumulate(demand + 1, demand + dim, 0.0);
  int cap_k = ceil(sum_demand / cap);
  num_vehicle = cap_k;
  cout << "LBNumVehicle= " << num_vehicle << endl;
}

void CVRP::augmentNonGuillotineRound(BbNode *node) {
  if (node->index) return;
  vector<yzzLong> old_mem = ng_mem4_vertex;
  for (int i = 1; i < dim; ++i) {
    set<int> tmp;
    set<int> tmp2;
    for (int j = 1; j < i; ++j) {
      if (arc_graph_revised[j][i] > TOLERANCE) {
        tmp.emplace(j);
      }
    }
    for (int j = i + 1; j < dim; ++j) {
      if (arc_graph_revised[i][j] > TOLERANCE) {
        tmp.emplace(j);
      }
    }
    for (int j = 1; j < dim; ++j) {
      if (ng_mem4_vertex[i][j]) tmp2.emplace(j);
    }
    set<int> tmp3;
    set_intersection(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end(), inserter(tmp3, tmp3.begin()));
    ng_mem4_vertex[i] = 0;
    ng_mem4_vertex[i].set(i);
    for (auto n : tmp3) {
      ng_mem4_vertex[i].set(n);
    }
  }
  auto last = last_max_time_labeling;
  int round = 0;
  bool if_empty;
  double old_val = node->value;
  double cost_time = 0;
  while (true) {
    ++round;
    cout << "DSSR: " << round << endl;
    findNonGuillotineMemorySets(node, if_empty);//do not test if empty
    auto beg = chrono::high_resolution_clock::now();
    force_not_rollback = true;
    solveLPInLabeling(node);
    force_not_rollback = false;
    auto end = chrono::high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    if (eps > cost_time) cost_time = eps;
    if (last_max_time_labeling > last and node->value + TOLERANCE < old_val) {
      cout << "DSSR: Reach the hard limit! We use ng rollback!" << endl;
      ng_mem4_vertex = old_mem;
      deleteColumnByNGMemory(node);
      force_not_rollback = true;
      solveLPInLabeling(node);
      force_not_rollback = false;
      return;
    }
    old_mem = ng_mem4_vertex;
    if (node->value + TOLERANCE > old_val) break;
  }
  round = 0;
  old_val = node->value;
  double std = TOLERANCE;
  while (true) {
    ++round;
    cout << "NGAugmentation: " << round << endl;
    findNonGuillotineMemorySets(node, if_empty);
    if (if_empty) break;
    auto beg = chrono::high_resolution_clock::now();
    solveLPInLabeling(node);
    auto end = chrono::high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    if (eps > Config::NGAugTimeHardThresholdFactor * cost_time) {
      rollback = 1;
    } else if (eps > Config::NGAugTimeSoftThresholdFactor * cost_time) {
      rollback = 3;
    }
    cout << "eps= " << eps << endl;
    cout << "cost_time= " << cost_time << endl;
    if (rollback == 3) break;
    else if (rollback == 1) {
      cout << "NGAugmentation: Reach the hard limit! We use ng rollback!" << endl;
      ng_mem4_vertex = old_mem;
      deleteColumnByNGMemory(node);
      force_not_rollback = true;
      solveLPInLabeling(node);
      force_not_rollback = false;
      return;
    }
    old_mem = ng_mem4_vertex;
    double tmp = abs(node->value - old_val) / node->value;
    if (tmp > std) std = tmp;
    else if (tmp < std * Config::NGAugTailOff) break;
    else {
      old_val = node->value;
    }
  }
}

void CVRP::findNonGuillotineMemorySets(BbNode *node, bool &if_empty) {
  if (node->index) {
    cerr << "findNonGuillotineMemorySets: node->index!=0" << endl;
    exit(1);
  }
  unordered_map<int, vector<pair<int, yzzLong>>> map_size_cycle;
  vector<size_t> tmp(dim);
  vector<pair<size_t, double>>
      OptCols(node->index_for_lp_solutions_in_column_pool.begin() + node->num_parent_cols_in_lp_solutions, node->index_for_lp_solutions_in_column_pool.end());
  sort(OptCols.begin(), OptCols.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
    return a.second > b.second;
  });
  vector<size_t> OptColIdx(OptCols.size());
  transform(OptCols.begin(),
            OptCols.end(),
            OptColIdx.begin(),
            [](const pair<int, double> &a) { return a.first; });

  int re_cols = 0;
  for (auto idx : OptColIdx) {
    memset(tmp.data(), 0, sizeof(size_t) * dim);
    bool if_re_col = false;
    for (auto j = idx + 1;; ++j) {
      int curr_node = col_pool4_pricing[j];
      if (!curr_node) break;
      if (tmp[curr_node]) {
        if_re_col = true;
        auto length = int(j - tmp[curr_node] - 1);
        map_size_cycle[length].emplace_back(curr_node, 0);
        auto &mem = map_size_cycle[length].back().second;
        for (auto k = tmp[curr_node] + 1; k < j; ++k) {
          mem.set(col_pool4_pricing[k]);
        }
      }
      tmp[curr_node] = j;//not in else!
    }
    if (if_re_col) {
      ++re_cols;
      if (re_cols > Config::MaxNumColsInNGAug) break;
    }
  }
  vector<pair<int, vector<pair<int, yzzLong>>>> size_cycle(map_size_cycle.begin(), map_size_cycle.end());
  sort(size_cycle.begin(), size_cycle.end(), [](const pair<int, vector<pair<int, yzzLong>>> &a,
                                                const pair<int, vector<pair<int, yzzLong>>> &b) {
    return a.first < b.first;
  });
  if (size_cycle.empty()) {
    cout << "size_cycle is empty!" << endl;
    if_empty = true;
    return;
  } else if_empty = false;
  for (auto &i : size_cycle[0].second) {
    int re = i.first;
    for (int j = 1; j < dim; ++j) {
      if (i.second[j]) {
        auto &mem = ng_mem4_vertex[j];
        mem.set(re);
        if (mem.count() == Config::MaxNGSize) break;
      }
    }
  }
  for (int i = 1; i < size_cycle.size(); ++i) {
    if (size_cycle[i].first > Config::CycleSize) break;
    for (auto &j : size_cycle[i].second) {
      int re = j.first;
      for (int k = 1; k < dim; ++k) {
        if (j.second[k]) {
          auto &mem = ng_mem4_vertex[k];
          mem.set(re);
          if (mem.count() == Config::MaxNGSize) break;
        }
      }
    }
  }
  deleteColumnByNGMemory(node);
}

void CVRP::deleteColumnByNGMemory(BbNode *node) {
  int len = 0, keep = node->num_parent_cols, current_node;
  bool if_break;
  yzzLong local_pi;
  vector<int> solver_ind(num_col);
  for (int i = keep; i < num_col; ++i) {
    if_break = false;
	local_pi = 0;
    for (auto j = node->index_columns[i] + 1;; ++j) {
      current_node = col_pool4_pricing[j];
      if (!current_node) break;
      if (local_pi[current_node]) {
        solver_ind[len++] = i;
        if_break = true;
        break;
      }
	  local_pi = local_pi & ng_mem4_vertex[current_node];
	  local_pi.set(current_node);
    }
    if (!if_break) {
      node->index_columns[keep++] = node->index_columns[i];
    }
  }

  if (len) {
    safe_solver(node->solver.delVars(len, solver_ind.data()))
    safe_solver(node->solver.reoptimize())
    safe_solver(node->solver.getNumCol(&num_col))
  }
}

double CVRP::calculateGapImprovement(double nowVal, double b4Val) const {
  double guess_UB = ub > 2 * b4Val ? b4Val * 1.008 : ub;
  return (nowVal - b4Val) / (guess_UB - b4Val);
}

void CVRP::deleteArcByFalseBranchConstraint(BbNode *node, std::pair<int, int> edge) const {
  int i = edge.first, j = edge.second;
  bool if_rep = false;
  AGAIN:
  for (int bin = 0; bin < num_buckets_per_vertex; ++bin) {
    auto if_find = std::find(node->all_forward_buckets[i][bin].bucket_arcs.begin(),
                             node->all_forward_buckets[i][bin].bucket_arcs.end(), j);
    if (if_find == node->all_forward_buckets[i][bin].bucket_arcs.end()) {
      auto iff = std::find_if(node->all_forward_buckets[i][bin].jump_arcs.begin(),
                              node->all_forward_buckets[i][bin].jump_arcs.end(),
                              [&](const pair<double, int> &p) { return p.second == j; });
      if (iff != node->all_forward_buckets[i][bin].jump_arcs.end()) {
        node->all_forward_buckets[i][bin].jump_arcs.erase(iff);
      }
    } else {
      node->all_forward_buckets[i][bin].bucket_arcs.erase(if_find);
    }
  }
  if (!if_rep) {
    if_rep = true;
    std::swap(i, j);
    goto AGAIN;
  }
}

void CVRP::checkIfCutsLegal(BbNode *node) const {
  vector<bool> if_takeup(num_row, false);
  for (int i = 0; i < dim; ++i)if_takeup[i] = true;
  for (auto &rcc : node->rccs) {
    if (if_takeup[rcc.idx_rcc]) throw std::runtime_error("Error in rccs");
    if_takeup[rcc.idx_rcc] = true;
  }
  for (auto &r1c : node->r1cs) {
    if (if_takeup[r1c.idx_r1c]) throw std::runtime_error("Error in r1cs");
    if_takeup[r1c.idx_r1c] = true;
  }
  for (auto &r1c_multi : node->r1cs_multi) {
    if (if_takeup[r1c_multi.idx_r1c]) throw std::runtime_error("Error in R1CMulti");
    if_takeup[r1c_multi.idx_r1c] = true;
  }
  for (auto &brc : node->brcs) {
    if (brc.idx_br_c != -1) {
      if (if_takeup[brc.idx_br_c])throw std::runtime_error("Error in brcs");
      if_takeup[brc.idx_br_c] = true;
    }
  }
  for (int i = 0; i < num_row; ++i) {
    if (!if_takeup[i]) {
      cout << "i: " << i << endl;
      throw std::runtime_error("Error in check_if_cuts_are_legal");
    }
  }

  safe_solver(node->solver.updateModel())
  for (int i = 0; i < num_row; i++) {
    double cof, rhs;
    safe_solver(node->solver.getCoeff(i, 0, cof))
    safe_solver(node->solver.getRhs(i, 1, &rhs))
    if (cof != rhs) {
      cout << "i: " << i << endl;
      for (auto &rcc : node->rccs) {
        if (i == rcc.idx_rcc) {
          cout << "rcc: " << rcc.rhs << endl;
          break;
        }
      }
      for (auto &r1c : node->r1cs) {
        if (i == r1c.idx_r1c) {
          cout << "r1c: " << r1c.rhs << endl;
          break;
        }
      }
      for (auto &r2c : node->r1cs_multi) {
        if (i == r2c.idx_r1c) {
          cout << "r2c: " << r2c.rhs << endl;
          break;
        }
      }
      cout << "cof: " << cof << " rhs: " << rhs << endl;
      throw std::runtime_error("Error in not equal");
    }
  }
  /**
   * more test about the coefficients
   */

  cout << "check_if_cuts_are_legal more!" << endl;
  cout << "tmp not check rcc..." << endl;

}

void CVRP::updateUpperBoundEdgeSolution() {
  edge_sol.clear();
  edge_sol.reserve(2 * dim);
  for (auto &r : ip_opt_sol) {
    int past_node = 0;
    for (int i = 1; i < r.size(); ++i) {
      int current_node = r[i];
      auto pr =
          past_node < current_node ? make_pair(past_node, current_node) : make_pair(current_node, past_node);
      edge_sol.emplace(pr);
      past_node = current_node;
    }
  }
}

void CVRP::readSolutionFile(bool if_force) {
  string fileName = "sol_used/" + file_name + ".sol";
  ifstream fin(fileName);
  if (!fin && if_force) {
    throw std::runtime_error("Error in readSolutionFile");
  }
  if (!fin) return;
  ip_opt_sol.clear();
  double optval = 0;
  string line, item;
  while (std::getline(fin, line)) {
    if (line.find("<Optval") != std::string::npos) {
      stringstream ss(line);
      string optval_str;
      ss >> optval_str >> optval;
    } else if (line[0] == '0' && line[1] == '-') {
      stringstream ss(line);
      vector<int> tmp;
      tmp.reserve(dim);
      while (getline(ss, item, '-')) {
        tmp.emplace_back(std::stoi(item));
      }
      ip_opt_sol.emplace_back(std::move(tmp));
    }
  }
  fin.close();
  if (ip_opt_sol.empty() && if_force) {
    throw std::runtime_error("Error in No solution!");
  }
  cout << "optval: " << optval << " ub: " << ub << endl;
  if (abs(optval - ub) < TOLERANCE || optval < ub) {
    cout << "ub is updated to: " << optval << endl;
  } else {
    cout << "we abandoned the old ub (though better, we want the solution) and update to " << optval << endl;
  }
  ub = optval;
  updateUpperBoundEdgeSolution();
}

void CVRP::updateLowerBound(double val) {
  lb = val;
  lb_transformed = ceilTransformedNumberRelated(lb - TOLERANCE);
}