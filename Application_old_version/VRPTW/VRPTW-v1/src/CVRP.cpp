

#include "CVRP.hpp"
#include "templateFunctors.hpp"
#include "heavierHeuristic.hpp"
#include "lighterHeuristic.hpp"

using namespace std;
using namespace chrono;

CVRP::CVRP(const InstanceData &instanceData) {
  dim = instanceData.dim;
  real_dim = dim - 1;
  /**
   * todo: read the information and add stronger cut about vehicles used!
   */
  max_num_vehicle= real_dim;
  num_vehicle = max(instanceData.k, 1);
  cap = instanceData.cap;
  info_vertex = instanceData.info_vertex;
  file_name = instanceData.name;
  safe_Hyperparameter(checkMaximumNumberCustomers())
}

void CVRP::lateProcessing() {
  max_num_route4_mip = Config::max_num_route4_mip;
  rank1_mem_size_limit = int(dim * Config::MaxCutMemFactor);
  ratio_gap_improved_vs_time_lb= Config::RatioGapImprovedVSTimeIncreased_LB;
  ratio_gap_improved_vs_time_ub= Config::RatioGapImprovedVSTimeIncreased_UB;
  soft_time= Config::CutGenTimeThresholdInPricingInitial;
  gap_tolerance4_arc_elimination_n_enumeration = Config::InitGapTolerance4ArcEliminationNEnumeration;
  gap_improved_4_arc_elimination_n_enumeration= Config::GapImproved4ArcEliminationNEnumeration;
  num_buckets_per_vertex = Config::InitialNumBuckets;
  arc_elimination_time= Config::HardTimeThresholdInArcEliminationLastHalf;
  cost_mat4_vertex.resize(dim, vector<double>(dim, 0));
  ng_mem4_vertex.resize(dim, 0);
#ifdef EXTRA_ARC_MEMORY
  neighborhood_indicator_4_extra_arc_mem.resize(dim, 0);
#endif
  ng_mem4_vertex_sort_dist.resize(dim, vector<int>());
  rank1_sep_heur_mem4_vertex.resize(dim, 0);
  prior_mem_sets.resize(dim, Config::InitialNGSize);
  chg_cost_mat4_vertex.resize(dim, vector<double>(dim));
  round_up_tolerance = -1. / pow_self(10, transformed_number) + TOLERANCE;
  Config::MaxHeurSepMem4RowRank1=min(Config::MaxHeurSepMem4RowRank1, dim-1);
  cg_v_cut_map.resize(dim);
  cg_v_v_use_states.resize(dim, vector<vector<int>>(dim));
  lp_v_cut_map.resize(dim);
  lp_v_v_use_states.resize(dim, vector<vector<int>>(dim));

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
  vector<vector<int>> tex(dim, vector<int>(Config::InitialNGSize, 0));
  for (int i = 0; i < dim; ++i) {
    for (int j = 1; j < dim; ++j) {
      cost[j].first = j;
      cost[j].second = cost_mat4_vertex[i][j];
    }
	cost[0].first = 0;
    cost[0].second = INFINITY;
    std::stable_sort(cost.begin(), cost.end(),
                     [](const pair<int, double> &a, const pair<int, double> &b) {
                       return a.second < b.second;
                     });
	for(int k=0; k<dim; ++k){
	  ng_mem4_vertex_sort_dist[i].emplace_back(cost[k].first);
	}
    yzzLong &vst = ng_mem4_vertex[i];
    for (int k = 0; k < Config::InitialNGSize; ++k) {
      vst.set(cost[k].first);
    }
    yzzLong &vst2 = rank1_sep_heur_mem4_vertex[i];
    for (int k = 0; k < Config::MaxHeurSepMem4RowRank1; ++k) {
      vst2.set(cost[k].first);
    }
#ifdef EXTRA_ARC_MEMORY
	yzzLong &vst3 = neighborhood_indicator_4_extra_arc_mem[i];
	for (int k = 0; k < Config::neighborhood_indicator_size; ++k) {
	  vst3.set(cost[k].first);
	}
#endif
  }
  demand = new double[dim];
  for (int i = 0; i < dim; ++i) {
    demand[i] = info_vertex[i][3];
  }
  getLowerBoundofMinimumNumberCars();
  generateOptimalMultiplierForR1CMulti();
  setResourceInBucketGraph();

  auto tmp= int(double(resource.first_res) / num_buckets_per_vertex/pow(2, MAX_NUM_REGENERATE_BUCKET));
  step_size =res_int(tmp* pow(2, MAX_NUM_REGENERATE_BUCKET));
  num_buckets_per_vertex = (int) floor(resource.first_res / step_size) + 1;


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

#ifdef MASTER_VALVE_ML
  ml.calculate_prerequisites();
#endif
}

void CVRP::initializeLabels() {
  all_label[0].end_vertex = 0;
  all_label[0].p_label = nullptr;
  auto &cut_map= all_label[0].r1c.cut_map;
  fill(cut_map, cut_map + MAX_NUM_R1CS_IN_PRICING, 0);
  for (int i = 1; i < dim; ++i) {
    all_label[i].pi.set(i);
    all_label[i].cost = cost_mat4_vertex[0][i];
    all_label[i].end_vertex = i;
    all_label[i].p_label = all_label;
    increaseMainResourceConsumption({}, all_label[i].res, 0, i);
  }
#ifdef SYMMETRY_PROHIBIT
  int max_num = 2 * dim - 1;
  for (int i = dim; i < max_num; ++i) {
    int point = i - dim + 1;
    all_label[i].pi.set(point);
    all_label[i].cost = cost_mat4_vertex[0][point];
    all_label[i].end_vertex = point;
    all_label[i].p_label = all_label;
    decreaseMainResourceConsumption(resource, all_label[i].res, 0, point);
  }
#endif
}

void CVRP::buildModel() {
  auto *node = new BbNode(500, this);
  initializeBucketGraphForNode(node);

  int tmp_cnt = 0;

  node->cols.resize(dim);
  int colCnt=0;
  auto &col = node->cols[colCnt++];
  double res=0;
  for (int i = 1; i < dim; ++i){
	col.col_seq.emplace_back(i);
	col.main_res.emplace_back(0);
  }
  col.forward_concatenate_pos=(int)col.col_seq.size()-1;

  for (int i = 1; i < dim; ++i) {
	auto &col_i = node->cols[colCnt++];
	col_i.col_seq.emplace_back(i);
	col_i.main_res.emplace_back(resource_across_arcs_in_forward_sense[0][i].first_res);
	col_i.forward_concatenate_pos=0;
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
  old_ub= ub;
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

  vector<double> rhs(dim);
  vector<char> sense(dim);
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
                                             sense.data(),
                                             rhs.data(),
                                             nullptr))
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getNumRow(&num_row))
  safe_solver(node->solver.getNumCol(&num_col))
  safe_solver(node->solver.optimize())

  idx_node = 0;
  node->tree_level = 0;
  node->value = 0;
  node->index = idx_node;


  bbt.push(node);
  lb = 0;
  lb_transformed = 0;

  num_arti_vars = dim;
#ifndef SYMMETRY_PROHIBIT
  seq_size_arti_vars = 2 * real_dim;
#else
  seq_size_arti_vars = 4 * real_dim;
#endif
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

#ifndef FIND_INDICATOR

  if (if_open_heur) {
    if (runColumnGenerationType(node, 1) || runColumnGenerationType(node, 2))
      goto QUIT;
  }

#endif

  if (if_open_exact) {
#ifdef MIX_STYLE
	if(!node->r1cs.empty()) {
	  if (runColumnGenerationType(node, 4)) {
		goto QUIT;
	  }
	}else {
	  if (runColumnGenerationType(node, 3)) {
		goto QUIT;
	  }
	}
#else
	if (runColumnGenerationType(node, 3)) {
		goto QUIT;
	  }
#endif

	if(node->is_terminated) goto QUIT;
    if (num_col > LP_COL_FINAL_LIMIT) cleanIndexColForNode(node, false);
    if_can_arc_elimination_by_exact_cg = true;
  }

  if (if_record_sol && rollback != ROLL_BACK_LONG_TIME) {
    recordOptimalColumn(node);
  }

  QUIT:
  if (if_changed) {
    safe_solver(node->solver.setEnvMethod(env_method))
  }
  if (node->is_terminated) return;
#ifndef FIND_INDICATOR
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
#endif
}


void CVRP::priceLabeling(BbNode *node,const std::vector<double> &pi_vector) {
  for (int i = 1; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      chg_cost_mat4_vertex[i][j] = cost_mat4_vertex[i][j] - 0.5 * (pi_vector[i - 1] + pi_vector[j - 1]);
    }
  }
  for (int i = 1; i < dim; ++i) {
    chg_cost_mat4_vertex[0][i] = cost_mat4_vertex[0][i] - 0.5 * (pi_vector[i - 1] + pi_vector[real_dim]);
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
      rc = pi_vector[rcc.idx_rcc];
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
      rc = pi_vector[rcc.idx_rcc];
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
      chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] -= pi_vector[brc.idx_br_c];
      chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] -= pi_vector[brc.idx_br_c];
    }
  }
  getRank1DualsInCG(node, pi_vector);
}



void CVRP::assignMemory() {
  route_in_pricing_assign = MAX_ROUTE_PRICING;
  aver_route_length = int(ceil(dim / num_vehicle)) + 5;
  mem4_pricing = aver_route_length * (route_in_pricing_assign);
  col_pool4_pricing = new int[mem4_pricing];

#ifndef READ_ENUMERATION_TREES
  label_assign = LABEL_ASSIGN/2;
  reallocateLabel();
#endif

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
  delete[]col_pool4_pricing;
  delete[]copy_col_pool4_pricing;
  delete[]demand;
  delete[]label_int_space;

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

res_int roundAndConvertResLong(double value) {
  double rounded = std::round(value * RESOURCE_FACTOR) / RESOURCE_FACTOR;
  if (abs(rounded - value) > TOLERANCE) throw std::runtime_error("RESOURCE_FACTOR is too small");
  rounded *= RESOURCE_FACTOR;
  if (rounded > std::numeric_limits<res_int>::max()) {
	throw std::overflow_error("Value exceeds the range of res_int");
  }
  return static_cast<res_int>(rounded);
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

  pair<vector<SequenceInfo>, vector<double>> local_sol;
  pair<vector<SequenceInfo>, vector<double>> local_all_sol;

  double route_size = 0;
  int cnt=0;
  for (int i = 0; i < num_col; ++i) {
	if (X[i] > TOLERANCE) {
	  local_all_sol.second.emplace_back(X[i]);
	  local_all_sol.first.emplace_back(node->cols[i]);
	  route_size+=local_all_sol.first.back().col_seq.size();
	  ++cnt;
	  if(X[i]<1-TOLERANCE){
		local_sol.second.emplace_back(X[i]);
		local_sol.first.emplace_back(node->cols[i]);
	  }
	}
  }

  route_size/=	cnt;
  route_size+=3;// two zeros and round up
  aver_route_length=max(aver_route_length, route_size);

  if (node->all_lp_sol == local_all_sol) {
#if VERBOSE_MODE == 1
	cout << "solution remains the same" << endl;
#endif
    if (if_force_rewrite) {
      goto REWRITE;
    } else {
      return;
    }
  }

  node->only_frac_sol = local_sol;
  node->all_lp_sol = local_all_sol;

  REWRITE:
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      arc_graph[i][j] = 0;
      arc_graph_revised[i][j] = 0;
    }
  }

  auto opt_size= (int)local_all_sol.first.size();
  for(int i=0;i<opt_size;++i){
	auto &seq= local_all_sol.first[i].col_seq;
	double val= local_all_sol.second[i];
	int b4=0;
	for(auto j:seq){
	  arc_graph[b4][j]+=val;
	  b4=j;
	}
	arc_graph[b4][0]+=val;
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
		auto &seq= node->cols[i.first].col_seq;
		if (seq.size()==1) arc_graph_revised[0][seq[0]] -= i.second;
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
#if SOLVER_VRPTW==1
  START_OVER:
#endif
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getNumCol(&num_col))
  safe_solver(node->solver.getObjVal(&lp_val))
  if (num_col > LP_COL_FINAL_LIMIT && abs(prior_value - lp_val) > TOLERANCE)cleanIndexColForNode(node, false);
  safe_solver(node->solver.getObjVal(&lp_val))
  vector<double> X(num_col);
  safe_solver(node->solver.getX(0, num_col, X.data()))

  pi4_labeling.resize(num_row);
#ifdef FIND_INDICATOR

  if(indicator_state==0){
	safe_solver(node->solver.getDual(0, num_row, pi4_labeling.data()))
  }else if(indicator_state==1){
	copyModel2GetDual(node, true);
  }else if(indicator_state==2){
	copyModel2GetDual(node, false);
  }

#else

    changeModelForBetterDual(node);

#endif

  node->is_integer = true;
  for (int i = 0; i < num_col; ++i) {
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      node->is_integer = false;
      break;
    }
  }
  if (node->is_integer) {
    if (ceilTransformedNumberRelated(lp_val - TOLERANCE) + TOLERANCE < ub) {
	  updateIPOptSol(node, X);
#if SOLVER_VRPTW ==1
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

void CVRP::cleanIndexColForNode(BbNode *node, bool if_only_rcfixing) {
  if (!if_change_col) return;
  int beg = if_in_enu_state ? 1 : dim;
  vector<int> col_idx;
  if (if_only_rcfixing) {
	safe_solver(node->solver.updateModel())
	opt_gap = calculateOptimalGap(node);
	size_t numnzP;
	safe_solver(node->solver.XgetVars(&numnzP, nullptr, nullptr, nullptr, 0, num_col))
	vector<size_t> vbeg(num_col + 1);
	vector<int> vind(numnzP);
	vector<double> vval(numnzP);
	safe_solver(node->solver.XgetVars(&numnzP, vbeg.data(), vind.data(), vval.data(), 0, num_col))
	vector<Eigen::Triplet<double>> triplet(numnzP);
	vbeg[num_col] = numnzP;
	numnzP = 0;
	for (int i = 0; i < num_col; ++i) {
	  for (int j = vbeg[i]; j < vbeg[i + 1]; ++j) {
		triplet[numnzP++] = Eigen::Triplet<double>(vind[j], i, vval[j]);
	  }
	}
	sparseColMatrixXd A(num_row, num_col);
	A.setFromTriplets(triplet.begin(), triplet.end());
	Eigen::Map<Eigen::RowVectorXd> pi(optimal_dual_vector.data(), num_row);
	vector<double> cost(num_col);
	safe_solver(node->solver.getObj(0, num_col, cost.data()))
	Eigen::Map<Eigen::RowVectorXd> c(cost.data(), num_col);
	RowVectorXd rc= c- pi * A;
	for (int i = beg; i < num_col; ++i) if (rc[i] > opt_gap) col_idx.emplace_back(i);
  } else {
	safe_solver(node->solver.optimize())
	vector<double> rc(num_col);
	safe_solver(node->solver.getRC(0, num_col, rc.data()))
	vector<double> rc_copy(rc.data() + beg, rc.data() + num_col);
	int n = int((num_col - beg) * COL_KEEP_FRAC);
	nth_element(rc_copy.begin(), rc_copy.begin() + n, rc_copy.end());
	double threshold = max(rc_copy[n], TOLERANCE);
	for (int i = beg; i < num_col; ++i) if (rc[i] > threshold) col_idx.emplace_back(i);
  }
  rmLPCols(node, col_idx);
  safe_solver(node->solver.reoptimize())
}


void CVRP::getNewConstraintCoefficientByEdge(BbNode *node,
											 const std::pair<int, int> &edge,
											 std::vector<int> &ind,
											 std::vector<double> &val) {
  int size= node->edge_col_map[edge].size();
  ind.resize(size);
  val.resize(size);
  int cnt=0;
  for(auto &col: node->edge_col_map[edge]){
	ind[cnt]=col.first;
	val[cnt++]=col.second;
  }
}

void CVRP::getCoefficientRCC(BbNode *node, Rcc &rcc, std::vector<int> &ind,std::vector<double> &val) {
  auto & map= node->edge_col_map;
  vector<double> sup(num_col, 0);
  if (rcc.form_rcc) {//normal case
    auto &customer_info = rcc.info_rcc_customer;
    for (auto i = customer_info.begin(); i != customer_info.end(); ++i) {
      auto j = i;
      ++j;
      for (; j != customer_info.end(); ++j) {
		auto pr= *i < *j? make_pair(*i, *j): make_pair(*j, *i);
		for(auto &col: map[pr]){
		  sup[col.first]+=col.second;
		}
      }
    }
  } else {//another case
    auto &customer_info = rcc.info_rcc_customer;
    auto &outside_customer_info = rcc.info_rcc_outside_customer;
    int cnt = 0;
    for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
      auto j = i;
      ++j;
      for (; j != outside_customer_info.end(); ++j) {
		auto pr= *i < *j? make_pair(*i, *j): make_pair(*j, *i);
		for(auto &col: map[pr]){
		  sup[col.first]+=col.second;
		}
      }
    }

    for (auto customer_it : outside_customer_info) {
	  for(auto &col: map[{0, customer_it}]){
		sup[col.first]+=0.5* col.second;
	  }
    }
    for (auto customer_it : customer_info) {
	  for(auto &col: map[{0, customer_it}]){
		sup[col.first]-=0.5* col.second;
	  }
    }
  }

  ind.clear();
  val.clear();
  sup[0] = rcc.rhs;
  for (int i = 0; i < num_col; ++i) {
	if (sup[0]!=0){
	  ind.emplace_back(i);
	  val.emplace_back(sup[i]);
	}
  }
}

void CVRP::reformulateIntegerProgramSolution(BbNode *node) {
  cout << "reformulateIntegerProgramSolution" << endl;
#if  SOLVER_VRPTW==1
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

	vector<double> x(num_col);
	safe_solver(node->solver.getX(0, num_col, x.data()))
	for(int i=0;i<num_col;++i){
	  if(abs(x[i]-1)<TOLERANCE) {
		ip_opt_sol.emplace_back(node->cols[i].col_seq);
	  }
	}

	for(auto &pr: node->only_frac_sol.first){//should use only_frac_sol
	  int b4=0;
	  for(auto i: pr.col_seq){
		if(b4<i) arc_graph_revised[b4][i]=0;
		else arc_graph_revised[i][b4]=0;
		b4=i;
	  }
	  arc_graph_revised[0][b4]=0;
	}

    for (int j = 1; j < dim; ++j) {
      if (arc_graph_revised[0][j] == 0) continue;
      vector<int> tmp;
      tmp.reserve(dim);
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
	  tmp.erase(tmp.end()-1);//remove 0
      ip_opt_sol.emplace_back(std::move(tmp));
    }
  }
}

void CVRP::cleanAllPointers(BbNode *const node, int mode, bool if_clear_concatenate) {
  if (mode == 1) {
    for (int i = 0; i < dim; ++i) {
      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        label_array_in_forward_sense[i][b].clear();
        if_exist_extra_labels_in_forward_sense[i][b].second = 0;
      }
    }
  }
#ifdef SYMMETRY_PROHIBIT
  else if (mode == 2) {
    for (int i = 0; i < dim; ++i) {
      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        label_array_in_backward_sense[i][b].clear();
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
	goto QUIT;
  } else if_can_arc_elimination_by_exact_cg = false;
  if (if_stop_arc_elimination) {
    final_decision_4_arc_elimination = false;
    if_stop_arc_elimination = false;
#if VERBOSE_MODE==1
    cout << "arc elimination is banned!" << endl;
#endif
    goto QUIT;
  }
  if (final_decision_4_arc_elimination) goto QUIT;
  if (abs(old_ub - ub) > TOLERANCE) {
#if VERBOSE_MODE==1
    cout << "ub is updated from " << old_ub << " to " << ub << " and perform arc elimination!" << endl;
#endif
    old_ub = ub;
    final_decision_4_arc_elimination = true;
    goto QUIT;
  }
  now_gap = (ub - node->value) / ub;
  if (node->last_gap / now_gap > gap_improved_4_arc_elimination_n_enumeration
      || now_gap < gap_tolerance4_arc_elimination_n_enumeration) {
    node->last_gap = now_gap;
    final_decision_4_arc_elimination = true;
    goto QUIT;
  }
  QUIT:
  return;
}

void CVRP::determineIfEnumeration(BbNode *node) {
  double gap;
  if (final_decision_4_enumeration) {
#if VERBOSE_MODE==1
    cout << "enumeration is forced!" << endl;
#endif
    if (!if_arc_elimination_succeed) {
#if VERBOSE_MODE==1
	  cout<<"arc elimination failed! enumeration is skipped!"<<endl;
#endif
      final_decision_4_enumeration = false;
    }
    goto QUIT;
  }
  if (!if_arc_elimination_succeed && !if_arc_elimination_tried_but_failed) {
    goto QUIT;
  }

  gap = (ub - node->value) / ub;
  if (gap > Config::MaxGap2TryEnumeration || gap > last_enumeration_fail_gap) {
    goto QUIT;
  }
  final_decision_4_enumeration = true;
  if (if_arc_elimination_tried_but_failed) {
	final_decision_4_enumeration = false;
  }
  QUIT:
  return;
}

void CVRP::rollbackEasyWay(BbNode *const node, int old_num) {
  cout << BIG_PHASE_SEPARATION;
  if (if_in_enu_state)throw runtime_error("error: rollback is not allowed in enumeration state!");
  cout << "roll back to the previous cutting state!" << endl;

  ++ counter_try_hard_pricing;
  ratio_gap_improved_vs_time_lb/= Config::Ratio_Adjusted_GapImprovedVSTimeIncreased;
  ratio_gap_improved_vs_time_ub*= Config::Ratio_Adjusted_GapImprovedVSTimeIncreased;

  vector<int> delete_cuts(num_row - old_num);
  iota(delete_cuts.begin(), delete_cuts.end(), old_num);
  deleteNonactiveCuts(node, delete_cuts);

  for (auto &i : reset_cut_mem) {
	auto &cut= node->r1cs[get<0>(i)];
	cut.mem = get<1>(i);
	cut.arc_mem = get<2>(i);
  }

  reset_cut_mem.clear();
  getVCutMapLP(node);

  if(rollback_solver.model){
	node->solver.freeModel();
	node->solver.model=rollback_solver.model;
	rollback_solver.model=nullptr;
	safe_solver(node->solver.optimize())
	safe_solver(node->solver.getNumCol(&num_col))
	safe_solver(node->solver.getNumRow(&num_row))
	node->cols=rollback_cols;
	recordOptimalColumn(node);
  }else{
	cout<<"re-run cg to get the complete model!"<<endl;
	cout<<"also we decrease the soft time!"<<endl;
	soft_time/=Config::SoftTimeDecreaseFactorInPricing;
	safe_solver(node->solver.optimize())
	rollback = ROLL_BACK_INIT;
	force_not_rollback=true;
	solveLPInLabeling(node);
	force_not_rollback=false;
  }
  cout << BIG_PHASE_SEPARATION;
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
  getTopologicalOrder(node);
  node->canLeaveDepot_forward=0;
  for (auto j : node->all_forward_buckets[0][0].bucket_arcs)  node->canLeaveDepot_forward.set(j);
#ifdef SYMMETRY_PROHIBIT
  node->canLeaveDepot_backward=0;
  for (auto j : node->all_backward_buckets[0][0].bucket_arcs)  node->canLeaveDepot_backward.set(j);
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

int CVRP::addBranchConstraint(
							std::vector<int> &cind,
							std::vector<double> &cval,
                          char sense,
                          double rhs,
                          const char *constrname,
                          Solver &local_solver) {
  int error = local_solver.addConstraint(cind.size(), cind.data(), cval.data(), sense, rhs, constrname);
  if (error) return error;
  error += local_solver.updateModel();
  if (error) return error;
  int idxconstrs;
  error += local_solver.getNumRow(&idxconstrs);
  if (error) return error;
  --idxconstrs;
  int vind = 0;
  error += local_solver.XchangeCoeffs(1, &idxconstrs, &vind, &rhs);
  return error;
}

void CVRP::changeBranchConstraint(const std::vector<int> &vind,
								const std::vector<double> &vval,
                           char sense,
                           double rhs,
						   int row_idx,
                           Solver &local_solver) {
  vector<int> ind(num_col);
  iota(ind.begin(), ind.end(), 0);
  vector<double> val(num_col,0);
  for(int i=0;i<vind.size();++i){
	val[vind[i]]=vval[i];
  }
  vector<int> cind(num_col, row_idx);
  safe_solver(local_solver.changeCoeffs(num_col, cind.data(), ind.data(), val.data()))
  safe_solver(local_solver.setRhs(row_idx, 1, &sense, &rhs))
  int vind2 = 0;
  safe_solver(local_solver.XchangeCoeffs(1, &row_idx, &vind2, &rhs))
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


bool CVRP::increaseMainResourceConsumption(const ResTuple &nowResource,
										   ResTuple &newResource, int start, int end) {
  newResource= nowResource + resource_across_arcs_in_forward_sense[start][end];
  if(tellResTupleRelations<'p'>(newResource, ub4_vertex[end])) return false;
  if(newResource.first_res<lb4_vertex[end].first_res) newResource.first_res=lb4_vertex[end].first_res;
#ifdef USE_TWO_RESOURCE
  if(newResource.second_res<lb4_vertex[end].second_res) newResource.second_res=lb4_vertex[end].second_res;
#endif
  if(tellResTupleRelations<'c'>(newResource, resource_across_arcs_in_forward_sense[end][0])) return false;
  return true;
}

bool CVRP::decreaseMainResourceConsumption(const ResTuple &nowResource,
										   ResTuple &newResource, int start, int end) {
  newResource= nowResource - resource_across_arcs_in_backward_sense[start][end];
  if(tellResTupleRelations<'p'>(lb4_vertex[end], newResource)) return false;
  if(newResource.first_res>ub4_vertex[end].first_res) newResource.first_res=ub4_vertex[end].first_res;
#ifdef USE_TWO_RESOURCE
  if(newResource.second_res>ub4_vertex[end].second_res) newResource.second_res=ub4_vertex[end].second_res;
#endif
  /**
   * not type c
   */
  if(tellResTupleRelations<'p'>({}, newResource - resource_across_arcs_in_backward_sense[end][0])) return false;
  return true;
}

bool CVRP::runColumnGenerationType(BbNode *node, int mode) {
  if (node->is_terminated) return true;
  int status = 0;
  bool switch_is_on = true;
#if VERBOSE_MODE==1
  switch (mode) {
	case 1:cout << "LighterHeur phase begin...\n";
	  break;
	case 2:cout << "HeavierHeur phase begin...\n";
	  break;
	case 3:cout << "Exact phase begin...\n";
	  break;
	case 4:cout << "Mixed phase begin...\n";
	  break;
	default: throw runtime_error("Error: None of these modes are used!");
  }
#endif
  int iter = 0, ccnt = 0, old_ncol = num_col, tag;
  double b4_node_val = node->value;
  auto beg = chrono::high_resolution_clock::now();
  auto end = beg;

  time_pricing4_iter=0;
  time_resolve_lp4_iter=0;

#ifdef FIND_INDICATOR

  lp_time=0;
  it_=0;
  cg_time=0;

#endif

#ifdef CHANGE_DUAL
  if(mode==3 || mode==4) reviseSubPricingModel(node);
#endif

  while (true) {
    beg = chrono::high_resolution_clock::now();

    switch (mode) {
      case 1:ccnt = generateColumnsByLighterHeuristic<
#ifdef SYMMETRY_PROHIBIT
			false
#else
			true
#endif
		>(node);
		++iter;
        break;
      case 2:ccnt = generateColumnsByHeavierHeuristic<
#ifdef SYMMETRY_PROHIBIT
			false
#else
			true
#endif
		>(node);
		++iter;
        break;
      case 3:
        ccnt = generateColsByBidir<
#ifdef SYMMETRY_PROHIBIT
            false
#else
            true
#endif
        >(node);
		++iter;
        if (rollback == ROLL_BACK_TAIL_OFF)  soft_time_reached=true;
        break;
	  case 4:
		ccnt = runMixedExactColumnGeneration(node);
		++iter;
		if (rollback == ROLL_BACK_TAIL_OFF)  soft_time_reached=true;
		break;
    }

    if(mode!=4){
	  end = chrono::high_resolution_clock::now();
	  time_pricing4_iter += chrono::duration<double>(end - beg).count();
	}

#ifdef FIND_INDICATOR
	cg_time+=time_pricing4_iter;
	it_=iter;
#endif

	if ((iter % PRINT_LABELING_STEP_SIZE)==0) {
#if VERBOSE_MODE==2
	  goto OUTSIDE;
#endif
      REPORT:
      glo_end = chrono::high_resolution_clock::now();
      glo_eps = chrono::duration<double>(glo_end - glo_beg).count();
      printInfoLabeling(iter, num_col - old_ncol, num_col, num_row, time_resolve_lp4_iter,
                        time_pricing4_iter, glo_eps,
                        lp_val, lb, ub);
      if (!switch_is_on) {
        goto QUIT;
      }
	  time_pricing4_iter = 0;
	  time_resolve_lp4_iter = 0;
      old_ncol = num_col;
    }
	OUTSIDE:

    if (ccnt==0) {
	  break;
    }
    if (node->is_terminated) {
      break;
    }

    beg = chrono::high_resolution_clock::now();
    tag = optimizeLPForOneIteration(node, b4_node_val);
    b4_node_val = lp_val;

    if(mode!=4){
	  end = chrono::high_resolution_clock::now();
	  time_resolve_lp4_iter += chrono::duration<double>(end - beg).count();
	}
#ifdef FIND_INDICATOR
	lp_time+=time_resolve_lp4_iter;
#endif
    if (tag) {
      status = 1;//QUIT
      goto QUIT;
    }
  }
  if (switch_is_on && (iter % PRINT_LABELING_STEP_SIZE)) {
    switch_is_on = false;
    goto REPORT;
  }
  QUIT:
  #ifdef CHANGE_DUAL
  if(sub_pricing_solver.model) sub_pricing_solver.freeModel();
  #endif
  #ifdef FIND_INDICATOR
  if(indicator_state==0) {
	normal_eps.emplace_back(cg_time);
	normal_lp.emplace_back(lp_time);
	normal_it.emplace_back(it_);
  }
  else if(indicator_state==1) {
	good_eps.emplace_back(cg_time);
	good_lp.emplace_back(lp_time);
	good_it.emplace_back(it_);
  }
  else if(indicator_state==2) {
	bad_eps.emplace_back(cg_time);
	bad_lp.emplace_back(lp_time);
	bad_it.emplace_back(it_);
  }
  #endif
  return status;
}

void CVRP::generateRCMatrixInPricing(BbNode *node,
									 sparseColMatrixXd &mat,
									 RowVectorXd &cost) {
  auto index = useful_label_vec.size();
  cost.resize((int)index);
  aux_label_vec = new AUX_LABEL[(int)index];
  for(int i=0;i<index;++i) {
	aux_label_vec[i].sparse_lp_states.resize(lp_r1c_denominator.size());
	aux_label_vec[i].states.resize(lp_r1c_denominator.size(),0);
  }
  index = 0;

  vector<Eigen::Triplet<double> > triplets;
  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> edge_map;

  for (auto &ki : useful_label_vec) {
	ki->aux_label = aux_label_vec+index;
	ki->aux_label->index4_rc_matrix = (int)index;
	generateRCMatrix4Label(triplets, ki, cost, edge_map);
	++index;
  }

  vector<vector<double>> rcc_mat(node->rccs.size(), vector<double>(index, 0));

  int cnt=-1;
  for (auto &rcc : node->rccs) {
	++cnt;
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
		  for (auto &col : edge_map[pr]) ++rcc_mat[cnt][col];
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
		  for (auto &it_map : edge_map[pr]) ++rcc_mat[cnt][it_map];
		}
	  }
	  for (int aj : outside_customer_info) {
		for (auto &it_map : edge_map[make_pair(0, aj)]) rcc_mat[cnt][it_map]+=0.5;
	  }
	  for (int aj : customer_info) {
		for (auto &it_map : edge_map[make_pair(0, aj)])  rcc_mat[cnt][it_map]-=0.5;
	  }
	}
  }

  vector<vector<double>> brc_mat(node->brcs.size(), vector<double>(index, 0));

  cnt=-1;
  for (auto &br : node->brcs) {
	++cnt;
	int ai = br.edge.first;
	int aj = br.edge.second;
	auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
	for (auto &it_map : edge_map[pr]) ++brc_mat[cnt][it_map];
  }

  auto n = triplets.size();

  triplets.resize(n + brc_mat.size() *index + index + rcc_mat.size()*index);

  for(int i=0;i<rcc_mat.size();++i){
	int idx=node->rccs[i].idx_rcc;
	for(int j=0;j<rcc_mat[i].size();++j) {
	  if(rcc_mat[i][j]!=0) triplets[n++] = {idx, j, rcc_mat[i][j]};
	}
  }

  for(int i=0;i<brc_mat.size();++i){
	int idx=node->brcs[i].idx_br_c;
	for(int j=0;j<brc_mat[i].size();++j) {
	  if (brc_mat[i][j] != 0) triplets[n++] = {idx, j, brc_mat[i][j]};
	}
  }

  for (int i = 0; i < index; ++i) {
	triplets[n++] = {real_dim, i, 0.5};//note here is 0.5
  }

  triplets.resize(n);

  mat.resize(num_row, (int)index);
  mat.setFromTriplets(triplets.begin(), triplets.end());
}


int CVRP::runMixedExactColumnGeneration(BbNode *node) {
  time_point<high_resolution_clock > beg4lp, beg4cg, end4lp, end4cg;
  beg4cg=high_resolution_clock::now();
  int ccnt = generateColsByBidir<
#ifdef SYMMETRY_PROHIBIT
	  false
#else
	  true
#endif
  >(node);
  end4cg=high_resolution_clock::now();
  time_pricing4_iter+=chrono::duration<double>(end4cg - beg4cg).count();
  if(ccnt==0) return 0;

  takeOutUsefulLabels<true,
#ifdef SYMMETRY_PROHIBIT
	  false
#else
	  true
#endif
  >();

  sparseColMatrixXd rc_mat;
  RowVectorXd cost;
  generateRCMatrixInPricing(node, rc_mat, cost);

  int heur_cnt=1;
  double b4_node_val = node->value;
  while(heur_cnt){
	beg4lp=high_resolution_clock::now();
	bool tag = optimizeLPForOneIteration(node, b4_node_val);
	end4lp=high_resolution_clock::now();
	time_resolve_lp4_iter+=chrono::duration<double>(end4lp - beg4lp).count();
	if(tag) break;
	b4_node_val = lp_val;
	beg4cg=high_resolution_clock::now();
	heur_cnt=retrieveLabeling<
#ifdef SYMMETRY_PROHIBIT
		false
#else
		true
#endif
	>(node, rc_mat, cost);
	end4cg=high_resolution_clock::now();
	time_pricing4_iter+=chrono::duration<double>(end4cg - beg4cg).count();
   cout<<"time for heur= "<<chrono::duration<double>(end4cg - beg4cg).count()<<endl;
	ccnt+=heur_cnt;
   cout<< " cnt= "<<heur_cnt<< " lp_val= " <<lp_val<<endl;
  }
  delete []aux_label_vec;
  aux_label_vec= nullptr;
  return ccnt;
}

double CVRP::transformCost(double x) {
  return std::floor(x + 0.5);
}

void CVRP::setTailOffStandardAndRollBackStandard() const {
  for (int b = 0; b < num_buckets_per_vertex; ++b) {
    int min_num_labels = MAX_INT;
    for (int i = 0; i < dim; ++i) {
      if (label_array_in_forward_sense[i][b].size() && min_num_labels > label_array_in_forward_sense[i][b].size()) {
        min_num_labels = label_array_in_forward_sense[i][b].size();
      }
    }
    if (min_num_labels == MAX_INT) {
      min_num_labels = 1;
    }
    for (int i = 0; i < dim; ++i) {
      int hard_max = max((int)label_array_in_forward_sense[i][b].size(), min_num_labels) * FACTOR_NUM_LABEL;
      hard_max = int(pow_self(2.0, int(ceil(log(hard_max) / log(2)) + TOLERANCE)));
      if_exist_extra_labels_in_forward_sense[i][b].first.resize(hard_max);
    }
#ifdef SYMMETRY_PROHIBIT
    min_num_labels = MAX_INT;
    for (int i = 0; i < dim; ++i) {
      if (label_array_in_backward_sense[i][b].size() && min_num_labels > label_array_in_backward_sense[i][b].size()) {
        min_num_labels = label_array_in_backward_sense[i][b].size();
      }
    }
    if (min_num_labels == MAX_INT) {
      min_num_labels = 1;
    }
    for (int i = 0; i < dim; ++i) {
      int hard_max = max((int)label_array_in_backward_sense[i][b].size(), min_num_labels) * FACTOR_NUM_LABEL;
      hard_max = (int) (pow_self(2.0, (int) ceil(log(hard_max) / log(2))) + TOLERANCE);
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


void CVRP::regenerateGraphBucket(BbNode *node) {
  if (step_size / 2 < 1) return;
  if(step_size%2) {
	cout<<"regenerateGraphBucket is banned since step_size is odd!"<<endl;
	return;
  }
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

  auto new_IfExistExtraLabelsInForwardSense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	new_IfExistExtraLabelsInForwardSense[i] = new VecLabel[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  new_IfExistExtraLabelsInForwardSense[i][j].first.resize(
		  if_exist_extra_labels_in_forward_sense[i][j / 2].first.size() / 2);
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

  auto new_IfExistExtraLabelsInBackwardSense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	new_IfExistExtraLabelsInBackwardSense[i] = new VecLabel[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  new_IfExistExtraLabelsInBackwardSense[i][j].first.resize(
		  if_exist_extra_labels_in_backward_sense[i][j / 2].first.size() / 2);
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
	label_array_in_forward_sense[i] = new ListLabel [num_buckets_per_vertex];
  }

  for (int i = 0; i < dim; ++i) {
	delete[]if_exist_extra_labels_in_forward_sense[i];
  }
  delete[]if_exist_extra_labels_in_forward_sense;
  if_exist_extra_labels_in_forward_sense = new_IfExistExtraLabelsInForwardSense;

#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < dim; ++i) {
	delete[]label_array_in_backward_sense[i];
	label_array_in_backward_sense[i] = new ListLabel [num_buckets_per_vertex];
  }
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

  getTopologicalOrder(node);
  node->canLeaveDepot_forward=0;
  for (auto j : node->all_forward_buckets[0][0].bucket_arcs)  node->canLeaveDepot_forward.set(j);
#ifdef SYMMETRY_PROHIBIT
  node->canLeaveDepot_backward=0;
  for (auto j : node->all_backward_buckets[0][0].bucket_arcs)  node->canLeaveDepot_backward.set(j);
#endif
#if VERBOSE_MODE==1
	cout << "new generated bucket graph: num_buckets_per_vertex= " << num_buckets_per_vertex << " step_size= "
		 << step_size
		 << endl;
	cout << "we cannot use arc elimination next!" << endl;
#endif
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
  double std;
  while (true) {
    ++round;
    cout << "DSSR: " << round << endl;
    findNGMemorySets(node, if_empty);
	if(if_empty) break;
    auto beg = chrono::high_resolution_clock::now();
    force_not_rollback = true;
    solveLPInLabeling(node);
    force_not_rollback = false;
    auto end = chrono::high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    if (eps > cost_time) cost_time = eps;
    if (last_max_time_labeling > last && node->value + TOLERANCE < old_val) {
#if VERBOSE_MODE==1
      cout << "DSSR: Reach the hard limit! We use ng rollback!" << endl;
#endif
      ng_mem4_vertex = old_mem;
      deleteColumnByNGMemory(node, dim);
      force_not_rollback = true;
      solveLPInLabeling(node);
      force_not_rollback = false;
      return;
    }
    old_mem = ng_mem4_vertex;
    if (node->value + TOLERANCE > old_val) break;
  }
  if(if_empty) goto QUIT;
  round = 0;
  old_val = node->value;
  std = TOLERANCE;
  while (true) {
    ++round;
    cout << "NGAugmentation: " << round << endl;
    findNGMemorySets(node, if_empty);
    if (if_empty) break;
    auto beg = chrono::high_resolution_clock::now();
    solveLPInLabeling(node);
    auto end = chrono::high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    if (eps > Config::NGAugTimeHardThresholdFactor * cost_time) {
      rollback = ROLL_BACK_LONG_TIME;
    } else if (eps > Config::NGAugTimeSoftThresholdFactor * cost_time) {
      rollback = ROLL_BACK_TAIL_OFF;
    }
#if VERBOSE_MODE==1
    cout << "eps= " << eps << endl;
    cout << "cost_time= " << cost_time << endl;
#endif
    if (rollback == ROLL_BACK_TAIL_OFF) break;
    else if (rollback == ROLL_BACK_LONG_TIME) {
      cout << "NGAugmentation: Reach the hard limit! We use ng rollback!" << endl;
      ng_mem4_vertex = old_mem;
      deleteColumnByNGMemory(node, dim);
      force_not_rollback = true;
      solveLPInLabeling(node);
      force_not_rollback = false;
      return;
    }
    old_mem = ng_mem4_vertex;
    double tmp = abs(node->value - old_val) / node->value;
    if (tmp > std) std = tmp;
    else if (tmp < std * Config::NGAugTailOff) break;
    else old_val = node->value;
  }
  QUIT:
  return;
}

void CVRP::findNGMemorySets(BbNode *node, bool &if_empty) {
  if (node->index) throw runtime_error("findNGMemorySets called when node->index!=0");
  unordered_map<int, vector<pair<int, yzzLong>>> map_size_cycle;
  vector<size_t> tmp(dim);
  vector<pair<int, double>> opt_col_idx(node->only_frac_sol.first.size());
  for (int i = 0; i < opt_col_idx.size(); ++i) {
	opt_col_idx[i] = make_pair(i, node->only_frac_sol.second[i]);
  }
  sort(opt_col_idx.begin(), opt_col_idx.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
    return a.second > b.second;
  });

  int re_cols = 0;
  for (auto &pr : opt_col_idx) {
	fill(tmp.begin(), tmp.end(), -1);
    bool if_re_col = false;
	auto &route= node->only_frac_sol.first[pr.first].col_seq;
	for(int i=0;i<route.size();++i){
	  int curr_node=route[i];
	  if (tmp[curr_node]!=-1) {
		if_re_col = true;
		auto length = int(i - tmp[curr_node] - 1);
		map_size_cycle[length].emplace_back(curr_node, 0);
		auto &mem = map_size_cycle[length].back().second;
		for (auto k = tmp[curr_node] + 1; k < i; ++k) {
		  mem.set(route[k]);
		}
	  }
	  tmp[curr_node] = i;//not in else!
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
#if VERBOSE_MODE==1
    cout << "no small cycles are found!" << endl;
#endif
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
  deleteColumnByNGMemory(node, dim);
}

void CVRP::deleteColumnByNGMemory(BbNode *node, int start, bool if_full_mem) {
  yzzLong local_pi;
  vector<int> col_idx;

  for(int i=start;i<num_col;++i){
	local_pi = 0;
	for(auto current_node: node->cols[i].col_seq){
	  if (local_pi[current_node]) {
		col_idx.emplace_back(i);
		break;
	  }
	  if(!if_full_mem)local_pi = local_pi & ng_mem4_vertex[current_node];
	  local_pi.set(current_node);
	}
  }
  rmLPCols(node, col_idx);
  safe_solver(node->solver.reoptimize())
  safe_solver(node->solver.getObjVal(&lp_val))
#if VERBOSE_MODE==1
  if(if_full_mem){
	cout << "after clean non-ele routes lpval= " << lp_val << endl;
  }else{
	cout<<"after clean ng routes lpval= "<<lp_val<<endl;
  }
#endif
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
		int current_node = stoi(item);
		if(!current_node) continue;
        tmp.emplace_back(current_node);
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
}

void CVRP::updateLowerBound(double val) {
  double tmp_lb;
  if(!bbt.empty()){
	tmp_lb= bbt.top()->value;
  }else{
	if(!sub_bbt.empty()) {
	  tmp_lb = sub_bbt.top()->value;
	}else{
	  tmp_lb=val;
	}
  }
  if(lb!=tmp_lb) {
	lb = tmp_lb;
	lb_transformed = ceilTransformedNumberRelated(lb - TOLERANCE);
  }
}

void CVRP::rmLPCols(BbNode*node,const std::vector<int> &col_idx){
  if(col_idx.empty()) return;
  vector<int> sort_col_idx=col_idx;
  if(!is_sorted(col_idx.begin(),col_idx.end())){
	sort(sort_col_idx.begin(),sort_col_idx.end());
  }

  int delta = 0;
  auto stop_sign = sort_col_idx.end() - 1;
  for (auto i = sort_col_idx.begin(); i < stop_sign; ++i) {
	++delta;
	for (int j = *i + 1; j < *(i + 1); ++j) node->cols[j - delta] = node->cols[j];
  }
  ++delta;
  for (int j = *stop_sign + 1; j < num_col; ++j) node->cols[j - delta] = node->cols[j];

  safe_solver(node->solver.delVars(sort_col_idx.size(), sort_col_idx.data()))
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getNumCol(&num_col))

  node->cols.resize(num_col);

#ifdef CHANGE_DUAL
  if (sub_pricing_solver.model) {
	cout << "revise sub_pricing_solver as well!" << endl;
	safe_solver(sub_pricing_solver.delVars(sort_col_idx.size(), sort_col_idx.data()))
	safe_solver(sub_pricing_solver.updateModel())
	int loacl_col;
	safe_solver(sub_pricing_solver.getNumCol(&loacl_col))
	if (loacl_col != num_col) {
	  cout << "loacl_col= " << loacl_col << " num_col= " << num_col << endl;
	  throw runtime_error("Error in sub_pricing_solver");
	}
  }
#endif
}

void CVRP::updateEdgeColMap(BbNode *node){
  /**
  * revise the edge_col_map
  */
  std::unordered_map<std::pair<int, int>, std::unordered_map<int, int>, PairHasher>
	  sup_map{};
  for(int i=0;i<num_col;++i) {
	int b4=0;
	for(auto current_node: node->cols[i].col_seq){
	  auto pair= b4<current_node? make_pair(b4,current_node):make_pair(current_node,b4);
	  ++sup_map[pair][i];
	  b4=current_node;
	}
	++sup_map[{0, b4}][i];
  }
  auto &map=node->edge_col_map;
  map.clear();
  for(auto &pr:sup_map){
	map[pr.first] = vector<pair<int, int>>(pr.second.begin(), pr.second.end());
  }
}


void CVRP::updateIPOptSol(BbNode *node, const std::vector<double> &X) {
  ip_opt_sol.clear();
  for(int i=0;i<num_col;++i){
	if (X[i]>TOLERANCE){
	  ip_opt_sol.emplace_back(node->cols[i].col_seq);
	}
  }
}