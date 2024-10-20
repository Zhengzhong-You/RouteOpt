#include "cvrp.hpp"
#include "template_functors.hpp"
#include "heavier_heuristic.hpp"
#include "lighter_heuristic.hpp"
#ifdef HEURISTIC
#include "heuristic.hpp"
#endif

#include "branching.hpp"
#if SOLUTION_TYPE == 1
#include "best_bound_first_branching.hpp"
#elif SOLUTION_TYPE == 2
#include "depth_first_branching.hpp"
#endif

#ifdef READ_ENUMERATION_TREES
#include "read_enumeration_tree.hpp"
#endif

#ifdef FASTER_DUAL
#include "faster_dual.hpp"
#endif

#ifdef DUAL_SMOOTHING
#include "dual_smoothing.hpp"
#endif

#ifdef COMBINE_FASTER_SMOOTHING
#include "combine_faster_smoothing.hpp"
#endif

#ifdef MASTER_VALVE_ML
#include "machine_learning.hpp"
#endif


#include "robust_control.hpp"

using namespace std;
using namespace chrono;

CVRP::CVRP(const InstanceData &instanceData) {
    dim = instanceData.dim;
    real_dim = dim - 1;
    /**
     * todo: read the information and add stronger cut about vehicles used!
     */
    max_num_vehicle = real_dim;
    num_vehicle = max(instanceData.k, 2);
    cap = instanceData.cap;
    info_vertex = instanceData.info_vertex;
    file_name = instanceData.name;
    safe_Hyperparameter(checkMaximumNumberCustomers())
    read_dual_sol_call(ReadDualSol::readDualSol(file_name))
}

void CVRP::initialProcessing() {
    MaxNumRoute4Mip = Config::MaxNumRoute4Mip;
    gap_tolerance4_arc_elimination_n_enumeration = Config::InitGapTolerance4ArcEliminationNEnumeration;
    gap_improved_4_arc_elimination_n_enumeration = Config::GapImproved4ArcEliminationNEnumeration;
    num_buckets_per_vertex = Config::InitialNumBuckets;
    arc_elimination_time = Config::HardTimeThresholdInArcEliminationLastHalf;
    cost_mat4_vertex.resize(dim, vector<double>(dim, 0));
    ng_mem4_vertex.resize(dim, 0);
    ng_mem4_vertex_sort_dist.resize(dim, vector<int>());
    prior_mem_sets.resize(dim, Config::InitialNGSize);
    chg_cost_mat4_vertex.resize(dim, vector<double>(dim));
    round_up_tolerance = -1. / pow_self(10, transformed_number) + TOLERANCE;
    Config::InitialNGSize = min(Config::InitialNGSize, dim - 1);
    cg_v_cut_map.resize(dim);
    cg_v_v_use_states.resize(dim, vector<vector<int> >(dim));
    lp_v_cut_map.resize(dim);
    lp_v_v_use_states.resize(dim, vector<vector<int> >(dim));

    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            cost_mat4_vertex[i][j] = transformCost(
                sqrt(float((info_vertex[i][1] - info_vertex[j][1]) * (info_vertex[i][1] - info_vertex[j][1]) +
                           (info_vertex[i][2] - info_vertex[j][2]) * (info_vertex[i][2] - info_vertex[j][2]))));
        }
    }
    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            cost_mat4_vertex[j][i] = cost_mat4_vertex[i][j];
        }
    }

    vector<pair<int, double> > cost(dim);
    vector<vector<int> > tex(dim, vector<int>(Config::InitialNGSize, 0));
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
        for (int k = 0; k < dim; ++k) {
            ng_mem4_vertex_sort_dist[i].emplace_back(cost[k].first);
        }
        yzzLong &vst = ng_mem4_vertex[i];
        for (int k = 0; k < Config::InitialNGSize; ++k) {
            vst.set(cost[k].first);
        }
    }
    demand = new double[dim];
    for (int i = 0; i < dim; ++i) {
        demand[i] = info_vertex[i][3];
    }
    getLowerBoundofMinimumNumberCars();
    Rank1CutsSeparator::setInitialInfo(Config::MaxRowRank1, Config::MaxNumR1C3PerRound, Config::MaxNumR1CPerRound, dim,
                                       solver, cost_mat4_vertex);
    setResourceInBucketGraph();

    auto tmp = int(double(resource.first_res) / num_buckets_per_vertex / pow(2, MAX_NUM_REGENERATE_BUCKET));
    step_size = res_int(max(tmp * pow(2, MAX_NUM_REGENERATE_BUCKET), 1.)); //at least be 1!
    num_buckets_per_vertex = (int) floor(resource.first_res / step_size) + 1;

    assignMemory();

    initializeBucketGraph();
    initializeLabels();

    int candi_size = dim * dim / 2;
    BaseBranching::branching_history.exact_improvement_up.reserve(candi_size);
    BaseBranching::branching_history.exact_improvement_down.reserve(candi_size);
    BaseBranching::branching_history.heuristic_improvement_up.reserve(candi_size);
    BaseBranching::branching_history.heuristic_improvement_down.reserve(candi_size);
    BaseBranching::branching_history.lp_testing_improvement_up.reserve(candi_size);
    BaseBranching::branching_history.lp_testing_improvement_down.reserve(candi_size);
}

void CVRP::initializeLabels() {
    all_label[0].end_vertex = 0;
    all_label[0].p_label = nullptr;
    auto &cut_map = all_label[0].r1c.cut_map;
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

void CVRP::setSolverEnv() {
    safe_solver(solver.loadEnv(nullptr))
    safe_solver(solver.setEnvThreads(NUM_THREADS_LP, true))
    safe_solver(solver.setEnvOutputFlag(0, true))
    safe_solver(solver.setEnvInfUnbdInfo(1, true))
    safe_solver(solver.setEnvMIPGap(MIP_GAP_TOLERANCE, true))
#if (SOLVER_TYPE == 1) && defined(VRPSOLVER_CPLEX_PARAMETER)
  {
	int status;
	auto &env = solver.env;
	status = CPXsetdblparam(env, CPXPARAM_Simplex_Tolerances_Optimality, 9.9999999999999995e-08);
	status = CPXsetdblparam(env, CPXPARAM_Simplex_Tolerances_Feasibility, 9.9999999999999995e-08);
	status = CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, 1.0000000000000001e-09);
	int int_value;
	double dbl_value;

	CPXgetdblparam(env, CPXPARAM_Simplex_Tolerances_Optimality, &dbl_value);
	printf("CPXPARAM_Simplex_Tolerances_Optimality: %.17g\n", dbl_value);

	CPXgetdblparam(env, CPXPARAM_Simplex_Tolerances_Feasibility, &dbl_value);
	printf("CPXPARAM_Simplex_Tolerances_Feasibility: %.17g\n", dbl_value);

	CPXgetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, &dbl_value);
	printf("CPXPARAM_MIP_Tolerances_MIPGap: %.17g\n", dbl_value);
  }
#endif
}

#ifdef FIND_ALL_SOLUTIONS
void addSol2Node(CVRP *cvrp, BbNode *node) {
  auto &ipopt = cvrp->ip_opt_sol;
  if (ipopt.empty()) return;
  auto &new_col = cvrp->new_cols;
  new_col.clear();
  for (auto &route : ipopt) {
	if (route.empty()) continue;
	new_col.emplace_back();
	auto &col = new_col.back();
	col.col_seq = route;
	col.main_res.emplace_back(0);
	col.forward_concatenate_pos = (int)col.col_seq.size() - 1;
  }
  int ccnt = (int)ipopt.size();
  cvrp->addColumns(node, ccnt, false);
  safe_solver(node->getSolver().reoptimize(SOLVER_PRIMAL_SIMPLEX))
  double lp_val;
  safe_solver(node->getSolver().getObjVal(&lp_val))
  cout << "lp_val= " << lp_val << endl;
  node->getCurrentNodeVal() = lp_val;
}
#endif

void CVRP::buildModel() {
    auto *node = new BbNode(500, this);
    initializeBucketGraphForNode(node);

    vector<int> solver_beg, solver_ind;
    vector<double> solver_val, solver_obj;

    solver_beg.emplace_back(0);
    node->getCols().emplace_back();
    auto &col = node->getCols()[0];
    double res = 0;
    for (int i = 1; i < dim; ++i) {
        col.col_seq.emplace_back(i);
        solver_ind.emplace_back(i - 1);
        solver_val.emplace_back(1);
    }
    solver_ind.emplace_back(real_dim);
    solver_val.emplace_back(num_vehicle);
    solver_obj.emplace_back(0);
    col.main_res.emplace_back(0);
    col.forward_concatenate_pos = (int) col.col_seq.size() - 1;

    vector<int> feasible_set(dim, 0);
    int std = 1;
    for (auto &i: node->all_forward_buckets[0][0].bucket_arcs) ++feasible_set[i];
    symmetry_prohibit_call(for (auto
            &i : node->all_backward_buckets[0][0].bucket_arcs)
        ++feasible_set[i];
        std = 2;)

    double obj_sum = 0;
    for (int i = 1; i < dim; ++i) {
        if (feasible_set[i] != std) {
            cout << i << " cannot be visited alone! model infeasible!" << endl;
            exit(0);
            continue;
        }
        node->getCols().emplace_back();
        auto &col_i = node->getCols().back();
        col_i.col_seq.emplace_back(i);
        col_i.main_res.emplace_back(resource_across_arcs_in_forward_sense[0][i].first_res);
        col_i.forward_concatenate_pos = 0;
        solver_beg.emplace_back((int) solver_ind.size());
        solver_ind.emplace_back(i - 1);
        solver_val.emplace_back(1);
        solver_ind.emplace_back(real_dim);
        solver_val.emplace_back(1);
        // remove_vehicle_constraint_call()
        solver_obj.emplace_back(2 * cost_mat4_vertex[0][i]);
        obj_sum += solver_obj.back();
    }

    solver_beg.emplace_back((int) solver_ind.size());
    BaseBranching::ub = Config::ub < obj_sum ? Config::ub : obj_sum;
    old_ub = BaseBranching::ub;
    cout << "ub= " << BaseBranching::ub << endl;

    double increase_UB = int(BaseBranching::ub * 10);
    cout << "now the increase_UB is " << increase_UB << endl;
    solver_obj[0] = min(increase_UB, obj_sum);

    setSolverEnv();
    node->getSolver().getEnv(&solver);
    rollback_solver.getEnv(&solver);
    const char *model_name = "CVRP.lp";

    cout << SMALL_PHASE_SEPARATION;
    cout << "<Instance  " << file_name << "  Capacity  " << cap << ">" << endl;

    vector<double> rhs(dim, 1);
    vector<char> sense(dim, SOLVER_EQUAL);

    rhs[real_dim] = num_vehicle;
    sense[real_dim] = SOLVER_GREATER_EQUAL;

    safe_solver(node->getSolver().newModel(model_name, 0, nullptr, nullptr, nullptr, nullptr, nullptr))
    safe_solver(node->getSolver().addConstraints(
            dim,
            0,
            nullptr,
            nullptr,
            nullptr,
            sense.data(),
            rhs.data(),
            nullptr)
    )
    safe_solver(node->getSolver().addVars((int)solver_beg.size() - 1,
        (int)solver_ind.size(),
        solver_beg.data(),
        solver_ind.data(),
        solver_val.data(),
        solver_obj.data(),
        nullptr,
        nullptr,
        nullptr,
        nullptr))
    safe_solver(node->getSolver().updateModel())
    safe_solver(node->getSolver().getNumRow(&num_row))
    safe_solver(node->getSolver().getNumCol(&num_col))
    safe_solver(node->getSolver().optimize())

#ifdef FIND_ALL_SOLUTIONS
  addSol2Node(this, node);
#endif

    idx_node = 0;
    node->getTreeLevel() = 0;
    node->getCurrentNodeVal() = 0;
    node->index = idx_node;

#if SOLUTION_TYPE == 1
    BestBoundFirstBranching::bbt.push(node);
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::addNodeIn(DepthFirstBranching::bbt, node);
#endif
    BaseBranching::lb = 0;
    BaseBranching::lb_transformed = 0;

    heuristic_call(Heuristic::heuristic_col_info = node->getCols();)
}

void CVRP::solveLPInLabeling(BbNode *node, bool if_open_heur, bool if_open_exact, bool if_record_sol) {
    dual_smoothing_call(DualSmoothing::buildLPMatrix())

    combine_faster_call(DualSmoothing::buildLPMatrix())

    if_roll_back = false;

    if (!node->index) {
        ratio_dominance_checks_non_dominant = {};
    }

    optimizeLPForOneIteration(node, node->getCurrentNodeVal() barrier_call(, SOLVER_BARRIER));

    if (if_open_heur) {
        if (runColumnGenerationType(node, 1) || runColumnGenerationType(node, 2))
            goto QUIT;
    }

    if (if_open_exact) {
        if_exact_cg_finished = true;
        if (
            false &&
            !node->r1cs.empty()) {
            if (runColumnGenerationType(node, 4)) {
                goto QUIT;
            }
        } else {
            if (runColumnGenerationType(node, 3)) {
                goto QUIT;
            }
        }
        if (num_col > LP_COL_FINAL_LIMIT) cleanIndexColForNode(node, false);
        if_can_arc_elimination_by_exact_cg = true;
    }

    if (if_record_sol && if_exact_cg_finished) {
        recordOptimalColumn(node);
        max_labeling_time_last_round_for_node_so_far = max_labeling_time_for_node_so_far;
    }

QUIT:dual_smoothing_call(DualSmoothing::freeLPMatrix())
    combine_faster_call(DualSmoothing::freeLPMatrix())
    if (node->getIfTerminated()) return;
    if (if_exact_cg_finished && !node->index && if_open_exact) {
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

void CVRP::priceLabeling(BbNode *node, const std::vector<double> &pi_vector) {
    for (int i = 1; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            chg_cost_mat4_vertex[i][j] = cost_mat4_vertex[i][j] - 0.5 * (pi_vector[i - 1] + pi_vector[j - 1]);
        }
    }
    for (int i = 1; i < dim; ++i) {
        chg_cost_mat4_vertex[0][i] =
                cost_mat4_vertex[0][i] - 0.5 * (pi_vector[i - 1] + pi_vector[real_dim]);
    }
    for (int i = 1; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            chg_cost_mat4_vertex[j][i] = chg_cost_mat4_vertex[i][j];
        }
    }
    for (int i = 1; i < dim; ++i) {
        chg_cost_mat4_vertex[i][0] = chg_cost_mat4_vertex[0][i];
    }

    priceRCC(node, pi_vector);
    priceBRC(node, pi_vector);
}

void CVRP::priceRCC(BbNode *node, const std::vector<double> &pi_vector) {
    double rc;
    for (auto &rcc: node->rccs) {
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
            for (auto it: outside_customer_info) {
                chg_cost_mat4_vertex[0][it] -= half_rc;
                chg_cost_mat4_vertex[it][0] -= half_rc;
            }
            for (auto it: customer_info) {
                chg_cost_mat4_vertex[0][it] += half_rc;
                chg_cost_mat4_vertex[it][0] += half_rc;
            }
        }
    }
}

void CVRP::priceBRC(BbNode *node, const std::vector<double> &pi_vector) {
    for (auto &brc: node->getBrCs()) {
        if (!brc.br_dir) {
            chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] = numeric_limits<float>::max();
            chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] =
                    numeric_limits<float>::max(); //do not use double since the number will overflow
        } else {
            chg_cost_mat4_vertex[brc.edge.first][brc.edge.second] -= pi_vector[brc.idx_br_c];
            chg_cost_mat4_vertex[brc.edge.second][brc.edge.first] -= pi_vector[brc.idx_br_c];
        }
    }
}

void CVRP::assignMemory() {
    route_in_pricing_assign = MAX_ROUTE_PRICING;
    aver_route_length = int(ceil(dim / num_vehicle)) + 5;
    mem4_pricing = aver_route_length * (route_in_pricing_assign);
    col_pool4_pricing = new int[mem4_pricing];

    label_assign = LABEL_ASSIGN / 2;
    reallocateLabel();

    arc_graph = new double *[dim];
    for (int i = 0; i < dim; ++i) {
        arc_graph[i] = new double[dim];
    }
    arc_graph_revised = new double *[dim]; //
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
    solver.freeEnv();
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
    if (node->getIfTerminated()) return;
    safe_solver(node->getSolver().updateModel())
    safe_solver(node->getSolver().getObjVal(&lp_val))
    node->getCurrentNodeVal() = lp_val;
    if (ceilTransformedNumberRelated(node->getCurrentNodeVal() - TOLERANCE) + TOLERANCE >= BaseBranching::ub) {
        node->getIfTerminated() = true;
        return;
    }
    vector<double> X(num_col);
    safe_solver(node->getSolver().getX(0, num_col, X.data()))

    pair<vector<SequenceInfo>, vector<double> > local_sol;
    pair<vector<SequenceInfo>, vector<double> > local_all_sol;

    double route_size = 0;
    int cnt = 0;
    for (int i = 0; i < num_col; ++i) {
        if (node->getCols()[i].col_seq.empty()) continue;
        if (X[i] > TOLERANCE) {
            local_all_sol.second.emplace_back(X[i]);
            local_all_sol.first.emplace_back(node->getCols()[i]);
            route_size += local_all_sol.first.back().col_seq.size();
            ++cnt;
            if (X[i] < 1 - TOLERANCE) {
                local_sol.second.emplace_back(X[i]);
                local_sol.first.emplace_back(node->getCols()[i]);
            }
        }
    }

    route_size /= cnt;
    route_size += 3; // two zeros and round up
    aver_route_length = max(aver_route_length, route_size);
    verbose_call(cout << "route_size= " << route_size << endl;)

    if (node->all_lp_sol == local_all_sol) {
        verbose_call(cout << "solution remains the same" << endl;)
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

    auto opt_size = (int) local_all_sol.first.size();
    for (int i = 0; i < opt_size; ++i) {
        auto &seq = local_all_sol.first[i].col_seq;
        double val = local_all_sol.second[i];
        int b4 = 0;
        for (auto j: seq) {
            arc_graph[b4][j] += val;
            b4 = j;
        }
        arc_graph[b4][0] += val;
    }

    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            arc_graph[i][j] += arc_graph[j][i];
            arc_graph_revised[i][j] = arc_graph[i][j];
        }
    }

    vector<pair<int, double> > tmp_first_dim;
    for (int i = 1; i < dim; ++i) if (X[i] > TOLERANCE) tmp_first_dim.emplace_back(i, X[i]);
    for (auto &i: tmp_first_dim) {
        auto &seq = node->getCols()[i.first].col_seq;
        if (seq.size() == 1) arc_graph_revised[0][seq[0]] -= i.second;
        else break;
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
                ml_call(MASTER_VALVE_ML != 0, MachineLearning::recordEdgeLongInfo(i, j))
            }
        }
    }
    RobustControl::recordSolutionArcs();
}

int CVRP::optimizeLPForOneIteration(BbNode *node, double prior_value, int lp_method) {
#if SOLVER_VRPTW == 1
  START_OVER:
#endif


    safe_solver(node->getSolver().reoptimize(lp_method))
    safe_solver(node->getSolver().getNumCol(&num_col))
    safe_solver(node->getSolver().getObjVal(&lp_val))


    double tol = max(TOLERANCE * lp_val, TOLERANCE);
    if (num_col > LP_COL_FINAL_LIMIT && abs(prior_value - lp_val) > tol)cleanIndexColForNode(node, false);
    vector<double> X(num_col);
    safe_solver(node->getSolver().getX(0, num_col, X.data()))

    pi4_labeling.resize(num_row);
    safe_solver(node->getSolver().getDual(0, num_row, pi4_labeling.data()))
    read_dual_sol_call(ReadDualSol::changeDualSol(pi4_labeling))
    combine_faster_call(CombineFasterSmoothing::updatePiLP1())
    combine_faster_call(FasterDual::changeModel4BetterDual())
    combine_faster_call(CombineFasterSmoothing::updatePiLP2())
    combine_faster_call(CombineFasterSmoothing::updateSmoothingState())

    faster_dual_call(FasterDual::changeModel4BetterDual())

    dual_smoothing_call(DualSmoothing::updatePiOut())
    dual_smoothing_call(DualSmoothing::updatePiPrice())
    dual_smoothing_call(
        DualSmoothing::
        calculateMostNegativeReducedCostInLP())
    dual_smoothing_call(DualSmoothing::recordLPValue())

    priceLabeling(node, pi4_labeling);

    if (tellIfInt(node, X)) {
        if (lp_val + TOLERANCE < BaseBranching::ub) {
            updateIPOptSol(node, X);
#if SOLVER_VRPTW == 1
	  bool if_feasible;
	  checkSolutionFeasibleByCapacity(if_feasible);
	  if (if_feasible) {
		BaseBranching::ub = lp_val;
	  } else {
		cout << "find an infeasible IP solution!" << endl;
		ip_opt_sol.clear();
		/**
		 * add rcc cuts and they are non-removable before enumeration!
		 */
		if (!BaseBranching::if_test_cg) {
		  verbose_call(cout << "add rcc cuts and they are non-removable before enumeration!" << endl;)
		  if_force_keep_rcc = true;
		  generateRCCs(node);
		  if_force_keep_rcc = false;
		  goto START_OVER;
		}
	  }
#else

            BaseBranching::ub = lp_val;
#endif
        }
    }
    return 0;
}

bool CVRP::tellIfInt(BbNode *node, const vector<double> &X) {
    bool is_integer = true;
    for (int i = 0; i < num_col; ++i) {
        if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
            is_integer = false;
            break;
        }
    }
    return is_integer;
}

void CVRP::cleanIndexColForNode(BbNode *node, bool if_only_rcfixing) {
    if (!if_allow_change_col) return;
    vector<int> col_idx;
    if (if_only_rcfixing) {
        int beg = 1; //rc fixing can safely delete the columns
        safe_solver(node->getSolver().updateModel())
        opt_gap = calculateOptimalGap(node);
        size_t numnzP;
        safe_solver(node->getSolver().XgetVars(&numnzP, nullptr, nullptr, nullptr, 0, num_col))
        vector<size_t> vbeg(num_col + 1);
        vector<int> vind(numnzP);
        vector<double> vval(numnzP);
        safe_solver(node->getSolver().XgetVars(&numnzP, vbeg.data(), vind.data(), vval.data(), 0, num_col))
        vector<Eigen::Triplet<double> > triplet(numnzP);
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
        safe_solver(node->getSolver().getObj(0, num_col, cost.data()))
        Eigen::Map<Eigen::RowVectorXd> c(cost.data(), num_col);
        RowVectorXd rc;
        safe_eigen(rc = c - pi * A;)
        for (int i = beg; i < num_col; ++i) if (rc[i] > opt_gap) col_idx.emplace_back(i);
    } else {
        int beg = if_in_enu_state ? 1 : dim;
        safe_solver(node->getSolver().reoptimize(SOLVER_PRIMAL_SIMPLEX))
        vector<double> rc(num_col);
        safe_solver(node->getSolver().getRC(0, num_col, rc.data()))
        vector<double> rc_copy(rc.data() + beg, rc.data() + num_col);
        int n = int((num_col - beg) * COL_KEEP_FRAC);
        nth_element(rc_copy.begin(), rc_copy.begin() + n, rc_copy.end());
        double threshold = max(rc_copy[n], TOLERANCE);
        for (int i = beg; i < num_col; ++i) if (rc[i] > threshold) col_idx.emplace_back(i);
    }
    rmLPCols(node, col_idx);
    safe_solver(node->getSolver().reoptimize(SOLVER_PRIMAL_SIMPLEX))
}

void CVRP::getNewConstraintCoefficientByEdge(BbNode *node,
                                             const std::pair<int, int> &edge,
                                             std::vector<int> &ind,
                                             std::vector<double> &val) {
    int size = node->edge_col_map[edge].size();
    ind.resize(size);
    val.resize(size);
    int cnt = 0;
    for (auto &col: node->edge_col_map[edge]) {
        ind[cnt] = col.first;
        val[cnt++] = col.second;
    }
}

void CVRP::getCoefficientRCC(BbNode *node, Rcc &rcc, std::vector<int> &ind, std::vector<double> &val) {
    auto &map = node->edge_col_map;
    vector<double> sup(num_col, 0);
    if (rcc.form_rcc) {
        auto &customer_info = rcc.info_rcc_customer;
        for (auto i = customer_info.begin(); i != customer_info.end(); ++i) {
            auto j = i;
            ++j;
            for (; j != customer_info.end(); ++j) {
                auto pr = *i < *j ? make_pair(*i, *j) : make_pair(*j, *i);
                for (auto &col: map[pr]) {
                    sup[col.first] += col.second;
                }
            }
        }
    } else {
        auto &customer_info = rcc.info_rcc_customer;
        auto &outside_customer_info = rcc.info_rcc_outside_customer;
        int cnt = 0;
        for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
            auto j = i;
            ++j;
            for (; j != outside_customer_info.end(); ++j) {
                auto pr = *i < *j ? make_pair(*i, *j) : make_pair(*j, *i);
                for (auto &col: map[pr]) {
                    sup[col.first] += col.second;
                }
            }
        }

        for (auto customer_it: outside_customer_info) {
            for (auto &col: map[{0, customer_it}]) {
                sup[col.first] += 0.5 * col.second;
            }
        }
        for (auto customer_it: customer_info) {
            for (auto &col: map[{0, customer_it}]) {
                sup[col.first] -= 0.5 * col.second;
            }
        }
    }

    ind.clear();
    val.clear();
    sup[0] = rcc.rhs;
    for (int i = 0; i < num_col; ++i) {
        if (sup[0] != 0) {
            ind.emplace_back(i);
            val.emplace_back(sup[i]);
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
    return (BaseBranching::ub - node->getCurrentNodeVal() + round_up_tolerance);
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
#if VERBOSE_MODE == 1
        cout << "arc elimination is banned!" << endl;
#endif
        goto QUIT;
    }
    if (final_decision_4_arc_elimination) goto QUIT;
    if (abs(old_ub - BaseBranching::ub) > TOLERANCE) {
#if VERBOSE_MODE == 1
        cout << "BaseBranching::ub is updated from " << old_ub << " to " << BaseBranching::ub
                << " and perform arc elimination!" << endl;
#endif
        old_ub = BaseBranching::ub;
        final_decision_4_arc_elimination = true;
        goto QUIT;
    }
    now_gap = (BaseBranching::ub - node->getCurrentNodeVal()) / BaseBranching::ub;
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
    if_force_enumeration_suc = false;
    if (final_decision_4_enumeration) {
#if VERBOSE_MODE == 1
        cout << "enumeration is forced!" << endl;
#endif
        if (!if_arc_elimination_succeed) {
#if VERBOSE_MODE == 1
            cout << "arc elimination failed! enumeration is skipped!" << endl;
#endif
            final_decision_4_enumeration = false;
        }
        goto QUIT;
    }
    if (!if_arc_elimination_succeed && !if_arc_elimination_tried_but_failed) {
        goto QUIT;
    }

    gap = (BaseBranching::ub - node->getCurrentNodeVal()) / BaseBranching::ub;
    if (gap > Config::MaxGap2TryEnumeration || gap > getGapStdTryEnumeration()) {
        goto QUIT;
    }
    if (gap < max_enumeration_success_gap) {
        final_decision_4_enumeration = true;
        if_force_enumeration_suc = true;
        goto QUIT;
    }
    if (node->num_forward_bucket_arcs < max_bucket_arc_suc_enumeration.first)
        if_force_enumeration_suc = true;
    symmetry_prohibit_call(
        if (node->num_backward_bucket_arcs < max_bucket_arc_suc_enumeration.second) if_force_enumeration_suc = true;)

    if (node->num_forward_bucket_arcs > min_bucket_arc_fail_enumeration.first)
        goto QUIT;
    symmetry_prohibit_call(if (node->num_backward_bucket_arcs > min_bucket_arc_fail_enumeration.second)
        goto QUIT;)
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
    RobustControl::if_ever_roll_back = true;

    vector<int> delete_cuts(num_row - old_num);
    iota(delete_cuts.begin(), delete_cuts.end(), old_num);
    deleteNonactiveCuts(node, delete_cuts);

    for (auto &i: reset_cut_mem) {
        auto &cut = node->r1cs[get<0>(i)];
        cut.arc_mem = get<2>(i);
    }

    reset_cut_mem.clear();
    getVCutMapLP(node);

    if (rollback_solver.model) {
        node->getSolver().freeModel();
        node->getSolver().model = rollback_solver.copyModel();
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        safe_solver(node->getSolver().getNumCol(&num_col))
        safe_solver(node->getSolver().getNumRow(&num_row))
        node->getCols() = rollback_cols;
        recordOptimalColumn(node);
    } else {
        cout << "re-run cg to get the complete model!" << endl;
        cout << "also we decrease the soft time!" << endl;
        Config::CutGenTimeThresholdInPricingInitial /= Config::SoftTimeDecreaseFactorInPricing;
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        if_roll_back = false;
        if_short_memory = false;
        force_not_rollback = true;
        solveLPInLabeling(node);
        force_not_rollback = false;
        verbose_call(cout << "copy model\n";)
        if (rollback_solver.model)rollback_solver.freeModel();
        rollback_solver.model = node->getSolver().copyModel();
        rollback_cols = node->getCols();
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
    int cnt = 0;
    for (int i = 1; i < dim; ++i) {
        ResTuple res{};
        if (!increaseMainResourceConsumption({}, res, 0, i)) continue;
        bucket.bucket_arcs[cnt++] = i;
    }
    bucket.bucket_arcs.resize(cnt);

    max_num_forward_graph_arc = num_buckets_per_vertex * (real_dim - 1) * real_dim;
    node->num_forward_bucket_arcs = max_num_forward_graph_arc;
#ifdef SYMMETRY_PROHIBIT
  for (int i = 1; i < dim; ++i) {
	for (int b = 0; b < num_buckets_per_vertex; ++b) {
	  node->all_backward_buckets[i][b] = node->all_forward_buckets[i][0];
	}
  }
  auto &bucket_back = node->all_backward_buckets[0][0];
  cnt = 0;
  bucket_back.bucket_arcs.resize(real_dim);
  for (int i = 1; i < dim; ++i) {
	ResTuple res{};
	if (!decreaseMainResourceConsumption(resource, res, 0, i)) continue;
	bucket_back.bucket_arcs[cnt++] = i;
  }
  bucket_back.bucket_arcs.resize(cnt);
  max_num_backward_graph_arc = num_buckets_per_vertex * (real_dim - 1) * real_dim;
  node->num_backward_bucket_arcs = max_num_backward_graph_arc;
#endif
    getTopologicalOrder(node);
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
    } else if (strList.size() == 4) {
        Config::ub = strtod(strList[1].c_str(), nullptr);
        Config::col_pool_path = strList[2].c_str();
        Config::tree_path = strList[3].c_str();
    }
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
        } else if (arg == "-c" && i + 1 < argc) {
            stringstream convert(argv[++i]);
            if (!(convert >> Config::col_pool_path)) {
                printf("Invalid string: %s\n", argv[i]);
                continue;
            }
            cout << Config::col_pool_path << endl;
        } else if (arg == "-t" && i + 1 < argc) {
            stringstream convert(argv[++i]);
            if (!(convert >> Config::tree_path)) {
                printf("Invalid string: %s\n", argv[i]);
                continue;
            }
            cout << Config::tree_path << endl;
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
    vector<double> val(num_col, 0);
    for (int i = 0; i < vind.size(); ++i) {
        val[vind[i]] = vval[i];
    }
    vector<int> cind(num_col, row_idx);
    safe_solver(local_solver.changeCoeffs(num_col, cind.data(), ind.data(), val.data()))
    safe_solver(local_solver.setRhs(row_idx, 1, &sense, &rhs))
    int vind2 = 0;
    safe_solver(local_solver.XchangeCoeffs(1, &row_idx, &vind2, &rhs))
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
    newResource = nowResource + resource_across_arcs_in_forward_sense[start][end];
    if (tellResTupleRelations<'p'>(newResource, ub4_vertex[end])) return false;
    if (newResource.first_res < lb4_vertex[end].first_res) newResource.first_res = lb4_vertex[end].first_res;
#ifdef USE_TWO_RESOURCE
  if (newResource.second_res < lb4_vertex[end].second_res) newResource.second_res = lb4_vertex[end].second_res;
#endif
    if (tellResTupleRelations<'c'>(newResource, resource_across_arcs_in_forward_sense[end][0])) return false;
    return true;
}

bool CVRP::decreaseMainResourceConsumption(const ResTuple &nowResource,
                                           ResTuple &newResource, int start, int end) {
    newResource = nowResource - resource_across_arcs_in_backward_sense[start][end];
    if (tellResTupleRelations<'p'>(lb4_vertex[end], newResource)) return false;
    if (newResource.first_res > ub4_vertex[end].first_res) newResource.first_res = ub4_vertex[end].first_res;
#ifdef USE_TWO_RESOURCE
  if (newResource.second_res > ub4_vertex[end].second_res) newResource.second_res = ub4_vertex[end].second_res;
#endif
    /**
     * not type c
     */
    if (tellResTupleRelations<'p'>({}, newResource - resource_across_arcs_in_backward_sense[end][0])) return false;
    return true;
}

void CVRP::lighterHeuristic(BbNode *node, int &num) {
    num = generateColumnsByLighterHeuristic<
#ifdef SYMMETRY_PROHIBIT
	  false
#else
        true
#endif
    >(node);
}

void CVRP::heavierHeuristic(BbNode *node, int &num) {
    num = generateColumnsByHeavierHeuristic<
#ifdef SYMMETRY_PROHIBIT
	  false
#else
        true
#endif
    >(node);
}

void CVRP::exactLabeling(BbNode *node, int &num) {
    num = generateColsByBidir<
#ifdef SYMMETRY_PROHIBIT
	  false
#else
        true
#endif
    >(node);
}

bool CVRP::runColumnGenerationType(BbNode *node, int mode) {
    if (node->getIfTerminated()) return true;

    combine_faster_call(DualSmoothing::reset())
    dual_smoothing_call(DualSmoothing::reset())
    read_dual_sol_call(ReadDualSol::resetCoeff())

    int status = 0;
    bool switch_is_on = true;
#if VERBOSE_MODE == 1
    switch (mode) {
        case 1: cout << "LighterHeur phase begin...\n";
            break;
        case 2: cout << "HeavierHeur phase begin...\n";
            break;
        case 3: cout << "Exact phase begin...\n";
            break;
        default: throw runtime_error("Error: None of these modes are used!");
    }
#endif
    int iter = 0, ccnt = 0, old_ncol = num_col, tag;
    double eps_CG = 0, eps_LP = 0, b4_node_val = node->getCurrentNodeVal();
    auto beg = chrono::high_resolution_clock::now();
    auto end = beg;
    bool old_force;


    while (true) {
        faster_dual_call(if (mode == 3) FasterDual::reviseSubPricingModel();)
        combine_faster_call(if (mode == 3) FasterDual::reviseSubPricingModel();)
        beg = chrono::high_resolution_clock::now();
        switch (mode) {
            case 1: lighterHeuristic(node, ccnt);
                ++iter;
                break;
            case 2: heavierHeuristic(node, ccnt);
                ++iter;
                break;
            case 3: exactLabeling(node, ccnt);
                ++iter;
                if (if_roll_back) goto QUIT;
                break;
        }

        end = chrono::high_resolution_clock::now();
        time_pricing4_iter = chrono::duration<double>(end - beg).count();
        eps_CG += time_pricing4_iter;

        read_dual_sol_call(ReadDualSol::checkIfMisPrice(ccnt))
        combine_faster_call(DualSmoothing::tellIfMisPrice(ccnt))
        dual_smoothing_call(DualSmoothing::tellIfMisPrice(ccnt))

        if ((iter % PRINT_LABELING_STEP_SIZE) == 0) {
#if VERBOSE_MODE == 2
	  goto OUTSIDE;
#endif
        REPORT:
            BaseBranching::glo_end = chrono::high_resolution_clock::now();
            BaseBranching::glo_eps = chrono::duration<double>(BaseBranching::glo_end - BaseBranching::glo_beg).count();
            printInfoLabeling(iter, num_col - old_ncol, num_col, num_row, eps_LP,
                              eps_CG, BaseBranching::glo_eps,
                              lp_val, BaseBranching::lb, BaseBranching::ub);
            if (!switch_is_on) {
                goto QUIT;
            }
            eps_CG = 0;
            eps_LP = 0;
            old_ncol = num_col;
        }

#if VERBOSE_MODE == 2
	OUTSIDE:
#endif
        if (ccnt == 0) {
            if (mode == 3) optimal_dual_vector = pi4_labeling;
            break;
        }
        if (node->getIfTerminated()) break;

        beg = chrono::high_resolution_clock::now();
        tag = optimizeLPForOneIteration(node, b4_node_val);
        b4_node_val = lp_val;
        end = chrono::high_resolution_clock::now();

        time_resolve_lp4_iter = chrono::duration<double>(end - beg).count();
        eps_LP += time_resolve_lp4_iter;

        if (tag) {
            status = 1; //QUIT
            goto QUIT;
        }

        if (BaseBranching::if_test_cg) {
            if (mode == 1 && double(num_col) / LP_COL_FINAL_LIMIT > HEURISTIC_LIGHT_TESTING_MAX_COLUMN_RATIO) {
                verbose_call(
                    cout << "LighterHeur testing is terminated since the number of columns is too large!" << endl;)
                break;
            } else if (mode == 2 && double(num_col) / LP_COL_FINAL_LIMIT > HEURISTIC_HEAVY_TESTING_MAX_COLUMN_RATIO) {
                verbose_call(
                    cout << "HeavierHeur testing is terminated since the number of columns is too large!" << endl;)
                break;
            }
        }
    }

    if (
        switch_is_on && (iter
                         % PRINT_LABELING_STEP_SIZE)) {
        switch_is_on = false;
        goto REPORT;
    }
QUIT:
    faster_dual_call(FasterDual::freeSubPricingModel())
    combine_faster_call(FasterDual::freeSubPricingModel())
    return status;
}

double CVRP::transformCost(double x) {
    return std::floor(x + 0.5);
}

void CVRP::assignInitialLabelingMemory() const {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
        int min_num_labels = numeric_limits<int>::max();
        for (int i = 0; i < dim; ++i) {
            if (label_array_in_forward_sense[i][b].size() && min_num_labels > label_array_in_forward_sense[i][b].
                size()) {
                min_num_labels = label_array_in_forward_sense[i][b].size();
            }
        }
        if (min_num_labels == numeric_limits<int>::max()) {
            min_num_labels = 1;
        }
        for (int i = 0; i < dim; ++i) {
            int hard_max = max((int) label_array_in_forward_sense[i][b].size(), min_num_labels) * FACTOR_NUM_LABEL;
            hard_max = int(pow_self(2.0, int(ceil(log(hard_max) / log(2)) + TOLERANCE)));
            if_exist_extra_labels_in_forward_sense[i][b].first.resize(hard_max);
        }
#ifdef SYMMETRY_PROHIBIT
	min_num_labels = numeric_limits<int>::max();
	for (int i = 0; i < dim; ++i) {
	  if (label_array_in_backward_sense[i][b].size() && min_num_labels > label_array_in_backward_sense[i][b].size()) {
		min_num_labels = label_array_in_backward_sense[i][b].size();
	  }
	}
	if (min_num_labels == numeric_limits<int>::max()) {
	  min_num_labels = 1;
	}
	for (int i = 0; i < dim; ++i) {
	  int hard_max = max((int)label_array_in_backward_sense[i][b].size(), min_num_labels) * FACTOR_NUM_LABEL;
	  hard_max = (int)(pow_self(2.0, (int)ceil(log(hard_max) / log(2))) + TOLERANCE);
	  if_exist_extra_labels_in_backward_sense[i][b].first.resize(hard_max);
	}
#endif
    }
}

void CVRP::getEdgeInfo(BbNode *node, bool if_br) const {
    double **local_arc_graph;
    if (if_br) local_arc_graph = arc_graph_revised;
    else local_arc_graph = arc_graph;
    int numEdge = num_edge + 1; //as required
    if (node->getEdgeHead().size() < numEdge) {
        node->allocateMem(numEdge);
    }
    node->getNumEdges() = 0;
    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            if (local_arc_graph[i][j] > TOLERANCE) {
                ++node->getNumEdges();
                node->getEdgeTail()[node->getNumEdges()] = i;
                node->getEdgeHead()[node->getNumEdges()] = j;
                node->getEdgeValue()[node->getNumEdges()] = local_arc_graph[i][j];
            }
        }
    }
}

void CVRP::regenerateGraphBucket(BbNode *node) {
    if (step_size / 2 < 1) return;
    if (step_size % 2) {
        cout << "regenerateGraphBucket is banned since step_size is odd!" << endl;
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
        label_array_in_forward_sense[i] = new ListLabel[num_buckets_per_vertex];
    }

    for (int i = 0; i < dim; ++i) {
        delete[]if_exist_extra_labels_in_forward_sense[i];
    }
    delete[]if_exist_extra_labels_in_forward_sense;
    if_exist_extra_labels_in_forward_sense = new_IfExistExtraLabelsInForwardSense;

#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < dim; ++i) {
	delete[]label_array_in_backward_sense[i];
	label_array_in_backward_sense[i] = new ListLabel[num_buckets_per_vertex];
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
#if VERBOSE_MODE == 1
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
        for (auto n: tmp3) {
            ng_mem4_vertex[i].set(n);
        }
    }
    auto last = max_labeling_time_for_node_so_far;
    int round = 0;
    bool if_empty;
    double old_val = node->getCurrentNodeVal();
    double cost_time = 0;
    double std;
    while (true) {
        ++round;
        cout << "DSSR: " << round << endl;
        findNGMemorySets(node, if_empty);
        if (if_empty) break;
        auto beg = chrono::high_resolution_clock::now();
        force_not_rollback = true;
        solveLPInLabeling(node);
        force_not_rollback = false;
        auto end = chrono::high_resolution_clock::now();
        auto eps = duration<double>(end - beg).count();
        if (eps > cost_time) cost_time = eps;
        if (max_labeling_time_for_node_so_far > last && node->getCurrentNodeVal() + TOLERANCE < old_val) {
#if VERBOSE_MODE == 1
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
        if (node->getCurrentNodeVal() + TOLERANCE > old_val) break;
    }
    if (round == 1 && if_empty) {
        ng_mem4_vertex = old_mem;
    }
    if (if_empty) goto QUIT;
    round = 0;
    old_val = node->getCurrentNodeVal();
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
            if_roll_back = true;
        } else if (eps > Config::NGAugTimeSoftThresholdFactor * cost_time) {
            if_tail_off = true;
        }
        verbose_call(cout << "eps= " << eps << endl;
            cout << "cost_time= " << cost_time << endl;)
        if (if_tail_off) break;
        else if (if_roll_back) {
            cout << "NGAugmentation: Reach the hard limit! We use ng rollback!" << endl;
            ng_mem4_vertex = old_mem;
            deleteColumnByNGMemory(node, dim);
            force_not_rollback = true;
            solveLPInLabeling(node);
            force_not_rollback = false;
            return;
        }
        old_mem = ng_mem4_vertex;
        double tmp = abs(node->getCurrentNodeVal() - old_val) / node->getCurrentNodeVal();
        if (tmp > std) std = tmp;
        else if (tmp < std * Config::NGAugTailOff) break;
        else old_val = node->getCurrentNodeVal();
    }
QUIT:
    return;
}

void CVRP::findNGMemorySets(BbNode *node, bool &if_empty) {
    if (node->index) throw runtime_error("findNGMemorySets called when node->index!=0");
    unordered_map<int, vector<pair<int, yzzLong> > > map_size_cycle;
    vector<size_t> tmp(dim);
    vector<pair<int, double> > opt_col_idx(node->only_frac_sol.first.size());
    for (int i = 0; i < opt_col_idx.size(); ++i) {
        opt_col_idx[i] = make_pair(i, node->only_frac_sol.second[i]);
    }
    sort(opt_col_idx.begin(), opt_col_idx.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
        return a.second > b.second;
    });

    int re_cols = 0;
    for (auto &pr: opt_col_idx) {
        fill(tmp.begin(), tmp.end(), -1);
        bool if_re_col = false;
        auto &route = node->only_frac_sol.first[pr.first].col_seq;
        for (int i = 0; i < route.size(); ++i) {
            int curr_node = route[i];
            if (tmp[curr_node] != -1) {
                if_re_col = true;
                auto length = int(i - tmp[curr_node] - 1);
                map_size_cycle[length].emplace_back(curr_node, 0);
                auto &mem = map_size_cycle[length].back().second;
                for (auto k = tmp[curr_node] + 1; k < i; ++k) {
                    mem.set(route[k]);
                }
            }
            tmp[curr_node] = i; //not in else!
        }
        if (if_re_col) {
            ++re_cols;
            if (re_cols > Config::MaxNumColsInNGAug) break;
        }
    }
    vector<pair<int, vector<pair<int, yzzLong> > > > size_cycle(map_size_cycle.begin(), map_size_cycle.end());
    sort(size_cycle.begin(), size_cycle.end(), [](const pair<int, vector<pair<int, yzzLong> > > &a,
                                                  const pair<int, vector<pair<int, yzzLong> > > &b) {
        return a.first < b.first;
    });
    if (size_cycle.empty()) {
        verbose_call(cout << "no small cycles are found!" << endl;)
        if_empty = true;
        return;
    } else if_empty = false;
    for (auto &i: size_cycle[0].second) {
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
        for (auto &j: size_cycle[i].second) {
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

    for (int i = start; i < num_col; ++i) {
        local_pi = 0;
        for (auto current_node: node->getCols()[i].col_seq) {
            if (local_pi[current_node]) {
                col_idx.emplace_back(i);
                break;
            }
            if (!if_full_mem)local_pi = local_pi & ng_mem4_vertex[current_node];
            local_pi.set(current_node);
        }
    }
    rmLPCols(node, col_idx);
    safe_solver(node->getSolver().reoptimize())
    safe_solver(node->getSolver().getObjVal(&lp_val))
#if VERBOSE_MODE == 1
    if (if_full_mem) {
        cout << "after clean non-ele routes lpval= " << lp_val << endl;
    } else {
        cout << "after clean ng routes lpval= " << lp_val << endl;
    }
#endif
}

double CVRP::calculateGapImprovement(double nowVal, double b4Val) const {
    double guess_UB = BaseBranching::ub > 2 * b4Val ? b4Val * 1.008 : BaseBranching::ub;
    if (guess_UB - b4Val < TOLERANCE) return 0;
    return (nowVal - b4Val) / (guess_UB - b4Val);
}

void CVRP::deleteArcByFalseBranchConstraint(Bucket **buckets, const std::pair<int, int> &edge) const {
    int i = edge.first, j = edge.second;
    bool if_rep = false;
AGAIN:
    for (int bin = 0; bin < num_buckets_per_vertex; ++bin) {
        auto if_find = std::find(buckets[i][bin].bucket_arcs.begin(),
                                 buckets[i][bin].bucket_arcs.end(), j);
        if (if_find == buckets[i][bin].bucket_arcs.end()) {
            auto iff = std::find_if(buckets[i][bin].jump_arcs.begin(),
                                    buckets[i][bin].jump_arcs.end(),
                                    [&](const pair<double, int> &p) { return p.second == j; });
            if (iff != buckets[i][bin].jump_arcs.end()) {
                buckets[i][bin].jump_arcs.erase(iff);
            }
        } else {
            buckets[i][bin].bucket_arcs.erase(if_find);
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
    for (auto &rcc: node->rccs) {
        if (if_takeup[rcc.idx_rcc]) throw std::runtime_error("Error in rccs");
        if_takeup[rcc.idx_rcc] = true;
    }
    for (auto &r1c: node->r1cs) {
        if (if_takeup[r1c.idx_r1c]) throw std::runtime_error("Error in r1cs");
        if_takeup[r1c.idx_r1c] = true;
    }
    for (auto &brc: node->getBrCs()) {
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

    safe_solver(node->getSolver().updateModel())
    for (int i = 0; i < num_row; i++) {
        double cof, rhs;
        safe_solver(node->getSolver().getCoeff(i, 0, cof))
        safe_solver(node->getSolver().getRhs(i, 1, &rhs))
        if (cof != rhs) {
            cout << "i: " << i << endl;
            for (auto &rcc: node->rccs) {
                if (i == rcc.idx_rcc) {
                    cout << "rcc: " << rcc.rhs << endl;
                    break;
                }
            }
            for (auto &r1c: node->r1cs) {
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
                if (!current_node) continue;
                tmp.emplace_back(current_node);
            }
            ip_opt_sol.emplace_back(std::move(tmp));
        }
    }
    fin.close();
    if (ip_opt_sol.empty() && if_force) {
        throw std::runtime_error("Error in No solution!");
    }
    cout << "optval: " << optval << " BaseBranching::ub: " << BaseBranching::ub << endl;
    if (abs(optval - BaseBranching::ub) < TOLERANCE || optval < BaseBranching::ub) {
        cout << "BaseBranching::ub is updated to: " << optval << endl;
    } else {
        cout << "we abandoned the old BaseBranching::ub (though better, we want the solution) and update to " << optval
                << endl;
    }
    BaseBranching::ub = optval;
}

void CVRP::rmLPCols(BbNode *node, const std::vector<int> &col_idx) {
    if (col_idx.empty()) return;
    vector<int> sort_col_idx = col_idx;
    if (!is_sorted(col_idx.begin(), col_idx.end())) {
        sort(sort_col_idx.begin(), sort_col_idx.end());
    }

    int delta = 0;
    auto stop_sign = sort_col_idx.end() - 1;
    for (auto i = sort_col_idx.begin(); i < stop_sign; ++i) {
        ++delta;
        for (int j = *i + 1; j < *(i + 1); ++j) node->getCols()[j - delta] = node->getCols()[j];
    }
    ++delta;
    for (int j = *stop_sign + 1; j < num_col; ++j) node->getCols()[j - delta] = node->getCols()[j];

    safe_solver(node->getSolver().delVars(sort_col_idx.size(), sort_col_idx.data()))
    safe_solver(node->getSolver().updateModel())
    safe_solver(node->getSolver().getNumCol(&num_col))

    node->getCols().resize(num_col);
    faster_dual_call(FasterDual::deleteCols4SubPricingModel(sort_col_idx);)
    dual_smoothing_call(DualSmoothing::deleteCol4LPMatrix(sort_col_idx))
    combine_faster_call({
        FasterDual::deleteCols4SubPricingModel(sort_col_idx);
        DualSmoothing::deleteCol4LPMatrix(sort_col_idx);
        })
}

void CVRP::updateEdgeColMap(BbNode *node, bool if_br) {
    /**
    * revise the edge_col_map
     * 1. the first column is not counted
     * 2. the 0-1-0, only counts one edge if it is in br
    */
    std::unordered_map<std::pair<int, int>, std::unordered_map<int, int>, PairHasher>
            sup_map{};
    for (int i = 1; i < num_col; ++i) {
        int b4 = 0;
        if (node->getCols()[i].col_seq.size() == 1) {
            b4 = node->getCols()[i].col_seq[0];
            if (!if_br) ++sup_map[{0, b4}][i];
        } else {
            for (auto current_node: node->getCols()[i].col_seq) {
                auto pair = b4 < current_node ? make_pair(b4, current_node) : make_pair(current_node, b4);
                ++sup_map[pair][i];
                b4 = current_node;
            }
        }
        ++sup_map[{0, b4}][i];
    }
    auto &map = node->edge_col_map;
    map.clear();
    for (auto &pr: sup_map) {
        map[pr.first] = vector<pair<int, int> >(pr.second.begin(), pr.second.end());
    }
}

void CVRP::updateIPOptSol(BbNode *node, const std::vector<double> &X) {
    ip_opt_sol.clear();
    for (int i = 0; i < num_col; ++i) {
        if (node->getCols()[i].col_seq.empty()) continue;
        if (X[i] > 0.5) {
            ip_opt_sol.emplace_back(node->getCols()[i].col_seq);
        }
    }
}

void CVRP::tellIfEnterMIP(BbNode *node) {
    if (node->getIfTerminated()) {
    HERE:
        BaseBranching::freeNode();
    } else if (num_col + node->valid_size <= MaxNumRoute4Mip) {
        terminateByMIP(node);
        if (node->getIfTerminated()) {
            goto HERE;
        }
    }
}

void CVRP::resizePoolWarning(size_t &pricing_warning) {
    if (pool_beg4_pricing >= pricing_warning) {
        cout << SMALL_PHASE_SEPARATION;
        cout << "Warning: the pricing pool is almost full!" << endl;
        cout << "pool_beg4_pricing=" << pool_beg4_pricing << endl;
        cout << "mem4_pricing=" << mem4_pricing << endl;
        cout << "we reallocate the pricing pool!" << endl;
        reallocatePricingPool();
        pricing_warning = (size_t) (0.9 * (double) mem4_pricing);
        cout << "the new mem4_pricing=" << mem4_pricing << endl;
        cout << SMALL_PHASE_SEPARATION;
    }
}

void self_mkdir(const std::string &path1) {
    namespace fs = std::filesystem;
    if (fs::exists(path1) && fs::is_directory(path1)) {
        cout << path1 + " already exists" << endl;
    } else {
        fs::create_directory(path1);
        cout << path1 + " created" << endl;
    }
}

void CVRP::prepareSolveEnuTree(BbNode *node) {
    if_tail_off = false;
    if_in_enu_state = true;
    node->if_just_enter_enu = true;
    for (auto &r1c: node->r1cs) {
        r1c.arc_mem.clear();
    }
    safe_solver(node->getSolver().updateModel());
    safe_solver(node->getSolver().getNumRow(&num_row));
    safe_solver(node->getSolver().getNumCol(&num_col));
    createBasicMatrix(node);
}


bool &CVRP::getIfInEnuState() {
    return if_in_enu_state;
}


void CVRP::startSolveNode(BbNode *node) {
    if (if_in_enu_state) {
        max_num_enu_col_pool = node->size_enumeration_col_pool;
        solveLPByInspection(node, false, false, true);
    } else {
        if (rollback_solver.model)rollback_solver.freeModel();
        max_labeling_time_for_node_so_far = 0;
        node->br_value_improved = 0; //need to be here just before enumeration!
        if_tail_off = false;
        force_not_rollback = true;
        solveLPInLabeling(node, true, true, true);
        force_not_rollback = false;
    }
}

void CVRP::postSolveNode(BbNode *&node) {
    if (if_in_enu_state) return;
    if (node->getTreeLevel() == 0) {
        final_decision_4_arc_elimination = true;
    }
#if SETTING == 1
    eliminateArcs(node);
    if (node->getTreeLevel() != 0) {
        enumerateMIP(node);
        if (!node) goto QUIT;
    }
    cleanIndexColForNode(node, true);
#elif SETTING == 2
    if (node->getTreeLevel() == 0) {
        eliminateArcs(node);
        cleanIndexColForNode(node, true);
    }
#endif
QUIT:;
}

void CVRP::finishSolveNode() {
    if_tail_off = false;
}

void CVRP::cuttingAfterEnu(BbNode *&node) {
#ifndef DELUXING_APPLIED
    tellIfEnterMIP(node);
    if (!node) goto QUIT;
#endif

    if (num_col + node->valid_size > MaxNumRoute4Mip) separateHybridCuts(node);

#ifndef DELUXING_APPLIED
    tellIfEnterMIP(node);
    if (!node) goto QUIT;
#else
    int ncol = num_col + node->valid_size;
    if (ncol <= 40000) {
        if (ncol > 1000) {
            cout << "apply DELUXING..." << endl;
            applyRCF(node, DELUXING_APPLIED, true);
            solveLPByInspection(node, false, false, true);
        }
        terminateByMIP(node);
    }
    if (node->getIfTerminated()) {
        BaseBranching::freeNode();
        goto QUIT;
    }
#endif
#ifdef WRITE_ENUMERATION_TREES
    WriteEnumerationTree::updateNode(node);
    WriteEnumerationTree::writeEnuTree();
    freeNode();
#endif
QUIT:;
}

void CVRP::enterCuttingPhase(BbNode *&node) {
    if (if_in_enu_state) {
        cuttingAfterEnu(node);
        if (!node) goto QUIT;
        regenerateEnumMat(node, nullptr, true);
        recordOptimalColumn(node);
        getEdgeInfo(node, true);
    } else {
        separateHybridCuts(node);
        if (!node) goto QUIT;
        if (node->getIfTerminated()) {
            BaseBranching::freeNode();
            goto QUIT;
        }
        if (num_col > LP_COL_FINAL_LIMIT) {
            cleanIndexColForNode(node, false);
        }
        recordOptimalColumn(node);
    }
QUIT:;
}


int &CVRP::getNumCol() {
    return num_col;
}

int &CVRP::getNumRow() {
    return num_row;
}


double &CVRP::getMaxEnumerationSuccessGap() {
    return max_enumeration_success_gap;
}

void CVRP::preprocess4TestCG(BbNode *node) {
    if (!getIfInEnuState()) {
        force_not_rollback = true;
        if_force_not_regenerate_bucket_graph = true;
    } else {
        sparseRowMatrixXd mat(1, node->size_enumeration_col_pool);
        node->matrix_in_enumeration.push_back(std::move(mat));
    }
}

void CVRP::postprocess4TestCG(BbNode *node) {
    if (!getIfInEnuState()) {
        force_not_rollback = false;
        if_force_not_regenerate_bucket_graph = false;
    } else {
        node->matrix_in_enumeration.pop_back();
    }
}

bool &CVRP::getIfAllowChangeCol() {
    return if_allow_change_col;
}
