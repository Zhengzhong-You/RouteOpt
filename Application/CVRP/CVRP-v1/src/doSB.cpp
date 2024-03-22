#include "CVRP.hpp"

using namespace std;
using namespace chrono;

bool CVRP::solveSBModel() {

  cout << "\n<solveSBModel>\n\n";

  lateProcessing();
  glo_beg = std::chrono::high_resolution_clock::now();

  buildModel();

  glo_end = chrono::high_resolution_clock::now();
  glo_eps = duration<double>(glo_end - glo_beg).count();
  cout << "Build model...  create art_var= " << dim << "  main_cstr= " << real_dim << "  lb= " << lb
       << "  ub= " << ub << "  et= " << glo_eps << endl;
  cout << BIG_PHASE_SEPARATION;

#ifdef WRITE_ENUMERATION_TREES
  enumeration_col_idx.reserve(size_t(2e7));
#endif

  while (!bbt.empty()) {
    glo_end = std::chrono::high_resolution_clock::now();
    glo_eps = duration<double>(glo_end - glo_beg).count();

	BbNode *node = bbt.top();
	updateLowerBound(node->value);
    bbt.pop();
    ++num_explored_nodes;

    if (lb_transformed >= ub) {
      lb_transformed = ub;
      cout << "Pruned! Optimality has been proven!  lb=ub= " << ub << endl;
      cout << SMALL_PHASE_SEPARATION;
      delete node;
      goto DELETE;
    }

    if (glo_eps > GLOBAL_TIME_LIMIT) {
      cout << SMALL_PHASE_SEPARATION;
      cout << "Ins= " << file_name << "  cap= " << cap
           << "  shut down due to reaching lm. gt.= " << GLOBAL_TIME_LIMIT << endl;
      cout << "lb= " << lb << "  ub= " << ub << "  gap(ub-lb/ub)= "
           << (ub - lb) / (ub) * 100 << "%  gt= "
           << glo_eps << " nd= "
           << num_explored_nodes << "  br= " << num_br << endl;
      delete node;
      goto DELETE;
    }

    safe_solver(node->solver.updateModel())
    safe_solver(node->solver.getNumRow(&num_row))
    safe_solver(node->solver.getNumCol(&num_col))

#if VERBOSE_MODE==1
    cout << "num_buckets_per_vertex= " << num_buckets_per_vertex << " step_size= " << step_size << endl;
#endif
    getVCutMapLP(node);

    cout << BIG_PHASE_SEPARATION;
    cout << "nd_ind= " << node->index << "  nd_col= " << num_col << "  nd_val= " << node->value << "  nd_dep= "
         << node->tree_level << "  et= " << glo_eps << "  lb= " << lb << "  ub= "
         << ub << "  nd_rmn= " << bbt.size() << endl;
	if(!node->brcs.empty()){
	  for (auto &brc : node->brcs) {
		cout << (brc.br_dir ? "true" : "false") << "(" << brc.edge.first << "," << brc.edge.second << ")" << " ";
	  }
	  cout << endl;
	}

    double old_val;
    if (node->index) {
      old_val = node->value;
    }
#if defined(TRAINING_DATA_TREE_LEVEL)
#ifdef TRAINING_DATA_TREE_LEVEL
    if (node->tree_level > TRAINING_DATA_TREE_LEVEL) {
      delete node;
      continue;
    }
#endif
#endif


#ifdef USE_M_DYNAMICS
    auto beg_c = std::chrono::high_resolution_clock::now();
#endif

	soft_time_reached= false;
	if(rollback_solver.model)rollback_solver.freeModel();
	last_max_time_labeling = 0;
    force_not_rollback = true;
    solveLPInLabeling(node, true, true, true);
    force_not_rollback = false;

    if (node->index) {
      auto &brc = node->brcs.back();
      if (brc.br_dir) {
        real_improvement_up[brc.edge].first += node->value - old_val;
        ++real_improvement_up[brc.edge].second;
      } else {
        real_improvement_down[brc.edge].first += node->value - old_val;
        ++real_improvement_down[brc.edge].second;
      }
	  if (rollback == ROLL_BACK_TAIL_OFF) {
		if (node->is_terminated) {
		  delete node;
		  continue;
		}
		soft_time_reached = true;
#if SETTING!=2
		eliminateArcs(node);
		enumerateMIP(node);
#endif
		if (!node) continue;
	  }
    }

    if (node->is_terminated) {
      delete node;
      continue;
    }

    if (node->index==0) {
	  updateLowerBound(node->value);
	  final_decision_4_arc_elimination = true;
	  eliminateArcs(node);
      if (!node) continue;
      cleanIndexColForNode(node, true);
    }

	if(!soft_time_reached || counter_try_hard_pricing<Config::CounterTryHardPricingButFailed){
	  separateHybridCuts(node);
	  if (!node) continue;
	  else if (node->is_terminated) {
		delete node;
		continue;
	  }
	}

    cout << BIG_PHASE_SEPARATION;

    if (num_col > LP_COL_FINAL_LIMIT) {
      cleanIndexColForNode(node, false);
#if VERBOSE_MODE==1
      cout << BIG_PHASE_SEPARATION;
      cout << "Run column reduction in memory... ncol= " << num_col << " are left!" << endl;
      cout << BIG_PHASE_SEPARATION;
#endif
    }

#ifdef USE_M_DYNAMICS
    if (node->index) {
      auto end_c = std::chrono::high_resolution_clock::now();
      auto eps = duration<double>(end_c - beg_c).count();
      BbNode::updateState(eps, node->c,(int) node->brcs.size() - 1);
      double new_r;
      node->calculateRStar(node->value - old_val, new_r, this);
      BbNode::updateState(new_r, node->geo_r_star, (int)node->brcs.size() - 1);
#if VERBOSE_MODE==1
	  cout << "node->c= " << node->c << " node->geo_r_star= " << node->geo_r_star << endl;
#endif
    }
#endif
    recordOptimalColumn(node);

    pair<int, int> info;
	doSB(node, info);
    addBranchCutToUnsolved(node, info);
    ++num_br;
  }

  glo_end=std::chrono::high_resolution_clock::now();
  glo_eps=duration<double>(glo_end-glo_beg).count();
  if(glo_eps<=GLOBAL_TIME_LIMIT) {
	lb_transformed = ub;
	cout << "All remaining nodes pruned! Optimality has been proven!  lb=ub= " << ub << endl;
	cout << BIG_PHASE_SEPARATION;
  }

  DELETE:
  while (!bbt.empty()) {
    BbNode *extra_node = bbt.top();
    bbt.pop();
    ++num_explored_nodes;
    delete extra_node;
  }
  solver.freeEnv();

  global_gap = (ub - lb_transformed) / ub;

  printOptIntSol();

  return false;
}

bool CVRP::addBranchCutToUnsolved(
    BbNode *const node,
    const pair<int, int> &info
) {
  if (node->is_terminated) return false;

  int numnz;
  int c = num_row;//old_num_cut
  Brc bf;
  bf.edge = info;
  bf.idx_br_c = c;
  ++node->tree_level;
  vector<int> solver_ind;
  vector<double> solver_val;
  getNewConstraintCoefficientByEdge(node, info, solver_ind, solver_val);

  bf.br_dir = true;
  auto node2 = new BbNode(node, idx_node + 2, bf, num_buckets_per_vertex);
  safe_solver(addBranchConstraint(solver_ind, solver_val, SOLVER_GREATER_EQUAL, 1, nullptr, node2->solver))
  safe_solver(addBranchConstraint(solver_ind, solver_val, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))

  bf.br_dir = false;
  node->brcs.emplace_back(bf);

  deleteArcByFalseBranchConstraint(node, info);
  node->index = ++idx_node;
  ++idx_node;
  ++num_row;

  bbt.push(node2);//solve true first
  bbt.push(node);

  return false;
}

