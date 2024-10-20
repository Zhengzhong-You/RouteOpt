
/**
 * SB_MemSave:
 * 1. solve the nodes first, solve the node with larger heuristic value, and then put them into the stack
 * TODO:
 * all the other things will not be implemented here!
 * this mode is only for the use of RCF!
 */


#ifdef BRANCH_FASHION_MEM_SAVING

#include "CVRP.hpp"

using namespace std;
using namespace chrono;


void CVRP::solveModel() {

  cout << "\n<solveModel>\n\n";

  lateProcessing();
  glo_beg = std::chrono::high_resolution_clock::now();

  buildModel();

#ifdef WRITE_ENUMERATION_TREES
  enumeration_col_idx.reserve(size_t(2e7));
#endif

  BbNode *node = bbt.top();
  bbt.pop();
  tackleUnsolvedNode(node);
  if (node) {
	bbt.push(node);
	supp_set.emplace(node->value);
	supp_queue.push(node->value);
  }

  SB_MemSave(bbt);

  lb_transformed = ub;
  cout << "All remaining nodes pruned! Optimality has been proven!  lb=ub= " << ub << endl;
  cout << BIG_PHASE_SEPARATION;

  global_gap = (ub - lb_transformed) / ub;
  Solver.freeEnv();

  printOptIntSol();

#ifdef WRITE_ENUMERATION_TREES
  writeEnuCols();
#endif
}

void CVRP::terminateNodeMemSave(BbNode *&root_node) {
  if_in_enu_state = true;
  int root_index = root_node->index;

  cout << "To terminate the current node, " <<
	   "chg node selection strategy by exploring the children nodes of the node first!\n";

  --num_explored_nodes;
  for (auto &r1c : root_node->r1cs) {
	r1c.mem.clear();
  }
  for (auto &r1c : root_node->r1cs_multi) {
	r1c.mem.clear();
  }

  tackleUnsolvedNode(root_node);
  if (root_node) {
	sub_bbt.push(root_node);
	supp_set.emplace(root_node->value);
	supp_queue.push(root_node->value);
  }

  SB_MemSave(sub_bbt);

  root_node = nullptr;;
  if_in_enu_state = false;
}

bool CVRP::checkTimeFail() {
  glo_end = high_resolution_clock::now();
  glo_eps = duration<double>(glo_end - glo_beg).count();
  if (glo_eps > GLOBAL_TIME_LIMIT) {
	cout << SMALL_PHASE_SEPARATION;
	cout << "Ins= " << file_name << "  cap= " << cap
		 << "  shut down due to reaching lm. gt.= " << GLOBAL_TIME_LIMIT << endl;
	cout << "lb= " << lb << "  ub= " << ub << "  gap(ub-lb/ub)= "
		 << (ub - lb) / (ub) * 100 << "%  gt= "
		 << glo_eps << " nd= "
		 << num_explored_nodes << "  br= " << num_br << endl;
	return true;
  }
  return false;
}

void CVRP::printInfo(BbNode *node) {
  cout << BIG_PHASE_SEPARATION;
  cout << "nd_ind= " << node->index << "  nd_col= " << num_col << "  nd_val= " << node->value << "  nd_dep= "
	   << node->tree_level << "  et= " << glo_eps << "  lb= " << lb << "  ub= "
	   << ub << "  nd_rmn= " << bbt.size() << "  sub_rmn= " << sub_bbt.size() << endl;
  int count = 0;
  for (auto &brc : node->brcs) {
	cout << (brc.br_dir ? "true" : "false") << "(" << brc.edge.first << "," << brc.edge.second << ")";
	count++;
	if (count % 8 == 0) cout << endl; else cout << " ";
  }
  cout << endl;
}

void CVRP::resetEnvironment(BbNode *node) {
  safe_solver(node->solver.updateModel())
  safe_solver(node->solver.getNumRow(&num_row))
  safe_solver(node->solver.getNumCol(&num_col))
  if (!if_in_enu_state) {
	convertVertexToR1CsInOneLP(node);
  }
}

void CVRP::SB_MemSave(StackTree &tree) {
  /**
   * take one -> SB -> cutting -> push into the stack
   */
  while (!tree.empty()) {

	if (checkTimeFail()) goto DELETE;

	BbNode *node = tree.top();
	tree.pop();

	supp_set.erase(node->value);
	if (abs(lb - supp_queue.top()) > TOLERANCE)updateLowerBound(supp_queue.top());
	if (supp_set.find(supp_queue.top()) == supp_set.end()) supp_queue.pop();

	resetEnvironment(node);

	printInfo(node);

	recordOptimalColumn(node, true);//the old data can be the other node, so force rewrite!
	if (if_in_enu_state) {
	  getEdgeInfo(node, true);
	  writeMapEdgeColIndexInEnum(node);
	}
	pair<int, int> info;
	doSB(node, info);
	BbNode *node2{};
	addBrCut(node, node2, info);
	++num_br;

	tackleUnsolvedNode(node);
	tackleUnsolvedNode(node2);

	if (node && node2) {
	  if (node->value < node2->value) swap(node, node2);
	}

	if (node) {
	  supp_set.emplace(node->value);
	  supp_queue.push(node->value);
	  tree.push(node);
	}
	if (node2) {
	  supp_set.emplace(node2->value);
	  supp_queue.push(node2->value);
	  tree.push(node2);
	}
  }
  DELETE:
  while (!tree.empty()) {
	BbNode *extra_node = tree.top();
	tree.pop();
	delete extra_node;
  }
}

void CVRP::addBrCut(
	BbNode *node,
	BbNode *&node2,
	const pair<int, int> &info
) {
  int numnz;
  Brc bf;
  bf.edge = info;
  bf.br_dir = true;
  ++node->tree_level;
  BidirectionalLinkedList *lnode{}, *rnode{};
  vector<int> solver_ind(num_col);
  vector<double> solver_val(num_col);
  if (!if_in_enu_state) {
	node->ptr->giveBirth(lnode, rnode);
	getNewConstraintCoefficientByEdge(node, info, solver_ind.data(), solver_val.data(), numnz);
	bf.idx_br_c = num_row;
	node2 = new BbNode(node, rnode, num_col, idx_node + 2, bf, num_buckets_per_vertex);
	safe_solver(addBranchConstraint(numnz, solver_ind.data(), solver_val.data(), SOLVER_GREATER_EQUAL, 1, nullptr, node2->solver))
	safe_solver(addBranchConstraint(numnz, solver_ind.data(), solver_val.data(), SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
  } else {
	bf.idx_br_c = -1;
	cstr_index.resize(num_row);
	iota(cstr_index.begin(), cstr_index.end(), 0);
	node2 = new BbNode(node, num_col, idx_node + 2, bf);
	reviseEnumColInfoByBrC(node, node2, bf);
  }

  bf.br_dir = false;
  node->brcs.emplace_back(bf);

  if (!if_in_enu_state) {
	node->num_parent_cols = num_col;
	deleteArcByFalseBranchConstraint(node, info);
	node->ptr = lnode;
	++num_row;
	safe_Hyperparameter(checkCSTLimit())
  } else {
	cstr_index.resize(num_row);
	iota(cstr_index.begin(), cstr_index.end(), 0);
	safe_solver(node->solver.getNumCol(&num_col))//cannot be deleted, recover the env is very important!
	reviseEnumColInfoByBrC(node, node, bf);
  }
  node->index = ++idx_node;
  ++idx_node;
}

void CVRP::tackleUnsolvedNode(BbNode *&node) {
  ++num_explored_nodes;
  int numParentCols;

  resetEnvironment(node);

  if (!if_in_enu_state) {
	cout << "num_buckets_per_vertex= " << num_buckets_per_vertex << " step_size= " << step_size << endl;
  }

  double old_val;
  if (node->index) old_val = node->value;

  bool if_sep_cuts = true;
  if (!if_in_enu_state) {
	pool_beg4_pricing = 0;
	last_max_time_labeling = 0;
	force_not_rollback = true;
	cout << "node->tree_level= " << node->tree_level << endl;
	solveLPInLabeling(node, true, true, true);
	force_not_rollback = false;
	if (rollback == 3) {
	  if_sep_cuts = false;
	  cout << "labels reach the soft limit!" << endl;
	}
  } else {
	max_num_enu_col_pool = node->size_enumeration_col_pool;
	solveLPByInspection(node, false, false, true);
  }

  if (node->index) {
	auto &brc = node->brcs.back();
	if (brc.br_dir) {
	  real_improvement_up[brc.edge].first += node->value - old_val;
	  ++real_improvement_up[brc.edge].second;
	} else {
	  real_improvement_down[brc.edge].first += node->value - old_val;
	  ++real_improvement_down[brc.edge].second;
	}
  }

  if (node->is_terminated) {
	delete node;
	node = nullptr;
	goto QUIT;
  } else if (if_in_enu_state && num_col + node->size_enumeration_col_pool <= max_num_route4_mip) {
	terminateByMIP(node);
	if (node->is_terminated) {
	  delete node;
	  node = nullptr;
	  goto QUIT;
	}
  }

  if (!if_in_enu_state && (!node->index || !if_sep_cuts)) {
	final_decision_4_arc_elimination = true;
	eliminateArcs(node);
	final_decision_4_enumeration = true;
	enumerateMIP(node);
	if (!node) goto QUIT;
	cleanIndexColForNode(node, node->num_parent_cols, true);
  }

  if (if_sep_cuts) {
	separateHybridCuts(node);
	if (!node) goto QUIT;
	else if (node->is_terminated) {
	  delete node;
	  node = nullptr;
	  goto QUIT;
	} else if (if_in_enu_state && num_col + node->size_enumeration_col_pool <= max_num_route4_mip) {
	  terminateByMIP(node);
	  if (node->is_terminated) {
		delete node;
		node = nullptr;
		goto QUIT;
	  }
	}
  } else {
	findNonActiveCuts(node);
	printCutsInformation(node);
  }

  if (!if_in_enu_state) {
	/**
	  * once mem is written, cannot delete any column before the branching selection and cg!
	  */
	numParentCols = node->num_parent_cols;

	writeColumnToMemory(node);

	if (num_col > LP_COL_FINAL_LIMIT) {
	  cleanIndexColForNode(node, dim);
	  cout << BIG_PHASE_SEPARATION;
	  cout << "Run column reduction in memory... ncol= " << num_col << " are left!" << endl;
	  cout << BIG_PHASE_SEPARATION;
	  node->ptr->becomeParent(node, this);
	} else {
	  constructMap(node, numParentCols);
	}
  } else {
#ifdef DELUXING_APPLIED
	cout << "Trying to apply RCF..." << endl;
	throw runtime_error("RCF is not seen in this version!");
#endif
#ifdef WRITE_ENUMERATION_TREES
	writeEnuTree(node);
	delete node;
	node = nullptr;
	goto QUIT;
#endif
	regenerateEnumMat(node, nullptr, true);
  }
  QUIT:
  return;
}

#endif
