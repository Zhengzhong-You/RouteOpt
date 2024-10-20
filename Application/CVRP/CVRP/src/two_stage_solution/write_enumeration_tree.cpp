

#include "write_enumeration_tree.hpp"
#include "branching.hpp"

/**
 * important: not every enumerated column should be written in the .tree file
 * since some column may have a larger cost than the recorded one
 * such column should be abandoned!
 */

using namespace std;

CVRP *WriteEnumerationTree::cvrp{};
BbNode *WriteEnumerationTree::node{};
std::unordered_map<yzzLong, std::tuple<std::vector<int>, double, size_t>> WriteEnumerationTree::enumeration_col_idx{};

void WriteEnumerationTree::writeEnuTree() {
  cout << "Must use after the matrix has been built!" << endl;
  /**
   * step 1: match the map
   */
  cvrp->regenerateEnumMat(node, nullptr, true);

  vector<pair<size_t, double> > col_infos(node->size_enumeration_col_pool);
  int cnt = 0;
  auto &ptr_col = node->index_columns_in_enumeration_column_pool;
  for (int i = 0; i < node->size_enumeration_col_pool; ++i) {
	col_infos[cnt++] = {ptr_col[i], node->cost_for_columns_in_enumeration_column_pool[i]};
  }
  vector<pair<size_t, double>> col_idx_cost;
  populateEnuCols(col_idx_cost, node->getCols(), col_infos);
  /**
   * step 2: write the .tree file
   * include: nd_ind, nd_col, nd_val, nd_dep, et, lb, ub, ColPool
   * and the col index, cuts index!!!
   */
  writeEnuCuts(col_idx_cost);
}

void WriteEnumerationTree::writeEnuCuts(
	const std::vector<std::pair<size_t, double>> &col_idx_cost
) {
  self_mkdir(TREE_FOLDER);
  string local_file_name = TREE_FOLDER + "/" + cvrp->file_name + "_" + to_string(node->index) + ".tree";
  ofstream file(local_file_name);
  if (!file.is_open()) throw runtime_error("Cannot open file:" + local_file_name);

  file << "nd_ind= " << node->index << "  nd_col= " << cvrp->getNumCol() << "  nd_cstr= " << cvrp->getNumRow() << "  nd_val= "
	   << node->getCurrentNodeVal()
	   << "  nd_dep= " << node->getTreeLevel() << "  lb= " << BaseBranching::lb << "  ub= "
	   << BaseBranching::ub << "  ColPool= " << node->size_enumeration_col_pool <<
	   endl;

  file << "ColPool: ";
  for (auto &i : col_idx_cost) {
	file << i.first << " " << i.second << ",";
  }
  file << endl;

  file << "rccs | form_rcc | rhs | idx_rcc | info_rcc_customer | info_rcc_outside_customer" <<
	   endl;
  for (auto &rcc : node->rccs) {
	file << rcc.form_rcc << " | " << rcc.rhs << " | " << rcc.idx_rcc << " | ";
	for (auto &i : rcc.info_rcc_customer) {
	  file << i << " ";
	}
	file << " | ";
	for (auto &i : rcc.info_rcc_outside_customer) {
	  file << i << " ";
	}
	file << endl;
  }

  file << "r1cs | rhs | idx_r1c | info_r1c | plan" << endl;
  for (auto &r1c : node->r1cs) {
	file << r1c.rhs << " | " << r1c.idx_r1c << " | ";
	for (auto &i : r1c.info_r1c.first) {
	  file << i << " ";
	}
	file << " | " << r1c.info_r1c.second;
	file << endl;
  }

  file << "Brc | ai | aj | Sense" <<
	   endl;
  for (auto &br : node->getBrCs()) {
	file << br.edge.first << " | " << br.edge.second << " | " << br.br_dir << endl;
  }
}

void WriteEnumerationTree::populateEnuCols(
	std::vector<std::pair<size_t, double>> &col_idx_cost,
	const std::vector<SequenceInfo> &lp_col_info,
	const std::vector<std::pair<size_t, double>> &enu_col_idx_cost) {
  col_idx_cost.resize(lp_col_info.size() + enu_col_idx_cost.size());
  int cnt = 0;
  auto col_pool4_pricing = cvrp->col_pool4_pricing;
  vector<double> lp_cost(cvrp->getNumCol());
  safe_solver(node->getSolver().getObj(0, cvrp->getNumCol(), lp_cost.data()));
  for (int i = 0; i < lp_col_info.size(); ++i) {
	auto &col = lp_col_info[i];
	yzzLong tmp{};
	for (auto &j : col.col_seq) {
	  tmp.set(j);
	}
	if (enumeration_col_idx.find(tmp) == enumeration_col_idx.end()) {
	  enumeration_col_idx[tmp] = make_tuple(col.col_seq, lp_cost[i], enumeration_col_idx.size());
	  col_idx_cost[cnt++] = {get<2>(enumeration_col_idx[tmp]), lp_cost[i]};
	} else {
	  if (abs(get<1>(enumeration_col_idx[tmp]) - lp_cost[i]) < TOLERANCE) {
		col_idx_cost[cnt++] = {get<2>(enumeration_col_idx[tmp]), lp_cost[i]};
	  } else if (get<1>(enumeration_col_idx[tmp]) > lp_cost[i]) {
		get<0>(enumeration_col_idx[tmp]) = col.col_seq;
		get<1>(enumeration_col_idx[tmp]) = lp_cost[i];
		col_idx_cost[cnt++] = {get<2>(enumeration_col_idx[tmp]), lp_cost[i]};
	  }//else: less than the recorded cost, abandon
	}
  }
  for (auto &i : enu_col_idx_cost) {
	yzzLong tmp{};
	for (auto j = i.first + 1;; ++j) {
	  int curr_node = col_pool4_pricing[j];
	  if (!curr_node) break;
	  tmp.set(curr_node);
	}
	if (enumeration_col_idx.find(tmp) == enumeration_col_idx.end()) {
	  vector<int> seq(tmp.count());
	  int idx = 0;
	  for (auto j = i.first + 1;; ++j) {
		int curr_node = col_pool4_pricing[j];
		if (!curr_node) break;
		seq[idx++] = curr_node;
	  }
	  enumeration_col_idx[tmp] = make_tuple(seq, i.second, enumeration_col_idx.size());
	  col_idx_cost[cnt++] = {get<2>(enumeration_col_idx[tmp]), i.second};
	} else {
	  if (abs(get<1>(enumeration_col_idx[tmp]) - i.second) < TOLERANCE) {
		col_idx_cost[cnt++] = {get<2>(enumeration_col_idx[tmp]), i.second};
	  } else if (get<1>(enumeration_col_idx[tmp]) > i.second) {
		vector<int> seq(tmp.count());
		int idx = 0;
		for (auto j = i.first + 1;; ++j) {
		  int curr_node = col_pool4_pricing[j];
		  if (!curr_node) break;
		  seq[idx++] = curr_node;
		}
		get<0>(enumeration_col_idx[tmp]) = seq;
		get<1>(enumeration_col_idx[tmp]) = i.second;
		col_idx_cost[cnt++] = {get<2>(enumeration_col_idx[tmp]), i.second};
	  }
	}
  }
  if (cnt < col_idx_cost.size()) {
	cout << "original size: " << col_idx_cost.size() << endl;
	cout << "reduced size: " << cnt << endl;
	cout << "eliminate: " << col_idx_cost.size() - cnt << endl;
  }
  col_idx_cost.resize(cnt);
}

void WriteEnumerationTree::writeEnuCols() {
  self_mkdir(COL_POOL_FOLDER);
  string local_file_name = COL_POOL_FOLDER + "/" + cvrp->file_name + ".pool";
  ofstream file(local_file_name);
  if (!file.is_open()) throw runtime_error("Cannot open file: " + local_file_name);
  vector<pair<yzzLong, size_t>> tmp(enumeration_col_idx.size());
  size_t cnt = 0;
  for (auto &i : enumeration_col_idx) {
	tmp[cnt++] = {i.first, get<2>(i.second)};
  }
  sort(tmp.begin(), tmp.end(), [](const pair<yzzLong, size_t> &a, const pair<yzzLong, size_t> &b) {
	return a.second < b.second;
  });
  for (auto &i : tmp) {
	auto &col_info = enumeration_col_idx[i.first];
	auto &col = get<0>(col_info);
	auto cost = get<1>(col_info);
	file << cost << " | ";
	for (auto &j : col) {
	  file << j << " ";
	}
	file << endl;
  }
}

void WriteEnumerationTree::init(CVRP *pr_cvrp) {
  cvrp = pr_cvrp;
  enumeration_col_idx.clear();
  enumeration_col_idx.reserve(size_t(2e7));
}

void WriteEnumerationTree::updateNode(BbNode *pr_node) {
  node = pr_node;
}