
#include "CVRP.hpp"

#ifdef WRITE_ENUMERATION_TREES

using namespace std;

void CVRP::writeEnuTree(BbNode *node) {
  cout << "Must use after the matrix has been built!" << endl;
  /**
   * step 1: match the map
   */
  if (node->valid_size != node->size_enumeration_col_pool) {
	regenerateEnumMat(node, nullptr, true);
  }

  int all_columns = num_col + node->size_enumeration_col_pool;
  vector<pair<size_t, double> > col_infos(all_columns);
  int cnt = 0;
  auto &ptr = node->index_columns;
  vector<double> solver_obj(num_col);
  safe_solver(node->solver.getObj(0, num_col, solver_obj.data()));
  for (int i = 0; i < num_col; ++i) {
	col_infos[cnt++] = {ptr[i], solver_obj[i]};
  }
  auto &ptr_col = node->index_columns_in_enumeration_column_pool;
  for (int i = 0; i < node->size_enumeration_col_pool; ++i) {
	col_infos[cnt++] = {ptr_col[i], node->cost_for_columns_in_enumeration_column_pool[i]};
  }
  vector<size_t> col_idx;
  populateEnuCols(col_idx, col_infos);

  /**
   * step 2: write the .tree file
   * include: nd_ind, nd_col, nd_val, nd_dep, et, lb, ub, ColPool
   * and the col index, cuts index!!!
   */

  writeEnuCuts(node, col_idx, col_infos);
}

void CVRP::writeEnuCuts(BbNode *node,
						const std::vector<size_t> &col_idx,
						const std::vector<std::pair<size_t, double>> &col_infos) {
  self_mkdir(TREE_FOLDER);
  string local_file_name = TREE_FOLDER + "/" + file_name + "_" + to_string(node->index) + ".tree";
  ofstream file(local_file_name);
  if (!file.is_open()) {
	cerr << "Cannot open file: " << local_file_name << endl;
	exit(1);
  }

  file << "nd_ind= " << node->index << "  nd_col= " << num_col << "  nd_cstr= " << num_row << "  nd_val= "
	   << node->value
	   << "  nd_dep= " << node->tree_level << "  lb= " << lb << "  ub= "
	   << ub << "  ColPool= " << node->size_enumeration_col_pool << endl;

  file << "ColPool: ";
  for (int i = 0; i < col_idx.size(); ++i) {
	file << col_idx[i] << " " << col_infos[i].second << ",";
  }
  file << endl;

  file << "rccs | form_rcc | rhs | idx_rcc | info_rcc_customer | info_rcc_outside_customer" << endl;
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

  file << "r1cs | rhs | idx_r1c | info_r1c" << endl;
  for (auto &r1c : node->r1cs) {
	file << r1c.rhs << " | " << r1c.idx_r1c << " | ";
	for (auto &i : r1c.info_r1c) {
	  file << i << " ";
	}
	file << endl;
  }

  file << "R1cMulti | rhs | idx_r1c | info_r1c" << endl;
  for (auto &r1c : node->r1cs_multi) {
	file << r1c.rhs << " | " << r1c.idx_r1c << " | ";
	for (auto &i : r1c.info_r1c.first) {
	  file << i << " ";
	}
	file << " | " << r1c.info_r1c.second << endl;
  }

  file << "Brc | ai | aj | Sense" << endl;
  for (auto &br : node->brcs) {
	file << br.edge.first << " | " << br.edge.second << " | " << br.br_dir << endl;
  }
}

void CVRP::populateEnuCols(std::vector<size_t> &col_idx, const std::vector<std::pair<size_t, double>> &col_infos) {
  col_idx.resize(col_infos.size());
  int cnt = 0;
  for (auto &i : col_infos) {
	yzzLong tmp;
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
	} else {
	  if (get<1>(enumeration_col_idx[tmp]) > i.second) {
		vector<int> seq(tmp.count());
		int idx = 0;
		for (auto j = i.first + 1;; ++j) {
		  int curr_node = col_pool4_pricing[j];
		  if (!curr_node) break;
		  seq[idx++] = curr_node;
		}
		get<0>(enumeration_col_idx[tmp]) = seq;
		get<1>(enumeration_col_idx[tmp]) = i.second;
	  }
	}
	col_idx[cnt++] = get<2>(enumeration_col_idx[tmp]);
  }
}

void CVRP::writeEnuCols() {
  self_mkdir(COL_POOL_FOLDER);
  string local_file_name = COL_POOL_FOLDER + "/" + file_name + ".pool";
  ofstream file(local_file_name);
  if (!file.is_open()) {
	cerr << "Cannot open file: " << local_file_name << endl;
	exit(1);
  }
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

#endif