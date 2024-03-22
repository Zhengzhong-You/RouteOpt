
#include "CVRP.hpp"

using namespace std;
using namespace Eigen;

void CVRP::rmColByBranchInEnuMatrix(BbNode *node,
									vector<bool> &deleted_columns_in_enumeration_pool,
									bool if_test_not_use,
									const vector<Brc> &brcs) const {
  auto &mat = node->matrix_in_enumeration.front();
  int size_pool = node->size_enumeration_col_pool;
  sparseRowMatrixXd tmp(1, size_pool);

  vector<int> must_use(brcs.size()), cannot_use(brcs.size());
  int cnt1 = 0, cnt2 = 0;
  for (int c = 0; c < brcs.size(); ++c) {
	if (brcs[c].br_dir) {
	  must_use[cnt1++] = c;
	} else cannot_use[cnt2++] = c;
  }
  must_use.resize(cnt1);
  cannot_use.resize(cnt2);

  for (auto i : must_use) {//must use: 1. use and only use one edge .2 use two but not next to each other
	auto &brc = brcs[i];
	int ai = brc.edge.first, aj = brc.edge.second;
	if (ai) {
	  tmp = mat.row(ai - 1) + mat.row(aj - 1);
	  for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
		if (it.value() > 0.5) {
		  if (deleted_columns_in_enumeration_pool[it.col()]) continue;
		  if (it.value() < 1.5) {//==1
			deleted_columns_in_enumeration_pool[it.col()] = true;
		  } else {//==2 further test
			for (auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;; ++j) {
			  int current_node = col_pool4_pricing[j];
			  if (!current_node) break;
			  if (current_node == ai) {
				if (col_pool4_pricing[j + 1] != aj && col_pool4_pricing[j - 1] != aj)
				  deleted_columns_in_enumeration_pool[it.col()] = true;
			  }
			}
		  }
		}
	  }
	} else {//now ai ==0
	  for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, aj - 1); it; ++it) {
		if (it.value() > 0.5 && !deleted_columns_in_enumeration_pool[it.col()]) {//further test
		  auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;
		  if (col_pool4_pricing[j] == aj) {
			continue;
		  }
		  for (;; ++j) {
			int current_node = col_pool4_pricing[j];
			if (!current_node) break;
		  }
		  if (col_pool4_pricing[j - 1] != aj) {
			deleted_columns_in_enumeration_pool[it.col()] = true;
		  }
		}
	  }
	}
  }

  if (if_test_not_use) {
	for (auto i : cannot_use) {//cannot use: use two and next to each other
	  auto &brc = brcs[i];
	  int ai = brc.edge.first, aj = brc.edge.second;
	  if (ai) {
		tmp = mat.row(ai - 1) + mat.row(aj - 1);
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
		  if (it.value() > 1.5) {
			if (deleted_columns_in_enumeration_pool[it.col()]) continue;
			for (auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;; ++j) {
			  int current_node = col_pool4_pricing[j];
			  if (!current_node) break;
			  if (current_node == ai) {
				if (col_pool4_pricing[j + 1] == aj || col_pool4_pricing[j - 1] == aj)
				  deleted_columns_in_enumeration_pool[it.col()] = true;
			  }
			}
		  }
		}
	  } else {//now ai ==0
		for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, aj - 1); it; ++it) {
		  if (it.value() > 0.5 && !deleted_columns_in_enumeration_pool[it.col()]) {//further test
			auto j = node->index_columns_in_enumeration_column_pool[it.col()] + 1;
			if (col_pool4_pricing[j] == aj) {
			  deleted_columns_in_enumeration_pool[it.col()] = true;
			  continue;
			}
			for (;; ++j) {
			  int current_node = col_pool4_pricing[j];
			  if (!current_node) break;
			}
			if (col_pool4_pricing[j - 1] == aj) {
			  deleted_columns_in_enumeration_pool[it.col()] = true;
			}
		  }
		}
	  }
	}
  }
}

void CVRP::buildRCCInEnuMatrix(BbNode *node,
							   vector<Eigen::Triplet<double>> &triplets, int old_num) const {
  if (node->rccs.empty()) return;
  vector<Eigen::Triplet<double>> another_triplets;
  another_triplets.reserve(node->size_enumeration_col_pool * node->rccs.size());
  unordered_map<int, int> new_idx_map;
  new_idx_map.reserve(node->rccs.size());
  int cnt = 0;
  for (int c = 0; c < node->rccs.size(); ++c) {
	auto &rcc = node->rccs[c];
	if (rcc.idx_rcc < old_num)continue;
	new_idx_map[cnt] = c;
	for (auto i : rcc.info_rcc_customer) {
	  another_triplets.emplace_back(cnt, i - 1, 1);
	}
	++cnt;
  }
  sparseRowMatrixXd tmp(cnt, real_dim);
  tmp.setFromTriplets(another_triplets.begin(), another_triplets.end());
  tmp = tmp * node->basic_matrix;
  for (int c = 0; c < cnt; ++c) {
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, c); it; ++it) {
	  triplets.emplace_back(node->rccs[new_idx_map[c]].idx_rcc - old_num, it.col(), 1);
	}
  }
}

template<typename R1CType>
void processRow(const R1CType &r1c,
				const sparseRowMatrixXd &tmp,
				int row,
				vector<Eigen::Triplet<double>> &triplets,
				int old_num) {
  for (sparseRowMatrixXd::InnerIterator it(tmp, row); it; ++it) {
	int val_ = int(it.value() + TOLERANCE);
	if (val_) {
	  triplets.emplace_back(r1c.idx_r1c - old_num, it.col(), val_);
	}
  }
}


void CVRP::buildAllR1CInEnuMatrix(BbNode *node,
								  vector<Eigen::Triplet<double>> &triplets, int old_num) {
  if (node->r1cs.empty()) return;
  auto &mat0 = node->matrix_in_enumeration.empty() ? row_basic_matrix : node->matrix_in_enumeration.front();
  sparseRowMatrixXd tmp(1, node->size_enumeration_col_pool);
  for (auto &r1c : node->r1cs) {
	if (r1c.idx_r1c < old_num)continue;
	tmp.setZero();
	auto &info = r1c.info_r1c;
	const auto &plan = map_rank1_multiplier[(int)info.first.size()][info.second];
	const auto &multi = get<0>(plan);
	double denominator = get<1>(plan);
	int count = 0;
	for (auto i : info.first) {
	  tmp += mat0.row(i - 1) * multi[count++];
	}
	tmp /= denominator;
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
	  int val_ = int(it.value() + TOLERANCE);
	  if (val_)triplets.emplace_back(r1c.idx_r1c - old_num, it.col(), val_);
	}
  }

  if (node->matrix_in_enumeration.empty()) row_basic_matrix.resize(0, 0);
}

void CVRP::createBasicMatrix(BbNode *node) const {
  node->basic_matrix = (node->matrix_in_enumeration.front()).block(0, 0, real_dim, node->size_enumeration_col_pool);
}

