//
// Created by You, Zhengzhong on 8/15/23.
//


#include "CVRP.hpp"

#ifdef writeEnumerationTrees

using namespace std;

void CVRP::writeEnuTree(BBNODE *node) {
  cout << "Must use after the matrix has been built!" << endl;
  /**
   * step 1: match the map
   */
  if (node->validSize != node->SizeEnuColPool) {
    regenerateEnuMat(node, nullptr, true);
  }

  int all_columns = NumCol + node->SizeEnuColPool;
  vector<pair<size_t, double> > col_infos(all_columns);
  int cnt = 0;
  auto &ptr = node->IdxCols;
  safe_solver(node->solver.SOLVERgetObj(0, NumCol, solver_obj))
  for (int i = 0; i < NumCol; ++i) {
    col_infos[cnt++] = {ptr[i], solver_obj[i]};
  }
  auto &ptr_col = node->IdxColsInEnuColPool;
  for (int i = 0; i < node->SizeEnuColPool; ++i) {
    col_infos[cnt++] = {ptr_col[i], node->Cost4ColsInEnuColPool[i]};
  }
  vector<size_t> col_idx;
  populateEnuCols(col_idx, col_infos);

  /**
   * step 2: write the .tree file
   * include: nd_ind, nd_col, nd_val, nd_dep, et, lb, ub, ColPool
   * and the col index, cuts index!!!
   */

  writeEnuCuts(node, col_idx);
}

void CVRP::writeEnuCuts(BBNODE *node,
                        const std::vector<size_t> &col_idx) {
  self_mkdir(tree_folder);
  string file_name = tree_folder + "/" + FileName + "_" + to_string(node->Idx) + ".tree";
  ofstream file(file_name);
  if (!file.is_open()) {
    cerr << "Cannot open file: " << file_name << endl;
    exit(1);
  }

  file << "nd_ind= " << node->Idx << "  nd_col= " << NumCol << "  nd_cstr= " << NumRow << "  nd_val= " << node->Val
       << "  nd_dep= " << node->TreeLevel << "  lb= " << LB << "  ub= "
       << UB << "  ColPool= " << node->SizeEnuColPool << endl;

  file << "ColPool: ";
  for (auto &i : col_idx) {
    file << i << " ";
  }
  file << endl;

  file << "RCCs | FormRCC | RHS | IdxRCC | InfoRCCCustomer | InfoRCCOutsideCustomer" << endl;
  for (auto &rcc : node->RCCs) {
    file << rcc.FormRCC << " | " << rcc.RHS << " | " << rcc.IdxRCC << " | ";
    for (auto &i : rcc.InfoRCCCustomer) {
      file << i << " ";
    }
    file << " | ";
    for (auto &i : rcc.InfoRCCOutsideCustomer) {
      file << i << " ";
    }
    file << endl;
  }

  file << "R1Cs | RHS | IdxR1C | InfoR1C" << endl;
  for (auto &r1c : node->R1Cs) {
    file << r1c.RHS << " | " << r1c.IdxR1C << " | ";
    for (auto &i : r1c.InfoR1C) {
      file << i << " ";
    }
    file << endl;
  }

  file << "R1C_multi | RHS | IdxR1C | InfoR1C" << endl;
  for (auto &r1c : node->R1Cs_multi) {
    file << r1c.RHS << " | " << r1c.IdxR1C << " | ";
    for (auto &i : r1c.InfoR1C.first) {
      file << i << " ";
    }
    file << " | " << r1c.InfoR1C.second << endl;
  }

  file << "BrC | ai | aj | Sense" << endl;
  for (auto &br : node->BrCs) {
    file << br.Edge.first << " | " << br.Edge.second << " | " << br.BrDir << endl;
  }
}

void CVRP::populateEnuCols(std::vector<size_t> &col_idx, const std::vector<std::pair<size_t, double>> &col_infos) {
  col_idx.resize(col_infos.size());
  int cnt = 0;
  for (auto &i : col_infos) {
    yzzLong tmp;
    for (auto j = i.first + 1;; ++j) {
      int curr_node = ColPool4Pricing[j];
      if (!curr_node) break;
      tmp.set(curr_node);
    }
    if (enumeration_col_idx.find(tmp) == enumeration_col_idx.end()) {
      //add new col
      vector<int> seq(tmp.count());
      int idx = 0;
      for (auto j = i.first + 1;; ++j) {
        int curr_node = ColPool4Pricing[j];
        if (!curr_node) break;
        seq[idx++] = curr_node;
      }
      enumeration_col_idx[tmp] = make_tuple(seq, i.second, enumeration_col_idx.size());
    } else {
      //update the cost
      if (get<1>(enumeration_col_idx[tmp]) > i.second) {
        vector<int> seq(tmp.count());
        int idx = 0;
        for (auto j = i.first + 1;; ++j) {
          int curr_node = ColPool4Pricing[j];
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
  self_mkdir(col_pool_folder);
  string file_name = col_pool_folder + "/" + FileName + ".pool";
  ofstream file(file_name);
  if (!file.is_open()) {
    cerr << "Cannot open file: " << file_name << endl;
    exit(1);
  }
  //sort by the index
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