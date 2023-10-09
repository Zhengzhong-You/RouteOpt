//
// Created by Zhengzhong You on 3/27/23.
//

#include <utility>

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::constructBranchingSet(BBNODE *node) {
  getInfoEdge(node, true);
  vector<tuple<int, int, double>> fracEdges;
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeVal[i] < 1 - TOLERANCE) {
      fracEdges.emplace_back(node->EdgeTail[i], node->EdgeHead[i], abs(node->EdgeVal[i] - 0.5));
    }
  }
  std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
    return get<2>(a) < get<2>(b);
  });
  Branch_pair.resize(int(fracEdges.size()));
  transform(fracEdges.begin(),
            fracEdges.end(),
            Branch_pair.begin(),
            [](const auto &a) {
              return make_pair(get<0>(a), get<1>(a));
            });
  cout << "constructBranchingSet Branch_pair size: " << Branch_pair.size() << endl;
}

void CVRP::InitialScreening(BBNODE *node,
                            bool if_record_source,
                            bool if_fill_candidates_gap,
                            int num,
                            double pseudo_frac) {
  getInfoEdge(node, true);
  vector<tuple<int, int, double>> fracEdges;
  vector<tuple<int, int, double>> OldBranch;
  int num_all_frac_edge = 0;
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeVal[i] < 1 - TOLERANCE) {
      ++num_all_frac_edge;
      int ai = node->EdgeTail[i];
      int aj = node->EdgeHead[i];
      double frac_down = node->EdgeVal[i];
      double frac_up = 1 - frac_down;
      pair<int, int> edge = {ai, aj};
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_up_ptr;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_down_ptr;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator up_iter;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator down_iter;
      int loop = 1;
      //test if in real br
      HERE:
      switch (loop) {
        case 1: {
          improvement_up_ptr = &RealImprovement_up;
          improvement_down_ptr = &RealImprovement_down;
          break;
        }
        case 2: {
          improvement_up_ptr = &HeuristicImprovement_up;
          improvement_down_ptr = &HeuristicImprovement_down;
          break;
        }
        case 3: {
          improvement_up_ptr = &LPTestingImprovement_up;
          improvement_down_ptr = &LPTestingImprovement_down;
          break;
        }
        default: {
          goto END;
        }
      }
      up_iter = improvement_up_ptr->find(edge);
      down_iter = improvement_down_ptr->find(edge);
      if (up_iter != improvement_up_ptr->end() && down_iter != improvement_down_ptr->end()) {
        double up = up_iter->second.first / up_iter->second.second;
        double down = down_iter->second.first / down_iter->second.second;
        OldBranch.emplace_back(ai, aj, up * frac_up * down * frac_down);
        continue;
      } else {
        ++loop;
        goto HERE;
      }
      END:
      fracEdges.emplace_back(ai, aj, abs(node->EdgeVal[i] - 0.5));
    }
  }
  //sort
  std::sort(OldBranch.begin(), OldBranch.end(), [](const auto &a, const auto &b) {
    return get<2>(a) > get<2>(b);
  });
  //reverse
  std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
    return get<2>(a) < get<2>(b);
  });
  int all_branch_phase1 = min(num, num_all_frac_edge);
  int sudo_cap = min(int(OldBranch.size()), int(all_branch_phase1 * pseudo_frac));
  int frac_cap = min(int(fracEdges.size()), all_branch_phase1 - sudo_cap);
  if (if_fill_candidates_gap) {
    if (frac_cap < all_branch_phase1 - sudo_cap) {
      sudo_cap = min(int(OldBranch.size()), all_branch_phase1 - frac_cap);
    }
  }
  all_branch_phase1 = sudo_cap + frac_cap;
  Branch_pair.resize(all_branch_phase1);
  transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, Branch_pair.begin(), [](const auto &a) {
    return make_pair(get<0>(a), get<1>(a));
  });
  transform(fracEdges.begin(), fracEdges.begin() + frac_cap, Branch_pair.begin() + sudo_cap, [](const auto &a) {
    return make_pair(get<0>(a), get<1>(a));
  });
  if (if_record_source) {
    Branch_pair_from_pseudo.resize(sudo_cap);
    transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, Branch_pair_from_pseudo.begin(), [](const auto &a) {
      return make_pair(get<0>(a), get<1>(a));
    });
    Branch_pair_from_fractional.resize(frac_cap);
    transform(fracEdges.begin(), fracEdges.begin() + frac_cap, Branch_pair_from_fractional.begin(), [](const auto &a) {
      return make_pair(get<0>(a), get<1>(a));
    });
    cout << "pseudo: " << sudo_cap << endl;
    cout << "frac: " << frac_cap << endl;
  }
  cout << "InitialScreening Branch_pair size: " << Branch_pair.size() << endl;
}

void CVRP::writeMap_Edge_ColIdx_in_Enu(BBNODE *node) {
  Map_Edge_ColIdx_in_Enu.clear();
  Map_Edge_ColIdx_in_Enu.reserve(Dim * Dim);
  pair<int, int> edge;
  for (int i = 1; i < NumCol; ++i) {// keep the first one
    int past_node = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int curr_node = ColPool4Pricing[j];
      if (past_node > curr_node) edge = {curr_node, past_node};
      else edge = {past_node, curr_node};
      Map_Edge_ColIdx_in_Enu[edge].emplace_back(i);
      if (!curr_node) break;
      past_node = curr_node;
    }
  }
}

double CVRP::calculateDif(double tmp_val, double prior_val, bool if_chg) const {
  if (if_chg) {
    throw std::logic_error("calculateDif: if_chg is true");
    if (ceil_transformed_number_related(tmp_val) >= UB) tmp_val = UB;
  }
  return max(tmp_val - prior_val, TOLERANCE);
}

void CVRP::randomSelection(BBNODE *node) {
  getInfoEdge(node, true);
  std::default_random_engine generator;
  generator.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(0, node->NumEdges - 1);

  int number = distribution(generator), num_all_frac_edge = 0;
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeVal[i] < 1 - TOLERANCE) {
      if (num_all_frac_edge == number) {
        Branch_pair = {make_pair(node->EdgeTail[i], node->EdgeHead[i])};
        break;
      }
      ++num_all_frac_edge;
    }
  }
  cout << "random selection! " << endl;
}

/**
 * this is testing code1
 * testing random pick lp
 */
void CVRP::randomPickLP(BBNODE *node, int num) {
  getInfoEdge(node, true);
  vector<int> fracEdges;
  fracEdges.resize(node->NumEdges);
  int cnt = 0;
  for (int i = 1; i <= node->NumEdges; ++i) {
    if (node->EdgeVal[i] < 1 - TOLERANCE) {
      fracEdges[cnt++] = i;
    }
  }
  fracEdges.resize(cnt);
  std::default_random_engine generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  num = min(num, int(fracEdges.size()));
  Branch_pair.resize(num);
  for (int i = 0; i < num; ++i) {
    std::uniform_int_distribution<int> distribution(i, (int) fracEdges.size() - 1);
    int index = distribution(generator);
    std::swap(fracEdges[i], fracEdges[index]);
    Branch_pair[i] = {node->EdgeTail[fracEdges[i]], node->EdgeHead[fracEdges[i]]};
  }
  cout << "randomPickLP!" << endl;
  LPTesting(node, 1, true);
}

/**
 * initial screening pick lp
 */

void CVRP::initialPickLP(BBNODE *node, int num) {
  InitialScreening(node, true, true, num, CONFIG::Frac4sudoCostBranPhase0);
  LPTesting(node, 1, true);
}
