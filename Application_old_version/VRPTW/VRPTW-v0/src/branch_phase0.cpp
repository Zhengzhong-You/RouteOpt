
#include <utility>

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::createBranchingSet(BbNode *node) {
  getEdgeInfo(node, true);
  vector<tuple<int, int, double>> fracEdges;
  for (int i = 1; i <= node->num_edges; ++i) {
    if (node->edge_value[i] < 1 - TOLERANCE) {
      fracEdges.emplace_back(node->edge_tail[i], node->edge_head[i], abs(node->edge_value[i] - 0.5));
    }
  }
  std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
    return get<2>(a) < get<2>(b);
  });
  branch_pair.resize(int(fracEdges.size()));
  transform(fracEdges.begin(),
            fracEdges.end(),
            branch_pair.begin(),
            [](const auto &a) {
              return make_pair(get<0>(a), get<1>(a));
            });
  cout << "createBranchingSet branch_pair size: " << branch_pair.size() << endl;
}

void CVRP::initialScreen(BbNode *node,
                            bool if_record_source,
                            bool if_fill_candidates_gap,
                            int num,
                            double pseudo_frac) {
  getEdgeInfo(node, true);
  vector<tuple<int, int, double>> fracEdges;
  vector<tuple<int, int, double>> OldBranch;
  int num_all_frac_edge = 0;
  for (int i = 1; i <= node->num_edges; ++i) {
    if (node->edge_value[i] < 1 - TOLERANCE) {
      ++num_all_frac_edge;
      int ai = node->edge_tail[i];
      int aj = node->edge_head[i];
      double frac_down = node->edge_value[i];
      double frac_up = 1 - frac_down;
      pair<int, int> edge = {ai, aj};
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_up_ptr;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_down_ptr;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator up_iter;
      std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator down_iter;
      int loop = 1;
      HERE:
      switch (loop) {
        case 1: {
          improvement_up_ptr = &real_improvement_up;
          improvement_down_ptr = &real_improvement_down;
          break;
        }
        case 2: {
          improvement_up_ptr = &heuristic_improvement_up;
          improvement_down_ptr = &heuristic_improvement_down;
          break;
        }
        case 3: {
          improvement_up_ptr = &lp_testing_improvement_up;
          improvement_down_ptr = &lp_testing_improvement_down;
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
      fracEdges.emplace_back(ai, aj, abs(node->edge_value[i] - 0.5));
    }
  }
  std::sort(OldBranch.begin(), OldBranch.end(), [](const auto &a, const auto &b) {
    return get<2>(a) > get<2>(b);
  });
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
  branch_pair.resize(all_branch_phase1);
  transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, branch_pair.begin(), [](const auto &a) {
    return make_pair(get<0>(a), get<1>(a));
  });
  transform(fracEdges.begin(), fracEdges.begin() + frac_cap, branch_pair.begin() + sudo_cap, [](const auto &a) {
    return make_pair(get<0>(a), get<1>(a));
  });
  if (if_record_source) {
    branch_pair_from_pseudo.resize(sudo_cap);
    transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, branch_pair_from_pseudo.begin(), [](const auto &a) {
      return make_pair(get<0>(a), get<1>(a));
    });
    branch_pair_from_fractional.resize(frac_cap);
    transform(fracEdges.begin(), fracEdges.begin() + frac_cap, branch_pair_from_fractional.begin(), [](const auto &a) {
      return make_pair(get<0>(a), get<1>(a));
    });
#if VERBOSE_MODE==1
    cout << "pseudo: " << sudo_cap << " | "<< "frac: " << frac_cap << endl;
#endif
  }
#if VERBOSE_MODE==1
  cout << "initialScreen branch_pair size: " << branch_pair.size() << endl;
#endif
}

void CVRP::writeMapEdgeColIndexInEnum(BbNode *node) {
  map_edge_col_idx_in_enu.clear();
  map_edge_col_idx_in_enu.reserve(dim * dim);
  pair<int, int> edge;
  for (int i = 1; i < num_col; ++i) {// keep the first one
    int past_node = 0;
    for (auto j = node->index_columns[i] + 1;; ++j) {
      int curr_node = col_pool4_pricing[j];
      if (past_node > curr_node) edge = {curr_node, past_node};
      else edge = {past_node, curr_node};
      map_edge_col_idx_in_enu[edge].emplace_back(i);
      if (!curr_node) break;
      past_node = curr_node;
    }
  }
}

void CVRP::randomSelection(BbNode *node) {
  getEdgeInfo(node, true);
  std::default_random_engine generator;
  generator.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(0, node->num_edges - 1);

  int number = distribution(generator), num_all_frac_edge = 0;
  for (int i = 1; i <= node->num_edges; ++i) {
    if (node->edge_value[i] < 1 - TOLERANCE) {
      if (num_all_frac_edge == number) {
        branch_pair = {make_pair(node->edge_tail[i], node->edge_head[i])};
        break;
      }
      ++num_all_frac_edge;
    }
  }
  cout << "random selection! " << endl;
}

void CVRP::initialPickLP(BbNode *node, int num) {
  initialScreen(node, true, true, num, Config::Frac4sudoCostBranPhase0);
  testLP(node, 1, true);
}

double CVRP::calculateDifference(double tmp_val, double prior_val, bool if_chg) const {
  if (if_chg) {
	throw std::logic_error("calculateDifference: if_chg is true");
	if (ceilTransformedNumberRelated(tmp_val) >= ub) tmp_val = ub;
  }
  return max(tmp_val - prior_val, TOLERANCE);
}