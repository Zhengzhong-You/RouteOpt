
#include "CVRP.hpp"

using namespace std;

template<typename T>
void CVRP::extendKernel4LightHeur(Label *&ki,
                                  int i,
                                  double res,
                                  const std::vector<T> &arc,
                                  const double *r1c_to_pi,
                                  const double *r1c_multi_to_pi,
                                  bool if_return_to_depot) {
  for (auto &pair : arc) {
    int j;
    if constexpr (std::is_same<T, int>::value) {
      j = pair;
    } else {
      j = pair.second;
      res = pair.first;
    }
    if (ki->pi[j]) continue;
    auto &tmp_mainResource = all_label[idx_glo].sum_main_resource;

    if (!increaseMainResourceConsumption(res, tmp_mainResource, i, j)) continue;

    int bj = int(tmp_mainResource / step_size);
    auto &labelList_j = label_array_in_forward_sense[j][bj].first;
    auto &valid_num = label_array_in_forward_sense[j][bj].second;
    auto &tmp_rc = all_label[idx_glo].rc;
    auto &tmp_PI = all_label[idx_glo].pi;
    auto &tmp_Rank1CutMem = all_label[idx_glo].rank1_cut_mem;
    tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
    tmp_Rank1CutMem = ki->rank1_cut_mem;
    auto &tmp_Rank1CutMem_multi = all_label[idx_glo].rank1_cut_mem_multi;
    copy(ki->rank1_cut_mem_multi, ki->rank1_cut_mem_multi + num_valid_r1c_multi_in_cg, tmp_Rank1CutMem_multi);

    for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
      if (tmp_Rank1CutMem[l]) {
        tmp_Rank1CutMem[l] = false;
        tmp_rc -= r1c_to_pi[l];
      } else tmp_Rank1CutMem[l] = true;
    }

    for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
      tmp_Rank1CutMem_multi[get<0>(l)] += get<1>(l);
      if (tmp_Rank1CutMem_multi[get<0>(l)] >= get<2>(l)) {
        tmp_rc -= r1c_multi_to_pi[get<0>(l)];
        tmp_Rank1CutMem_multi[get<0>(l)] -= get<2>(l);
      }
    }

    if (valid_num) {
      auto kj = labelList_j[0];
      if (tmp_rc < kj->rc) {
        kj->is_extended = true;
      } else {
        continue;
      }
    }

    tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);

    for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;

    tmp_PI = (ki->pi) & (ng_mem4_vertex[j]);
    tmp_PI.set(j);

    labelList_j[0] = all_label + idx_glo;

    valid_num = 1;

    all_label[idx_glo].p_label = ki;
    all_label[idx_glo].end_vertex = j;
    all_label[idx_glo].is_extended = false;
    auto &bucket = if_exist_extra_labels_in_forward_sense[j][bj];
    bucket.first[0] = all_label + idx_glo;
    bucket.second = 1;
    double path_rc = tmp_rc + chg_cost_mat4_vertex[j][0];
    if (if_return_to_depot) {
      addPathByRC(path_rc, all_label + idx_glo, nullptr, Config::MaxNumRoutesInLighterHeur);
    }
    ++idx_glo;
  }
}

int CVRP::runLighterHeuristicLabeling(BbNode *const node) {

  int index = num_col - 1;

  PtrAllR1CS ptrAllR1Cs(node, this);

  prior_pool_beg4_pricing = pool_beg4_pricing;

  initializeLabels(node, 1, true, {true, 1, true});

  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  Label *ki, *b4_ki;
  bool if_break;

  unordered_set<int> depot_set;
  depot_set.reserve(dim);
#ifdef SYMMETRY_PROHIBIT
  for (auto j : node->all_backward_buckets[0][0].bucket_arcs) depot_set.emplace(j);
#else
  for (auto j : node->all_forward_buckets[0][0].bucket_arcs) depot_set.emplace(j);
#endif

  for (int b = 0; b < num_buckets_per_vertex; ++b) {
    int i = 1;
    STILL_EXIST:
    for (; i < dim; ++i) {
      auto &un_extend_num = if_exist_extra_labels_in_forward_sense[i][b].second;
      if (!un_extend_num) continue;
      auto &label_array = if_exist_extra_labels_in_forward_sense[i][b].first;
      ki = label_array[0];
      if (ki->is_extended) {
        un_extend_num = 0;
        continue;
      }
      if_break = false;
      for (int b4_b = b - 1; b4_b >= 0; --b4_b) {
        if (rc2_till_this_bin_in_forward_sense[i][b4_b] > ki->rc) break;
        auto &b4_valid_num = label_array_in_forward_sense[i][b4_b].second;
        if (!b4_valid_num) continue;
        auto &b4_label_list = label_array_in_forward_sense[i][b4_b].first;
        b4_ki = b4_label_list[0];
        if (ki->rc > b4_ki->rc) {
          if_break = true;
          goto b4_break;
        }
      }
      b4_break:
      ki->is_extended = true;
      if (if_break) continue;
      bool if_return_to_depot = depot_set.find(i) != depot_set.end();
      extendKernel4LightHeur(ki,
                             i,
                             ki->sum_main_resource,
                             node->all_forward_buckets[i][b].bucket_arcs,
                             r1c_to_pi,
                             r1c_multi_to_pi, if_return_to_depot);
      extendKernel4LightHeur(ki,
                             i,
                             0,
                             node->all_forward_buckets[i][b].jump_arcs,
                             r1c_to_pi,
                             r1c_multi_to_pi,
                             if_return_to_depot);
      un_extend_num = 0;
    }
    for (i = 1; i < dim; ++i) {
      if (if_exist_extra_labels_in_forward_sense[i][b].second)
        goto STILL_EXIST;
    }
    for (i = 1; i < dim; ++i) {
      std::stable_sort(label_array_in_forward_sense[i][b].first.begin(),
                       label_array_in_forward_sense[i][b].first.begin() + label_array_in_forward_sense[i][b].second,
                       CmpLabelRCLess);
    }
    if (b) {
      for (i = 1; i < dim; ++i) {
        if (label_array_in_forward_sense[i][b].second) {
          rc2_till_this_bin_in_forward_sense[i][b] = min(rc2_till_this_bin_in_forward_sense[i][b - 1],
                                                   label_array_in_forward_sense[i][b].first[0]->rc);
        } else {
          rc2_till_this_bin_in_forward_sense[i][b] = rc2_till_this_bin_in_forward_sense[i][b - 1];
        }
      }
    } else {
      for (i = 1; i < dim; ++i) {
        if (label_array_in_forward_sense[i][b].second) {
          rc2_till_this_bin_in_forward_sense[i][b] = label_array_in_forward_sense[i][b].first[0]->rc;
        } else {
          rc2_till_this_bin_in_forward_sense[i][b] = LARGE_FLOAT;
        }
      }
    }
  }

  writeColumnsInPricingPool(node, index);

  return index - num_col + 1;
}

int CVRP::generateColumnsByLighterHeuristic(BbNode *const node) {

  int ccnt = runLighterHeuristicLabeling(node);

  if (!ccnt) return 0;

  addColumns(node, ccnt);

  return ccnt;
}

