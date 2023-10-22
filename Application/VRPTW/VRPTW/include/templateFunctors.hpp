
#include "CVRP.hpp"

template<typename T, bool dir, bool if_last_half, bool if_symmetry>
int CVRP::extendKernel4Exact(Label *&ki,
                             int i,
                             double res,
                             const std::vector<T> &arc,
                             const double *r1c_to_pi,
                             const double *r1c_multi_to_pi,
                             int min_sorted_b) {
  int state = 0;
  int bj;
  bool if_suc;
  /**
   *
   */
#ifdef test_time
  double test_time1 = 0, test_time2 = 0;
#endif
  for (auto &pair : arc) {
    int j;
    if constexpr (std::is_same<T, int>::value) {
      j = pair;
    } else {
      j = pair.second;
      res = pair.first;
    }
#ifdef test_time
    auto beg = std::chrono::high_resolution_clock::now();
#endif
    updateLabel<dir, if_last_half, if_symmetry, false, false>(res, ki, i, j, bj, r1c_to_pi, r1c_multi_to_pi,
                                                              min_sorted_b,
                                                              if_suc);
#ifdef test_time
    auto end = std::chrono::high_resolution_clock::now();
    if constexpr (if_last_half)
      test_time1 += std::chrono::duration<double>(end - beg).count();
#endif
    if (!if_suc) continue;

#ifdef test_time
    beg = std::chrono::high_resolution_clock::now();
#endif
    doDominance<dir>(ki, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
#ifdef test_time
    end = std::chrono::high_resolution_clock::now();
    if constexpr (if_last_half)
      test_time2 += std::chrono::duration<double>(end - beg).count();
#endif
    if (!if_suc) continue;

    if constexpr (!if_last_half) {
      addPathByRC(all_label[idx_glo].rc + chg_cost_mat4_vertex[j][0],
                  all_label + idx_glo,
                  nullptr,
                  Config::MaxNumRoutesInExact);
    }

    ++idx_glo;
    if (idx_glo == label_assign) {
      rollback = 2;
      state = 2;//QUIT
      goto QUIT;
    }
  }
  QUIT:
#ifdef test_time
  if constexpr (if_last_half) {
    test_time_map["updateLabel"] += test_time1;
    test_time_map["doDominance"] += test_time2;
  }
#endif
  return state;
}

template<bool dir, bool if_last_half, bool if_symmetry, bool if_std_optgap, bool if_res_updated>
void CVRP::updateLabel(double res, Label *&ki, int i, int j, int &bj,
                       const double *r1c_to_pi,
                       const double *r1c_multi_to_pi,
                       int min_sorted_b,
                       bool &if_suc) {
  if_suc = false;
  if (ki->pi[j]) return;
  double &which_rc = (if_std_optgap ? opt_gap : rc_std);
  auto &tmp_mainResource = all_label[idx_glo].sum_main_resource;
  auto &tmp_rc = all_label[idx_glo].rc;
  tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];//real rc
  if constexpr (!if_res_updated) {
    if constexpr (dir) {
      if (!increaseMainResourceConsumption(res, tmp_mainResource, i, j)) return;
      if constexpr (!if_last_half) {
        if (tmp_mainResource > meet_point_resource_in_bi_dir) {
          concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_mainResource);
          return;
        }
      }
    } else {
      if (!decreaseMainResourceConsumption(res, tmp_mainResource, i, j)) return;
      if constexpr (!if_last_half) {
        if (tmp_mainResource < meet_point_resource_in_bi_dir) {
          concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_mainResource);
          return;
        }
      }
    }
  }
  if constexpr (if_last_half) {
    /**
     * here we only test the RC2TillThisBinInXXXwardSense, we do not get into details of the labels in the bins
     */
    int concate_bj;
    int if_state;
    if constexpr (if_symmetry) {
      concate_bj = int((max_main_resource - tmp_mainResource) / step_size);
    } else {
      concate_bj = int(tmp_mainResource / step_size);
    }
    constexpr bool if_dif = dir ^ if_symmetry;
    concatenateOneLabelWithOtherLabels<dir, if_symmetry, true, true>(ki,
                                                                     j,
                                                                     concate_bj,
                                                                     tmp_rc,
                                                                     tmp_mainResource,
                                                                     r1c_to_pi,
                                                                     r1c_multi_to_pi,
                                                                     if_state);
    if (if_state == -2) return;
    else if (if_state >= 0) goto outside;
    for (if_dif ? ++concate_bj : --concate_bj; if_dif ? concate_bj < num_buckets_per_vertex : concate_bj >= 0;
         if_dif ? ++concate_bj : --concate_bj) {
      concatenateOneLabelWithOtherLabels<dir, if_symmetry, false, true>(ki,
                                                                        j,
                                                                        concate_bj,
                                                                        tmp_rc,
                                                                        tmp_mainResource,
                                                                        r1c_to_pi,
                                                                        r1c_multi_to_pi,
                                                                        if_state);
      if (if_state == -2) break;
      else if (if_state >= 0) goto outside;
    }
    outside:
    if (if_state < 0) {
      return;
    }
  }

  if_suc = true;
  bj = int(tmp_mainResource / step_size);
  auto &tmp_PI = all_label[idx_glo].pi;
  auto &tmp_Rank1CutMem = all_label[idx_glo].rank1_cut_mem;
  auto &tmp_num_valid_rank1_cut = all_label[idx_glo].num_valid_rank1_cut;
  auto &tmp_valid_rank1_cut = all_label[idx_glo].valid_rank1_cut;
  auto &tmp_num_valid_rank1_cut_multi = all_label[idx_glo].num_valid_rank1_cut_multi;
  auto &tmp_valid_rank1_cut_multi = all_label[idx_glo].valid_rank1_cut_multi;

  tmp_PI = (ki->pi) & (ng_mem4_vertex[j]);
  tmp_PI.set(j);
  tmp_Rank1CutMem = ki->rank1_cut_mem;

  auto &tmp_Rank1CutMem_multi = all_label[idx_glo].rank1_cut_mem_multi;
  std::copy(ki->rank1_cut_mem_multi, ki->rank1_cut_mem_multi + num_valid_r1c_multi_in_cg, tmp_Rank1CutMem_multi);

  for (auto l : std::get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
    if (tmp_Rank1CutMem[l]) {
      tmp_Rank1CutMem[l] = false;
      tmp_rc -= r1c_to_pi[l];
    } else tmp_Rank1CutMem[l] = true;
  }
  tmp_Rank1CutMem &= std::get<1>(Vertex2ActiveInOnePricingR1Cs[j]);
  tmp_num_valid_rank1_cut = 0;
  for (auto l : std::get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
    if (tmp_Rank1CutMem[l]) {
      tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
    }
  }

  for (auto &l : std::get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
    int tmp_cut = std::get<0>(l);
    tmp_Rank1CutMem_multi[tmp_cut] += std::get<1>(l);
    if (tmp_Rank1CutMem_multi[tmp_cut] >= std::get<2>(l)) {
      tmp_rc -= r1c_multi_to_pi[tmp_cut];
      tmp_Rank1CutMem_multi[tmp_cut] -= std::get<2>(l);
    }
  }

  for (auto l : std::get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;

  tmp_num_valid_rank1_cut_multi = 0;
  for (auto l : std::get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
    if (tmp_Rank1CutMem_multi[l]) {
      tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
    }
  }
}

template<bool dir>
void CVRP::doDominance(Label *&ki, int j, int bj, const double *r1c_to_pi,
                       const double *r1c_multi_to_pi, bool &if_suc) {
  auto &labelList_j = dir ? label_array_in_forward_sense[j][bj].first : label_array_in_backward_sense[j][bj].first;
  auto &valid_num_j = dir ? label_array_in_forward_sense[j][bj].second : label_array_in_backward_sense[j][bj].second;
  auto &tmp_mainResource = all_label[idx_glo].sum_main_resource;
  auto &tmp_rc = all_label[idx_glo].rc;
  auto &tmp_PI = all_label[idx_glo].pi;
  auto &tmp_Rank1CutMem = all_label[idx_glo].rank1_cut_mem;
  auto &tmp_num_valid_rank1_cut = all_label[idx_glo].num_valid_rank1_cut;
  auto &tmp_valid_rank1_cut = all_label[idx_glo].valid_rank1_cut;
  auto &tmp_Rank1CutMem_multi = all_label[idx_glo].rank1_cut_mem_multi;
  auto &tmp_num_valid_rank1_cut_multi = all_label[idx_glo].num_valid_rank1_cut_multi;
  auto &tmp_valid_rank1_cut_multi = all_label[idx_glo].valid_rank1_cut_multi;

  double tmp_rc_add = tmp_rc + RC_TOLERANCE, tmp_rc_sub = tmp_rc - RC_TOLERANCE;
  if_suc = true;
  double dif;
  for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
    auto &kj = labelList_j[vec_index_j];
    ++num_dominance_checks;
    if (define_dir_resource<dir>(kj->sum_main_resource, tmp_mainResource)) {
      if (kj->rc < tmp_rc_sub) {
        if (((tmp_PI & kj->pi) ^ (kj->pi)).none()) {
          dif = tmp_rc_sub + gap_between_last_smallest_rc_and_rc_threshold;
          for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
            if (!tmp_Rank1CutMem[kj->valid_rank1_cut[l]]) {
              dif += r1c_to_pi[kj->valid_rank1_cut[l]];
              if (dif < kj->rc) goto here;
            }
          }

          for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
            int tmp = kj->valid_rank1_cut_multi[l];
            if (kj->rank1_cut_mem_multi[tmp]
                > tmp_Rank1CutMem_multi[tmp]) {
              dif += r1c_multi_to_pi[tmp];
              if (dif < kj->rc) goto here;
            }
          }
          if_suc = false;
          break;//taken
          here:;
        }
      }
      ++vec_index_j;
    } else if (define_dir_resource<dir>(tmp_mainResource, kj->sum_main_resource)) {
      if (tmp_rc_add < kj->rc) {
        if (((tmp_PI & kj->pi) ^ (tmp_PI)).none()) {
          dif = tmp_rc_add - gap_between_last_smallest_rc_and_rc_threshold;
          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
            if (!kj->rank1_cut_mem[tmp_valid_rank1_cut[l]]) {
              dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
              if (dif > kj->rc) goto there;
            }
          }

          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
            int tmp = tmp_valid_rank1_cut_multi[l];
            if (tmp_Rank1CutMem_multi[tmp]
                > kj->rank1_cut_mem_multi[tmp]) {
              dif -= r1c_multi_to_pi[tmp];
              if (dif > kj->rc) goto there;
            }
          }
          kj->is_extended = true;
          kj = labelList_j[--valid_num_j];
          goto _there;
          there:
          ++vec_index_j;
          _there:;
        } else ++vec_index_j;
      } else ++vec_index_j;
    } else {
      if (kj->rc < tmp_rc_add) {
        if (((tmp_PI & kj->pi) ^ (kj->pi)).none()) {
          dif = tmp_rc_sub + gap_between_last_smallest_rc_and_rc_threshold;
          for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
            if (!tmp_Rank1CutMem[kj->valid_rank1_cut[l]]) {
              dif += r1c_to_pi[kj->valid_rank1_cut[l]];
              if (dif < kj->rc) goto here1;
            }
          }

          for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
            int tmp = kj->valid_rank1_cut_multi[l];
            if (kj->rank1_cut_mem_multi[tmp]
                > tmp_Rank1CutMem_multi[tmp]) {
              dif += r1c_multi_to_pi[tmp];
              if (dif < kj->rc) goto here1;
            }
          }
          if_suc = false;
          break;//taken
          here1:;
        }
        ++vec_index_j;
      } else if (tmp_rc_sub < kj->rc) {
        if (((tmp_PI & kj->pi) ^ (tmp_PI)).none()) {
          dif = tmp_rc_add - gap_between_last_smallest_rc_and_rc_threshold;
          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
            if (!kj->rank1_cut_mem[tmp_valid_rank1_cut[l]]) {
              dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
              if (dif > kj->rc) goto there1;
            }
          }

          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
            int tmp = tmp_valid_rank1_cut_multi[l];
            if (tmp_Rank1CutMem_multi[tmp]
                > kj->rank1_cut_mem_multi[tmp]) {
              dif -= r1c_multi_to_pi[tmp];
              if (dif > kj->rc) goto there1;
            }
          }
          kj->is_extended = true;
          kj = labelList_j[--valid_num_j];
          goto _there1;
          there1:
          ++vec_index_j;
          _there1:;
        } else ++vec_index_j;
      } else {//q & rc all equal now?
        if ((tmp_PI ^ kj->pi).none()) {//all equal
          std::bitset<2> who_win = 1;//1==tmp_win,2==kj win,0==all keep
          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
            if (!kj->rank1_cut_mem[tmp_valid_rank1_cut[l]]) {
              who_win = 2;
              goto next_test1;
            }
          }
          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
            int tmp = tmp_valid_rank1_cut_multi[l];
            if (tmp_Rank1CutMem_multi[tmp]
                > kj->rank1_cut_mem_multi[tmp]) {
              who_win = 2;
              goto next_test1;
            }
          }
          next_test1:
          if (who_win == 1) {
            kj->is_extended = true;
            kj = labelList_j[--valid_num_j];
          } else {
            for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
              if (!tmp_Rank1CutMem[kj->valid_rank1_cut[l]]) {
                who_win = 0;
                goto next_test2;
              }
            }
            for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
              int tmp = kj->valid_rank1_cut_multi[l];
              if (kj->rank1_cut_mem_multi[tmp]
                  > tmp_Rank1CutMem_multi[tmp]) {
                who_win = 0;
                goto next_test2;
              }
            }
            next_test2:
            if (who_win == 2) {
              if_suc = false;
              break;//taken
            } else ++vec_index_j;
          }
        } else {
          yzzLong tmp = tmp_PI & kj->pi;
          if ((tmp ^ tmp_PI).none()) {
            dif = tmp_rc_add - gap_between_last_smallest_rc_and_rc_threshold;
            for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
              if (!kj->rank1_cut_mem[tmp_valid_rank1_cut[l]]) {
                dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
                if (dif > kj->rc) goto there2;
              }
            }

            for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
              int tmp2 = tmp_valid_rank1_cut_multi[l];
              if (tmp_Rank1CutMem_multi[tmp2]
                  > kj->rank1_cut_mem_multi[tmp2]) {
                dif -= r1c_multi_to_pi[tmp2];
                if (dif > kj->rc) goto there2;
              }
            }

            kj->is_extended = true;
            kj = labelList_j[--valid_num_j];
            goto _there2;

            there2:
            ++vec_index_j;
            _there2:;

          } else if ((tmp ^ kj->pi).none()) {
            dif = tmp_rc_sub + gap_between_last_smallest_rc_and_rc_threshold;
            for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
              if (!tmp_Rank1CutMem[kj->valid_rank1_cut[l]]) {
                dif += r1c_to_pi[kj->valid_rank1_cut[l]];
                if (dif < kj->rc) goto here2;
              }
            }

            for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
              int tmp2 = kj->valid_rank1_cut_multi[l];
              if (kj->rank1_cut_mem_multi[tmp2]
                  > tmp_Rank1CutMem_multi[tmp2]) {
                dif += r1c_multi_to_pi[tmp2];
                if (dif < kj->rc) goto here2;
              }
            }
            if_suc = false;
            break;//taken
            here2:
            ++vec_index_j;
          } else ++vec_index_j;
        }
      }
    }
  }
  if (if_suc) {
    labelList_j[valid_num_j++] = all_label + idx_glo;
    if (valid_num_j == labelList_j.size()) {
      labelList_j.resize(labelList_j.size() * 2);
    }

    all_label[idx_glo].p_label = ki;
    all_label[idx_glo].end_vertex = j;
    all_label[idx_glo].is_extended = false;
    auto &bucket = dir ? if_exist_extra_labels_in_forward_sense[j][bj] : if_exist_extra_labels_in_backward_sense[j][bj];
    bucket.first[bucket.second++] = all_label + idx_glo;
    if (bucket.second == bucket.first.size()) {
      bucket.first.resize(bucket.first.size() * 2);
    }
  }
}
template<bool dir>
void CVRP::checkIfDominated(Label *&ki, int i, int b, const double *r1c_to_pi,
                            const double *r1c_multi_to_pi,
                            bool &if_suc) {
  if_suc = true;
  double dif;
  double tmp_ki_rc_sub = ki->rc - RC_TOLERANCE;
  for (int b4_b = (dir ? b - 1 : b + 1); dir ? b4_b >= 0 : b4_b < num_buckets_per_vertex; dir ? --b4_b : ++b4_b) {
    auto &b4_label_list = dir ? label_array_in_forward_sense[i][b4_b].first : label_array_in_backward_sense[i][b4_b].first;
    auto &b4_valid_num = dir ? label_array_in_forward_sense[i][b4_b].second : label_array_in_backward_sense[i][b4_b].second;
    if ((dir ? rc2_till_this_bin_in_forward_sense[i][b4_b] : rc2_till_this_bin_in_backward_sense[i][b4_b]) > tmp_ki_rc_sub) break;
    for (int vec_b4 = 0; vec_b4 < b4_valid_num; ++vec_b4) {
      auto &b4_ki = b4_label_list[vec_b4];
      if (b4_ki->rc > tmp_ki_rc_sub) break;
      if (((ki->pi & b4_ki->pi) ^ (b4_ki->pi)).none()) {
        dif = tmp_ki_rc_sub + gap_between_last_smallest_rc_and_rc_threshold;
        for (int l = 0; l < b4_ki->num_valid_rank1_cut; ++l) {
          if (!ki->rank1_cut_mem[b4_ki->valid_rank1_cut[l]]) {
            dif += r1c_to_pi[b4_ki->valid_rank1_cut[l]];
            if (dif < b4_ki->rc) goto here;
          }
        }

        for (int l = 0; l < b4_ki->num_valid_rank1_cut_multi; ++l) {
          int tmp2 = b4_ki->valid_rank1_cut_multi[l];
          if (b4_ki->rank1_cut_mem_multi[tmp2] > ki->rank1_cut_mem_multi[tmp2]) {
            dif += r1c_multi_to_pi[tmp2];
            if (dif < b4_ki->rc) goto here;
          }
        }
        if_suc = false;
        return;
        here:;
      }
    }
  }
}
template<bool dir>
void CVRP::checkIfNoLabelsLeft(int &i, int b, int &min_sorted_b, bool &if_suc) {
  for (i = 1; i < dim; ++i) {
    if (dir ? if_exist_extra_labels_in_forward_sense[i][b].second : if_exist_extra_labels_in_backward_sense[i][b].second) {
      if_suc = false;
      return;
    }
  }
  if_suc = true;
  for (i = 1; i < dim; ++i) {
    auto &label_array = (dir ? label_array_in_forward_sense : label_array_in_backward_sense);
    auto &b4_label_list = label_array[i][b].first;
    auto &b4_valid_num = label_array[i][b].second;
    std::sort(b4_label_list.begin(), b4_label_list.begin() + b4_valid_num, CmpLabelRCLess);
    min_sorted_b = b;
  }
  for (i = 1; i < dim; ++i) {
    auto &label_bin = (dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b]);
    auto &rc2_till_this_bin = (dir ? rc2_till_this_bin_in_forward_sense[i][b] : rc2_till_this_bin_in_backward_sense[i][b]);
    if constexpr (dir) {
      if (label_bin.second) {
        rc2_till_this_bin = (b ? std::min(rc2_till_this_bin_in_forward_sense[i][b - 1], label_bin.first[0]->rc)
                               : label_bin.first[0]->rc);
      } else {
        rc2_till_this_bin = (b ? rc2_till_this_bin_in_forward_sense[i][b - 1] : LARGE_FLOAT);
      }
    } else {
      if (label_bin.second) {
        rc2_till_this_bin = (b < num_buckets_per_vertex - 1 ? std::min(rc2_till_this_bin_in_backward_sense[i][b + 1],
                                                                    label_bin.first[0]->rc)
                                                         : label_bin.first[0]->rc);
      } else {
        rc2_till_this_bin = (b < num_buckets_per_vertex - 1 ? rc2_till_this_bin_in_backward_sense[i][b + 1] : LARGE_FLOAT);
      }
    }
  }
}
template<bool dir, bool if_symmetry, bool if_check_res, bool if_std_optgap>
void CVRP::concatenateOneLabelWithOtherLabels(Label *ki, int j, int arr_bj, double tmp_rc, double tmp_mainResource,
                                              const double *r1c_to_pi,
                                              const double *r1c_multi_to_pi,
                                              int &if_state) {
  double path_rc;
  double &which_rc = if_std_optgap ? opt_gap : rc_std;
  auto ptr_rc_till_this_bin = &rc2_till_this_bin_in_forward_sense[j][arr_bj];
  auto ptr_rc_bin = &rc2_bin_in_forward_sense[j][arr_bj];
  if constexpr ((dir && !if_symmetry) || (!dir && if_symmetry)) {
    ptr_rc_till_this_bin = &rc2_till_this_bin_in_backward_sense[j][arr_bj];
    ptr_rc_bin = &rc2_bin_in_backward_sense[j][arr_bj];
  }

  if (*ptr_rc_till_this_bin + tmp_rc > which_rc) {
    if_state = -2;//no need to continue
    return;
  }
  if_state = -1;//do some operations
  if (*ptr_rc_bin + tmp_rc < which_rc) {//most_negative_rc_in_this_bin
    auto &label_arr =
        ((!dir && !if_symmetry) || (dir && if_symmetry)) ? label_array_in_forward_sense[j][arr_bj].first
                                                         : label_array_in_backward_sense[j][arr_bj].first;
    auto &label_valid_num =
        ((!dir && !if_symmetry) || (dir && if_symmetry)) ? label_array_in_forward_sense[j][arr_bj].second
                                                         : label_array_in_backward_sense[j][arr_bj].second;
    for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
      auto &kj = label_arr[vec_index];
      path_rc = kj->rc + tmp_rc;
      if (path_rc > which_rc) break;
      if constexpr (if_check_res) {
        if constexpr (if_symmetry) {
          if (tmp_mainResource + kj->sum_main_resource > max_main_resource) continue;
        } else {
          if (dir ? tmp_mainResource > kj->sum_main_resource : tmp_mainResource < kj->sum_main_resource) continue;
        }
      }
      if ((ki->pi & kj->pi).any()) continue;

      if (ki->num_valid_rank1_cut < kj->num_valid_rank1_cut) {
        for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
          if (kj->rank1_cut_mem[ki->valid_rank1_cut[l]]) {
            path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
            if (path_rc > which_rc) goto here;
          }
        }
      } else {
        for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
          if (ki->rank1_cut_mem[kj->valid_rank1_cut[l]]) {
            path_rc -= r1c_to_pi[kj->valid_rank1_cut[l]];
            if (path_rc > which_rc) goto here;
          }
        }
      }

      if (ki->num_valid_rank1_cut_multi < kj->num_valid_rank1_cut_multi) {
        for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
          int tmp_cut = ki->valid_rank1_cut_multi[l];
          if (kj->rank1_cut_mem_multi[tmp_cut] +
              ki->rank1_cut_mem_multi[tmp_cut]
              >= r1c_multi_denominator_in_cg[tmp_cut]
              ) {
            path_rc -= r1c_multi_to_pi[tmp_cut];
            if (path_rc > which_rc) goto here;
          }
        }
      } else {
        for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
          int tmp_cut = kj->valid_rank1_cut_multi[l];
          if (ki->rank1_cut_mem_multi[tmp_cut] +
              kj->rank1_cut_mem_multi[tmp_cut]
              >= r1c_multi_denominator_in_cg[tmp_cut]
              ) {
            path_rc -= r1c_multi_to_pi[tmp_cut];
            if (path_rc > which_rc) goto here;
          }
        }
      }
      if constexpr (!if_std_optgap) {
        addPathByRC(path_rc, ki, kj, Config::MaxNumRoutesInExact);
      } else {
        if (path_rc < opt_gap) {
          if_state = arr_bj;//what bin to concatenate
          break;
        }
      }
      here:;
    }
  }
}
template<bool dir, bool if_last_half, bool if_symmetry>
void CVRP::runLabeling(BbNode *node, const PtrAllR1CS &ptrAllR1Cs) {
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  bool if_suc;

  if constexpr (!if_last_half) {
    if constexpr (dir) {
      initializeLabels(node, 1, true, {true, 1, true});
    } else {
      initializeLabels(node, 2, false, {true, 2, false});
    }
  }
  auto beg = std::chrono::high_resolution_clock::now();
  auto end = beg;
#ifdef test_time
  auto new_beg = beg;
#endif
  double eps;
  int min_sorted_b = dir ? -1 : num_buckets_per_vertex;
#ifdef test_time
  double test_time1 = 0, test_time2 = 0, test_time3 = 0, test_time4 = 0;
#endif
  for (int b = (dir ? 0 : num_buckets_per_vertex - 1); (dir ? b < num_buckets_per_vertex : b >= 0); (dir ? ++b : --b)) {
    int i = 1;
    STILL_EXIST:
    for (; i < dim; ++i) {
      auto &valid_num =
          dir ? if_exist_extra_labels_in_forward_sense[i][b].second : if_exist_extra_labels_in_backward_sense[i][b].second;
      if (!valid_num) continue;
      auto &label_array =
          dir ? if_exist_extra_labels_in_forward_sense[i][b].first : if_exist_extra_labels_in_backward_sense[i][b].first;
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->is_extended) continue;
#ifdef test_time
        new_beg = std::chrono::high_resolution_clock::now();
#endif

        checkIfDominated<dir>(ki, i, b, r1c_to_pi, r1c_multi_to_pi, if_suc);
#ifdef test_time
        end = std::chrono::high_resolution_clock::now();
        test_time1 += std::chrono::duration<double>(end - new_beg).count();
#endif
        ki->is_extended = true;
        if (!if_suc)continue;
#ifdef test_time
        new_beg = std::chrono::high_resolution_clock::now();
#endif
        auto sig = extendKernel4Exact<int,
                                      dir, if_last_half, if_symmetry>(ki, i, ki->sum_main_resource,
                                                                      dir ? node->all_forward_buckets[i][b].bucket_arcs
                                                                          : node->all_backward_buckets[i][b].bucket_arcs,
                                                                      r1c_to_pi,
                                                                      r1c_multi_to_pi,
                                                                      min_sorted_b);
#ifdef test_time
        end = std::chrono::high_resolution_clock::now();
        test_time2 += std::chrono::duration<double>(end - new_beg).count();
#endif
        if (sig == 2) goto populateBin;
#ifdef test_time
        new_beg = std::chrono::high_resolution_clock::now();
#endif
        sig = extendKernel4Exact<std::pair<double, int>,
                                 dir, if_last_half, if_symmetry>(ki, i, 0,
                                                                 dir ? node->all_forward_buckets[i][b].jump_arcs
                                                                     : node->all_backward_buckets[i][b].jump_arcs,
                                                                 r1c_to_pi,
                                                                 r1c_multi_to_pi, min_sorted_b);
#ifdef test_time
        end = std::chrono::high_resolution_clock::now();
        test_time3 += std::chrono::duration<double>(end - new_beg).count();
#endif
        if (sig == 2) goto populateBin;
      }
      valid_num = 0;
    }
    end = std::chrono::high_resolution_clock::now();
    eps = std::chrono::duration<double>(end - beg).count();
    if constexpr (if_last_half) {
      if (eps > Config::HardTimeThresholdInArcEliminationLastHalf) {
        rollback = 1;
        goto QUIT;
      }
    } else {
      if (eps > cut_gen_time_threshold_in_pricing) {
        rollback = 3;
        if (eps > Config::HardTimeThresholdInPricing) {
          if (!force_not_rollback) {
            rollback = 1;
            goto QUIT;
          }
        }
      }
    }

    checkIfNoLabelsLeft<dir>(i, b, min_sorted_b, if_suc);
    if (!if_suc) goto STILL_EXIST;
  }

  populateBin:
  if constexpr (!if_last_half) {
    end = std::chrono::high_resolution_clock::now();
    eps = std::chrono::duration<double>(end - beg).count();//overall time
    last_max_time_labeling = std::max(last_max_time_labeling, eps);
    if (rollback == 0) {//check if reach the soft time
      if (last_max_time_labeling * hard_rollback_factor > Config::HardTimeThresholdInPricing) {
        rollback = 3;
      }
    }
    if (dir ? min_sorted_b < num_buckets_per_vertex : min_sorted_b >= 0) {
      for (int i = 1; i < dim; ++i) {
        for (int b = dir ? min_sorted_b + 1 : min_sorted_b - 1; dir ? b < num_buckets_per_vertex : b >= 0;
             dir ? ++b : --b) {
          auto &label_pr = dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b];
          auto &valid_num = label_pr.second;
          if (!valid_num) continue;
          auto &label_bin = label_pr.first;
          std::sort(label_bin.begin(), label_bin.begin() + valid_num, CmpLabelRCLess);
        }
      }
    }
  }
#ifdef test_time
  else {
    std::cout << "if_dominated= " << test_time1 << std::endl;
    std::cout << "bucket_extend= " << test_time2 << std::endl;
    std::cout << "jump_extend= " << test_time3 << std::endl;
    for (auto &it : test_time_map) {
      std::cout << it.first << " " << it.second << std::endl;
    }
    test_time_map.clear();
  }
#endif
  populateRC2TillThisBinNRC2Bin<dir>(node);

  QUIT:
#ifdef DETAILED_EXACT_PRINT_INFO
  cout << "这里是测试！" << endl;
  {
    if (gap_between_last_smallest_rc_and_rc_threshold < TOLERANCE) {
      cout << "we find if there exists extra labels that can be dominated but saved for some reasons!" << endl;
      safe_solver(node->solver.getDual(0, num_row, pi))
      int all_num = int(node->r1cs.size() + node->r1cs_multi.size());
      vector<int> denominator(all_num);
      vector<double> pi(all_num);
      vector<vector<std::pair<int, int>>> cuts(dim);//cut_idx, cut_augment
      vector<unordered_set<int>> mem(dim);
      vector<vector<int>> no_mem(dim);
      int cnt = 0;
      for (auto &r1c : node->r1cs) {
        for (auto &i : r1c.info_r1c) {
          cuts[i].emplace_back(cnt, 1);
          mem[i].emplace(cnt);
        }
        for (auto &i : r1c.mem) {
          mem[i].emplace(cnt);
        }
        denominator[cnt] = 2;
        pi[cnt] = pi[r1c.idx_r1c];
        ++cnt;
      }

      for (auto &r1c : node->r1cs_multi) {
        auto &plan = map_rank1_multiplier[(int) r1c.info_r1c.first.size()][r1c.info_r1c.second];
        auto &multi = std::get<0>(plan);
        int deno = std::get<1>(plan);
        int count = 0;
        for (auto &i : r1c.info_r1c.first) {
          cuts[i].emplace_back(cnt, multi[count]);
          mem[i].emplace(cnt);
          ++count;
        }
        for (auto &v : r1c.mem) {
          mem[v].emplace(cnt);
        }
        denominator[cnt] = deno;
        pi[cnt] = pi[r1c.idx_r1c];
        ++cnt;
      }

      for (int i = 1; i < dim; ++i) {
        for (int j = 0; j < all_num; ++j) {
          if (mem[i].find(j) == mem[i].end()) {
            no_mem[i].emplace_back(j);
          }
        }
      }

      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        for (int i = 1; i < dim; ++i) {
          auto &label_list = label_array_in_forward_sense[i][b].first;
          auto &valid_num = label_array_in_forward_sense[i][b].second;
          int label_count = 0;
          HERE:
          if (valid_num <= label_count) continue;
          auto ki = label_list[label_count++];
          auto ng = ki->pi;
          for (int j = label_count; j < valid_num; ++j) {
            auto kj = label_list[j];
            if (ki->sum_main_resource < kj->sum_main_resource) {
              if (((kj->pi & ng) ^ (ng)).none()) {
                double rc_dif = kj->rc - ki->rc;
                if (rc_dif > TOLERANCE) {
                  vector<int> seq_ki;
                  vector<int> seq_kj;
                  auto tmp_ki = ki;
                  while (tmp_ki->p_label) {
                    seq_ki.emplace_back(tmp_ki->end_vertex);
                    tmp_ki = tmp_ki->p_label;
                  }
                  std::reverse(seq_ki.begin(), seq_ki.end());
                  auto tmp_kj = kj;
                  while (tmp_kj->p_label) {
                    seq_kj.emplace_back(tmp_kj->end_vertex);
                    tmp_kj = tmp_kj->p_label;
                  }
                  std::reverse(seq_kj.begin(), seq_kj.end());
                  vector<int> state_ki(all_num, 0);
                  vector<int> state_kj(all_num, 0);
                  for (auto &v : seq_ki) {
                    if (v == 0) continue;
                    for (auto &cut : cuts[v]) {
                      state_ki[cut.first] += cut.second;
                      if (state_ki[cut.first] >= denominator[cut.first]) {
                        state_ki[cut.first] -= denominator[cut.first];
                      }
                    }
                    for (auto &cut : no_mem[v]) {
                      state_ki[cut] = 0;
                    }
                  }
                  for (auto &v : seq_kj) {
                    if (v == 0) continue;
                    for (auto &cut : cuts[v]) {
                      state_kj[cut.first] += cut.second;
                      if (state_kj[cut.first] >= denominator[cut.first]) {
                        state_kj[cut.first] -= denominator[cut.first];
                      }
                    }
                    for (auto &cut : no_mem[v]) {
                      state_kj[cut] = 0;
                    }
                  }

                  double diff = 0;
                  for (int k = 0; k < all_num; ++k) {
                    if (state_ki[k] > state_kj[k]) diff -= pi[k];
                  }
                  cout << "diff: " << diff << " rc_dif: " << rc_dif << endl;
                  auto min_it = min_element(pi.begin(), pi.end());
                  vector<std::pair<double, int>> pi_copy;
                  for (int k = 0; k < pi.size(); ++k) {
                    if (k >= node->r1cs.size()) {
                      int num = k - node->r1cs.size();
                      pi_copy.emplace_back(pi[k], node->r1cs_multi[num].info_r1c.first.size());
                    } else {
                      pi_copy.emplace_back(pi[k], node->r1cs[k].info_r1c.size());
                    }
                  }
                  sort(pi_copy.begin(), pi_copy.end(), [](const std::pair<double, bool> &a, const std::pair<double, bool> &b) {
                    return a.first < b.first;
                  });
                  for (int k = 0; k < pi_copy.size(); ++k) {
                    cout << "(" << pi_copy[k].first << ", " << pi_copy[k].second << ") ";
                  }
                  cout << endl;
                  cout << "min pi: " << *min_it << endl;
                  int cut = min_it - pi.begin();
                  if (cut > node->r1cs.size()) {
                    cut -= node->r1cs.size();
                    for (int i = 0; i < cut; ++i) {
                      auto &r1c = node->r1cs_multi[i];
                      auto &plan = map_rank1_multiplier[(int) r1c.info_r1c.first.size()][r1c.info_r1c.second];
                      auto &multi = std::get<0>(plan);
                      int deno = std::get<1>(plan);
                      int count = 0;
                      for (auto &i : r1c.info_r1c.first) {
                        cout << i << " ";
                      }
                      cout << endl;
                      cout << "multi: ";
                      for (auto &i : multi) {
                        cout << i << " ";
                      }
                      cout << endl;
                      cout << "deno: " << deno << endl;
                    }
                  } else {
                    auto &r1c = node->r1cs[cut];
                    for (auto &i : r1c.info_r1c) {
                      cout << i << " ";
                    }
                    cout << endl;
                  }

                  if (diff < rc_dif) {
                    cout << "diff: " << diff << " rc_dif: " << rc_dif << endl;
                    cout << "seq_ki: ";
                    for (auto &v : seq_ki) {
                      cout << v << " ";
                    }
                    cout << endl;
                    cout << "seq_kj: ";
                    for (auto &v : seq_kj) {
                      cout << v << " ";
                    }
                    cout << endl;
                    cout << "state_ki: ";
                    for (auto &v : state_ki) {
                      cout << v << " ";
                    }
                    cout << endl;
                    cout << "state_kj: ";
                    for (auto &v : state_kj) {
                      cout << v << " ";
                    }
                    cout << endl;
                    cout << "pi: ";
                    for (auto &v : pi) {
                      cout << v << " ";
                    }
                    cout << endl;
                    cout << "denominator: ";
                    for (auto &v : denominator) {
                      cout << v << " ";
                    }
                    cout << endl;
                    cout << "cuts: " << endl;
                    for (auto &v : cuts) {
                      for (auto &vv : v) {
                        cout << vv.first << " " << vv.second << " ";
                      }
                      cout << endl;
                    }
                    cout << "no_mem: " << endl;
                    for (auto &v : no_mem) {
                      for (auto &vv : v) {
                        cout << vv << " ";
                      }
                      cout << endl;
                    }
                    cout << "mem: " << endl;
                    for (auto &v : mem) {
                      for (auto &vv : v) {
                        cout << vv << " ";
                      }
                      cout << endl;
                    }
                  }
                }
              } else goto HERE;
            }

          }
        }
      }
    }
  }
  vector<int> tmp_forwad_soft;
  int times = 0;
  int hard_times = 0;
  for (int i = 1; i < dim; ++i) {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      tmp_forwad_soft.emplace_back(label_array_in_forward_sense[i][b].second);
    }
  }
  std::sort(tmp_forwad_soft.begin(), tmp_forwad_soft.end(), greater<>());
  int top_0_5 = (int) ((double) tmp_forwad_soft.size() * 0.005);
  int top_1 = (int) ((double) tmp_forwad_soft.size() * 0.01);
  int top_5 = (int) ((double) tmp_forwad_soft.size() * 0.05);
  cout << "for forward exact///max= " << tmp_forwad_soft[0] << "  0.5%=" << tmp_forwad_soft[top_0_5]
       << "  1%=" << tmp_forwad_soft[top_1] << "  5%=" << tmp_forwad_soft[top_5] << endl;

  vector<std::pair<int, double>> record;
  record.reserve(num_row);

  for (int i = 0; i < node->r1cs.size(); ++i) {
    int row_idx = node->r1cs[i].idx_r1c;
    record.emplace_back(i, abs(pi[row_idx]) * pow(node->r1cs[i].mem.size(), Config::MemFactor));
  }

  for (int i = 0; i < node->r1cs_multi.size(); ++i) {
    int row_idx = node->r1cs_multi[i].idx_r1c;
    record.emplace_back(i + node->r1cs.size(),
                        abs(pi[row_idx]) * pow(node->r1cs_multi[i].mem.size(), Config::MemFactor));
  }

  double
      sum = accumulate(record.begin(), record.end(), 0.0, [](double a, std::pair<int, double> b) { return a + b.second; });

  cout << "sum= " << sum << endl;
  end = std::chrono::high_resolution_clock::now();
  eps = std::chrono::duration<double>(end - beg).count();
  cout << "forward exact time= " << eps << endl;
#endif
  if (rollback == 1) {
    std::cout << "rollback to original states!" << std::endl;
  } else if (rollback == 2) {
    std::cout << "rollback with larger mem!" << std::endl;
  }
}

template<bool dir, bool if_symmetry>
void CVRP::concatenatePhaseInArcElimination(const double *r1c_to_pi,
                                            const double *r1c_multi_to_pi) {
  bool if_find;
  int concate_bj;
  int if_state;
  int bj;
  bool if_suc;
  auto beg = std::chrono::high_resolution_clock::now();
  auto end = beg;
  double eps;
  constexpr bool if_dif = dir ^ if_symmetry;
  for (auto &label_list : dir ? concatenate_labels_in_forward_cg : concatenate_labels_in_backward_cg) {
    int i = label_list.first.first;
    int j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      auto &tmp_mainResource = all_label[idx_glo].sum_main_resource;
      tmp_mainResource = pr.second;

      updateLabel<dir, true, if_symmetry, true, true>(tmp_mainResource, ki, i, j, bj, r1c_to_pi, r1c_multi_to_pi,
                                                      dir ? num_buckets_per_vertex - 1 : 0, if_suc);
      if (!if_suc) continue;

      doDominance<dir>(ki, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
      if (!if_suc) continue;

      ++idx_glo;
      if (idx_glo == label_assign) {
        rollback = 2;
        goto QUIT;
      }
    }
  }
  end = std::chrono::high_resolution_clock::now();
  eps = std::chrono::duration<double>(end - beg).count();
  if (eps > Config::HardTimeThresholdInArcEliminationMidConcatenate) {
    rollback = 1;
    std::cout << "concatenatePhaseInArcElimination time= " << eps << " failed!" << std::endl;
  }
  QUIT:
  return;
}

template<bool dir, bool if_symmetry>
void CVRP::eliminateBucketArcs(BbNode *node, const double *r1c_to_pi,
                              const double *r1c_multi_to_pi,
                              int dim_sq,
                              bool *stateBetween2Buckets,
                              int *latest_bucket) {
  memset(stateBetween2Buckets, 0, dim_sq * num_buckets_per_vertex * sizeof(bool));
  constexpr bool if_dif = dir ^ if_symmetry;
  if constexpr (if_dif) {
    memset(latest_bucket, -1, dim_sq * sizeof(int));
  } else {
    std::fill_n(latest_bucket, dim_sq, num_buckets_per_vertex);
  }
  double tmp_rc, path_rc, tmp_mainResource;
  int num_bucket_arcs = 0;
  for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
    for (int i = 1; i < dim; ++i) {
      concatenateTestKernelInArcElimination<std::pair<double, int>, dir, if_symmetry>(i,
                                                                                      b,
                                                                                      (dir
                                                                                       ? node->all_forward_buckets[i][b].jump_arcs
                                                                                       : node->all_backward_buckets[i][b].jump_arcs),
                                                                                      dim_sq,
                                                                                      stateBetween2Buckets,
                                                                                      latest_bucket,
                                                                                      r1c_to_pi,
                                                                                      r1c_multi_to_pi);
      concatenateTestKernelInArcElimination<int, dir, if_symmetry>(i,
                                                                   b,
                                                                   (dir
                                                                    ? node->all_forward_buckets[i][b].bucket_arcs
                                                                    : node->all_backward_buckets[i][b].bucket_arcs),
                                                                   dim_sq,
                                                                   stateBetween2Buckets,
                                                                   latest_bucket,
                                                                   r1c_to_pi,
                                                                   r1c_multi_to_pi);
    }
  }
#ifdef find_lost_arcs
  int i = 37, b = 10, j = 32, map = i * dim + j;
  int latest = -1;
  auto &label_array = label_array_in_forward_sense[i][b].first;
  auto &label_num = label_array_in_forward_sense[i][b].second;
  for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
    auto *ki = label_array[vec_index_i];
    if (!increaseMainResourceConsumption(ki->sum_main_resource, tmp_mainResource, i, j)) continue;
    tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
#ifdef SYMMETRY_PROHIBIT
    int arr_bj = max(int((tmp_mainResource) / step_size), latest + 1);
          for (int bj = num_buckets_per_vertex - 1; bj >= arr_bj; --bj) {
            if (tmp_rc + rc2_bin_in_backward_sense[j][bj] > opt_gap)continue;
            auto &label_arr = label_array_in_backward_sense[j][bj].first;
            auto &label_valid_num = label_array_in_backward_sense[j][bj].second;
#else
    int arr_bj = int((max_main_resource - tmp_mainResource) / step_size);
    for (int bj = 0; bj <= arr_bj; ++bj) {
      if (tmp_rc + rc2_bin_in_forward_sense[j][bj] > opt_gap)continue;
      auto &label_arr = label_array_in_forward_sense[j][bj].first;
      auto &label_valid_num = label_array_in_forward_sense[j][bj].second;
#endif
      for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
        auto &kj = label_arr[vec_index];
        if ((kj->pi & ki->pi).any()) continue;
#ifdef SYMMETRY_PROHIBIT
        if (tmp_mainResource > kj->sum_main_resource) continue;
#else
        if (tmp_mainResource + kj->sum_main_resource > max_main_resource) continue;
#endif
        path_rc = tmp_rc + kj->rc;
        if (path_rc > opt_gap) break;

        if (ki->num_valid_rank1_cut < kj->num_valid_rank1_cut) {
          for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
            if (kj->rank1_cut_mem[ki->valid_rank1_cut[l]])path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
          }
        } else {
          for (int l = 0; l < kj->num_valid_rank1_cut; ++l) {
            if (ki->rank1_cut_mem[kj->valid_rank1_cut[l]])path_rc -= r1c_to_pi[kj->valid_rank1_cut[l]];
          }
        }

        if (ki->num_valid_rank1_cut_multi < kj->num_valid_rank1_cut_multi) {
          for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
            int tmp_cut = ki->valid_rank1_cut_multi[l];
            if (kj->rank1_cut_mem_multi[tmp_cut] +
                ki->rank1_cut_mem_multi[tmp_cut]
                >= r1c_multi_denominator_in_cg[tmp_cut]
                )
              path_rc -= r1c_multi_to_pi[tmp_cut];
          }
        } else {
          for (int l = 0; l < kj->num_valid_rank1_cut_multi; ++l) {
            int tmp_cut = kj->valid_rank1_cut_multi[l];
            if (ki->rank1_cut_mem_multi[tmp_cut] +
                kj->rank1_cut_mem_multi[tmp_cut]
                >= r1c_multi_denominator_in_cg[tmp_cut]
                )
              path_rc -= r1c_multi_to_pi[tmp_cut];
          }
        }

        if (path_rc < opt_gap) {
          std::cout << "arc should exists!" <<  std::endl;
          std::cout << "bin= " << bj <<  std::endl;
          latest = bj;
          break;
        }
      }
    }
  }
  if (latest != -1) {
    int chg_latest = int((max_main_resource - latest * step_size) / step_size);
    std::cout << "tell= " << tell_which_bin4_arc_elimination_in_forward_sense[map + chg_latest * dim_sq] <<  std::endl;
  }
#endif

  std::vector<int> tmp_vec;
  tmp_vec.reserve(dim);
  int map1 = 0;
  for (int b = 0; b < num_buckets_per_vertex; ++b, map1 += dim_sq) {
    int map2 = map1 + dim;
    for (int i = 1; i < dim; ++i, map2 += dim) {
      tmp_vec.clear();
      for (int j : dir ? node->all_forward_buckets[i][b].bucket_arcs : node->all_backward_buckets[i][b].bucket_arcs) {
        if (stateBetween2Buckets[map2 + j]) {
          tmp_vec.emplace_back(j);
          ++num_bucket_arcs;
        }
      }
      (dir ? node->all_forward_buckets[i][b].bucket_arcs : node->all_backward_buckets[i][b].bucket_arcs) = tmp_vec;
    }
  }

  eliminatebuketArc4Depot<dir, if_symmetry>(node);

  std::cout << "Num of " << (dir ? "Forward" : "Backward") << " bucket_arcs= " << num_bucket_arcs << " prev.= "
            << double(num_bucket_arcs) / (dir ? node->num_forward_bucket_arcs : node->num_backward_bucket_arcs) * 100 << "%"
            << " max.= " << double(num_bucket_arcs) / max_num_forward_graph_arc * 100 << "%"
            << std::endl;
  (dir ? node->num_forward_bucket_arcs : node->num_backward_bucket_arcs) = num_bucket_arcs;
}


template<bool dir, bool if_symmetry>
void CVRP::eliminatebuketArc4Depot(BbNode *node) {
  bool constexpr if_dif = dir ^ if_symmetry;
  auto &allBuckets = dir ? node->all_forward_buckets : node->all_backward_buckets;
  auto &RC2TillThisBin = if_dif ? rc2_till_this_bin_in_backward_sense : rc2_till_this_bin_in_forward_sense;
  auto &RC2Bin = if_dif ? rc2_bin_in_backward_sense : rc2_bin_in_forward_sense;
  auto &labelArray = if_dif ? label_array_in_backward_sense : label_array_in_forward_sense;

  std::unordered_set<int> depot_set;
  for (auto &arc : allBuckets[0][0].bucket_arcs) depot_set.emplace(arc);

  for (int i = 1; i < dim; ++i) {
    if (depot_set.find(i) == depot_set.end()) continue;
    bool if_delete = true;

    for (int b = (if_dif ? 0 : (num_buckets_per_vertex - 1)); if_dif ? (b < num_buckets_per_vertex) : (b >= 0);
         if_dif ? ++b : --b) {
      if (RC2TillThisBin[i][b] + chg_cost_mat4_vertex[i][0] > opt_gap) break;
      if (RC2Bin[i][b] + chg_cost_mat4_vertex[i][0] > opt_gap) continue;

      auto &num = labelArray[i][b].second;
      if (num) {
        auto &labels = labelArray[i][b].first;
        for (int j = 0; j < num; ++j) {
          if (labels[j]->rc + chg_cost_mat4_vertex[i][0] > opt_gap) break;
          double tmp_res = labels[j]->sum_main_resource;
          if (if_dif ? (decreaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) : (
              increaseMainResourceConsumption(tmp_res, tmp_res, i, 0))) {
            if_delete = false;
            break;
          }
        }
      }
    }

    if (if_delete) {
      depot_set.erase(i);
    }
  }

  if (depot_set.size() != allBuckets[0][0].bucket_arcs.size()) {
    allBuckets[0][0].bucket_arcs.clear();
    allBuckets[0][0].bucket_arcs.assign(depot_set.begin(), depot_set.end());
    std::sort(allBuckets[0][0].bucket_arcs.begin(), allBuckets[0][0].bucket_arcs.end());
  }
}

template<typename T, bool dir, bool if_symmetry>
void CVRP::concatenateTestKernelInArcElimination(int i,
                                                 int b,
                                                 const std::vector<T> &arc,
                                                 int dim_sq,
                                                 bool *stateBetween2Buckets,
                                                 int *latest_bucket,
                                                 const double *r1c_to_pi,
                                                 const double *r1c_multi_to_pi) {
  double tmp_mainResource, tmp_rc;
  int if_state, arr_bj;
  for (auto &pr : arc) {
    int j;
    if constexpr (std::is_same<T, int>::value) {
      j = pr;
    } else {
      j = pr.second;
      if (dir ? !increaseMainResourceConsumption(pr.first, tmp_mainResource, i, j) :
          !decreaseMainResourceConsumption(pr.first, tmp_mainResource, i, j))
        continue;
    }
    int map = i * dim + j;
    int &latest = latest_bucket[map];
    int old_latest = latest;
    auto &label_array = dir ? label_array_in_forward_sense[i][b].first : label_array_in_backward_sense[i][b].first;
    auto &label_num = dir ? label_array_in_forward_sense[i][b].second : label_array_in_backward_sense[i][b].second;
    constexpr bool if_dif = dir ^ if_symmetry;
    for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
      auto *ki = label_array[vec_index_i];
      if constexpr (std::is_same<T, int>::value) {
        if (dir ? !increaseMainResourceConsumption(ki->sum_main_resource, tmp_mainResource, i, j) :
            !decreaseMainResourceConsumption(ki->sum_main_resource, tmp_mainResource, i, j))
          continue;
      }
      tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
      if constexpr (dir) {
        arr_bj =
            if_symmetry ? std::min(int((max_main_resource - tmp_mainResource) / step_size), latest - 1) : std::max(int(
                (tmp_mainResource) / step_size), latest + 1);
      } else {
        arr_bj =
            if_symmetry ? std::max(int((max_main_resource - tmp_mainResource) / step_size), latest + 1) : std::min(int(
                (tmp_mainResource) / step_size), latest - 1);
      }

      for (int bj = if_dif ? num_buckets_per_vertex - 1 : 0; if_dif ? bj >= arr_bj : bj <= arr_bj;
           if_dif ? --bj : ++bj) {
        concatenateOneLabelWithOtherLabels<dir, if_symmetry, true, true>(ki,
                                                                         j,
                                                                         bj,
                                                                         tmp_rc,
                                                                         tmp_mainResource,
                                                                         r1c_to_pi,
                                                                         r1c_multi_to_pi,
                                                                         if_state);
        if (if_state >= 0) {
          latest = if_state;
          goto outside;
        }
      }
      outside:
      continue;
    }
    if (latest != old_latest) {
      int bi, chg_latest;
      if constexpr (if_symmetry) {
        chg_latest = int((max_main_resource - latest * step_size) / step_size);
        if (chg_latest < 0) chg_latest = 0;
        else if (chg_latest >= num_buckets_per_vertex) chg_latest = num_buckets_per_vertex - 1;
      } else chg_latest = latest;
      bi = (dir ? tell_which_bin4_arc_elimination_in_forward_sense[map + chg_latest * dim_sq]
                : tell_which_bin4_arc_elimination_in_backward_sense[map + chg_latest * dim_sq]);
      int map2 = map + b * dim_sq;
      for (int k = b; dir ? k <= bi : k >= bi; dir ? ++k : --k, dir ? map2 += dim_sq : map2 -= dim_sq)
        stateBetween2Buckets[map2] = true;
    }
  }
}

template<bool dir>
void CVRP::populateRC2TillThisBinNRC2Bin(BbNode *const node) const {
  double rc, min_rc_per_bin;
  for (int i = 1; i < dim; ++i) {
    rc = LARGE_FLOAT;
    for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
      auto &label_list = dir ? label_array_in_forward_sense[i][b].first : label_array_in_backward_sense[i][b].first;
      auto &num_labels = dir ? label_array_in_forward_sense[i][b].second : label_array_in_backward_sense[i][b].second;
      if (num_labels) min_rc_per_bin = label_list[0]->rc;
      else min_rc_per_bin = LARGE_FLOAT;
      if constexpr (dir) {
        rc2_bin_in_forward_sense[i][b] = min_rc_per_bin;
        rc = std::min(rc, min_rc_per_bin);
        rc2_till_this_bin_in_forward_sense[i][b] = rc;
      } else {
        rc2_bin_in_backward_sense[i][b] = min_rc_per_bin;
        rc = std::min(rc, min_rc_per_bin);
        rc2_till_this_bin_in_backward_sense[i][b] = rc;
      }
    }
  }
}
template<bool if_symmetry>
int CVRP::generateColsByBidir(BbNode *node) {

  PtrAllR1CS ptrr1cs(node, this);

  RE_TRY:
  rollback = 0;
  runLabeling<true, false, if_symmetry>(node, ptrr1cs);

  if (rollback == 2) {
    reallocateLabel();
    goto RE_TRY;
  } else if (rollback == 1) {
    return 0;
  }

  if (!if_symmetry) {
    RE_TRY2:
    rollback = 0;
    runLabeling<false, false, if_symmetry>(node, ptrr1cs);//ccnt can only be applied for forward case

    if (rollback == 2) {
      reallocateLabel();
      goto RE_TRY2;
    } else if (rollback == 1) {
      return 0;
    }
  }


  int ccnt = concatenateCols_prior_forward<if_symmetry>(node, ptrr1cs);

  if (abs(gap_between_last_smallest_rc_and_rc_threshold) < TOLERANCE) {
    double NumExistedLabels = 0;
    double NumExistedLabel_back = 0;
    for (int i = 1; i < dim; ++i) {
      for (int b = 0; b < num_buckets_per_vertex; ++b) {
        NumExistedLabels += label_array_in_forward_sense[i][b].second;
        if constexpr (!if_symmetry) {
          NumExistedLabel_back += label_array_in_backward_sense[i][b].second;
        }
      }
    }
    ratio_dominance_checks_non_dominant.first +=
        if_symmetry ? num_dominance_checks / NumExistedLabels : num_dominance_checks
            / (NumExistedLabels + NumExistedLabel_back);
    ++ratio_dominance_checks_non_dominant.second;
    if constexpr (!if_symmetry) {
      double dif = abs(NumExistedLabels - NumExistedLabel_back);
      double over = dif / std::min(NumExistedLabels, NumExistedLabel_back);
#ifdef DETAILED_EXACT_PRINT_INFO
      std::cout << "over= " << over << std::endl;
#endif
      if (over > Config::NumberOfOverLabelsInMeetPoint) {
#ifdef DETAILED_EXACT_PRINT_INFO
        std::cout << "we adjust the meetpoint!" << std::endl;
#endif
        if (NumExistedLabels > NumExistedLabel_back) {
          meet_point_resource_in_bi_dir *= (1 - Config::MeetPointFactor);
        } else {
          meet_point_resource_in_bi_dir *= (1 + Config::MeetPointFactor);
        }
#ifdef DETAILED_EXACT_PRINT_INFO
        std::cout << "meet_point_resource_in_bi_dir= " << meet_point_resource_in_bi_dir << std::endl;
#endif
      }
    }
  }

  if (ccnt) {
    smallest_rc = std::get<2>(negative_rc_label_tuple[0]);
    if (abs(gap_between_last_smallest_rc_and_rc_threshold) > TOLERANCE) {
      gap_between_last_smallest_rc_and_rc_threshold = std::get<2>(negative_rc_label_tuple.back()) - smallest_rc;
      if (gap_between_last_smallest_rc_and_rc_threshold < TOLERANCE)gap_between_last_smallest_rc_and_rc_threshold = 0;
    } else {
      if (ceilTransformedNumberRelated(smallest_rc * num_vehicle + lp_val + RC_TOLERANCE) + TOLERANCE >= ub) {
        node->value = lp_val;
        node->is_terminated = true;
        std::cout << TERMINATED_MESSAGE_PROMISING_VEHICLES;
        return 0;
      }
    }
#ifdef DETAILED_EXACT_PRINT_INFO
    std::cout << "Smallest= " << smallest_rc << endl;
#endif
  } else {
    smallest_rc = 0;
    gap_between_last_smallest_rc_and_rc_threshold = 0;
  }
#ifdef DETAILED_EXACT_PRINT_INFO
  std::cout << "ccnt= " << ccnt << endl;
#endif
  if (!ccnt) return 0;

  addColumns(node, ccnt);

  return ccnt;
}
template<bool if_symmetry>
int CVRP::concatenateCols_prior_forward(BbNode *node, const PtrAllR1CS &ptrAllR1Cs) {
  prior_pool_beg4_pricing = pool_beg4_pricing;
  if (checkPricingPool()) reallocatePricingPool();
  int index = num_col - 1;
  double tmp_rc, tmp_mainResource;
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  int *tmp_seq = new int[int(max_main_resource) + 3];
  int size_tmp_seq;
  int i, j, arr_bj;
  Label *p;
  int if_state;//==0 means work

  for (auto &label_list : concatenate_labels_in_forward_cg) {
    i = label_list.first.first;
    j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      tmp_mainResource = pr.second;
      arr_bj =
          if_symmetry ? int((max_main_resource - tmp_mainResource) / step_size) : (int) (tmp_mainResource / step_size);
      tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];

      concatenateOneLabelWithOtherLabels<true, if_symmetry, true, false>(ki,
                                                                         j,
                                                                         arr_bj,
                                                                         tmp_rc,
                                                                         tmp_mainResource,
                                                                         r1c_to_pi,
                                                                         r1c_multi_to_pi,
                                                                         if_state);
      if (if_state == -2)continue;
      for (if_symmetry ? --arr_bj : ++arr_bj; if_symmetry ? arr_bj >= 0 : arr_bj < num_buckets_per_vertex;
           if_symmetry ? --arr_bj : ++arr_bj) {
        concatenateOneLabelWithOtherLabels<true, if_symmetry, false, false>(ki,
                                                                            j,
                                                                            arr_bj,
                                                                            tmp_rc,
                                                                            tmp_mainResource,
                                                                            r1c_to_pi,
                                                                            r1c_multi_to_pi,
                                                                            if_state);
        if (if_state == -2)break;
      }
    }
  }

  writeColumnsInPricingPool(node, index);

  delete[]tmp_seq;
  return index - num_col + 1;
}

template<bool dir>
void CVRP::populateTellWhichBin4ArcElimination() {
  /**
   * one time calculation except for regenerate the bucket graph
   */
  if constexpr (dir) {
    if (!tell_which_bin4_arc_elimination_in_forward_sense.empty()) return;
  } else {
    if (!tell_which_bin4_arc_elimination_in_backward_sense.empty()) return;
  }
  std::cout << "populateTellWhichBin4ArcElimination" << std::endl;
  size_t size = dim * dim * num_buckets_per_vertex;
  int dim_sq = dim * dim;
  double tmp_mainResource;
  if constexpr (dir) tell_which_bin4_arc_elimination_in_forward_sense.reserve(size);
  else tell_which_bin4_arc_elimination_in_backward_sense.reserve(size);
  std::unordered_map<int, int> map_bj_B;
  std::vector<std::pair<int, int>> vec_bj_B(num_buckets_per_vertex);
  map_bj_B.reserve(num_buckets_per_vertex + 1);
  for (int b = 0; b <= num_buckets_per_vertex; ++b) map_bj_B[b] = dir ? -1 : num_buckets_per_vertex;
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      int base = i * dim + j;
      for (auto &key_value : map_bj_B) key_value.second = dir ? -1 : num_buckets_per_vertex;
      for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
        double res = dir ? b * step_size : std::min((b + 1) * step_size - TOLERANCE, max_main_resource);
        if (dir ? !increaseMainResourceConsumption(res, tmp_mainResource, i, j) :
            !decreaseMainResourceConsumption(res, tmp_mainResource, i, j))
          break;
        int bj = int(tmp_mainResource / step_size);
        if (dir ? map_bj_B[bj] < b : map_bj_B[bj] > b) map_bj_B[bj] = b;
      }
      for (int bj = dir ? num_buckets_per_vertex - 1 : 0; dir ? bj >= 0 : bj < num_buckets_per_vertex; dir ? --bj : ++bj) {
        if (map_bj_B[bj] != (dir ? -1 : num_buckets_per_vertex)) {
          for (int bj2 = dir ? bj + 1 : bj - 1; dir ? bj2 < num_buckets_per_vertex : bj2 >= 0; dir ? ++bj2 : --bj2)
            map_bj_B[bj2] = map_bj_B[bj];
          break;
        }
      }
      for (int bj = 0; bj < num_buckets_per_vertex; ++bj) {
        (dir ? tell_which_bin4_arc_elimination_in_forward_sense[base + bj * dim_sq] :
		 tell_which_bin4_arc_elimination_in_backward_sense[base + bj * dim_sq]) = map_bj_B[bj];
      }
    }
  }
}

template<bool dir>
void CVRP::obtainjumpArcs(BbNode *node, std::bitset<2> **bitMap) const {
  int num_jump_arcs = 0;
  bool if_used;

  for (int i = 1; i < dim; ++i) {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      if constexpr (dir) {
        node->all_forward_buckets[i][b].jump_arcs.clear();
      } else {
        node->all_backward_buckets[i][b].jump_arcs.clear();
      }
      for (int j = 1; j < dim; ++j) bitMap[j][b] = 2;
      bitMap[i][b] = 1;
      for (int j : (dir ? node->all_forward_buckets[i][b].bucket_arcs :
                    node->all_backward_buckets[i][b].bucket_arcs))
        bitMap[j][b] = 0;
    }
    for (int b = (dir ? 0 : num_buckets_per_vertex - 1); dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
      for (int j = 1; j < dim; ++j) {
        if (bitMap[j][b] == 2) {//需要有jump arc
          if_used = false;
          for (int b4_i = (dir ? b + 1 : b - 1); dir ? b4_i < num_buckets_per_vertex : b4_i >= 0;
               dir ? ++b4_i : --b4_i) {
            if (bitMap[j][b4_i] == 0) {//发现有去的arc了
              std::pair<int, int> map = (dir ? std::make_pair(b4_i * step_size, j) :
                                         std::make_pair((b4_i + 1) * step_size - TOLERANCE, j));
              for (int tmp_b = b; dir ? tmp_b < b4_i : tmp_b > b4_i; dir ? ++tmp_b : --tmp_b) {//不取等号< qi
                bitMap[j][tmp_b] = 1;
                if constexpr (dir) {
                  node->all_forward_buckets[i][tmp_b].jump_arcs.emplace_back(map);
                } else {
                  node->all_backward_buckets[i][tmp_b].jump_arcs.emplace_back(map);
                }
                ++num_jump_arcs;
              }
              if_used = true;
              break;
            }
          }
          if (!if_used) {
            for (int tmp_b = b; dir ? tmp_b < num_buckets_per_vertex : tmp_b >= 0; dir ? ++tmp_b : --tmp_b)
              bitMap[j][tmp_b] = 1;
          }
        }
      }
    }
  }

  (dir ? node->num_forward_jump_arcs : node->num_backward_jump_arcs) = num_jump_arcs;

  std::cout << "Obtain" << (dir ? "Forward" : "Backward") << " Jump Arcs= " << num_jump_arcs << std::endl;
}

template<bool dir, bool if_symmetry>
int CVRP::enumerateHalfwardRoutes(BbNode *node,
                                  const double *r1c_to_pi,
                                  const double *r1c_multi_to_pi,
                                  std::unordered_map<yzzLong, std::tuple<Label *, Label *, double>> &Tags,
                                  std::vector<Label *> **copy_bucket,
                                  int &num_routes_now) {
  int status = 0;
  int edgemap;
  (dir ? num_forward_labels_in_enu : num_backward_labels_in_enu) = 0;
  int Max_routes_phase1 = Config::MaxNumRouteInEnumeration_half;
  bool if_keep, if_break;
  double path_rc, path_cost;
  constexpr bool if_dif = dir ^ if_symmetry;
  if constexpr (dir) {
    initializeLabels(node, 1, false, {true, 1, true});
  } else {
    initializeLabels(node, 2, false, {true, 2, false});
  }

  auto beg = std::chrono::high_resolution_clock::now();
  auto end = beg;
  auto b4_end = beg;
  auto af_end = beg;
  double eps;
  double eps2;
  double left_time = Config::HardTimeThresholdInAllEnumeration;

  for (int b = dir ? 0 : num_buckets_per_vertex - 1; dir ? b < num_buckets_per_vertex : b >= 0; dir ? ++b : --b) {
    int i = 1;
    STILL_EXIST:
    for (; i < dim; ++i) {
      end = std::chrono::high_resolution_clock::now();
      eps = std::chrono::duration<double>(end - b4_end).count();
      if (eps > left_time) {
        status = 2;
        goto outside;
      }
      auto &valid_num = (dir ? if_exist_extra_labels_in_forward_sense[i][b].second :
                         if_exist_extra_labels_in_backward_sense[i][b].second);
      if (!valid_num) continue;
      auto &label_array = (dir ? if_exist_extra_labels_in_forward_sense[i][b].first :
                           if_exist_extra_labels_in_backward_sense[i][b].first);
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->is_extended) continue;
        ki->is_extended = true;
        for (int j : (dir ? node->all_forward_buckets[i][b].bucket_arcs :
                      node->all_backward_buckets[i][b].bucket_arcs)) {
          if (ki->pi[j]) continue;
          auto &tmp_mainResource = all_label[idx_glo].sum_main_resource;
          if (dir ? !increaseMainResourceConsumption(ki->sum_main_resource, tmp_mainResource, i, j) :
              !decreaseMainResourceConsumption(ki->sum_main_resource, tmp_mainResource, i, j))
            continue;
          auto &tmp_rc = all_label[idx_glo].rc;
          tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];//real rc
          if_keep = false;
          int arr_bj = (if_symmetry ? int((max_main_resource - tmp_mainResource) / step_size) :
                        int(tmp_mainResource / step_size));
          if constexpr (dir) {
            if (tmp_mainResource > meet_point_resource_in_bi_dir_enu) {
              if constexpr (if_symmetry) {
                if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap)
                  concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_mainResource);
              } else {
                if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap)
                  concatenate_labels_in_forward_cg[{i, j}].emplace_back(ki, tmp_mainResource);
              }
              continue;
            }
          } else {
            if (tmp_mainResource < meet_point_resource_in_bi_dir_enu) {
              if constexpr (if_symmetry) {
                if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap)
                  concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_mainResource);
              } else {
                if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap)
                  concatenate_labels_in_backward_cg[{i, j}].emplace_back(ki, tmp_mainResource);
              }
              continue;
            }
          }
          if (tmp_rc + (if_dif ? rc2_till_this_bin_in_backward_sense[j][arr_bj] : rc2_till_this_bin_in_forward_sense[j][arr_bj])
              > opt_gap)
            continue;
          if (tmp_rc + (if_dif ? rc2_bin_in_backward_sense[j][arr_bj] : rc2_bin_in_forward_sense[j][arr_bj])
              < opt_gap) {
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if constexpr (if_symmetry) {
                if (tmp_mainResource + kkj->sum_main_resource > max_main_resource) continue;
              } else {
                if constexpr (dir) {
                  if (tmp_mainResource > kkj->sum_main_resource) continue;
                } else {
                  if (tmp_mainResource < kkj->sum_main_resource) continue;
                }
              }
              if ((ki->pi & kkj->pi).any()) continue;
              path_rc = tmp_rc + kkj->rc;
              if (path_rc > opt_gap) break;
              if (ki->num_valid_rank1_cut < kkj->num_valid_rank1_cut) {
                for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
                  if (kkj->rank1_cut_mem[ki->valid_rank1_cut[l]]) {
                    path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
                    if (path_rc > opt_gap)goto here;
                  }
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut; ++l) {
                  if (ki->rank1_cut_mem[kkj->valid_rank1_cut[l]]) {
                    path_rc -= r1c_to_pi[kkj->valid_rank1_cut[l]];
                    if (path_rc > opt_gap)goto here;
                  }
                }
              }

              if (ki->num_valid_rank1_cut_multi < kkj->num_valid_rank1_cut_multi) {
                for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = ki->valid_rank1_cut_multi[l];
                  if (kkj->rank1_cut_mem_multi[tmp_cut] +
                      ki->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > opt_gap)goto here;
                  }
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = kkj->valid_rank1_cut_multi[l];
                  if (ki->rank1_cut_mem_multi[tmp_cut] +
                      kkj->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > opt_gap)goto here;
                  }
                }
              }
              if_keep = true;
              goto outside1;
              here:;
            }
          }
          for (if_dif ? ++arr_bj : --arr_bj; if_dif ? arr_bj < num_buckets_per_vertex : arr_bj >= 0;
               if_dif ? ++arr_bj : --arr_bj) {
            if (tmp_rc + (if_dif ? rc2_till_this_bin_in_backward_sense[j][arr_bj] : rc2_till_this_bin_in_forward_sense[j][arr_bj])
                > opt_gap)
              break;
            if (tmp_rc + (if_dif ? rc2_bin_in_backward_sense[j][arr_bj] : rc2_bin_in_forward_sense[j][arr_bj])
                > opt_gap)
              continue;
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if ((ki->pi & kkj->pi).any()) continue;
              path_rc = tmp_rc + kkj->rc;
              if (path_rc > opt_gap) break;
              if (ki->num_valid_rank1_cut < kkj->num_valid_rank1_cut) {
                for (int l = 0; l < ki->num_valid_rank1_cut; ++l) {
                  if (kkj->rank1_cut_mem[ki->valid_rank1_cut[l]]) {
                    path_rc -= r1c_to_pi[ki->valid_rank1_cut[l]];
                    if (path_rc > opt_gap)goto here1;
                  }
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut; ++l) {
                  if (ki->rank1_cut_mem[kkj->valid_rank1_cut[l]]) {
                    path_rc -= r1c_to_pi[kkj->valid_rank1_cut[l]];
                    if (path_rc > opt_gap)goto here1;
                  }
                }
              }

              if (ki->num_valid_rank1_cut_multi < kkj->num_valid_rank1_cut_multi) {
                for (int l = 0; l < ki->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = ki->valid_rank1_cut_multi[l];
                  if (kkj->rank1_cut_mem_multi[tmp_cut] +
                      ki->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > opt_gap)goto here1;
                  }
                }
              } else {
                for (int l = 0; l < kkj->num_valid_rank1_cut_multi; ++l) {
                  int tmp_cut = kkj->valid_rank1_cut_multi[l];
                  if (ki->rank1_cut_mem_multi[tmp_cut] +
                      kkj->rank1_cut_mem_multi[tmp_cut]
                      >= r1c_multi_denominator_in_cg[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > opt_gap)goto here1;
                  }
                }
              }

              if_keep = true;
              goto outside1;
              here1:;
            }
          }
          outside1:
          if (!if_keep) continue;
          int bj = int(tmp_mainResource / step_size);
          auto &labelList_j = dir ? label_array_in_forward_sense[j][bj].first : label_array_in_backward_sense[j][bj].first;
          auto &valid_num_j = dir ? label_array_in_forward_sense[j][bj].second : label_array_in_backward_sense[j][bj].second;
          auto &tmp_PI = all_label[idx_glo].pi;
          auto &tmp_Cost = all_label[idx_glo].cost;
          auto &tmp_Rank1CutMem = all_label[idx_glo].rank1_cut_mem;
          auto &tmp_num_valid_rank1_cut = all_label[idx_glo].num_valid_rank1_cut;
          auto &tmp_valid_rank1_cut = all_label[idx_glo].valid_rank1_cut;
          tmp_PI = ki->pi;
          tmp_PI.set(j);
          tmp_Cost = ki->cost + cost_mat4_vertex[i][j];
          tmp_Rank1CutMem = ki->rank1_cut_mem;
          for (auto l : std::get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_Rank1CutMem[l] = false;
              tmp_rc -= r1c_to_pi[l];
            } else tmp_Rank1CutMem[l] = true;
          }
          tmp_Rank1CutMem &= std::get<1>(Vertex2ActiveInOnePricingR1Cs[j]);

          auto &tmp_Rank1CutMem_multi = all_label[idx_glo].rank1_cut_mem_multi;
          auto &tmp_num_valid_rank1_cut_multi = all_label[idx_glo].num_valid_rank1_cut_multi;
          auto &tmp_valid_rank1_cut_multi = all_label[idx_glo].valid_rank1_cut_multi;
          std::copy(ki->rank1_cut_mem_multi, ki->rank1_cut_mem_multi + num_valid_r1c_multi_in_cg, tmp_Rank1CutMem_multi);
          for (auto &l : std::get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
            int tmp_cut = std::get<0>(l);
            tmp_Rank1CutMem_multi[tmp_cut] += std::get<1>(l);
            if (tmp_Rank1CutMem_multi[tmp_cut] >= std::get<2>(l)) {
              tmp_rc -= r1c_multi_to_pi[tmp_cut];
              tmp_Rank1CutMem_multi[tmp_cut] -= std::get<2>(l);
            }
          }
          for (auto l : std::get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;

          if_break = false;

          for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
            auto &kj = labelList_j[vec_index_j];
#ifdef CAPACITY_AS_MAIN_RESOURCE
            if (abs(kj->sum_main_resource - tmp_mainResource) > TOLERANCE) {
              ++vec_index_j;
              continue;
            }
            if ((kj->pi ^ tmp_PI).none()) {
              if (kj->cost > tmp_Cost) {
                kj->is_extended = true;
                kj = labelList_j[--valid_num_j];
                dir ? --num_forward_labels_in_enu : --num_backward_labels_in_enu;
              } else {
                if_break = true;
                break;
              }
            } else ++vec_index_j;
#else
            if (kj->cost > tmp_Cost) {
              if (dir ? kj->sum_main_resource > tmp_mainResource : kj->sum_main_resource < tmp_mainResource) {
                if ((kj->pi ^ tmp_PI).none()) {
                  kj->is_extended = true;
                  kj = labelList_j[--valid_num_j];
                  dir ? --num_forward_labels_in_enu : --num_backward_labels_in_enu;
                } else ++vec_index_j;
              } else ++vec_index_j;
            } else {
              if (dir ? kj->sum_main_resource < tmp_mainResource : kj->sum_main_resource > tmp_mainResource) {
                if ((kj->pi ^ tmp_PI).none()) {
                  if_break = true;
                  break;
                } else ++vec_index_j;
              } else ++vec_index_j;
            }
#endif
          }

          if (if_break) continue;

          tmp_num_valid_rank1_cut = 0;
          for (auto l : std::get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
            }
          }

          tmp_num_valid_rank1_cut_multi = 0;
          for (auto l : std::get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
            if (tmp_Rank1CutMem_multi[l]) {
              tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
            }
          }

          labelList_j[valid_num_j++] = all_label + idx_glo;
          if (valid_num_j == labelList_j.size()) {
            labelList_j.resize(labelList_j.size() * 2);
          }

          all_label[idx_glo].p_label= ki;
          all_label[idx_glo].end_vertex = j;
          all_label[idx_glo].is_extended = false;
          auto &bucket = (dir ? if_exist_extra_labels_in_forward_sense[j][bj] : if_exist_extra_labels_in_backward_sense[j][bj]);
          bucket.first[bucket.second++] = all_label + idx_glo;
          if (bucket.second == bucket.first.size()) {
            bucket.first.resize(bucket.first.size() * 2);
          }

          if constexpr (dir) {
            if (tmp_rc + chg_cost_mat4_vertex[j][0] < opt_gap) {
              path_cost = tmp_Cost + cost_mat4_vertex[j][0];

              if (Tags.find(tmp_PI) == Tags.end()) {
                Tags[tmp_PI] = {all_label + idx_glo, nullptr, path_cost};
                ++num_routes_now;
                if (num_routes_now > Max_routes_phase1) {
                  status = 2;//routes limit
                  goto outside;
                }
              } else if (std::get<2>(Tags[tmp_PI]) > path_cost) {
                Tags[tmp_PI] = {all_label + idx_glo, nullptr, path_cost};
              }
            }
          }

          if ((dir ? ++num_forward_labels_in_enu : ++num_backward_labels_in_enu) > Config::MaxNumLabelInEnumeration) {
            status = 3;//all labels limit
            goto outside;
          }
          ++idx_glo;//can be put here, because once go outside, the function will end
          if (idx_glo == label_assign) {
            rollback = 2;
            goto outside;
          }
        }
      }
      valid_num = 0;
    }
    for (i = 1; i < dim; ++i) {
      if (dir ? if_exist_extra_labels_in_forward_sense[i][b].second : if_exist_extra_labels_in_backward_sense[i][b].second)
        goto STILL_EXIST;
    }
    af_end = std::chrono::high_resolution_clock::now();
    eps2 = std::chrono::duration<double>(af_end - b4_end).count();
    eps = std::chrono::duration<double>(af_end - beg).count();
    left_time = (Config::HardTimeThresholdInAllEnumeration - eps) / (dir ? (num_buckets_per_vertex - b) : (b + 1));
    if (eps2 > left_time) {
      status = 2;
      goto outside;
    }
    b4_end = af_end;
  }
  outside:
  for (int i = 0; i < dim; ++i) {
    delete[]copy_bucket[i];
  }
  delete[] copy_bucket;
  std::cout << "Half" << (dir ? "Forward" : "Backward") << " labeling: num_labels= "
            << (dir ? num_forward_labels_in_enu : num_backward_labels_in_enu)
            << " num_routes= " << num_routes_now <<
            std::endl;
  if (status)return status;

  for (int i = 1; i < dim; ++i) {
    for (int b = 0; b < num_buckets_per_vertex; ++b) {
      auto &label_array = dir ? label_array_in_forward_sense[i][b] : label_array_in_backward_sense[i][b];
      std::sort(label_array.first.begin(), label_array.first.begin() + label_array.second, CmpLabelRCLess);
    }
  }
  return 0;
}