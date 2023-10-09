//
// Created by Zhengzhong You on 4/5/23.
//

#include "CVRP.hpp"

template<typename T, bool dir, bool if_last_half, bool if_symmetry>
int CVRP::extendKernel4Exact(LABEL *&ki,
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
      addPathByRC(AllLabel[IdxGlo].RC + ChgCostMat4Vertex[j][0],
                  AllLabel + IdxGlo,
                  nullptr,
                  CONFIG::MaxNumRoutesInExact);
    }

    ++IdxGlo;
    if (IdxGlo == LabelAssign) {
      Rollback = 2;
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
void CVRP::updateLabel(double res, LABEL *&ki, int i, int j, int &bj,
                       const double *r1c_to_pi,
                       const double *r1c_multi_to_pi,
                       int min_sorted_b,
                       bool &if_suc) {
  if_suc = false;
  if (ki->PI[j]) return;
  double &which_rc = (if_std_optgap ? OptGap : rc_std);
  auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
  auto &tmp_rc = AllLabel[IdxGlo].RC;
  tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
  if constexpr (!if_res_updated) {
    if constexpr (dir) {
      if (!increaseMainResourceConsumption(res, tmp_mainResource, i, j)) return;
      if constexpr (!if_last_half) {
        if (tmp_mainResource > MeetPointResourceInBiDir) {
//          int arr_bj = (if_symmetry ? int((MaxMainResource - tmp_mainResource) / StepSize) :
//                        int(tmp_mainResource / StepSize));
//          if constexpr (if_symmetry) {
//            if ((min_sorted_b < arr_bj) || (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc < which_rc)) {
//              concatenateLabelsInForwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
//            }
//          } else {
//            if ((min_sorted_b > arr_bj) || (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc < which_rc)) {
//              concatenateLabelsInForwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
//            }
//          }//cannot use this, because the labels might be used in arc elimination!
          concatenateLabelsInForwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
          return;
        }
      }
    } else {
      if (!decreaseMainResourceConsumption(res, tmp_mainResource, i, j)) return;
      if constexpr (!if_last_half) {
        if (tmp_mainResource < MeetPointResourceInBiDir) {
//          int arr_bj = if_symmetry ? int((MaxMainResource - tmp_mainResource) / StepSize) :
//                       int(tmp_mainResource / StepSize);
//          if constexpr (if_symmetry) {
//            if ((min_sorted_b > arr_bj) || (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc < which_rc)) {
//              concatenateLabelsInBackwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
//            }
//          } else {
//            if ((min_sorted_b < arr_bj) || (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc < which_rc)) {
//              concatenateLabelsInBackwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
//            }
//          }
          concatenateLabelsInBackwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
          return;
        }
      }
    }
  }
  if constexpr (if_last_half) {
    //check if the label can pass the completion bounds test
    /**
     * here we only test the RC2TillThisBinInXXXwardSense, we do not get into details of the labels in the bins
     */
//    if constexpr (if_symmetry) {
//      int concate_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//      if (tmp_rc + (dir ? RC2TillThisBinInForwardSense[j][concate_bj] : RC2TillThisBinInBackwardSense[j][concate_bj])
//          > OptGap)
//        return;
//    } else {
//      int concate_bj = int(tmp_mainResource / StepSize);
//      if (tmp_rc + (dir ? RC2TillThisBinInBackwardSense[j][concate_bj] : RC2TillThisBinInForwardSense[j][concate_bj])
//          > OptGap)
//        return;
//    }
    int concate_bj;
    int if_state;
    if constexpr (if_symmetry) {
      concate_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
    } else {
      concate_bj = int(tmp_mainResource / StepSize);
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
    for (if_dif ? ++concate_bj : --concate_bj; if_dif ? concate_bj < NumBucketsPerVertex : concate_bj >= 0;
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
  bj = int(tmp_mainResource / StepSize);
  auto &tmp_PI = AllLabel[IdxGlo].PI;
  auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
  auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
  auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
  auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
  auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;

  tmp_PI = (ki->PI) & (NGMem4Vertex[j]);
  tmp_PI.set(j);
  tmp_Rank1CutMem = ki->Rank1CutMem;

  auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
  std::copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);

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
void CVRP::doDominance(LABEL *&ki, int j, int bj, const double *r1c_to_pi,
                       const double *r1c_multi_to_pi, bool &if_suc) {
  auto &labelList_j = dir ? LabelArrayInForwardSense[j][bj].first : LabelArrayInBackwardSense[j][bj].first;
  auto &valid_num_j = dir ? LabelArrayInForwardSense[j][bj].second : LabelArrayInBackwardSense[j][bj].second;
  auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
  auto &tmp_rc = AllLabel[IdxGlo].RC;
  auto &tmp_PI = AllLabel[IdxGlo].PI;
  auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
  auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
  auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
  auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
  auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
  auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;

//  num_dominance_test = 0;
  double tmp_rc_add = tmp_rc + RC_TOLERANCE, tmp_rc_sub = tmp_rc - RC_TOLERANCE;
  if_suc = true;
  double dif;
  for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
    auto &kj = labelList_j[vec_index_j];
//    ++num_dominance_test;
    ++NumDominanceChecks;
    if (define_dir_resource<dir>(kj->Sum_MainResource, tmp_mainResource)) {
      if (kj->RC < tmp_rc_sub) {
        if (((tmp_PI & kj->PI) ^ (kj->PI)).none()) {
          dif = tmp_rc_sub + GapBetweenLastSmallestRCAndRCThreshold;
          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
            if (!tmp_Rank1CutMem[kj->validRank1Cut[l]]) {
              dif += r1c_to_pi[kj->validRank1Cut[l]];
              if (dif < kj->RC) goto here;
            }
          }

          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
            int tmp = kj->validRank1Cut_multi[l];
            if (kj->Rank1CutMem_multi[tmp]
                > tmp_Rank1CutMem_multi[tmp]) {
              dif += r1c_multi_to_pi[tmp];
              if (dif < kj->RC) goto here;
            }
          }
          if_suc = false;
          break;//taken
          here:;
        }
      }
      ++vec_index_j;
    } else if (define_dir_resource<dir>(tmp_mainResource, kj->Sum_MainResource)) {
      if (tmp_rc_add < kj->RC) {
        if (((tmp_PI & kj->PI) ^ (tmp_PI)).none()) {
          dif = tmp_rc_add - GapBetweenLastSmallestRCAndRCThreshold;
          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]]) {
              dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
              if (dif > kj->RC) goto there;
            }
          }

          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
            int tmp = tmp_valid_rank1_cut_multi[l];
            if (tmp_Rank1CutMem_multi[tmp]
                > kj->Rank1CutMem_multi[tmp]) {
              dif -= r1c_multi_to_pi[tmp];
              if (dif > kj->RC) goto there;
            }
          }
          kj->if_extended = true;
          kj = labelList_j[--valid_num_j];
          goto _there;
          there:
          ++vec_index_j;
          _there:;
        } else ++vec_index_j;
      } else ++vec_index_j;
    } else {
      if (kj->RC < tmp_rc_add) {
        if (((tmp_PI & kj->PI) ^ (kj->PI)).none()) {
          dif = tmp_rc_sub + GapBetweenLastSmallestRCAndRCThreshold;
          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
            if (!tmp_Rank1CutMem[kj->validRank1Cut[l]]) {
              dif += r1c_to_pi[kj->validRank1Cut[l]];
              if (dif < kj->RC) goto here1;
            }
          }

          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
            int tmp = kj->validRank1Cut_multi[l];
            if (kj->Rank1CutMem_multi[tmp]
                > tmp_Rank1CutMem_multi[tmp]) {
              dif += r1c_multi_to_pi[tmp];
              if (dif < kj->RC) goto here1;
            }
          }
          if_suc = false;
          break;//taken
          here1:;
        }
        ++vec_index_j;
      } else if (tmp_rc_sub < kj->RC) {
        if (((tmp_PI & kj->PI) ^ (tmp_PI)).none()) {
          dif = tmp_rc_add - GapBetweenLastSmallestRCAndRCThreshold;
          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]]) {
              dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
              if (dif > kj->RC) goto there1;
            }
          }

          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
            int tmp = tmp_valid_rank1_cut_multi[l];
            if (tmp_Rank1CutMem_multi[tmp]
                > kj->Rank1CutMem_multi[tmp]) {
              dif -= r1c_multi_to_pi[tmp];
              if (dif > kj->RC) goto there1;
            }
          }
          kj->if_extended = true;
          kj = labelList_j[--valid_num_j];
          goto _there1;
          there1:
          ++vec_index_j;
          _there1:;
        } else ++vec_index_j;
      } else {//q & rc all equal now?
        if ((tmp_PI ^ kj->PI).none()) {//all equal
          std::bitset<2> who_win = 1;//1==tmp_win,2==kj win,0==all keep
          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]]) {
              who_win = 2;
              goto next_test1;
            }
          }
          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
            int tmp = tmp_valid_rank1_cut_multi[l];
            if (tmp_Rank1CutMem_multi[tmp]
                > kj->Rank1CutMem_multi[tmp]) {
              who_win = 2;
              goto next_test1;
            }
          }
          next_test1:
          if (who_win == 1) {
            kj->if_extended = true;
            kj = labelList_j[--valid_num_j];
          } else {
            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
              if (!tmp_Rank1CutMem[kj->validRank1Cut[l]]) {
                who_win = 0;
                goto next_test2;
              }
            }
            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
              int tmp = kj->validRank1Cut_multi[l];
              if (kj->Rank1CutMem_multi[tmp]
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
          yzzLong tmp = tmp_PI & kj->PI;
          if ((tmp ^ tmp_PI).none()) {
            dif = tmp_rc_add - GapBetweenLastSmallestRCAndRCThreshold;
            for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
              if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]]) {
                dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
                if (dif > kj->RC) goto there2;
              }
            }

            for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
              int tmp2 = tmp_valid_rank1_cut_multi[l];
              if (tmp_Rank1CutMem_multi[tmp2]
                  > kj->Rank1CutMem_multi[tmp2]) {
                dif -= r1c_multi_to_pi[tmp2];
                if (dif > kj->RC) goto there2;
              }
            }

            kj->if_extended = true;
            kj = labelList_j[--valid_num_j];
            goto _there2;

            there2:
            ++vec_index_j;
            _there2:;

          } else if ((tmp ^ kj->PI).none()) {
            dif = tmp_rc_sub + GapBetweenLastSmallestRCAndRCThreshold;
            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
              if (!tmp_Rank1CutMem[kj->validRank1Cut[l]]) {
                dif += r1c_to_pi[kj->validRank1Cut[l]];
                if (dif < kj->RC) goto here2;
              }
            }

            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
              int tmp2 = kj->validRank1Cut_multi[l];
              if (kj->Rank1CutMem_multi[tmp2]
                  > tmp_Rank1CutMem_multi[tmp2]) {
                dif += r1c_multi_to_pi[tmp2];
                if (dif < kj->RC) goto here2;
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
//    succ_ratio.first += num_dominance_test;
//    succ_ratio.second++;
    labelList_j[valid_num_j++] = AllLabel + IdxGlo;
    if (valid_num_j == labelList_j.size()) {
      labelList_j.resize(labelList_j.size() * 2);
    }

    //Seq
    AllLabel[IdxGlo].PLabel = ki;
    //EndVertex
    AllLabel[IdxGlo].EndVertex = j;
    //if_extended
    AllLabel[IdxGlo].if_extended = false;
    //bucket
    auto &bucket = dir ? IfExistExtraLabelsInForwardSense[j][bj] : IfExistExtraLabelsInBackwardSense[j][bj];
    bucket.first[bucket.second++] = AllLabel + IdxGlo;
    if (bucket.second == bucket.first.size()) {
      bucket.first.resize(bucket.first.size() * 2);
    }
  }
//  else {
//    fail_ratio.first += num_dominance_test;
//    fail_ratio.second++;
//  }
}
template<bool dir>
void CVRP::checkIfDominated(LABEL *&ki, int i, int b, const double *r1c_to_pi,
                            const double *r1c_multi_to_pi,
                            bool &if_suc) {
  if_suc = true;
  double dif;
  double tmp_ki_rc_sub = ki->RC - RC_TOLERANCE;
  for (int b4_b = (dir ? b - 1 : b + 1); dir ? b4_b >= 0 : b4_b < NumBucketsPerVertex; dir ? --b4_b : ++b4_b) {
    auto &b4_label_list = dir ? LabelArrayInForwardSense[i][b4_b].first : LabelArrayInBackwardSense[i][b4_b].first;
    auto &b4_valid_num = dir ? LabelArrayInForwardSense[i][b4_b].second : LabelArrayInBackwardSense[i][b4_b].second;
    if ((dir ? RC2TillThisBinInForwardSense[i][b4_b] : RC2TillThisBinInBackwardSense[i][b4_b]) > tmp_ki_rc_sub) break;
    for (int vec_b4 = 0; vec_b4 < b4_valid_num; ++vec_b4) {
      auto &b4_ki = b4_label_list[vec_b4];
      if (b4_ki->RC > tmp_ki_rc_sub) break;
      if (((ki->PI & b4_ki->PI) ^ (b4_ki->PI)).none()) {
        dif = tmp_ki_rc_sub + GapBetweenLastSmallestRCAndRCThreshold;
        for (int l = 0; l < b4_ki->numValidRank1Cut; ++l) {
          if (!ki->Rank1CutMem[b4_ki->validRank1Cut[l]]) {
            dif += r1c_to_pi[b4_ki->validRank1Cut[l]];
            if (dif < b4_ki->RC) goto here;
          }
        }

        for (int l = 0; l < b4_ki->numValidRank1Cut_multi; ++l) {
          int tmp2 = b4_ki->validRank1Cut_multi[l];
          if (b4_ki->Rank1CutMem_multi[tmp2] > ki->Rank1CutMem_multi[tmp2]) {
            dif += r1c_multi_to_pi[tmp2];
            if (dif < b4_ki->RC) goto here;
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
  for (i = 1; i < Dim; ++i) {
    if (dir ? IfExistExtraLabelsInForwardSense[i][b].second : IfExistExtraLabelsInBackwardSense[i][b].second) {
      if_suc = false;
      return;
    }
  }
  if_suc = true;
  //write data
  for (i = 1; i < Dim; ++i) {
    auto &label_array = (dir ? LabelArrayInForwardSense : LabelArrayInBackwardSense);
    auto &b4_label_list = label_array[i][b].first;
    auto &b4_valid_num = label_array[i][b].second;
    std::sort(b4_label_list.begin(), b4_label_list.begin() + b4_valid_num, CmpLabelRCLess);
    min_sorted_b = b;
  }
  for (i = 1; i < Dim; ++i) {
    auto &label_bin = (dir ? LabelArrayInForwardSense[i][b] : LabelArrayInBackwardSense[i][b]);
    auto &rc2_till_this_bin = (dir ? RC2TillThisBinInForwardSense[i][b] : RC2TillThisBinInBackwardSense[i][b]);
    if constexpr (dir) {
      if (label_bin.second) {
        rc2_till_this_bin = (b ? std::min(RC2TillThisBinInForwardSense[i][b - 1], label_bin.first[0]->RC)
                               : label_bin.first[0]->RC);
      } else {
        rc2_till_this_bin = (b ? RC2TillThisBinInForwardSense[i][b - 1] : LARGEFLOAT);
      }
    } else {
      if (label_bin.second) {
        rc2_till_this_bin = (b < NumBucketsPerVertex - 1 ? std::min(RC2TillThisBinInBackwardSense[i][b + 1],
                                                                    label_bin.first[0]->RC)
                                                         : label_bin.first[0]->RC);
      } else {
        rc2_till_this_bin = (b < NumBucketsPerVertex - 1 ? RC2TillThisBinInBackwardSense[i][b + 1] : LARGEFLOAT);
      }
    }
  }
}
template<bool dir, bool if_symmetry, bool if_check_res, bool if_std_optgap>
void CVRP::concatenateOneLabelWithOtherLabels(LABEL *ki, int j, int arr_bj, double tmp_rc, double tmp_mainResource,
                                              const double *r1c_to_pi,
                                              const double *r1c_multi_to_pi,
                                              int &if_state) {
  double path_rc;
  double &which_rc = if_std_optgap ? OptGap : rc_std;
  auto ptr_rc_till_this_bin = &RC2TillThisBinInForwardSense[j][arr_bj];
  auto ptr_rc_bin = &RC2BinInForwardSense[j][arr_bj];
  if constexpr ((dir && !if_symmetry) || (!dir && if_symmetry)) {
    ptr_rc_till_this_bin = &RC2TillThisBinInBackwardSense[j][arr_bj];
    ptr_rc_bin = &RC2BinInBackwardSense[j][arr_bj];
  }

  if (*ptr_rc_till_this_bin + tmp_rc > which_rc) {
    if_state = -2;//no need to continue
    return;
  }
  if_state = -1;//do some operations
  if (*ptr_rc_bin + tmp_rc < which_rc) {//most_negative_rc_in_this_bin
    //add one more condition for testing capacity
    auto &label_arr =
        ((!dir && !if_symmetry) || (dir && if_symmetry)) ? LabelArrayInForwardSense[j][arr_bj].first
                                                         : LabelArrayInBackwardSense[j][arr_bj].first;
    auto &label_valid_num =
        ((!dir && !if_symmetry) || (dir && if_symmetry)) ? LabelArrayInForwardSense[j][arr_bj].second
                                                         : LabelArrayInBackwardSense[j][arr_bj].second;
    for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
      auto &kj = label_arr[vec_index];
      path_rc = kj->RC + tmp_rc;
      if (path_rc > which_rc) break;
      if constexpr (if_check_res) {
        if constexpr (if_symmetry) {
          if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
        } else {
          if (dir ? tmp_mainResource > kj->Sum_MainResource : tmp_mainResource < kj->Sum_MainResource) continue;
        }
      }
      if ((ki->PI & kj->PI).any()) continue;

      if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
        for (int l = 0; l < ki->numValidRank1Cut; ++l) {
          if (kj->Rank1CutMem[ki->validRank1Cut[l]]) {
            path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
            if (path_rc > which_rc) goto here;
          }
        }
      } else {
        for (int l = 0; l < kj->numValidRank1Cut; ++l) {
          if (ki->Rank1CutMem[kj->validRank1Cut[l]]) {
            path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
            if (path_rc > which_rc) goto here;
          }
        }
      }

      if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
        for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
          int tmp_cut = ki->validRank1Cut_multi[l];
          if (kj->Rank1CutMem_multi[tmp_cut] +
              ki->Rank1CutMem_multi[tmp_cut]
              >= R1C_multi_denominator_InCG[tmp_cut]
              ) {
            path_rc -= r1c_multi_to_pi[tmp_cut];
            if (path_rc > which_rc) goto here;
          }
        }
      } else {
        for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
          int tmp_cut = kj->validRank1Cut_multi[l];
          if (ki->Rank1CutMem_multi[tmp_cut] +
              kj->Rank1CutMem_multi[tmp_cut]
              >= R1C_multi_denominator_InCG[tmp_cut]
              ) {
            path_rc -= r1c_multi_to_pi[tmp_cut];
            if (path_rc > which_rc) goto here;
          }
        }
      }
      if constexpr (!if_std_optgap) {
        addPathByRC(path_rc, ki, kj, CONFIG::MaxNumRoutesInExact);
      } else {
        if (path_rc < OptGap) {
          if_state = arr_bj;//what bin to concatenate
          break;
        }
      }
      here:;
    }
  }
}
template<bool dir, bool if_last_half, bool if_symmetry>
void CVRP::runLabeling(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs) {
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
  int min_sorted_b = dir ? -1 : NumBucketsPerVertex;
#ifdef test_time
  double test_time1 = 0, test_time2 = 0, test_time3 = 0, test_time4 = 0;
#endif
  for (int b = (dir ? 0 : NumBucketsPerVertex - 1); (dir ? b < NumBucketsPerVertex : b >= 0); (dir ? ++b : --b)) {
    int i = 1;
    STILL_EXIST:
    for (; i < Dim; ++i) {
      auto &valid_num =
          dir ? IfExistExtraLabelsInForwardSense[i][b].second : IfExistExtraLabelsInBackwardSense[i][b].second;
      if (!valid_num) continue;
      auto &label_array =
          dir ? IfExistExtraLabelsInForwardSense[i][b].first : IfExistExtraLabelsInBackwardSense[i][b].first;
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->if_extended) continue;
#ifdef test_time
        new_beg = std::chrono::high_resolution_clock::now();
#endif

        checkIfDominated<dir>(ki, i, b, r1c_to_pi, r1c_multi_to_pi, if_suc);
#ifdef test_time
        end = std::chrono::high_resolution_clock::now();
        test_time1 += std::chrono::duration<double>(end - new_beg).count();
#endif
        ki->if_extended = true;
        if (!if_suc)continue;
#ifdef test_time
        new_beg = std::chrono::high_resolution_clock::now();
#endif
        auto sig = extendKernel4Exact<int,
                                      dir, if_last_half, if_symmetry>(ki, i, ki->Sum_MainResource,
                                                                      dir ? node->AllForwardBuckets[i][b].BucketArcs
                                                                          : node->AllBackwardBuckets[i][b].BucketArcs,
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
                                                                 dir ? node->AllForwardBuckets[i][b].JumpArcs
                                                                     : node->AllBackwardBuckets[i][b].JumpArcs,
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
    //test if all labels are extended
    end = std::chrono::high_resolution_clock::now();
    eps = std::chrono::duration<double>(end - beg).count();
    if constexpr (if_last_half) {
      if (eps > CONFIG::HardTimeThresholdInArcElimination_last_half) {
        Rollback = 1;
        goto QUIT;
      }
    } else {
      if (eps > CutGenTimeThresholdInPricing) {
        Rollback = 3;
        if (eps > CONFIG::HardTimeThresholdInPricing) {
          if (!ForceNotRollback) {
            Rollback = 1;
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
    LastMaxTimeLabeling = std::max(LastMaxTimeLabeling, eps);
    if (Rollback == 0) {//check if reach the soft time
      if (LastMaxTimeLabeling * HardRollBackFactor > CONFIG::HardTimeThresholdInPricing) {
        Rollback = 3;
      }
    }
    if (dir ? min_sorted_b < NumBucketsPerVertex : min_sorted_b >= 0) {
      for (int i = 1; i < Dim; ++i) {
        for (int b = dir ? min_sorted_b + 1 : min_sorted_b - 1; dir ? b < NumBucketsPerVertex : b >= 0;
             dir ? ++b : --b) {
          auto &label_pr = dir ? LabelArrayInForwardSense[i][b] : LabelArrayInBackwardSense[i][b];
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
    //print the time map
    for (auto &it : test_time_map) {
      std::cout << it.first << " " << it.second << std::endl;
    }
    test_time_map.clear();
  }
#endif
  //changing to arc architecture does not affect the following code
  populateRC2TillThisBinNRC2Bin<dir>(node);

  QUIT:
#ifdef DETAILED_EXACT_PRINT_INFO
  cout << "这里是测试！" << endl;
  {
    if (GapBetweenLastSmallestRCAndRCThreshold < TOLERANCE) {
      cout << "we find if there exists extra labels that can be dominated but saved for some reasons!" << endl;
      //we here generate the cuts
      //rank1 cuts
      safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
      int all_num = int(node->R1Cs.size() + node->R1Cs_multi.size());
      vector<int> denominator(all_num);
      vector<double> pi(all_num);
      vector<vector<std::pair<int, int>>> cuts(Dim);//cut_idx, cut_augment
      vector<unordered_set<int>> mem(Dim);
      vector<vector<int>> no_mem(Dim);
      int cnt = 0;
      for (auto &r1c : node->R1Cs) {
        for (auto &i : r1c.InfoR1C) {
          cuts[i].emplace_back(cnt, 1);
          mem[i].emplace(cnt);
        }
        for (auto &i : r1c.Mem) {
          mem[i].emplace(cnt);
        }
        denominator[cnt] = 2;
        pi[cnt] = Pi[r1c.IdxR1C];
        ++cnt;
      }

      for (auto &r1c : node->R1Cs_multi) {
        auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
        auto &multi = std::get<0>(plan);
        int deno = std::get<1>(plan);
        int count = 0;
        for (auto &i : r1c.InfoR1C.first) {
          cuts[i].emplace_back(cnt, multi[count]);
          mem[i].emplace(cnt);
          ++count;
        }
        for (auto &v : r1c.Mem) {
          mem[v].emplace(cnt);
        }
        denominator[cnt] = deno;
        pi[cnt] = Pi[r1c.IdxR1C];
        ++cnt;
      }

      for (int i = 1; i < Dim; ++i) {
        for (int j = 0; j < all_num; ++j) {
          if (mem[i].find(j) == mem[i].end()) {
            no_mem[i].emplace_back(j);
          }
        }
      }

      for (int b = 0; b < NumBucketsPerVertex; ++b) {
        for (int i = 1; i < Dim; ++i) {
          auto &label_list = LabelArrayInForwardSense[i][b].first;
          auto &valid_num = LabelArrayInForwardSense[i][b].second;
          int label_count = 0;
          HERE:
          if (valid_num <= label_count) continue;
          auto ki = label_list[label_count++];
          auto ng = ki->PI;
          for (int j = label_count; j < valid_num; ++j) {
            auto kj = label_list[j];
            if (ki->Sum_MainResource < kj->Sum_MainResource) {
//              cout << "pass sum main resource!" << endl;
              if (((kj->PI & ng) ^ (ng)).none()) {
//                cout << "now this label pass ng test and rc!" << endl;
//                cout << "now we find out why this label could not be dominate!" << endl;
                double rc_dif = kj->RC - ki->RC;
                //if rank1_dif can bring more rc_dif, then this label indeed cannot be dominated!
                if (rc_dif > TOLERANCE) {
//                  cout << "we find out why this label could not be dominate!" << endl;
                  //find the sequence of ki and kj
                  vector<int> seq_ki;
                  vector<int> seq_kj;
                  auto tmp_ki = ki;
                  while (tmp_ki->PLabel) {
                    seq_ki.emplace_back(tmp_ki->EndVertex);
                    tmp_ki = tmp_ki->PLabel;
                  }
                  std::reverse(seq_ki.begin(), seq_ki.end());
                  auto tmp_kj = kj;
                  while (tmp_kj->PLabel) {
                    seq_kj.emplace_back(tmp_kj->EndVertex);
                    tmp_kj = tmp_kj->PLabel;
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
                    if (k >= node->R1Cs.size()) {
                      int num = k - node->R1Cs.size();
                      pi_copy.emplace_back(pi[k], node->R1Cs_multi[num].InfoR1C.first.size());
                    } else {
                      pi_copy.emplace_back(pi[k], node->R1Cs[k].InfoR1C.size());
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
                  if (cut > node->R1Cs.size()) {
                    cut -= node->R1Cs.size();
                    for (int i = 0; i < cut; ++i) {
                      auto &r1c = node->R1Cs_multi[i];
                      auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
                      auto &multi = std::get<0>(plan);
                      int deno = std::get<1>(plan);
                      int count = 0;
                      // print info
                      for (auto &i : r1c.InfoR1C.first) {
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
                    auto &r1c = node->R1Cs[cut];
                    // print info
                    for (auto &i : r1c.InfoR1C) {
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
  for (int i = 1; i < Dim; ++i) {
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      tmp_forwad_soft.emplace_back(LabelArrayInForwardSense[i][b].second);
    }
  }
  std::sort(tmp_forwad_soft.begin(), tmp_forwad_soft.end(), greater<>());
  int top_0_5 = (int) ((double) tmp_forwad_soft.size() * 0.005);
  int top_1 = (int) ((double) tmp_forwad_soft.size() * 0.01);
  int top_5 = (int) ((double) tmp_forwad_soft.size() * 0.05);
  cout << "for forward exact///max= " << tmp_forwad_soft[0] << "  0.5%=" << tmp_forwad_soft[top_0_5]
       << "  1%=" << tmp_forwad_soft[top_1] << "  5%=" << tmp_forwad_soft[top_5] << endl;

  vector<std::pair<int, double>> record;
  record.reserve(NumRow);

  for (int i = 0; i < node->R1Cs.size(); ++i) {
    int row_idx = node->R1Cs[i].IdxR1C;
    record.emplace_back(i, abs(Pi[row_idx]) * pow(node->R1Cs[i].Mem.size(), CONFIG::MemFactor));
  }

  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
    int row_idx = node->R1Cs_multi[i].IdxR1C;
    record.emplace_back(i + node->R1Cs.size(),
                        abs(Pi[row_idx]) * pow(node->R1Cs_multi[i].Mem.size(), CONFIG::MemFactor));
  }

  double
      sum = accumulate(record.begin(), record.end(), 0.0, [](double a, std::pair<int, double> b) { return a + b.second; });

  cout << "sum= " << sum << endl;
  end = std::chrono::high_resolution_clock::now();
  eps = std::chrono::duration<double>(end - beg).count();
  cout << "forward exact time= " << eps << endl;
#endif
  if (Rollback == 1) {
    std::cout << "Rollback to original states!" << std::endl;
  } else if (Rollback == 2) {
    std::cout << "Rollback with larger Mem!" << std::endl;
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
//  int num_labels = 0;
  constexpr bool if_dif = dir ^ if_symmetry;
  for (auto &label_list : dir ? concatenateLabelsInForwardCG : concatenateLabelsInBackwardCG) {
    int i = label_list.first.first;
    int j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
      tmp_mainResource = pr.second;
//      auto &tmp_rc = AllLabel[IdxGlo].RC;
//      tmp_rc = ki->RC + ChgCostMat4Vertex[i][j]; //no need

//      if constexpr (if_symmetry) {
//        concate_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//      } else {
//        concate_bj = int(tmp_mainResource / StepSize);
//      }
//      concatenateOneLabelWithOtherLabels<dir, if_symmetry, true, true>(ki,
//                                                                       j,
//                                                                       concate_bj,
//                                                                       tmp_rc,
//                                                                       tmp_mainResource,
//
//                                                                       r1c_to_pi,
//                                                                       r1c_multi_to_pi,
//                                                                       if_state);
//      if (if_state == -2) continue;
//      else if (if_state >= 0) goto outside;
//      for (if_dif ? ++concate_bj : --concate_bj; if_dif ? concate_bj < NumBucketsPerVertex : concate_bj >= 0;
//           if_dif ? ++concate_bj : --concate_bj) {
//        concatenateOneLabelWithOtherLabels<dir, if_symmetry, false, true>(ki,
//                                                                          j,
//                                                                          concate_bj,
//                                                                          tmp_rc,
//                                                                          tmp_mainResource,
//
//                                                                          r1c_to_pi,
//                                                                          r1c_multi_to_pi,
//                                                                          if_state);
//        if (if_state == -2) break;
//        else if (if_state >= 0) goto outside;
//      }
//      outside:
//      if (if_state < 0) {
//        continue;
//      }
//      ++num_labels;
      //begin to extend by this direction
      updateLabel<dir, true, if_symmetry, true, true>(tmp_mainResource, ki, i, j, bj, r1c_to_pi, r1c_multi_to_pi,
                                                      dir ? NumBucketsPerVertex - 1 : 0, if_suc);
      if (!if_suc) continue;

      doDominance<dir>(ki, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
      if (!if_suc) continue;

      ++IdxGlo;
      if (IdxGlo == LabelAssign) {
        Rollback = 2;
        goto QUIT;
      }
    }
  }
  end = std::chrono::high_resolution_clock::now();
  eps = std::chrono::duration<double>(end - beg).count();
  if (eps > CONFIG::HardTimeThresholdInArcElimination_Mid_concatenate) {
    Rollback = 1;
    std::cout << "concatenatePhaseInArcElimination time= " << eps << " failed!" << std::endl;
  }
  QUIT:
//  std::cout << "num_labels= " << num_labels << std::endl;
  return;
}

template<bool dir, bool if_symmetry>
void CVRP::eliminatebuketArcs(BBNODE *node, const double *r1c_to_pi,
                              const double *r1c_multi_to_pi,
                              int dim_sq,
                              bool *stateBetween2Buckets,
                              int *latest_bucket) {
  memset(stateBetween2Buckets, 0, dim_sq * NumBucketsPerVertex * sizeof(bool));
  constexpr bool if_dif = dir ^ if_symmetry;
  if constexpr (if_dif) {
    memset(latest_bucket, -1, dim_sq * sizeof(int));
  } else {
    std::fill_n(latest_bucket, dim_sq, NumBucketsPerVertex);
  }
  double tmp_rc, path_rc, tmp_mainResource;
  int num_bucket_arcs = 0;
  for (int b = dir ? 0 : NumBucketsPerVertex - 1; dir ? b < NumBucketsPerVertex : b >= 0; dir ? ++b : --b) {
    for (int i = 1; i < Dim; ++i) {
      concatenateTestKernelInArcElimination<std::pair<double, int>, dir, if_symmetry>(i,
                                                                                      b,
                                                                                      (dir
                                                                                       ? node->AllForwardBuckets[i][b].JumpArcs
                                                                                       : node->AllBackwardBuckets[i][b].JumpArcs),
                                                                                      dim_sq,
                                                                                      stateBetween2Buckets,
                                                                                      latest_bucket,
                                                                                      r1c_to_pi,
                                                                                      r1c_multi_to_pi);
      concatenateTestKernelInArcElimination<int, dir, if_symmetry>(i,
                                                                   b,
                                                                   (dir
                                                                    ? node->AllForwardBuckets[i][b].BucketArcs
                                                                    : node->AllBackwardBuckets[i][b].BucketArcs),
                                                                   dim_sq,
                                                                   stateBetween2Buckets,
                                                                   latest_bucket,
                                                                   r1c_to_pi,
                                                                   r1c_multi_to_pi);
    }
  }
#ifdef find_lost_arcs
  int i = 37, b = 10, j = 32, map = i * Dim + j;
  int latest = -1;
  auto &label_array = LabelArrayInForwardSense[i][b].first;
  auto &label_num = LabelArrayInForwardSense[i][b].second;
  for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
    auto *ki = label_array[vec_index_i];
    if (!increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
    tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
    //keep updating the latest_bucket
#ifdef SYMMETRY_PROHIBIT
    int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
            //first test
            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
            //real test
            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
#else
    int arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
    for (int bj = 0; bj <= arr_bj; ++bj) {
      //first test
      if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
      //real test
      auto &label_arr = LabelArrayInForwardSense[j][bj].first;
      auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
#endif
      for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
        auto &kj = label_arr[vec_index];
        if ((kj->PI & ki->PI).any()) continue;
#ifdef SYMMETRY_PROHIBIT
        if (tmp_mainResource > kj->Sum_MainResource) continue;
#else
        if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
#endif
        path_rc = tmp_rc + kj->RC;
        if (path_rc > OptGap) break;

        if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
          for (int l = 0; l < ki->numValidRank1Cut; ++l) {
            if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
          }
        } else {
          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
            if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
          }
        }

        if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
          for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
            int tmp_cut = ki->validRank1Cut_multi[l];
            if (kj->Rank1CutMem_multi[tmp_cut] +
                ki->Rank1CutMem_multi[tmp_cut]
                >= R1C_multi_denominator_InCG[tmp_cut]
                )
              path_rc -= r1c_multi_to_pi[tmp_cut];
          }
        } else {
          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
            int tmp_cut = kj->validRank1Cut_multi[l];
            if (ki->Rank1CutMem_multi[tmp_cut] +
                kj->Rank1CutMem_multi[tmp_cut]
                >= R1C_multi_denominator_InCG[tmp_cut]
                )
              path_rc -= r1c_multi_to_pi[tmp_cut];
          }
        }

        if (path_rc < OptGap) {
          std::cout << "arc should exists!" <<  std::endl;
          std::cout << "bin= " << bj <<  std::endl;
          latest = bj;
          break;
        }
      }
    }
  }
  if (latest != -1) {
    int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
    std::cout << "tell= " << TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq] <<  std::endl;
  }
#endif

  std::vector<int> tmp_vec;
  tmp_vec.reserve(Dim);
  int map1 = 0;
  for (int b = 0; b < NumBucketsPerVertex; ++b, map1 += dim_sq) {
    int map2 = map1 + Dim;
    for (int i = 1; i < Dim; ++i, map2 += Dim) {
      tmp_vec.clear();
      for (int j : dir ? node->AllForwardBuckets[i][b].BucketArcs : node->AllBackwardBuckets[i][b].BucketArcs) {
        if (stateBetween2Buckets[map2 + j]) {
          tmp_vec.emplace_back(j);
          ++num_bucket_arcs;
        }
      }
      (dir ? node->AllForwardBuckets[i][b].BucketArcs : node->AllBackwardBuckets[i][b].BucketArcs) = tmp_vec;
    }
  }

  eliminatebuketArc4Depot<dir, if_symmetry>(node);

  std::cout << "Num of " << (dir ? "Forward" : "Backward") << " BucketArcs= " << num_bucket_arcs << " prev.= "
            << double(num_bucket_arcs) / (dir ? node->NumForwardBucketArcs : node->NumBackwardBucketArcs) * 100 << "%"
            << " max.= " << double(num_bucket_arcs) / MaxNumForwardGraphArc * 100 << "%"
            << std::endl;
  (dir ? node->NumForwardBucketArcs : node->NumBackwardBucketArcs) = num_bucket_arcs;
}

//template<bool dir, bool if_symmetry>
//void CVRP::eliminatebuketArc4Depot(BBNODE *node) {
//  if constexpr (if_symmetry) {
//    if constexpr (dir) {
//      std::unordered_set<int> depot_set;
//      for (auto &arc : node->AllForwardBuckets[0][0].BucketArcs) depot_set.emplace(arc);
//      for (int i = 1; i < Dim; ++i) {
//        if (depot_set.find(i) == depot_set.end())continue;
//        bool if_delete = true;
//        for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
//          if (RC2TillThisBinInForwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap)break;
//          if (RC2BinInForwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap) {
//            continue;
//          } else {
//            auto &num = LabelArrayInForwardSense[i][b].second;
//            if (num) {
//              auto &labels = LabelArrayInForwardSense[i][b].first;
//              for (int j = 0; j < num; ++j) {
//                if (labels[j]->RC + ChgCostMat4Vertex[i][0] > OptGap)break;
//                double tmp_res = labels[j]->Sum_MainResource;
//                if (increaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) {
//                  if_delete = false;
//                  break;
//                }
//              }
//            }
//          }
//        }
//        if (if_delete) {
//          depot_set.erase(i);
//        }
//      }
//      if (depot_set.size() != node->AllForwardBuckets[0][0].BucketArcs.size()) {
//        node->AllForwardBuckets[0][0].BucketArcs.clear();
//        node->AllForwardBuckets[0][0].BucketArcs.assign(depot_set.begin(), depot_set.end());
//        std::sort(node->AllForwardBuckets[0][0].BucketArcs.begin(), node->AllForwardBuckets[0][0].BucketArcs.end());
//      }
//    } else {
//      std::unordered_set<int> depot_set;
//      for (auto &arc : node->AllBackwardBuckets[0][0].BucketArcs) depot_set.emplace(arc);
//      for (int i = 1; i < Dim; ++i) {
//        if (depot_set.find(i) == depot_set.end())continue;
//        bool if_delete = true;
//        for (int b = 0; b < NumBucketsPerVertex; ++b) {
//          if (RC2TillThisBinInBackwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap)break;
//          if (RC2BinInBackwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap) {
//            continue;
//          } else {
//            auto &num = LabelArrayInBackwardSense[i][b].second;
//            if (num) {
//              auto &labels = LabelArrayInBackwardSense[i][b].first;
//              for (int j = 0; j < num; ++j) {
//                if (labels[j]->RC + ChgCostMat4Vertex[i][0] > OptGap)break;
//                double tmp_res = labels[j]->Sum_MainResource;
//                if (decreaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) {
//                  if_delete = false;
//                  break;
//                }
//              }
//            }
//          }
//        }
//        if (if_delete) {
//          depot_set.erase(i);
//        }
//      }
//      if (depot_set.size() != node->AllBackwardBuckets[0][0].BucketArcs.size()) {
//        node->AllBackwardBuckets[0][0].BucketArcs.clear();
//        node->AllBackwardBuckets[0][0].BucketArcs.assign(depot_set.begin(), depot_set.end());
//        std::sort(node->AllBackwardBuckets[0][0].BucketArcs.begin(), node->AllBackwardBuckets[0][0].BucketArcs.end());
//      }
//    }
//  } else {
//    if constexpr (dir) {
//      std::unordered_set<int> depot_set;
//      for (auto &arc : node->AllForwardBuckets[0][0].BucketArcs) depot_set.emplace(arc);
//      for (int i = 1; i < Dim; ++i) {
//        if (depot_set.find(i) == depot_set.end())continue;
//        bool if_delete = true;
//        for (int b = 0; b < NumBucketsPerVertex; ++b) {
//          if (RC2TillThisBinInBackwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap)break;
//          if (RC2BinInBackwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap) {
//            continue;
//          } else {
//            auto &num = LabelArrayInBackwardSense[i][b].second;
//            if (num) {
//              auto &labels = LabelArrayInBackwardSense[i][b].first;
//              for (int j = 0; j < num; ++j) {
//                if (labels[j]->RC + ChgCostMat4Vertex[i][0] > OptGap)break;
//                double tmp_res = labels[j]->Sum_MainResource;
//                if (decreaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) {
//                  if_delete = false;
//                  break;
//                }
//              }
//            }
//          }
//        }
//        if (if_delete) {
//          depot_set.erase(i);
//        }
//      }
//      if (depot_set.size() != node->AllForwardBuckets[0][0].BucketArcs.size()) {
//        node->AllForwardBuckets[0][0].BucketArcs.clear();
//        node->AllForwardBuckets[0][0].BucketArcs.assign(depot_set.begin(), depot_set.end());
//        std::sort(node->AllForwardBuckets[0][0].BucketArcs.begin(), node->AllForwardBuckets[0][0].BucketArcs.end());
//      }
//    } else {
//      std::unordered_set<int> depot_set;
//      for (auto &arc : node->AllBackwardBuckets[0][0].BucketArcs) depot_set.emplace(arc);
//      for (int i = 1; i < Dim; ++i) {
//        if (depot_set.find(i) == depot_set.end())continue;
//        bool if_delete = true;
//        for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
//          if (RC2TillThisBinInForwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap)break;
//          if (RC2BinInForwardSense[i][b] + ChgCostMat4Vertex[i][0] > OptGap) {
//            continue;
//          } else {
//            auto &num = LabelArrayInForwardSense[i][b].second;
//            if (num) {
//              auto &labels = LabelArrayInForwardSense[i][b].first;
//              for (int j = 0; j < num; ++j) {
//                if (labels[j]->RC + ChgCostMat4Vertex[i][0] > OptGap)break;
//                double tmp_res = labels[j]->Sum_MainResource;
//                if (increaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) {
//                  if_delete = false;
//                  break;
//                }
//              }
//            }
//          }
//        }
//        if (if_delete) {
//          depot_set.erase(i);
//        }
//      }
//      if (depot_set.size() != node->AllBackwardBuckets[0][0].BucketArcs.size()) {
//        node->AllBackwardBuckets[0][0].BucketArcs.clear();
//        node->AllBackwardBuckets[0][0].BucketArcs.assign(depot_set.begin(), depot_set.end());
//        std::sort(node->AllBackwardBuckets[0][0].BucketArcs.begin(), node->AllBackwardBuckets[0][0].BucketArcs.end());
//      }
//    }
//  }
//}

template<bool dir, bool if_symmetry>
void CVRP::eliminatebuketArc4Depot(BBNODE *node) {
  bool constexpr if_dif = dir ^ if_symmetry;
  auto &allBuckets = dir ? node->AllForwardBuckets : node->AllBackwardBuckets;
  auto &RC2TillThisBin = if_dif ? RC2TillThisBinInBackwardSense : RC2TillThisBinInForwardSense;
  auto &RC2Bin = if_dif ? RC2BinInBackwardSense : RC2BinInForwardSense;
  auto &labelArray = if_dif ? LabelArrayInBackwardSense : LabelArrayInForwardSense;

  std::unordered_set<int> depot_set;
  for (auto &arc : allBuckets[0][0].BucketArcs) depot_set.emplace(arc);

  for (int i = 1; i < Dim; ++i) {
    if (depot_set.find(i) == depot_set.end()) continue;
    bool if_delete = true;

    for (int b = (if_dif ? 0 : (NumBucketsPerVertex - 1)); if_dif ? (b < NumBucketsPerVertex) : (b >= 0);
         if_dif ? ++b : --b) {
      if (RC2TillThisBin[i][b] + ChgCostMat4Vertex[i][0] > OptGap) break;
      if (RC2Bin[i][b] + ChgCostMat4Vertex[i][0] > OptGap) continue;

      auto &num = labelArray[i][b].second;
      if (num) {
        auto &labels = labelArray[i][b].first;
        for (int j = 0; j < num; ++j) {
          if (labels[j]->RC + ChgCostMat4Vertex[i][0] > OptGap) break;
          double tmp_res = labels[j]->Sum_MainResource;
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

  if (depot_set.size() != allBuckets[0][0].BucketArcs.size()) {
    allBuckets[0][0].BucketArcs.clear();
    allBuckets[0][0].BucketArcs.assign(depot_set.begin(), depot_set.end());
    std::sort(allBuckets[0][0].BucketArcs.begin(), allBuckets[0][0].BucketArcs.end());
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
    int map = i * Dim + j;
    int &latest = latest_bucket[map];
    int old_latest = latest;
    auto &label_array = dir ? LabelArrayInForwardSense[i][b].first : LabelArrayInBackwardSense[i][b].first;
    auto &label_num = dir ? LabelArrayInForwardSense[i][b].second : LabelArrayInBackwardSense[i][b].second;
    constexpr bool if_dif = dir ^ if_symmetry;
    for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
      auto *ki = label_array[vec_index_i];
      if constexpr (std::is_same<T, int>::value) {
        if (dir ? !increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j) :
            !decreaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j))
          continue;
      }
      tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
      //keep updating the latest_bucket
      if constexpr (dir) {
        arr_bj =
            if_symmetry ? std::min(int((MaxMainResource - tmp_mainResource) / StepSize), latest - 1) : std::max(int(
                (tmp_mainResource) / StepSize), latest + 1);
      } else {
        arr_bj =
            if_symmetry ? std::max(int((MaxMainResource - tmp_mainResource) / StepSize), latest + 1) : std::min(int(
                (tmp_mainResource) / StepSize), latest - 1);
      }

      for (int bj = if_dif ? NumBucketsPerVertex - 1 : 0; if_dif ? bj >= arr_bj : bj <= arr_bj;
           if_dif ? --bj : ++bj) {
        concatenateOneLabelWithOtherLabels<dir, if_symmetry, true, true>(ki,
                                                                         j,
                                                                         bj,
                                                                         tmp_rc,
                                                                         tmp_mainResource,
                                                                         r1c_to_pi,
                                                                         r1c_multi_to_pi,
                                                                         if_state);
//        if (if_state == -2) break;// since it was starting from 0
        if (if_state >= 0) {
          latest = if_state;
          goto outside;
        }
      }
      outside:
      continue;
    }
    //update the bins
    if (latest != old_latest) {
      int bi, chg_latest;
      if constexpr (if_symmetry) {
        chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
        if (chg_latest < 0) chg_latest = 0;
        else if (chg_latest >= NumBucketsPerVertex) chg_latest = NumBucketsPerVertex - 1;
      } else chg_latest = latest;
      bi = (dir ? TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq]
                : TellWhichBin4ArcEliminationInBackwardSense[map + chg_latest * dim_sq]);
      int map2 = map + b * dim_sq;
      for (int k = b; dir ? k <= bi : k >= bi; dir ? ++k : --k, dir ? map2 += dim_sq : map2 -= dim_sq)
        stateBetween2Buckets[map2] = true;
    }
  }
}

template<bool dir>
void CVRP::populateRC2TillThisBinNRC2Bin(BBNODE *const node) const {
  double rc, min_rc_per_bin;
  for (int i = 1; i < Dim; ++i) {
    rc = LARGEFLOAT;
    for (int b = dir ? 0 : NumBucketsPerVertex - 1; dir ? b < NumBucketsPerVertex : b >= 0; dir ? ++b : --b) {
      auto &label_list = dir ? LabelArrayInForwardSense[i][b].first : LabelArrayInBackwardSense[i][b].first;
      auto &num_labels = dir ? LabelArrayInForwardSense[i][b].second : LabelArrayInBackwardSense[i][b].second;
      if (num_labels) min_rc_per_bin = label_list[0]->RC;
      else min_rc_per_bin = LARGEFLOAT;
      if constexpr (dir) {
        RC2BinInForwardSense[i][b] = min_rc_per_bin;
        rc = std::min(rc, min_rc_per_bin);
        RC2TillThisBinInForwardSense[i][b] = rc;
      } else {
        RC2BinInBackwardSense[i][b] = min_rc_per_bin;
        rc = std::min(rc, min_rc_per_bin);
        RC2TillThisBinInBackwardSense[i][b] = rc;
      }
    }
  }
}
template<bool if_symmetry>
int CVRP::generateColsByBidir(BBNODE *node) {

  PtrAllR1Cs ptrr1cs(node, this);

  RE_TRY:
  Rollback = 0;
  runLabeling<true, false, if_symmetry>(node, ptrr1cs);

  if (Rollback == 2) {
    reallocateLabel();
    goto RE_TRY;
  } else if (Rollback == 1) {
    return 0;
  }

  if (!if_symmetry) {
    RE_TRY2:
    Rollback = 0;
    runLabeling<false, false, if_symmetry>(node, ptrr1cs);//ccnt can only be applied for forward case

    if (Rollback == 2) {
      reallocateLabel();
      goto RE_TRY2;
    } else if (Rollback == 1) {
      return 0;
    }
  }

//  int cnt = 0;
//  for (int i = 1; i < Dim; ++i) {
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      cnt += LabelArrayInForwardSense[i][b].second;
//    }
//  }
//
//  std::cout << "cnt= " << cnt << " all= " << IdxGlo << " ratio= " << double(cnt) / IdxGlo << std::endl;

  int ccnt = concatenateCols_prior_forward<if_symmetry>(node, ptrr1cs);
  //return -1

  if (abs(GapBetweenLastSmallestRCAndRCThreshold) < TOLERANCE) {
    double NumExistedLabels = 0;
    double NumExistedLabel_back = 0;
    for (int i = 1; i < Dim; ++i) {
      for (int b = 0; b < NumBucketsPerVertex; ++b) {
        NumExistedLabels += LabelArrayInForwardSense[i][b].second;
        if constexpr (!if_symmetry) {
          NumExistedLabel_back += LabelArrayInBackwardSense[i][b].second;
        }
      }
    }
    Ratio_DominanceChecks_NonDominant.first +=
        if_symmetry ? NumDominanceChecks / NumExistedLabels : NumDominanceChecks
            / (NumExistedLabels + NumExistedLabel_back);
    ++Ratio_DominanceChecks_NonDominant.second;
    if constexpr (!if_symmetry) {
      double dif = abs(NumExistedLabels - NumExistedLabel_back);
      double over = dif / std::min(NumExistedLabels, NumExistedLabel_back);
#ifdef DETAILED_EXACT_PRINT_INFO
      std::cout << "over= " << over << std::endl;
#endif
      if (over > CONFIG::NumberOfOverLabelsInMeetPoint) {
#ifdef DETAILED_EXACT_PRINT_INFO
        std::cout << "we adjust the meetpoint!" << std::endl;
#endif
        if (NumExistedLabels > NumExistedLabel_back) {
          MeetPointResourceInBiDir *= (1 - CONFIG::MeetPointFactor);
        } else {
          MeetPointResourceInBiDir *= (1 + CONFIG::MeetPointFactor);
        }
#ifdef DETAILED_EXACT_PRINT_INFO
        std::cout << "MeetPointResourceInBiDir= " << MeetPointResourceInBiDir << std::endl;
#endif
      }
    }
  }

  if (ccnt) {
    SmallestRC = std::get<2>(NegativeRCLabelTuple[0]);
    if (abs(GapBetweenLastSmallestRCAndRCThreshold) > TOLERANCE) {
      GapBetweenLastSmallestRCAndRCThreshold = std::get<2>(NegativeRCLabelTuple.back()) - SmallestRC;
      if (GapBetweenLastSmallestRCAndRCThreshold < TOLERANCE)GapBetweenLastSmallestRCAndRCThreshold = 0;
    } else {
      if (ceil_transformed_number_related(SmallestRC * K + LPVal + RC_TOLERANCE) + TOLERANCE >= UB) {
        node->Val = LPVal;
        node->if_terminated = true;
        std::cout << TERMINATED_MESSAGE_PROMISING_VEHICLES;
        return 0;
      }
    }
#ifdef DETAILED_EXACT_PRINT_INFO
    std::cout << "Smallest= " << SmallestRC << endl;
#endif
  } else {
    SmallestRC = 0;
    GapBetweenLastSmallestRCAndRCThreshold = 0;
  }
#ifdef DETAILED_EXACT_PRINT_INFO
  std::cout << "ccnt= " << ccnt << endl;
#endif
  //return 0
  if (!ccnt) return 0;

  addCols(node, ccnt);

  return ccnt;
}
template<bool if_symmetry>
int CVRP::concatenateCols_prior_forward(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs) {
  //update most_negative_rc
  //write in PoolBeg4Pricing
  PriorPoolBeg4Pricing = PoolBeg4Pricing;
  if (checkPricingPool()) reallocatePricingPool();
  int index = NumCol - 1;
  double tmp_rc, tmp_mainResource;
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  int *tmp_seq = new int[int(MaxMainResource) + 3];
  int size_tmp_seq;
  int i, j, arr_bj;
  LABEL *p;
  int if_state;//==0 means work

  for (auto &label_list : concatenateLabelsInForwardCG) {
    i = label_list.first.first;
    j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      tmp_mainResource = pr.second;
      arr_bj =
          if_symmetry ? int((MaxMainResource - tmp_mainResource) / StepSize) : (int) (tmp_mainResource / StepSize);
      tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];

      //arr_bj
      concatenateOneLabelWithOtherLabels<true, if_symmetry, true, false>(ki,
                                                                         j,
                                                                         arr_bj,
                                                                         tmp_rc,
                                                                         tmp_mainResource,
                                                                         r1c_to_pi,
                                                                         r1c_multi_to_pi,
                                                                         if_state);
      if (if_state == -2)continue;
      //arr_bj-1
      for (if_symmetry ? --arr_bj : ++arr_bj; if_symmetry ? arr_bj >= 0 : arr_bj < NumBucketsPerVertex;
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

  writeColsInPricingPool(node, index);

  delete[]tmp_seq;
  return index - NumCol + 1;
}

template<bool dir>
void CVRP::populateTellWhichBin4ArcElimination() {
  /**
   * one time calculation except for regenerate the bucket graph
   */
  if constexpr (dir) {
    if (!TellWhichBin4ArcEliminationInForwardSense.empty()) return;
  } else {
    if (!TellWhichBin4ArcEliminationInBackwardSense.empty()) return;
  }
  std::cout << "populateTellWhichBin4ArcElimination" << std::endl;
  size_t size = Dim * Dim * NumBucketsPerVertex;
  int dim_sq = Dim * Dim;
  double tmp_mainResource;
  if constexpr (dir) TellWhichBin4ArcEliminationInForwardSense.reserve(size);
  else TellWhichBin4ArcEliminationInBackwardSense.reserve(size);
  //calculate according the bucket Graph, no! since every node is different!
  std::unordered_map<int, int> map_bj_B;
  std::vector<std::pair<int, int>> vec_bj_B(NumBucketsPerVertex);
  map_bj_B.reserve(NumBucketsPerVertex + 1);
  for (int b = 0; b <= NumBucketsPerVertex; ++b) map_bj_B[b] = dir ? -1 : NumBucketsPerVertex;
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
      int base = i * Dim + j;
      for (auto &key_value : map_bj_B) key_value.second = dir ? -1 : NumBucketsPerVertex;
      for (int b = dir ? 0 : NumBucketsPerVertex - 1; dir ? b < NumBucketsPerVertex : b >= 0; dir ? ++b : --b) {
        double res = dir ? b * StepSize : std::min((b + 1) * StepSize - TOLERANCE, MaxMainResource);
        if (dir ? !increaseMainResourceConsumption(res, tmp_mainResource, i, j) :
            !decreaseMainResourceConsumption(res, tmp_mainResource, i, j))
          break;
        int bj = int(tmp_mainResource / StepSize);
        if (dir ? map_bj_B[bj] < b : map_bj_B[bj] > b) map_bj_B[bj] = b;
      }
      for (int bj = dir ? NumBucketsPerVertex - 1 : 0; dir ? bj >= 0 : bj < NumBucketsPerVertex; dir ? --bj : ++bj) {
        if (map_bj_B[bj] != (dir ? -1 : NumBucketsPerVertex)) {
          for (int bj2 = dir ? bj + 1 : bj - 1; dir ? bj2 < NumBucketsPerVertex : bj2 >= 0; dir ? ++bj2 : --bj2)
            map_bj_B[bj2] = map_bj_B[bj];
          break;
        }
      }
      for (int bj = 0; bj < NumBucketsPerVertex; ++bj) {
        (dir ? TellWhichBin4ArcEliminationInForwardSense[base + bj * dim_sq] :
         TellWhichBin4ArcEliminationInBackwardSense[base + bj * dim_sq]) = map_bj_B[bj];
        //should be very careful about expression like this
      }
    }
  }
}

template<bool dir>
void CVRP::obtainjumpArcs(BBNODE *node, std::bitset<2> **bitMap) const {
  int num_jump_arcs = 0;
  bool if_used;

  for (int i = 1; i < Dim; ++i) {
    //start to collect data
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      if constexpr (dir) {
        node->AllForwardBuckets[i][b].JumpArcs.clear();
      } else {
        node->AllBackwardBuckets[i][b].JumpArcs.clear();
      }
      //2 represents the bin is not checked
      //1 represents the bin is checked
      //0 represents the bin has arc itself
      for (int j = 1; j < Dim; ++j) bitMap[j][b] = 2;
      bitMap[i][b] = 1;
      for (int j : (dir ? node->AllForwardBuckets[i][b].BucketArcs :
                    node->AllBackwardBuckets[i][b].BucketArcs))
        bitMap[j][b] = 0;
    }
    for (int b = (dir ? 0 : NumBucketsPerVertex - 1); dir ? b < NumBucketsPerVertex : b >= 0; dir ? ++b : --b) {
      for (int j = 1; j < Dim; ++j) {
        if (bitMap[j][b] == 2) {//需要有jump arc
          //开始找最近的step个有去j这个地方的bin
          if_used = false;
          for (int b4_i = (dir ? b + 1 : b - 1); dir ? b4_i < NumBucketsPerVertex : b4_i >= 0;
               dir ? ++b4_i : --b4_i) {
            if (bitMap[j][b4_i] == 0) {//发现有去的arc了
              std::pair<int, int> map = (dir ? std::make_pair(b4_i * StepSize, j) :
                                         std::make_pair((b4_i + 1) * StepSize - TOLERANCE, j));
              for (int tmp_b = b; dir ? tmp_b < b4_i : tmp_b > b4_i; dir ? ++tmp_b : --tmp_b) {//不取等号< qi
                bitMap[j][tmp_b] = 1;
                if constexpr (dir) {
                  node->AllForwardBuckets[i][tmp_b].JumpArcs.emplace_back(map);
                } else {
                  node->AllBackwardBuckets[i][tmp_b].JumpArcs.emplace_back(map);
                }
                ++num_jump_arcs;
              }
              if_used = true;
              break;
            }
          }
          if (!if_used) {
            for (int tmp_b = b; dir ? tmp_b < NumBucketsPerVertex : tmp_b >= 0; dir ? ++tmp_b : --tmp_b)
              bitMap[j][tmp_b] = 1;
          }
        }
      }
    }
  }

  (dir ? node->NumForwardJumpArcs : node->NumBackwardJumpArcs) = num_jump_arcs;

  std::cout << "Obtain" << (dir ? "Forward" : "Backward") << " Jump Arcs= " << num_jump_arcs << std::endl;
}

template<bool dir, bool if_symmetry>
int CVRP::enumerateHalfwardRoutes(BBNODE *node,
                                  const double *r1c_to_pi,
                                  const double *r1c_multi_to_pi,
                                  std::unordered_map<yzzLong, std::tuple<LABEL *, LABEL *, double>> &Tags,
                                  std::vector<LABEL *> **copy_bucket,
                                  int &num_routes_now) {
  //no sort in this case
  int status = 0;
  int edgemap;
  (dir ? NumForwardLabelsInEnu : NumBackwardLabelsInEnu) = 0;
  int Max_routes_phase1 = CONFIG::MaxNumRouteInEnumeration_half;
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
  double left_time = CONFIG::HardTimeThresholdInAllEnumeration;

  for (int b = dir ? 0 : NumBucketsPerVertex - 1; dir ? b < NumBucketsPerVertex : b >= 0; dir ? ++b : --b) {
    int i = 1;
    STILL_EXIST:
    for (; i < Dim; ++i) {
      end = std::chrono::high_resolution_clock::now();
      eps = std::chrono::duration<double>(end - b4_end).count();
      if (eps > left_time) {
        status = 2;
        goto outside;
      }
      auto &valid_num = (dir ? IfExistExtraLabelsInForwardSense[i][b].second :
                         IfExistExtraLabelsInBackwardSense[i][b].second);
      if (!valid_num) continue;
      auto &label_array = (dir ? IfExistExtraLabelsInForwardSense[i][b].first :
                           IfExistExtraLabelsInBackwardSense[i][b].first);
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->if_extended) continue;
        ki->if_extended = true;
        for (int j : (dir ? node->AllForwardBuckets[i][b].BucketArcs :
                      node->AllBackwardBuckets[i][b].BucketArcs)) {
          if (ki->PI[j]) continue;
          auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
          if (dir ? !increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j) :
              !decreaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j))
            continue;
          auto &tmp_rc = AllLabel[IdxGlo].RC;
          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
          if_keep = false;
          int arr_bj = (if_symmetry ? int((MaxMainResource - tmp_mainResource) / StepSize) :
                        int(tmp_mainResource / StepSize));
          if constexpr (dir) {
            if (tmp_mainResource > MeetPointResourceInBiDirEnu) {
              if constexpr (if_symmetry) {
                if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc < OptGap)
                  concatenateLabelsInForwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
              } else {
                if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc < OptGap)
                  concatenateLabelsInForwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
              }
              continue;
            }
          } else {
            if (tmp_mainResource < MeetPointResourceInBiDirEnu) {
              if constexpr (if_symmetry) {
                if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc < OptGap)
                  concatenateLabelsInBackwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
              } else {
                if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc < OptGap)
                  concatenateLabelsInBackwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
              }
              continue;
            }
          }
          if (tmp_rc + (if_dif ? RC2TillThisBinInBackwardSense[j][arr_bj] : RC2TillThisBinInForwardSense[j][arr_bj])
              > OptGap)
            continue;
          if (tmp_rc + (if_dif ? RC2BinInBackwardSense[j][arr_bj] : RC2BinInForwardSense[j][arr_bj])
              < OptGap) {
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if constexpr (if_symmetry) {
                if (tmp_mainResource + kkj->Sum_MainResource > MaxMainResource) continue;
              } else {
                if constexpr (dir) {
                  if (tmp_mainResource > kkj->Sum_MainResource) continue;
                } else {
                  if (tmp_mainResource < kkj->Sum_MainResource) continue;
                }
              }
              if ((ki->PI & kkj->PI).any()) continue;
              path_rc = tmp_rc + kkj->RC;
              if (path_rc > OptGap) break;
              if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
                  if (kkj->Rank1CutMem[ki->validRank1Cut[l]]) {
                    path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
                    if (path_rc > OptGap)goto here;
                  }
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
                  if (ki->Rank1CutMem[kkj->validRank1Cut[l]]) {
                    path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
                    if (path_rc > OptGap)goto here;
                  }
                }
              }

              if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = ki->validRank1Cut_multi[l];
                  if (kkj->Rank1CutMem_multi[tmp_cut] +
                      ki->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > OptGap)goto here;
                  }
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = kkj->validRank1Cut_multi[l];
                  if (ki->Rank1CutMem_multi[tmp_cut] +
                      kkj->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > OptGap)goto here;
                  }
                }
              }
              if_keep = true;
              goto outside1;
              here:;
            }
          }
          for (if_dif ? ++arr_bj : --arr_bj; if_dif ? arr_bj < NumBucketsPerVertex : arr_bj >= 0;
               if_dif ? ++arr_bj : --arr_bj) {
            //first test
            if (tmp_rc + (if_dif ? RC2TillThisBinInBackwardSense[j][arr_bj] : RC2TillThisBinInForwardSense[j][arr_bj])
                > OptGap)
              break;
            //second test
            if (tmp_rc + (if_dif ? RC2BinInBackwardSense[j][arr_bj] : RC2BinInForwardSense[j][arr_bj])
                > OptGap)
              continue;
            //real test
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if ((ki->PI & kkj->PI).any()) continue;
              path_rc = tmp_rc + kkj->RC;
              if (path_rc > OptGap) break;
              if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
                  if (kkj->Rank1CutMem[ki->validRank1Cut[l]]) {
                    path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
                    if (path_rc > OptGap)goto here1;
                  }
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
                  if (ki->Rank1CutMem[kkj->validRank1Cut[l]]) {
                    path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
                    if (path_rc > OptGap)goto here1;
                  }
                }
              }

              if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = ki->validRank1Cut_multi[l];
                  if (kkj->Rank1CutMem_multi[tmp_cut] +
                      ki->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > OptGap)goto here1;
                  }
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = kkj->validRank1Cut_multi[l];
                  if (ki->Rank1CutMem_multi[tmp_cut] +
                      kkj->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      ) {
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                    if (path_rc > OptGap)goto here1;
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
          int bj = int(tmp_mainResource / StepSize);
          auto &labelList_j = dir ? LabelArrayInForwardSense[j][bj].first : LabelArrayInBackwardSense[j][bj].first;
          auto &valid_num_j = dir ? LabelArrayInForwardSense[j][bj].second : LabelArrayInBackwardSense[j][bj].second;
          auto &tmp_PI = AllLabel[IdxGlo].PI;
          auto &tmp_Cost = AllLabel[IdxGlo].Cost;
          auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
          auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
          auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
          tmp_PI = ki->PI;
          tmp_PI.set(j);
          tmp_Cost = ki->Cost + CostMat4Vertex[i][j];
          tmp_Rank1CutMem = ki->Rank1CutMem;
          //tmp_num_valid_rank1_cut do not have to copy
          for (auto l : std::get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_Rank1CutMem[l] = false;
              tmp_rc -= r1c_to_pi[l];
            } else tmp_Rank1CutMem[l] = true;
          }
          tmp_Rank1CutMem &= std::get<1>(Vertex2ActiveInOnePricingR1Cs[j]);

          auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
          auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
          auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;
          std::copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);
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
            if (abs(kj->Sum_MainResource - tmp_mainResource) > TOLERANCE) {
              ++vec_index_j;
              continue;
            }
            if ((kj->PI ^ tmp_PI).none()) {
              if (kj->Cost > tmp_Cost) {
                kj->if_extended = true;
                kj = labelList_j[--valid_num_j];
                dir ? --NumForwardLabelsInEnu : --NumBackwardLabelsInEnu;
              } else {
                if_break = true;
                break;
              }
            } else ++vec_index_j;
#else
            if (kj->Cost > tmp_Cost) {
              if (dir ? kj->Sum_MainResource > tmp_mainResource : kj->Sum_MainResource < tmp_mainResource) {
                if ((kj->PI ^ tmp_PI).none()) {
                  kj->if_extended = true;
                  kj = labelList_j[--valid_num_j];
                  dir ? --NumForwardLabelsInEnu : --NumBackwardLabelsInEnu;
                } else ++vec_index_j;
              } else ++vec_index_j;
            } else {
              if (dir ? kj->Sum_MainResource < tmp_mainResource : kj->Sum_MainResource > tmp_mainResource) {
                if ((kj->PI ^ tmp_PI).none()) {
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

          labelList_j[valid_num_j++] = AllLabel + IdxGlo;
          if (valid_num_j == labelList_j.size()) {
            labelList_j.resize(labelList_j.size() * 2);
          }

          //Seq
          AllLabel[IdxGlo].Seq = AllSeq + SeqBeg;
          std::copy(ki->Seq, ki->Seq + ki->IdxEndSeq + 1, AllLabel[IdxGlo].Seq);
          AllLabel[IdxGlo].IdxEndSeq = ki->IdxEndSeq + 1;
          *(AllLabel[IdxGlo].Seq + AllLabel[IdxGlo].IdxEndSeq) = j;
          SeqBeg += AllLabel[IdxGlo].IdxEndSeq + 1;
          //EndVertex
          AllLabel[IdxGlo].EndVertex = j;
          //if_extended
          AllLabel[IdxGlo].if_extended = false;
          //bucket
          auto &bucket = (dir ? IfExistExtraLabelsInForwardSense[j][bj] : IfExistExtraLabelsInBackwardSense[j][bj]);
          bucket.first[bucket.second++] = AllLabel + IdxGlo;
          if (bucket.second == bucket.first.size()) {
            bucket.first.resize(bucket.first.size() * 2);
          }

          if constexpr (dir) {
            if (tmp_rc + ChgCostMat4Vertex[j][0] < OptGap) {
              path_cost = tmp_Cost + CostMat4Vertex[j][0];

              if (Tags.find(tmp_PI) == Tags.end()) {
                Tags[tmp_PI] = {AllLabel + IdxGlo, nullptr, path_cost};
                ++num_routes_now;
                if (num_routes_now > Max_routes_phase1) {
                  status = 2;//routes limit
                  goto outside;
                }
              } else if (std::get<2>(Tags[tmp_PI]) > path_cost) {
                Tags[tmp_PI] = {AllLabel + IdxGlo, nullptr, path_cost};
              }
            }
          }

          if ((dir ? ++NumForwardLabelsInEnu : ++NumBackwardLabelsInEnu) > CONFIG::MaxNumLabelInEnumeration) {
            status = 3;//all labels limit
            goto outside;
          }
          ++IdxGlo;//can be put here, because once go to outside, the function will end
          if (IdxGlo == LabelAssign) {
            Rollback = 2;
            goto outside;
          }
        }
      }
      valid_num = 0;
    }
//test if all labels are extended
    for (i = 1; i < Dim; ++i) {
      if (dir ? IfExistExtraLabelsInForwardSense[i][b].second : IfExistExtraLabelsInBackwardSense[i][b].second)
        goto STILL_EXIST;
    }
    af_end = std::chrono::high_resolution_clock::now();
    eps2 = std::chrono::duration<double>(af_end - b4_end).count();
    eps = std::chrono::duration<double>(af_end - beg).count();
    left_time = (CONFIG::HardTimeThresholdInAllEnumeration - eps) / (dir ? (NumBucketsPerVertex - b) : (b + 1));
    if (eps2 > left_time) {
      status = 2;
      goto outside;
    }
    b4_end = af_end;
  }
  outside:
  for (int i = 0; i < Dim; ++i) {
    delete[]copy_bucket[i];
  }
  delete[] copy_bucket;
  std::cout << "Half" << (dir ? "Forward" : "Backward") << " labeling: num_labels= "
            << (dir ? NumForwardLabelsInEnu : NumBackwardLabelsInEnu)
            << " num_routes= " << num_routes_now <<
            std::endl;
  if (status)return status;
//there is one situation that rc of one col could be smaller than optGap
//but could be dominated by cols whose rc is larger than optGap
//this could lead to a discrepancy between the number of routes

//for reasons above, we revise the number of routes and sort the routes
  for (int i = 1; i < Dim; ++i) {
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      auto &label_array = dir ? LabelArrayInForwardSense[i][b] : LabelArrayInBackwardSense[i][b];
      std::sort(label_array.first.begin(), label_array.first.begin() + label_array.second, CmpLabelRCLess);
    }
  }
//  populateRC2TillThisBinNRC2Bin(node, 1);//never use this function in phase1
  return 0;
}

template<typename T, bool if_enu>
void CVRP::addSelectedR1C_N_multiCuts(BBNODE *node,
                                      std::vector<std::tuple<int, std::set<int>, double, int, int>> &cut_info_set,
    //cut_idx, mem, vio, lp_index(if max, then new cut!)
                                      T &cut_info) {
  int num_cuts_can_be_added = 0;
  if constexpr (std::is_same<T, std::vector<std::pair<std::vector<int>, double>>>::value) {
    num_cuts_can_be_added = MaxNum_R1Cs - (int) node->R1Cs.size() - 1;
  } else {
    num_cuts_can_be_added = MaxNum_R1C_multi - (int) node->R1Cs_multi.size() - 1;
  }
  if (num_cuts_can_be_added <= 0) {
    return;
  }
  //revise vio
  std::vector<double> vio_vec(cut_info_set.size());
  if constexpr (std::is_same<T, std::vector<std::pair<std::vector<int>, double>>>::value) {
    fill_n(vio_vec.begin(), cut_info_set.size(), 0.5);
  } else {
    int cnt = 0;
    for (auto &i : cut_info_set) {
      auto &cut = cut_info[std::get<0>(i)];
      auto deno = std::get<1>(map_rank1_multiplier[(int) std::get<0>(cut).size()][std::get<1>(cut)]);
      vio_vec[cnt++] = 1. - 1. / deno;
    }
  }

  //transform (vio_vec[i] - get<2>(cut_info_set[i])
  transform(vio_vec.begin(), vio_vec.end(), cut_info_set.begin(), cut_info_set.begin(),
            [](double &a, std::tuple<int, std::set<int>, double, int, int> &b) {
              std::get<2>(b) = a - std::get<2>(b);
              return b;
            });

//sort by vio first and then by sie of mem, when sort by vio, the tolerance is 1e-6
  sort(cut_info_set.begin(), cut_info_set.end(), [](const std::tuple<int, std::set<int>, double, int, int> &a,
                                                    const std::tuple<int, std::set<int>, double, int, int> &b) {
    if (abs(std::get<2>(a) - std::get<2>(b)) < 1e-6) {
      return std::get<4>(a) < std::get<4>(b);
    } else return std::get<2>(a) < std::get<2>(b);
  });

  std::unordered_map<int, double> min_gap;
  if constexpr (!if_enu) {
    //for each mem size, record the max vio
    for (auto &cut : cut_info_set) {
      int mem_size = std::get<4>(cut);
      if (min_gap.find(mem_size) == min_gap.end()) {
        min_gap[mem_size] = std::get<2>(cut);
      } else {
        min_gap[mem_size] = std::min(min_gap[mem_size], std::get<2>(cut));
      }
    }
    //subtract by 1e-6 for each element
    transform(cut_info_set.begin(), cut_info_set.end(), cut_info_set.begin(),
              [](std::tuple<int, std::set<int>, double, int, int> &a) {
                std::get<2>(a) += 1e-6;
                return a;
              });
  }


  //we check if the cut's mem is greater than corresponding max vio,
  //it should be greater than the max vio, when faced with less mem size
  std::unordered_map<int, int> num_r1c;//size & num
  std::unordered_map<int, int> r1c_limit;//size & limit
  for (int i = 1; i <= CONFIG::MaxRowRank1; ++i) r1c_limit[i] = CONFIG::MaxNumR1CPerRoundSecondSelect;
  r1c_limit[3] = CONFIG::MaxNumR1C3PerRoundSecondSelect;
#ifdef NewWay2AddCuts
  std::cout << "No selection!" << std::endl;
#endif
  for (auto &cut : cut_info_set) {
    int mem_size;
    bool if_should_add = true;
#ifndef NewWay2AddCuts
    if constexpr (!if_enu) {
      mem_size = std::get<4>(cut);
      for (int i = mem_size - 1; i >= 0; --i) {
        if (min_gap.find(i) != min_gap.end() && std::get<2>(cut) > min_gap[i]) {
          if_should_add = false;
          break;
        }
      }
    }
#endif
    if (if_should_add) {
      //if cutinfo is vector<pair<vector<int>, double>>, we use addR1C
      //else we use addR1C_multi
      int size = (int) std::get<0>(cut_info[std::get<0>(cut)]).size();
      if (++num_r1c[size] > r1c_limit[size]) {
        continue;
      }
      if constexpr (std::is_same<T, std::vector<std::pair<std::vector<int>, double>>>::value) {
        if constexpr (if_enu) {
          addR1CInEnu(node, cut_info[std::get<0>(cut)].first);
        } else
          addR1C(node, cut_info[std::get<0>(cut)].first, std::get<1>(cut), std::get<3>(cut));
      } else {
        if constexpr (if_enu) {
          addR1C_multiInEnu(node,
                            make_pair(std::get<0>(cut_info[std::get<0>(cut)]),
                                      std::get<1>(cut_info[std::get<0>(cut)])));
        } else
          addR1C_multi(node,
                       make_pair(std::get<0>(cut_info[std::get<0>(cut)]), std::get<1>(cut_info[std::get<0>(cut)])),
                       std::get<1>(cut),
                       std::get<3>(cut));
      }
      --num_cuts_can_be_added;
      if (num_cuts_can_be_added <= 0) break;
    }
  }
}
