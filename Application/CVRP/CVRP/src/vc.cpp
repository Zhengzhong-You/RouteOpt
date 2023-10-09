////
//// Created by Zhengzhong You on 4/5/23.
////
//#include "CVRP.hpp"
//#include "templateFunctors.hpp"
//
//using namespace std;
//using namespace Eigen;
//using namespace chrono;
//
//void CVRP::runHalfForwardLabeling(BBNODE *const node, const PtrAllR1Cs &ptrAllR1Cs) {
//  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
//  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
//  bool if_suc;
//
//  initializeLabels(node, 1, true, {true, 1, true});
//
//  auto beg = high_resolution_clock::now();
//  auto end = beg;
//  double eps;
//  int min_sorted_b = 0;
//  for (int b = 0; b < NumBucketsPerVertex; ++b) {
//    int i = 1;
//    STILL_EXIST:
//    for (; i < Dim; ++i) {
//      auto &valid_num = IfExistExtraLabelsInForwardSense[i][b].second;
//      if (!valid_num) continue;
//      auto &label_array = IfExistExtraLabelsInForwardSense[i][b].first;
//      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
//        auto &ki = label_array[vec_index];
//        if (ki->if_extended) continue;
//        checkIfDominated<true>(ki, i, b, r1c_to_pi, r1c_multi_to_pi, if_suc);
//        ki->if_extended = true;
//        if (!if_suc)continue;
//        auto sig = extendKernel4Exact<int, true, false, true>(ki,
//                                                              i,
//                                                              ki->Sum_MainResource,
//                                                              node->AllForwardBuckets[i][b].BucketArcs,
//                                                              r1c_to_pi,
//                                                              r1c_multi_to_pi);
//        if (sig == 3) goto populateBin;
//        sig = extendKernel4Exact<pair<double, int>, true, false, true>(ki,
//                                                                       i,
//                                                                       0,
//                                                                       node->AllForwardBuckets[i][b].JumpArcs,
//                                                                       r1c_to_pi,
//                                                                       r1c_multi_to_pi);
//        if (sig == 3) goto populateBin;
//      }
//      valid_num = 0;
//    }
//    //test if all labels are extended
//    end = high_resolution_clock::now();
//    eps = duration<double>(end - beg).count();
//    if (eps > CutGenTimeThresholdInPricing) {
//      Rollback = 3;
//      if (eps > CONFIG::HardTimeThresholdInPricing) {
//        if (!ForceNotRollback) {
//          Rollback = 1;
//          goto QUIT;
//        }
//      }
//    }
//    checkIfNoLabelsLeft<true>(i, b, min_sorted_b, if_suc);
//    if (!if_suc) goto STILL_EXIST;
//  }
//
//  populateBin:
//  if (min_sorted_b < NumBucketsPerVertex) {
//    for (int i = 1; i < Dim; ++i) {
//      for (int b = min_sorted_b; b < NumBucketsPerVertex; ++b) {
//        std::stable_sort(LabelArrayInForwardSense[i][b].first.begin(),
//                         LabelArrayInForwardSense[i][b].first.begin() + LabelArrayInForwardSense[i][b].second,
//                         CmpLabelRCLess);
//      }
//    }
//  }
//  //changing to arc architecture does not affect the following code
//  populateRC2TillThisBinNRC2Bin(node);
//  end = high_resolution_clock::now();
//  eps = duration<double>(end - beg).count();
//  LastMaxTimeLabeling = max(LastMaxTimeLabeling, eps);
//  if (Rollback != 3) {
//    if (eps * HardRollBackFactor > CONFIG::HardTimeThresholdInPricing) {
//      Rollback = 3;
//    }
//  }
//  QUIT:
//#ifdef DETAILED_EXACT_PRINT_INFO
//  cout << "这里是测试！" << endl;
//  {
//    if (GapBetweenLastSmallestRCAndRCThreshold < TOLERANCE) {
//      cout << "we find if there exists extra labels that can be dominated but saved for some reasons!" << endl;
//      //we here generate the cuts
//      //rank1 cuts
//      safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//      int all_num = int(node->R1Cs.size() + node->R1Cs_multi.size());
//      vector<int> denominator(all_num);
//      vector<double> pi(all_num);
//      vector<vector<pair<int, int>>> cuts(Dim);//cut_idx, cut_augment
//      vector<unordered_set<int>> mem(Dim);
//      vector<vector<int>> no_mem(Dim);
//      int cnt = 0;
//      for (auto &r1c : node->R1Cs) {
//        for (auto &i : r1c.InfoR1C) {
//          cuts[i].emplace_back(cnt, 1);
//          mem[i].emplace(cnt);
//        }
//        for (auto &i : r1c.Mem) {
//          mem[i].emplace(cnt);
//        }
//        denominator[cnt] = 2;
//        pi[cnt] = Pi[r1c.IdxR1C];
//        ++cnt;
//      }
//
//      for (auto &r1c : node->R1Cs_multi) {
//        auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
//        auto &multi = get<0>(plan);
//        int deno = get<1>(plan);
//        int count = 0;
//        for (auto &i : r1c.InfoR1C.first) {
//          cuts[i].emplace_back(cnt, multi[count]);
//          mem[i].emplace(cnt);
//          ++count;
//        }
//        for (auto &v : r1c.Mem) {
//          mem[v].emplace(cnt);
//        }
//        denominator[cnt] = deno;
//        pi[cnt] = Pi[r1c.IdxR1C];
//        ++cnt;
//      }
//
//      for (int i = 1; i < Dim; ++i) {
//        for (int j = 0; j < all_num; ++j) {
//          if (mem[i].find(j) == mem[i].end()) {
//            no_mem[i].emplace_back(j);
//          }
//        }
//      }
//
//      for (int b = 0; b < NumBucketsPerVertex; ++b) {
//        for (int i = 1; i < Dim; ++i) {
//          auto &label_list = LabelArrayInForwardSense[i][b].first;
//          auto &valid_num = LabelArrayInForwardSense[i][b].second;
//          int label_count = 0;
//          HERE:
//          if (valid_num <= label_count) continue;
//          auto ki = label_list[label_count++];
//          auto ng = ki->PI;
//          for (int j = label_count; j < valid_num; ++j) {
//            auto kj = label_list[j];
//            if (ki->Sum_MainResource < kj->Sum_MainResource) {
////              cout << "pass sum main resource!" << endl;
//              if (((kj->PI & ng) ^ (ng)).none()) {
////                cout << "now this label pass ng test and rc!" << endl;
////                cout << "now we find out why this label could not be dominate!" << endl;
//                double rc_dif = kj->RC - ki->RC;
//                //if rank1_dif can bring more rc_dif, then this label indeed cannot be dominated!
//                if (rc_dif > TOLERANCE) {
////                  cout << "we find out why this label could not be dominate!" << endl;
//                  //find the sequence of ki and kj
//                  vector<int> seq_ki;
//                  vector<int> seq_kj;
//                  auto tmp_ki = ki;
//                  while (tmp_ki->PLabel) {
//                    seq_ki.emplace_back(tmp_ki->EndVertex);
//                    tmp_ki = tmp_ki->PLabel;
//                  }
//                  std::reverse(seq_ki.begin(), seq_ki.end());
//                  auto tmp_kj = kj;
//                  while (tmp_kj->PLabel) {
//                    seq_kj.emplace_back(tmp_kj->EndVertex);
//                    tmp_kj = tmp_kj->PLabel;
//                  }
//                  std::reverse(seq_kj.begin(), seq_kj.end());
//                  vector<int> state_ki(all_num, 0);
//                  vector<int> state_kj(all_num, 0);
//                  for (auto &v : seq_ki) {
//                    if (v == 0) continue;
//                    for (auto &cut : cuts[v]) {
//                      state_ki[cut.first] += cut.second;
//                      if (state_ki[cut.first] >= denominator[cut.first]) {
//                        state_ki[cut.first] -= denominator[cut.first];
//                      }
//                    }
//                    for (auto &cut : no_mem[v]) {
//                      state_ki[cut] = 0;
//                    }
//                  }
//                  for (auto &v : seq_kj) {
//                    if (v == 0) continue;
//                    for (auto &cut : cuts[v]) {
//                      state_kj[cut.first] += cut.second;
//                      if (state_kj[cut.first] >= denominator[cut.first]) {
//                        state_kj[cut.first] -= denominator[cut.first];
//                      }
//                    }
//                    for (auto &cut : no_mem[v]) {
//                      state_kj[cut] = 0;
//                    }
//                  }
//
//                  double diff = 0;
//                  for (int k = 0; k < all_num; ++k) {
//                    if (state_ki[k] > state_kj[k]) diff -= pi[k];
//                  }
//                  cout << "diff: " << diff << " rc_dif: " << rc_dif << endl;
//                  auto min_it = min_element(pi.begin(), pi.end());
//                  vector<pair<double, int>> pi_copy;
//                  for (int k = 0; k < pi.size(); ++k) {
//                    if (k >= node->R1Cs.size()) {
//                      int num = k - node->R1Cs.size();
//                      pi_copy.emplace_back(pi[k], node->R1Cs_multi[num].InfoR1C.first.size());
//                    } else {
//                      pi_copy.emplace_back(pi[k], node->R1Cs[k].InfoR1C.size());
//                    }
//                  }
//                  sort(pi_copy.begin(), pi_copy.end(), [](const pair<double, bool> &a, const pair<double, bool> &b) {
//                    return a.first < b.first;
//                  });
//                  for (int k = 0; k < pi_copy.size(); ++k) {
//                    cout << "(" << pi_copy[k].first << ", " << pi_copy[k].second << ") ";
//                  }
//                  cout << endl;
//                  cout << "min pi: " << *min_it << endl;
//                  int cut = min_it - pi.begin();
//                  if (cut > node->R1Cs.size()) {
//                    cut -= node->R1Cs.size();
//                    for (int i = 0; i < cut; ++i) {
//                      auto &r1c = node->R1Cs_multi[i];
//                      auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
//                      auto &multi = get<0>(plan);
//                      int deno = get<1>(plan);
//                      int count = 0;
//                      // print info
//                      for (auto &i : r1c.InfoR1C.first) {
//                        cout << i << " ";
//                      }
//                      cout << endl;
//                      cout << "multi: ";
//                      for (auto &i : multi) {
//                        cout << i << " ";
//                      }
//                      cout << endl;
//                      cout << "deno: " << deno << endl;
//                    }
//                  } else {
//                    auto &r1c = node->R1Cs[cut];
//                    // print info
//                    for (auto &i : r1c.InfoR1C) {
//                      cout << i << " ";
//                    }
//                    cout << endl;
//                  }
//
//                  if (diff < rc_dif) {
//                    cout << "diff: " << diff << " rc_dif: " << rc_dif << endl;
//                    cout << "seq_ki: ";
//                    for (auto &v : seq_ki) {
//                      cout << v << " ";
//                    }
//                    cout << endl;
//                    cout << "seq_kj: ";
//                    for (auto &v : seq_kj) {
//                      cout << v << " ";
//                    }
//                    cout << endl;
//                    cout << "state_ki: ";
//                    for (auto &v : state_ki) {
//                      cout << v << " ";
//                    }
//                    cout << endl;
//                    cout << "state_kj: ";
//                    for (auto &v : state_kj) {
//                      cout << v << " ";
//                    }
//                    cout << endl;
//                    cout << "pi: ";
//                    for (auto &v : pi) {
//                      cout << v << " ";
//                    }
//                    cout << endl;
//                    cout << "denominator: ";
//                    for (auto &v : denominator) {
//                      cout << v << " ";
//                    }
//                    cout << endl;
//                    cout << "cuts: " << endl;
//                    for (auto &v : cuts) {
//                      for (auto &vv : v) {
//                        cout << vv.first << " " << vv.second << " ";
//                      }
//                      cout << endl;
//                    }
//                    cout << "no_mem: " << endl;
//                    for (auto &v : no_mem) {
//                      for (auto &vv : v) {
//                        cout << vv << " ";
//                      }
//                      cout << endl;
//                    }
//                    cout << "mem: " << endl;
//                    for (auto &v : mem) {
//                      for (auto &vv : v) {
//                        cout << vv << " ";
//                      }
//                      cout << endl;
//                    }
//                  }
//                }
//              } else goto HERE;
//            }
//
//          }
//        }
//      }
//    }
//  }
//  vector<int> tmp_forwad_soft;
//  int times = 0;
//  int hard_times = 0;
//  for (int i = 1; i < Dim; ++i) {
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      tmp_forwad_soft.emplace_back(LabelArrayInForwardSense[i][b].second);
//    }
//  }
//  std::sort(tmp_forwad_soft.begin(), tmp_forwad_soft.end(), greater<>());
//  int top_0_5 = (int) ((double) tmp_forwad_soft.size() * 0.005);
//  int top_1 = (int) ((double) tmp_forwad_soft.size() * 0.01);
//  int top_5 = (int) ((double) tmp_forwad_soft.size() * 0.05);
//  cout << "for forward exact///max= " << tmp_forwad_soft[0] << "  0.5%=" << tmp_forwad_soft[top_0_5]
//       << "  1%=" << tmp_forwad_soft[top_1] << "  5%=" << tmp_forwad_soft[top_5] << endl;
//
//  vector<pair<int, double>> record;
//  record.reserve(NumRow);
//
//  for (int i = 0; i < node->R1Cs.size(); ++i) {
//    int row_idx = node->R1Cs[i].IdxR1C;
//    record.emplace_back(i, abs(Pi[row_idx]) * pow(node->R1Cs[i].Mem.size(), CONFIG::MemFactor));
//  }
//
//  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
//    int row_idx = node->R1Cs_multi[i].IdxR1C;
//    record.emplace_back(i + node->R1Cs.size(),
//                        abs(Pi[row_idx]) * pow(node->R1Cs_multi[i].Mem.size(), CONFIG::MemFactor));
//  }
//
//  double
//      sum = accumulate(record.begin(), record.end(), 0.0, [](double a, pair<int, double> b) { return a + b.second; });
//
//  cout << "sum= " << sum << endl;
//  end = high_resolution_clock::now();
//  eps = duration<double>(end - beg).count();
//  cout << "forward exact time= " << eps << endl;
//#endif
//  if (Rollback == 1) {
//    cout << "Rollback to original states!" << endl;
//  } else if (Rollback == 2) {
//    cout << "Rollback with larger Mem!" << endl;
//  }
//}
//
//#ifdef SYMMETRY_PROHIBIT
//void CVRP::runHalfBackwardLabeling(BBNODE *const node, const PtrAllR1Cs &ptrAllR1Cs) {
//  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
//  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
//  double dif;
//  bool if_suc;
//
//  initializeLabels(node, 2, false, {true, 2, false});
//
//  auto beg = high_resolution_clock::now();
//  auto end = beg;
//  double eps;
//
//  int min_sorted_b = NumBucketsPerVertex - 1;
//  for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
//    int i = 1;
//    STILL_EXIST:
//    for (; i < Dim; ++i) {
//      auto &valid_num = IfExistExtraLabelsInBackwardSense[i][b].second;
//      if (!valid_num) continue;
//      auto &label_array = IfExistExtraLabelsInBackwardSense[i][b].first;
//      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
//        auto &ki = label_array[vec_index];
//        if (ki->if_extended) continue;
//        checkIfDominated<false>(ki, i, b, r1c_to_pi, r1c_multi_to_pi, if_suc);
//        ki->if_extended = true;
//        if (!if_suc) continue;//don't have to ++vec_index
//        auto sig = extendKernel4Exact_Backward(ki,
//                                               i,
//                                               ki->Sum_MainResource,
//                                               node->AllBackwardBuckets[i][b].BucketArcs,
//                                               r1c_to_pi,
//                                               r1c_multi_to_pi);
//        if (sig == 3) goto populateBin;
//        sig = extendKernel4Exact_Backward(ki, i, 0, node->AllBackwardBuckets[i][b].JumpArcs,
//                                          r1c_to_pi, r1c_multi_to_pi);
//        if (sig == 3) goto populateBin;
//      }
//      valid_num = 0;
//    }
//    //test if all labels are extended
//    end = high_resolution_clock::now();
//    eps = duration<double>(end - beg).count();
//    if (eps > CutGenTimeThresholdInPricing) {
//      Rollback = 3;
//      if (eps > CONFIG::HardTimeThresholdInPricing) {
//        if (!ForceNotRollback) {
//          Rollback = 1;
//          goto QUIT;
//        }
//      }
//    }
//    checkIfNoLabelsLeft<false>(i, b, min_sorted_b, if_suc);
//    if (!if_suc) goto STILL_EXIST;
//  }
//
//  populateBin:
//  if (min_sorted_b >= 0) {
//    for (int i = 1; i < Dim; ++i) {
//      for (int b = min_sorted_b; b >= 0; --b) {
//        std::stable_sort(LabelArrayInBackwardSense[i][b].first.begin(),
//                         LabelArrayInBackwardSense[i][b].first.begin() + LabelArrayInBackwardSense[i][b].second,
//                         CmpLabelRCLess);
//      }
//    }
//  }
//  //changing to arc architecture does not affect the following code
//  populateRC2TillThisBinNRC2Bin(node, 2);
//
//  QUIT:
//  if (Rollback == 1) {
//    cout << "Rollback to original states!" << endl;
////    return -1;
//  } else if (Rollback == 2) {
//    cout << "Rollback with larger Mem!" << endl;
////    return -2;
//  }
////  return ccnt;
//}
//
//template<typename T>
//int CVRP::extendKernel4Exact_Backward(LABEL *&ki,
//                                      int i,
//                                      double res,
//                                      const vector<T> &arc,
//                                      const double *r1c_to_pi,
//                                      const double *r1c_multi_to_pi) {
//  /*
//   */
//  int state = 0;
//  for (auto &pair : arc) {
//    int j;
//    if constexpr (std::is_same<T, int>::value) {
//      j = pair;
//    } else {
//      j = pair.second;
//      res = pair.first;
//    }
//    if (ki->PI[j]) continue;
//
//    int bj;
//    bool if_suc;
//    updateLabel<false, false>(res, ki, i, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
//    if (!if_suc) continue;
//
//    doDominance<false>(ki, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
//    if (!if_suc) continue;
//
//    addPathByRC(AllLabel[IdxGlo].RC + ChgCostMat4Vertex[j][0], AllLabel + IdxGlo, nullptr);
//
//    ++IdxGlo;
//    if (IdxGlo == LabelAssign) {
//      Rollback = 2;
//      state = 2;//QUIT
//      goto QUIT;
//    }
//  }
//  QUIT:
//  return state;
//}
//
//int CVRP::generateColsByBidir_NONE_Symmetry(BBNODE *node) {
//
//  PtrAllR1Cs ptrr1cs(node, this);
//
//  RE_TRY:
//  Rollback = 0;
//  runHalfLabeling<true>(node, ptrr1cs);//ccnt can only be applied for forward case
//
//  if (Rollback == 2) {
//    reallocateLabel();
//    goto RE_TRY;
//  } else if (Rollback == 1) {
//    return 0;
//  }
//
//  RE_TRY2:
//  Rollback = 0;
//  runHalfLabeling<false>(node, ptrr1cs);//ccnt can only be applied for forward case
//
//  if (Rollback == 2) {
//    reallocateLabel();
//    goto RE_TRY2;
//  } else if (Rollback == 1) {
//    return 0;
//  }
//
//  int ccnt = concatCols_NONE_Symmetry(node, ptrr1cs);
//
//  if (abs(GapBetweenLastSmallestRCAndRCThreshold) < TOLERANCE) {
//    double NumExistedLabels_Forward = 0;
//    double NumExistedLabels_Backward = 0;
//    for (int i = 1; i < Dim; ++i) {
//      for (int b = 0; b < NumBucketsPerVertex; ++b) {
//        NumExistedLabels_Forward += LabelArrayInForwardSense[i][b].second;
//        NumExistedLabels_Backward += LabelArrayInBackwardSense[i][b].second;
//      }
//    }
//    Ratio_DominanceChecks_NonDominant.first +=
//        NumDominanceChecks / (NumExistedLabels_Forward + NumExistedLabels_Backward);
//    ++Ratio_DominanceChecks_NonDominant.second;
//    double dif = abs(NumExistedLabels_Forward - NumExistedLabels_Backward);
//    double over = dif / min(NumExistedLabels_Forward, NumExistedLabels_Backward);
//#ifdef DETAILED_EXACT_PRINT_INFO
//    cout << "over= " << over << endl;
//#endif
//    if (over > CONFIG::NumberOfOverLabelsInMeetPoint) {
//#ifdef DETAILED_EXACT_PRINT_INFO
//      cout << "we adjust the meetpoint!" << endl;
//#endif
//      if (NumExistedLabels_Forward > NumExistedLabels_Backward) {
//        MeetPointResourceInBiDir *= (1 - CONFIG::MeetPointFactor);
//      } else {
//        MeetPointResourceInBiDir *= (1 + CONFIG::MeetPointFactor);
//      }
//#ifdef DETAILED_EXACT_PRINT_INFO
//      cout << "MeetPointResourceInBiDir= " << MeetPointResourceInBiDir << endl;
//#endif
//    }
//  }
//
//  if (ccnt) {
//    SmallestRC = get<2>(NegativeRCLabelTuple[0]);
//    if (abs(GapBetweenLastSmallestRCAndRCThreshold) > TOLERANCE) {
//      GapBetweenLastSmallestRCAndRCThreshold = get<2>(NegativeRCLabelTuple.back()) - SmallestRC;
//      if (GapBetweenLastSmallestRCAndRCThreshold < TOLERANCE)GapBetweenLastSmallestRCAndRCThreshold = 0;
//    } else {
//      if (ceil_transformed_number_related(SmallestRC * K + LPVal + RC_TOLERANCE) + TOLERANCE >= UB) {
//        node->Val = LPVal;
//        node->if_terminated = true;
//        cout << TERMINATED_MESSAGE_PROMISING_VEHICLES;
//        return 0;
//      }
//    }
//#ifdef DETAILED_EXACT_PRINT_INFO
//    cout << "Smallest= " << SmallestRC << endl;
//#endif
//  } else {
//    SmallestRC = 0;
//    GapBetweenLastSmallestRCAndRCThreshold = 0;
//  }
//
//#ifdef DETAILED_EXACT_PRINT_INFO
//  cout << "ccnt: " << ccnt << endl;
//#endif
//  //return 0
//  if (!ccnt) return 0;
//
//  addCols(node, ccnt);
//
//  return ccnt;
//}
//
//int CVRP::concatCols_NONE_Symmetry(BBNODE *node, const PtrAllR1Cs &ptrAllR1Cs) {
//  //update most_negative_rc
//  //write in PoolBeg4Pricing
//  PriorPoolBeg4Pricing = PoolBeg4Pricing;
//  if (checkPricingPool()) reallocatePricingPool();
//  int index = NumCol - 1;
//  double path_rc, tmp_rc, tmp_mainResource;
//  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
//  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
//  int *tmp_seq = new int[int(MaxMainResource) + 3];
//  int size_tmp_seq;
//  int i, j, arr_bj;
//  LABEL *p;
//  bool if_suc;
//
//  for (auto &label_list : concatenateLabelsInForwardCG) {
//    i = label_list.first.first;
//    j = label_list.first.second;
//    auto &label_vec = label_list.second;
//    for (auto &pr : label_vec) {
//      auto &ki = pr.first;
//      tmp_mainResource = pr.second;
//      arr_bj = int(tmp_mainResource / StepSize);
//      tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//
//      concatenateOneLabelWithOtherLabels<false, true>(ki, j, arr_bj, tmp_rc, tmp_mainResource, path_rc,
//                                                      r1c_to_pi, r1c_multi_to_pi, if_suc);
//      if (!if_suc)continue;
//
//
//      //bj-1
//      for (++arr_bj; arr_bj < NumBucketsPerVertex; ++arr_bj) {
//        concatenateOneLabelWithOtherLabels<false, false>(ki, j, arr_bj, tmp_rc, tmp_mainResource, path_rc,
//                                                         r1c_to_pi, r1c_multi_to_pi, if_suc);
//        if (!if_suc)break;
//      }
//    }
//  }
//
//  writeColsInPricingPool(node, index);
//  delete[]tmp_seq;
//  return index - NumCol + 1;
//}
//#endif
//
//int CVRP::forwardConcatenateInArcElimination(const double *r1c_to_pi,
//                                             const double *r1c_multi_to_pi) {
//  int status = 0;
//  bool if_find;
//  double path_rc;
//  int arr_bj;
//  auto beg = high_resolution_clock::now();
//  auto end = beg;
//  double eps;
//  for (auto &label_list : concatenateLabelsInForwardCG) {
//    int i = label_list.first.first;
//    int j = label_list.first.second;
//    auto &label_vec = label_list.second;
//    for (auto &pr : label_vec) {
//      auto &ki = pr.first;
//      auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
//      tmp_mainResource = pr.second;
//      auto &tmp_rc = AllLabel[IdxGlo].RC;
//      tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//      if_find = false;
//#ifdef SYMMETRY_PROHIBIT
//      arr_bj = int(tmp_mainResource / StepSize);
//      if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_till_this_bin
//        continue;
//
//      if (RC2BinInBackwardSense[j][arr_bj] + tmp_rc < OptGap) {//most_negative_rc_in_this_bin
//        //add one more condition for testing capacity
//        auto &label_arr = LabelArrayInBackwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInBackwardSense[j][arr_bj].second;
//#else
//      arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//      //bj
//      if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_till_this_bin
//        continue;
//
//      if (RC2BinInForwardSense[j][arr_bj] + tmp_rc < OptGap) {//most_negative_rc_in_this_bin
//        //add one more condition for testing capacity
//        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
//#endif
//        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//          auto &kj = label_arr[vec_index];
//          path_rc = kj->RC + tmp_rc;
//          if (path_rc > OptGap) break;
//
//          if ((ki->PI & kj->PI).any()) continue;
//
//#ifdef SYMMETRY_PROHIBIT
//          if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//          if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
//
//          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//              if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//            }
//          }
//
//          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = ki->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut] +
//                  ki->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (ki->Rank1CutMem_multi[tmp_cut] +
//                  kj->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          }
//
//          if (path_rc < OptGap) {
//            if_find = true;
//            goto outside;
//          }
//        }
//      }
//#ifdef SYMMETRY_PROHIBIT
//      //bj-1
//      for (++arr_bj; arr_bj < NumBucketsPerVertex; ++arr_bj) {
//        if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_till_this_bin
//          break;
//
//        if (RC2BinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_in_this_bin
//          continue;
//
//        auto &label_arr = LabelArrayInBackwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInBackwardSense[j][arr_bj].second;
//#else
//      //bj-1
//      for (--arr_bj; arr_bj >= 0; --arr_bj) {
//        if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_till_this_bin
//          break;
//
//        if (RC2BinInForwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_in_this_bin
//          continue;
//
//        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
//#endif
//        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//          auto &kj = label_arr[vec_index];
//          path_rc = kj->RC + tmp_rc;
//          if (path_rc > OptGap) break;
//
//          if ((ki->PI & kj->PI).any()) continue;
//
//          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//              if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//            }
//          }
//
//          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = ki->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut] +
//                  ki->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (ki->Rank1CutMem_multi[tmp_cut] +
//                  kj->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          }
//
//          if (path_rc < OptGap) {
//            if_find = true;
//            goto outside;
//          }
//        }
//      }
//
//      outside:
//      if (!if_find) {
//        continue;
//      }
//      //begin to extend by this direction
//      auto sig = extendKernel4ArcElimination_inner_Forward(ki,
//                                                           i,
//                                                           j,
//                                                           tmp_mainResource,
//                                                           r1c_to_pi,
//                                                           r1c_multi_to_pi);
//      if (sig == 1) continue;
//      if (sig == 2) {
//        status = 2;
//        goto QUIT;
//      }
//    }
//  }
//  end = chrono::high_resolution_clock::now();
//  eps = chrono::duration<double>(end - beg).count();
//  if (eps > CONFIG::HardTimeThresholdInArcElimination_Mid_concatenate) {
//    Rollback = 1;
//    status = 2;
//  }
//  QUIT:
//  return status;
//}
//
//template<typename T>
//int CVRP::extendKernel4ArcElimination_last_half_Forward(LABEL *&ki,
//                                                        int i,
//                                                        double res,
//                                                        const vector<T> &arc,
//                                                        const double *r1c_to_pi,
//                                                        const double *r1c_multi_to_pi) {
//  int state = 0;
//  int bj;
//  bool if_suc;
//  for (auto &pair : arc) {
//    int j;
//    if constexpr (std::is_same<T, int>::value) {
//      j = pair;
//    } else {
//      j = pair.second;
//      res = pair.first;
//    }
//
//    updateLabel<true, true, true, false>(res, ki, i, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
//    if (!if_suc) continue;
//
//    doDominance<true>(ki, j, bj, r1c_to_pi, r1c_multi_to_pi, if_suc);
//    if (!if_suc) continue;
//
//    ++IdxGlo;
//    if (IdxGlo == LabelAssign) {
//      Rollback = 2;
//      state = 2;//QUIT
//      goto QUIT;
//    }
//  }
//  QUIT:
//  return state;
//}
//
////int CVRP::firstHalfForwardInArcElimination(BBNODE *node, const double *r1c_to_pi,
////                                           const double *r1c_multi_to_pi) {
////  initializeLabels(node, 1, true, {true, 1, true});
////  bool if_break;
////  double dif;
////  int sig;
////
////  auto beg = high_resolution_clock::now();
////  auto end = beg;
////  double eps;
////
////
////  //first we extend to half way and obtain the info of rc_bin
////  for (int b = 0; b < NumBucketsPerVertex; ++b) {
////    int i = 1;
////    STILL_EXIST:
////    for (; i < Dim; ++i) {
////      if (!IfExistExtraLabelsInForwardSense[i][b]) continue;
////      auto &label_array = LabelArrayInForwardSense[i][b];
////      auto &valid_num = ValidNumLabelsInForwardSense[i][b];
////      for (int vec_index = 0; vec_index < valid_num;) {
////        auto &ki = label_array[vec_index];
////        if (ki->if_extended) {
////          ++vec_index;
////          continue;
////        }
////        if_break = false;
////        double tmp_ki_rc_sub = ki->RC - RC_TOLERANCE;
////        for (int b4_b = b - 1; b4_b >= 0; --b4_b) {
////          auto &b4_label_list = LabelArrayInForwardSense[i][b4_b];
////          auto &b4_valid_num = ValidNumLabelsInForwardSense[i][b4_b];
////          if (RC2TillThisBinInForwardSense[i][b4_b] > tmp_ki_rc_sub) break;
////          for (int vec_b4 = 0; vec_b4 < b4_valid_num; ++vec_b4) {
////            auto &b4_ki = b4_label_list[vec_b4];
////            if (b4_ki->RC > tmp_ki_rc_sub) break;
////            if (((ki->PI & b4_ki->PI) ^ (b4_ki->PI)).none()) {
////              dif = tmp_ki_rc_sub;
////              for (int l = 0; l < b4_ki->numValidRank1Cut; ++l) {
////                if (!ki->Rank1CutMem[b4_ki->validRank1Cut[l]]) dif += r1c_to_pi[b4_ki->validRank1Cut[l]];
////              }
////              for (int l = 0; l < b4_ki->numValidRank1Cut_multi; ++l) {
////                int tmp_cut = b4_ki->validRank1Cut_multi[l];
////                if (b4_ki->Rank1CutMem_multi[tmp_cut]
////                    > ki->Rank1CutMem_multi[tmp_cut]) {
////                  dif += r1c_multi_to_pi[tmp_cut];
////                }
////              }
////              if (dif > b4_ki->RC) {
////                ki = label_array[--valid_num];//-- should be first
////                if_break = true;
////                goto b4_break;
////              }
////            }
////          }
////        }
////        b4_break:
////        if (if_break) continue;//don't have to ++vec_index
////        ki->if_extended = true;
////        vector<pair<double, int>> bucketarcs(node->AllForwardBuckets[i][b].BucketArcs.size());
////        for (int j = 0; j < node->AllForwardBuckets[i][b].BucketArcs.size(); ++j) {
////          bucketarcs[j] = {
////              ki->Sum_MainResource,
////              node->AllForwardBuckets[i][b].BucketArcs[j]};
////        }
////        sig = extendKernel4ArcElimination_first_half_Forward(ki, i, bucketarcs, r1c_to_pi, r1c_multi_to_pi);
////        if (sig == 2) goto QUIT;
////        sig = extendKernel4ArcElimination_first_half_Forward(ki,
////                                                             i,
////                                                             node->AllForwardBuckets[i][b].JumpArcs,
////                                                             r1c_to_pi,
////                                                             r1c_multi_to_pi);
////        if (sig == 2) goto QUIT;
////        ++vec_index;
////      }
////      IfExistExtraLabelsInForwardSense[i][b] = false;
////    }
////    //test if all labels are extended
////    end = high_resolution_clock::now();
////    eps = duration<double>(end - beg).count();
////    if (eps > CONFIG::HardTimeThresholdInArcElimination_first_half) {
////      Rollback = 1;
////      sig = 2;
////      goto QUIT;
////    }
////    for (i = 1; i < Dim; ++i) {
////      if (IfExistExtraLabelsInForwardSense[i][b])
////        goto STILL_EXIST;
////    }
////    //write data
////    for (i = 1; i < Dim; ++i) {
////      std::stable_sort(LabelArrayInForwardSense[i][b].begin(),
////                       LabelArrayInForwardSense[i][b].begin() + ValidNumLabelsInForwardSense[i][b],
////                       CmpLabelRCLess);
////    }
////    if (b) {
////      for (i = 1; i < Dim; ++i) {
////        if (ValidNumLabelsInForwardSense[i][b]) {
////          RC2TillThisBinInForwardSense[i][b] = min(RC2TillThisBinInForwardSense[i][b - 1],
////                                                   LabelArrayInForwardSense[i][b][0]->RC);
////        } else {
////          RC2TillThisBinInForwardSense[i][b] = RC2TillThisBinInForwardSense[i][b - 1];
////        }
////      }
////    } else {
////      for (i = 1; i < Dim; ++i) {
////        if (ValidNumLabelsInForwardSense[i][b]) {
////          RC2TillThisBinInForwardSense[i][b] = LabelArrayInForwardSense[i][b][0]->RC;
////        } else {
////          RC2TillThisBinInForwardSense[i][b] = LARGEFLOAT;
////        }
////      }
////    }
////  }
////  populateRC2TillThisBinNRC2Bin(node);
////  QUIT:
////  return sig;
////}
//
//int CVRP::lastHalfForwardInArcElimination(BBNODE *node, const double *r1c_to_pi,
//                                          const double *r1c_multi_to_pi) {
//  bool if_break;
//  double dif;
//  int sig;
//  int status = 0;
//  auto beg = high_resolution_clock::now();
//  auto end = beg;
//  double eps;
//  for (int b = 0; b < NumBucketsPerVertex; ++b) {
//    int i = 1;
//    int cnt;
//    double sum_label;
//    STILL_EXIST1:
//    for (; i < Dim; ++i) {
//      auto &valid_num = IfExistExtraLabelsInForwardSense[i][b].second;
//      if (!valid_num) continue;
//      auto &label_array = IfExistExtraLabelsInForwardSense[i][b].first;
//      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
//        auto &ki = label_array[vec_index];
//        if (ki->if_extended) continue;
//        if_break = false;
//        double tmp_ki_rc_sub = ki->RC - RC_TOLERANCE;
//        for (int b4_b = b - 1; b4_b >= 0; --b4_b) {
//          auto &b4_label_list = LabelArrayInForwardSense[i][b4_b].first;
//          auto &b4_valid_num = LabelArrayInForwardSense[i][b4_b].second;
//          if (RC2TillThisBinInForwardSense[i][b4_b] > tmp_ki_rc_sub) break;
//          for (int vec_b4 = 0; vec_b4 < b4_valid_num; ++vec_b4) {
//            auto &b4_ki = b4_label_list[vec_b4];
//            if (b4_ki->RC > tmp_ki_rc_sub) break;
//            if (((ki->PI & b4_ki->PI) ^ (b4_ki->PI)).none()) {
//              dif = tmp_ki_rc_sub;
//              for (int l = 0; l < b4_ki->numValidRank1Cut; ++l) {
//                if (!ki->Rank1CutMem[b4_ki->validRank1Cut[l]]) dif += r1c_to_pi[b4_ki->validRank1Cut[l]];
//              }
//              for (int l = 0; l < b4_ki->numValidRank1Cut_multi; ++l) {
//                int tmp_cut = b4_ki->validRank1Cut_multi[l];
//                if (b4_ki->Rank1CutMem_multi[tmp_cut]
//                    > ki->Rank1CutMem_multi[tmp_cut]) {
//                  dif += r1c_multi_to_pi[tmp_cut];
//                }
//              }
//              if (dif > b4_ki->RC) {
//                if_break = true;
//                goto b4_break1;
//              }
//            }
//          }
//        }
//        b4_break1:
//        ki->if_extended = true;
//        if (if_break) continue;//don't have to ++vec_index
//        sig = extendKernel4ArcElimination_last_half_Forward(ki,
//                                                            i,
//                                                            ki->Sum_MainResource,
//                                                            node->AllForwardBuckets[i][b].BucketArcs,
//                                                            r1c_to_pi,
//                                                            r1c_multi_to_pi);
//        if (sig == 2) {
//          status = 2;
//          goto QUIT;
//        }
//        sig = extendKernel4ArcElimination_last_half_Forward(ki,
//                                                            i,
//                                                            0,
//                                                            node->AllForwardBuckets[i][b].JumpArcs,
//                                                            r1c_to_pi, r1c_multi_to_pi);
//        if (sig == 2) {
//          status = 2;
//          goto QUIT;
//        }
//      }
//      valid_num = 0;
//    }
//    //test if all labels are extended
//    end = high_resolution_clock::now();
//    eps = duration<double>(end - beg).count();
//    if (eps > CONFIG::HardTimeThresholdInArcElimination_last_half) {
//      status = 1;
//      Rollback = 1;
//      goto QUIT;
//    }
//    for (i = 1; i < Dim; ++i) {
//      if (IfExistExtraLabelsInForwardSense[i][b].second)
//        goto STILL_EXIST1;
//    }
//    //write data
//    for (i = 1; i < Dim; ++i) {
//      std::stable_sort(LabelArrayInForwardSense[i][b].first.begin(),
//                       LabelArrayInForwardSense[i][b].first.begin() + LabelArrayInForwardSense[i][b].second,
//                       CmpLabelRCLess);
//    }
//    if (b) {
//      for (i = 1; i < Dim; ++i) {
//        if (LabelArrayInForwardSense[i][b].second) {
//          RC2TillThisBinInForwardSense[i][b] = min(RC2TillThisBinInForwardSense[i][b - 1],
//                                                   LabelArrayInForwardSense[i][b].first[0]->RC);
//        } else {
//          RC2TillThisBinInForwardSense[i][b] = RC2TillThisBinInForwardSense[i][b - 1];
//        }
//      }
//    } else {
//      for (i = 1; i < Dim; ++i) {
//        if (LabelArrayInForwardSense[i][b].second) {
//          RC2TillThisBinInForwardSense[i][b] = LabelArrayInForwardSense[i][b].first[0]->RC;
//        } else {
//          RC2TillThisBinInForwardSense[i][b] = LARGEFLOAT;
//        }
//      }
//    }
//  }
//  populateRC2TillThisBinNRC2Bin(node);
//  QUIT:
//  return status;
//}
//
//
////int CVRP::extendKernel4ArcElimination_first_half_Forward(LABEL *&ki,
////                                                         int i,
////                                                         const vector<pair<double, int>> &arc,
////                                                         const double *r1c_to_pi,
////                                                         const double *r1c_multi_to_pi) {
////  int state = 0;
////  for (auto &pair : arc) {
////    int j = pair.second;
////    if (ki->PI[j]) continue;
////
////    auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
////    if (!increaseMainResourceConsumption(pair.first, tmp_mainResource, i, j)) continue;
////
////    if (tmp_mainResource > MeetPointResourceInBiDir) {
////      int edgemap;
////      if (i < j) edgemap = i * Dim + j;
////      else edgemap = j * Dim + i;
////      concatenateLabels4ArcEliminationInForwardSense[edgemap].emplace_back(ki);
////      continue;
////    }
////
////    auto sig = extendKernel4ArcElimination_inner_Forward(ki, i, j, tmp_mainResource, r1c_to_pi, r1c_multi_to_pi);
////    if (sig == 1) continue;
////    if (sig == 2) {
////      state = 2;
////      goto QUIT;
////    }
////  }
////  QUIT:
////  return state;
////}
//
//
//
//int CVRP::extendKernel4ArcElimination_inner_Forward(LABEL *&ki,
//                                                    int i, int j, double &tmp_mainResource,
//                                                    const double *r1c_to_pi,
//                                                    const double *r1c_multi_to_pi) {
//  int state = 0;
//  int bj = int(tmp_mainResource / StepSize);
//  auto &labelList_j = LabelArrayInForwardSense[j][bj].first;
//  auto &valid_num_j = LabelArrayInForwardSense[j][bj].second;
//  auto &tmp_rc = AllLabel[IdxGlo].RC;
//  auto &tmp_PI = AllLabel[IdxGlo].PI;
//  auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
//  auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
//  auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
//  tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
//  tmp_PI = (ki->PI) & (NGMem4Vertex[j]);
//  tmp_PI.set(j);
//  tmp_Rank1CutMem = ki->Rank1CutMem;
//  //tmp_num_valid_rank1_cut do not have to copy
//  for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
//    if (tmp_Rank1CutMem[l]) {
//      tmp_Rank1CutMem[l] = false;
//      tmp_rc -= r1c_to_pi[l];
//    } else tmp_Rank1CutMem[l] = true;
//  }
//  tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);
//  tmp_num_valid_rank1_cut = 0;
//  for (auto l : get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
//    if (tmp_Rank1CutMem[l]) {
//      tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
//    }
//  }
//
//  auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
//  auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
//  auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;
//  copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);
//
//  for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
//    int tmp_cut = get<0>(l);
//    tmp_Rank1CutMem_multi[tmp_cut] += get<1>(l);
//    if (tmp_Rank1CutMem_multi[tmp_cut] >= get<2>(l)) {
//      tmp_rc -= r1c_multi_to_pi[tmp_cut];
//      tmp_Rank1CutMem_multi[tmp_cut] -= get<2>(l);
//    }
//  }
//
//  for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;
//
//  tmp_num_valid_rank1_cut_multi = 0;
//  for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
//    if (tmp_Rank1CutMem_multi[l]) {
//      tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
//    }
//  }
//
//  double tmp_rc_add = tmp_rc + RC_TOLERANCE, tmp_rc_sub = tmp_rc - RC_TOLERANCE;
//
//  bool if_break = false;
//  double dif;
//  for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
//    auto &kj = labelList_j[vec_index_j];
//    if (kj->Sum_MainResource < tmp_mainResource) {
//      if (kj->RC < tmp_rc_sub) {
//        if (((tmp_PI & kj->PI) ^ (kj->PI)).none()) {
//          dif = tmp_rc_sub;
//          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//            if (!tmp_Rank1CutMem[kj->validRank1Cut[l]])dif += r1c_to_pi[kj->validRank1Cut[l]];
//          }
//          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//            int tmp = kj->validRank1Cut_multi[l];
//            if (kj->Rank1CutMem_multi[tmp]
//                > tmp_Rank1CutMem_multi[tmp]) {
//              dif += r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif > kj->RC) {
//            if_break = true;
//            break;//taken
//          }
//        }
//      }
//      ++vec_index_j;
//    } else if (tmp_mainResource < kj->Sum_MainResource) {
//      if (tmp_rc_add < kj->RC) {
//        if (((tmp_PI & kj->PI) ^ (tmp_PI)).none()) {
//          dif = tmp_rc_add;
//          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]])dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
//          }
//          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//            int tmp = tmp_valid_rank1_cut_multi[l];
//            if (tmp_Rank1CutMem_multi[tmp]
//                > kj->Rank1CutMem_multi[tmp]) {
//              dif -= r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif < kj->RC) {
//            kj->if_extended = true;
//            kj = labelList_j[--valid_num_j];
//          } else ++vec_index_j;
//        } else ++vec_index_j;
//      } else ++vec_index_j;
//    } else {
//      if (kj->RC < tmp_rc_add) {
//        if (((tmp_PI & kj->PI) ^ (kj->PI)).none()) {
//          dif = tmp_rc_sub;
//          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//            if (!tmp_Rank1CutMem[kj->validRank1Cut[l]])dif += r1c_to_pi[kj->validRank1Cut[l]];
//          }
//          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//            int tmp = kj->validRank1Cut_multi[l];
//            if (kj->Rank1CutMem_multi[tmp]
//                > tmp_Rank1CutMem_multi[tmp]) {
//              dif += r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif > kj->RC) {
//            if_break = true;
//            break;//taken
//          }
//        }
//        ++vec_index_j;
//      } else if (tmp_rc_sub < kj->RC) {
//        if (((tmp_PI & kj->PI) ^ (tmp_PI)).none()) {
//          dif = tmp_rc_add;
//          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]])dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
//          }
//          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//            int tmp = tmp_valid_rank1_cut_multi[l];
//            if (tmp_Rank1CutMem_multi[tmp]
//                > kj->Rank1CutMem_multi[tmp]) {
//              dif -= r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif < kj->RC) {
//            kj->if_extended = true;
//            kj = labelList_j[--valid_num_j];
//          } else ++vec_index_j;
//        } else ++vec_index_j;
//      } else {//q & rc all equal now?
//        if ((tmp_PI ^ kj->PI).none()) {//all equal
//          bitset<2> who_win = 1;//1==tmp_win,2==kj win,0==all keep
//          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]]) {
//              who_win = 2;
//              goto next_test3;
//            }
//          }
//          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//            int tmp = tmp_valid_rank1_cut_multi[l];
//            if (tmp_Rank1CutMem_multi[tmp]
//                > kj->Rank1CutMem_multi[tmp]) {
//              who_win = 2;
//              goto next_test3;
//            }
//          }
//          next_test3:
//          if (who_win == 1) {
//            kj->if_extended = true;
//            kj = labelList_j[--valid_num_j];
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (!tmp_Rank1CutMem[kj->validRank1Cut[l]]) {
//                who_win = 0;
//                goto next_test4;
//              }
//            }
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp = kj->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp]
//                  > tmp_Rank1CutMem_multi[tmp]) {
//                who_win = 0;
//                goto next_test4;
//              }
//            }
//            next_test4:
//            if (who_win == 2) {
//              if_break = true;
//              break;//taken
//            } else ++vec_index_j;
//          }
//        } else {
//          yzzLong tmp = tmp_PI & kj->PI;
//          if ((tmp ^ tmp_PI).none()) {
//            dif = tmp_rc_add;
//            for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//              if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]])dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
//            }
//            for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//              int tmp_cut = tmp_valid_rank1_cut_multi[l];
//              if (tmp_Rank1CutMem_multi[tmp_cut]
//                  > kj->Rank1CutMem_multi[tmp_cut]) {
//                dif -= r1c_multi_to_pi[tmp_cut];
//              }
//            }
//            if (dif < kj->RC) {
//              kj->if_extended = true;
//              kj = labelList_j[--valid_num_j];
//            } else ++vec_index_j;
//          } else if ((tmp ^ kj->PI).none()) {
//            dif = tmp_rc_sub;
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (!tmp_Rank1CutMem[kj->validRank1Cut[l]])dif += r1c_to_pi[kj->validRank1Cut[l]];
//            }
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut]
//                  > tmp_Rank1CutMem_multi[tmp_cut]) {
//                dif += r1c_multi_to_pi[tmp_cut];
//              }
//            }
//            if (dif > kj->RC) {
//              if_break = true;
//              break;//taken
//            }
//            ++vec_index_j;
//          } else ++vec_index_j;
//        }
//      }
//    }
//  }
//
//  auto &bucket = IfExistExtraLabelsInForwardSense[j][bj];
//  if (if_break) {
//    state = 1;
//    goto QUIT;
//  }
//
//  labelList_j[valid_num_j++] = AllLabel + IdxGlo;
//  if (valid_num_j == (int) labelList_j.size()) {
//    labelList_j.resize(labelList_j.size() * 2);
//  }
//
//  //EndVertex
//  AllLabel[IdxGlo].EndVertex = j;
//  //if_extended
//  AllLabel[IdxGlo].if_extended = false;
//  //bucket
//  bucket.first[bucket.second++] = AllLabel + IdxGlo;
//  if (bucket.second == bucket.first.size()) {
//    bucket.first.resize(bucket.first.size() * 2);
//  }
//
//  //state_b4
//  ++IdxGlo;
//
//  if (IdxGlo == LabelAssign) {
//    Rollback = 2;
//    state = 2;
//    goto QUIT;
//  }
//  QUIT:
//  return state;
//}
//
//#ifdef SYMMETRY_PROHIBIT
////int CVRP::extendKernel4ArcElimination_first_half_Backward(LABEL *&ki,
////                                                          int i,
////                                                          const std::vector<std::pair<double, int>> &arc,
////                                                          const double *r1c_to_pi,
////                                                          const double *r1c_multi_to_pi) {
////  int state = 0;
////  for (auto &pair : arc) {
////    int j = pair.second;
////    if (ki->PI[j]) continue;
////
////    auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
////    if (!decreaseMainResourceConsumption(pair.first, tmp_mainResource, i, j)) continue;
////
////    if (tmp_mainResource < MeetPointResourceInBiDir) {
////      int edgemap;
////      if (i < j) edgemap = i * Dim + j;
////      else edgemap = j * Dim + i;
////      concatenateLabels4ArcEliminationInBackwardSense[edgemap].emplace_back(ki);
////      continue;
////    }
////
////    auto sig = extendKernel4ArcElimination_inner_Backward(ki, i, j, tmp_mainResource, r1c_to_pi, r1c_multi_to_pi);
////    if (sig == 1) continue;
////    if (sig == 2) {
////      state = 2;
////      goto QUIT;
////    }
////  }
////  QUIT:
////  return state;
////}
//template<typename T>
//int CVRP::extendKernel4ArcElimination_last_half_Backward(LABEL *&ki,
//                                                         int i,
//                                                         double res,
//                                                         const std::vector<T> &arc,
//                                                         const double *r1c_to_pi,
//                                                         const double *r1c_multi_to_pi) {
//  int state = 0;
//  for (auto &pair : arc) {
//    int j;
//    if constexpr (std::is_same<T, int>::value) {
//      j = pair;
//    } else {
//      j = pair.second;
//      res = pair.first;
//    }
//    if (ki->PI[j]) continue;
//
//    auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
//    if (!decreaseMainResourceConsumption(res, tmp_mainResource, i, j)) continue;
//
//    auto &tmp_rc = AllLabel[IdxGlo].RC;
//    tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
//
//    bool if_keep = false;
//    double path_rc;
//    int arr_bj = int(tmp_mainResource / StepSize);
//    //first test
//    if (tmp_rc + RC2TillThisBinInForwardSense[j][arr_bj] > OptGap) continue;
//    //second test
//    if (tmp_rc + RC2BinInForwardSense[j][arr_bj] < OptGap) {
//      auto &label_list_arr_bj = LabelArrayInForwardSense[j][arr_bj].first;
//      auto &val_num_arr_bj = LabelArrayInForwardSense[j][arr_bj].second;
//      for (int v_bj = 0; v_bj < val_num_arr_bj; ++v_bj) {
//        auto &kkj = label_list_arr_bj[v_bj];
//        if (tmp_mainResource < kkj->Sum_MainResource) continue;
//        if ((ki->PI & kkj->PI).any()) continue;
//        path_rc = tmp_rc + kkj->RC;
//        if (path_rc > OptGap) break;
//
//        if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
//          for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//            if (kkj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//          }
//        } else {
//          for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
//            if (ki->Rank1CutMem[kkj->validRank1Cut[l]])path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
//          }
//        }
//
//        if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
//          for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = ki->validRank1Cut_multi[l];
//            if (kkj->Rank1CutMem_multi[tmp_cut] +
//                ki->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        } else {
//          for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = kkj->validRank1Cut_multi[l];
//            if (ki->Rank1CutMem_multi[tmp_cut] +
//                kkj->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        }
//
//        if (path_rc < OptGap) {
//          if_keep = true;
//          goto outside2;
//        }
//      }
//    }
//    //real test
//    for (--arr_bj; arr_bj >= 0; --arr_bj) {
//      //first test
//      if (tmp_rc + RC2TillThisBinInForwardSense[j][arr_bj] > OptGap) break;
//      //second test
//      if (tmp_rc + RC2BinInForwardSense[j][arr_bj] > OptGap) continue;
//      //real test
//      auto &label_list_arr_bj = LabelArrayInForwardSense[j][arr_bj].first;
//      auto &val_num_arr_bj = LabelArrayInForwardSense[j][arr_bj].second;
//      for (int v_bj = 0; v_bj < val_num_arr_bj; ++v_bj) {
//        auto &kkj = label_list_arr_bj[v_bj];
//        if ((ki->PI & kkj->PI).any()) continue;
//        path_rc = tmp_rc + kkj->RC;
//        if (path_rc > OptGap) break;
//
//        if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
//          for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//            if (kkj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//          }
//        } else {
//          for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
//            if (ki->Rank1CutMem[kkj->validRank1Cut[l]])path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
//          }
//        }
//
//        if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
//          for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = ki->validRank1Cut_multi[l];
//            if (kkj->Rank1CutMem_multi[tmp_cut] +
//                ki->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        } else {
//          for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = kkj->validRank1Cut_multi[l];
//            if (ki->Rank1CutMem_multi[tmp_cut] +
//                kkj->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        }
//
//        if (path_rc < OptGap) {
//          if_keep = true;
//          goto outside2;
//        }
//      }
//    }
//    outside2:
//    if (!if_keep) continue;
//    auto sig = extendKernel4ArcElimination_inner_Backward(ki, i, j, tmp_mainResource, r1c_to_pi, r1c_multi_to_pi);
//    if (sig == 1) continue;
//    if (sig == 2) {
//      state = 2;
//      goto QUIT;
//    }
//  }
//  QUIT:
//  return state;
//}
//
//int CVRP::extendKernel4ArcElimination_inner_Backward(LABEL *&ki,
//                                                     int i, int j, double &tmp_mainResource,
//                                                     const double *r1c_to_pi,
//                                                     const double *r1c_multi_to_pi) {
//  int state = 0;
//  int bj = int(tmp_mainResource / StepSize);
//  auto &labelList_j = LabelArrayInBackwardSense[j][bj].first;
//  auto &valid_num_j = LabelArrayInBackwardSense[j][bj].second;
//  auto &tmp_rc = AllLabel[IdxGlo].RC;
//  auto &tmp_PI = AllLabel[IdxGlo].PI;
//  auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
//  auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
//  auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
//  tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
//  tmp_PI = (ki->PI) & (NGMem4Vertex[j]);
//  tmp_PI.set(j);
//  tmp_Rank1CutMem = ki->Rank1CutMem;
//  //tmp_num_valid_rank1_cut do not have to copy
//  for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
//    if (tmp_Rank1CutMem[l]) {
//      tmp_Rank1CutMem[l] = false;
//      tmp_rc -= r1c_to_pi[l];
//    } else tmp_Rank1CutMem[l] = true;
//  }
//  tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);
//  tmp_num_valid_rank1_cut = 0;
//  for (auto l : get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
//    if (tmp_Rank1CutMem[l]) {
//      tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
//    }
//  }
//
//  auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
//  auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
//  auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;
//  copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);
//
//  for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
//    int tmp_cut = get<0>(l);
//    tmp_Rank1CutMem_multi[tmp_cut] += get<1>(l);
//    if (tmp_Rank1CutMem_multi[tmp_cut] >= get<2>(l)) {
//      tmp_rc -= r1c_multi_to_pi[tmp_cut];
//      tmp_Rank1CutMem_multi[tmp_cut] -= get<2>(l);
//    }
//  }
//
//  for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;
//
//  tmp_num_valid_rank1_cut_multi = 0;
//  for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
//    if (tmp_Rank1CutMem_multi[l]) {
//      tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
//    }
//  }
//
//  double tmp_rc_add = tmp_rc + RC_TOLERANCE, tmp_rc_sub = tmp_rc - RC_TOLERANCE;
//
//  bool if_break = false;
//  double dif;
//  for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
//    auto &kj = labelList_j[vec_index_j];
//    if (kj->Sum_MainResource > tmp_mainResource) {
//      if (kj->RC < tmp_rc_sub) {
//        if (((tmp_PI & kj->PI) ^ (kj->PI)).none()) {
//          dif = tmp_rc_sub;
//          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//            if (!tmp_Rank1CutMem[kj->validRank1Cut[l]])dif += r1c_to_pi[kj->validRank1Cut[l]];
//          }
//          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//            int tmp = kj->validRank1Cut_multi[l];
//            if (kj->Rank1CutMem_multi[tmp]
//                > tmp_Rank1CutMem_multi[tmp]) {
//              dif += r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif > kj->RC) {
//            if_break = true;
//            break;//taken
//          }
//        }
//      }
//      ++vec_index_j;
//    } else if (tmp_mainResource > kj->Sum_MainResource) {
//      if (tmp_rc_add < kj->RC) {
//        if (((tmp_PI & kj->PI) ^ (tmp_PI)).none()) {
//          dif = tmp_rc_add;
//          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]])dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
//          }
//          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//            int tmp = tmp_valid_rank1_cut_multi[l];
//            if (tmp_Rank1CutMem_multi[tmp]
//                > kj->Rank1CutMem_multi[tmp]) {
//              dif -= r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif < kj->RC) {
//            kj->if_extended = true;
//            kj = labelList_j[--valid_num_j];
//          } else ++vec_index_j;
//        } else ++vec_index_j;
//      } else ++vec_index_j;
//    } else {
//      if (kj->RC < tmp_rc_add) {
//        if (((tmp_PI & kj->PI) ^ (kj->PI)).none()) {
//          dif = tmp_rc_sub;
//          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//            if (!tmp_Rank1CutMem[kj->validRank1Cut[l]])dif += r1c_to_pi[kj->validRank1Cut[l]];
//          }
//          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//            int tmp = kj->validRank1Cut_multi[l];
//            if (kj->Rank1CutMem_multi[tmp]
//                > tmp_Rank1CutMem_multi[tmp]) {
//              dif += r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif > kj->RC) {
//            if_break = true;
//            break;//taken
//          }
//        }
//        ++vec_index_j;
//      } else if (tmp_rc_sub < kj->RC) {
//        if (((tmp_PI & kj->PI) ^ (tmp_PI)).none()) {
//          dif = tmp_rc_add;
//          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]])dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
//          }
//          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//            int tmp = tmp_valid_rank1_cut_multi[l];
//            if (tmp_Rank1CutMem_multi[tmp]
//                > kj->Rank1CutMem_multi[tmp]) {
//              dif -= r1c_multi_to_pi[tmp];
//            }
//          }
//          if (dif < kj->RC) {
//            kj->if_extended = true;
//            kj = labelList_j[--valid_num_j];
//          } else ++vec_index_j;
//        } else ++vec_index_j;
//      } else {//q & rc all equal now?
//        if ((tmp_PI ^ kj->PI).none()) {//all equal
//          bitset<2> who_win = 1;//1==tmp_win,2==kj win,0==all keep
//          for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//            if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]]) {
//              who_win = 2;
//              goto next_test3;
//            }
//          }
//          for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//            int tmp = tmp_valid_rank1_cut_multi[l];
//            if (tmp_Rank1CutMem_multi[tmp]
//                > kj->Rank1CutMem_multi[tmp]) {
//              who_win = 2;
//              goto next_test3;
//            }
//          }
//          next_test3:
//          if (who_win == 1) {
//            kj->if_extended = true;
//            kj = labelList_j[--valid_num_j];
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (!tmp_Rank1CutMem[kj->validRank1Cut[l]]) {
//                who_win = 0;
//                goto next_test4;
//              }
//            }
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp = kj->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp]
//                  > tmp_Rank1CutMem_multi[tmp]) {
//                who_win = 0;
//                goto next_test4;
//              }
//            }
//            next_test4:
//            if (who_win == 2) {
//              if_break = true;
//              break;//taken
//            } else ++vec_index_j;
//          }
//        } else {
//          yzzLong tmp = tmp_PI & kj->PI;
//          if ((tmp ^ tmp_PI).none()) {
//            dif = tmp_rc_add;
//            for (int l = 0; l < tmp_num_valid_rank1_cut; ++l) {
//              if (!kj->Rank1CutMem[tmp_valid_rank1_cut[l]])dif -= r1c_to_pi[tmp_valid_rank1_cut[l]];
//            }
//            for (int l = 0; l < tmp_num_valid_rank1_cut_multi; ++l) {
//              int tmp_cut = tmp_valid_rank1_cut_multi[l];
//              if (tmp_Rank1CutMem_multi[tmp_cut]
//                  > kj->Rank1CutMem_multi[tmp_cut]) {
//                dif -= r1c_multi_to_pi[tmp_cut];
//              }
//            }
//            if (dif < kj->RC) {
//              kj->if_extended = true;
//              kj = labelList_j[--valid_num_j];
//            } else ++vec_index_j;
//          } else if ((tmp ^ kj->PI).none()) {
//            dif = tmp_rc_sub;
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (!tmp_Rank1CutMem[kj->validRank1Cut[l]])dif += r1c_to_pi[kj->validRank1Cut[l]];
//            }
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut]
//                  > tmp_Rank1CutMem_multi[tmp_cut]) {
//                dif += r1c_multi_to_pi[tmp_cut];
//              }
//            }
//            if (dif > kj->RC) {
//              if_break = true;
//              break;//taken
//            }
//            ++vec_index_j;
//          } else ++vec_index_j;
//        }
//      }
//    }
//  }
//
//  auto &bucket = IfExistExtraLabelsInBackwardSense[j][bj];
//  if (if_break) {
//    state = 1;
//    goto QUIT;
//  }
//
//  labelList_j[valid_num_j++] = AllLabel + IdxGlo;
//  if (valid_num_j == (int) labelList_j.size()) {
//    labelList_j.resize(labelList_j.size() * 2);
//  }
//
//  //EndVertex
//  AllLabel[IdxGlo].EndVertex = j;
//  //if_extended
//  AllLabel[IdxGlo].if_extended = false;
//  //bucket
//  bucket.first[bucket.second++] = AllLabel + IdxGlo;
//  if (bucket.second == bucket.first.size()) {
//    bucket.first.resize(bucket.first.size() * 2);
//  }
//  //state_b4
//  ++IdxGlo;
//
//  if (IdxGlo == LabelAssign) {
//    Rollback = 2;
//    state = 2;
//    goto QUIT;
//  }
//  QUIT:
//  return state;
//}
//
////int CVRP::firstHalfBackwardInArcElimination(BBNODE *node, const double *r1c_to_pi, const double *r1c_multi_to_pi) {
////  initializeLabels(node, 2, false, {true, 2, false});
////  bool if_break;
////  double dif;
////  int sig;
////
////  auto beg = high_resolution_clock::now();
////  auto end = beg;
////  double eps;
////
////  //first we extend to half way and obtain the info of rc_bin
////  for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
////    int i = 1;
////    STILL_EXIST:
////    for (; i < Dim; ++i) {
////      if (!IfExistExtraLabelsInBackwardSense[i][b]) continue;
////      auto &label_array = LabelArrayInBackwardSense[i][b];
////      auto &valid_num = ValidNumLabelsInBackwardSense[i][b];
////      for (int vec_index = 0; vec_index < valid_num;) {
////        auto &ki = label_array[vec_index];
////        if (ki->if_extended) {
////          ++vec_index;
////          continue;
////        }
////        if_break = false;
////        double tmp_ki_rc_sub = ki->RC - RC_TOLERANCE;
////        for (int b4_b = b + 1; b4_b < NumBucketsPerVertex; ++b4_b) {
////          auto &b4_label_list = LabelArrayInBackwardSense[i][b4_b];
////          auto &b4_valid_num = ValidNumLabelsInBackwardSense[i][b4_b];
////          if (RC2TillThisBinInBackwardSense[i][b4_b] > tmp_ki_rc_sub) break;
////          for (int vec_b4 = 0; vec_b4 < b4_valid_num; ++vec_b4) {
////            auto &b4_ki = b4_label_list[vec_b4];
////            if (b4_ki->RC > tmp_ki_rc_sub) break;
////            if (((ki->PI & b4_ki->PI) ^ (b4_ki->PI)).none()) {
////              dif = tmp_ki_rc_sub;
////              for (int l = 0; l < b4_ki->numValidRank1Cut; ++l) {
////                if (!ki->Rank1CutMem[b4_ki->validRank1Cut[l]]) dif += r1c_to_pi[b4_ki->validRank1Cut[l]];
////              }
////              for (int l = 0; l < b4_ki->numValidRank1Cut_multi; ++l) {
////                int tmp_cut = b4_ki->validRank1Cut_multi[l];
////                if (b4_ki->Rank1CutMem_multi[tmp_cut]
////                    > ki->Rank1CutMem_multi[tmp_cut]) {
////                  dif += r1c_multi_to_pi[tmp_cut];
////                }
////              }
////
////              if (dif > b4_ki->RC) {
////                ki = label_array[--valid_num];//-- should be first
////                if_break = true;
////                goto b4_break;
////              }
////            }
////          }
////        }
////        b4_break:
////        if (if_break) continue;//don't have to ++vec_index
////        ki->if_extended = true;
////        vector<pair<double, int>> bucketarcs(node->AllBackwardBuckets[i][b].BucketArcs.size());
////        for (int j = 0; j < node->AllBackwardBuckets[i][b].BucketArcs.size(); ++j) {
////          bucketarcs[j] = {
////              ki->Sum_MainResource,
////              node->AllBackwardBuckets[i][b].BucketArcs[j]};
////        }
////        sig = extendKernel4ArcElimination_first_half_Backward(ki, i, bucketarcs, r1c_to_pi, r1c_multi_to_pi);
////        if (sig == 2) goto QUIT;
////        sig =
////            extendKernel4ArcElimination_first_half_Backward(ki,
////                                                            i,
////                                                            node->AllBackwardBuckets[i][b].JumpArcs,
////                                                            r1c_to_pi,
////                                                            r1c_multi_to_pi);
////        if (sig == 2) goto QUIT;
////        ++vec_index;
////      }
////      IfExistExtraLabelsInBackwardSense[i][b] = false;
////    }
////    //test if all labels are extended
////    end = high_resolution_clock::now();
////    eps = duration<double>(end - beg).count();
////    if (eps > CONFIG::HardTimeThresholdInArcElimination_first_half) {
////      sig = 2;
////      goto QUIT;
////    }
////    for (i = 1; i < Dim; ++i) {
////      if (IfExistExtraLabelsInBackwardSense[i][b])
////        goto STILL_EXIST;
////    }
////    //write data
////    for (i = 1; i < Dim; ++i) {
////      std::stable_sort(LabelArrayInBackwardSense[i][b].begin(),
////                       LabelArrayInBackwardSense[i][b].begin() + ValidNumLabelsInBackwardSense[i][b],
////                       CmpLabelRCLess);
////    }
////    if (b != NumBucketsPerVertex - 1) {
////      for (i = 1; i < Dim; ++i) {
////        if (ValidNumLabelsInBackwardSense[i][b]) {
////          RC2TillThisBinInBackwardSense[i][b] = min(RC2TillThisBinInBackwardSense[i][b + 1],
////                                                    LabelArrayInBackwardSense[i][b][0]->RC);
////        } else {
////          RC2TillThisBinInBackwardSense[i][b] = RC2TillThisBinInBackwardSense[i][b + 1];
////        }
////      }
////    } else {
////      for (i = 1; i < Dim; ++i) {
////        if (ValidNumLabelsInBackwardSense[i][b]) {
////          RC2TillThisBinInBackwardSense[i][b] = LabelArrayInBackwardSense[i][b][0]->RC;
////        } else {
////          RC2TillThisBinInBackwardSense[i][b] = LARGEFLOAT;
////        }
////      }
////    }
////  }
////  populateRC2TillThisBinNRC2Bin(node, 2);
////  QUIT:
////  return sig;
////}
//
//int CVRP::lastHalfBackwardInArcElimination(BBNODE *node,
//                                           const double *r1c_to_pi,
//                                           const double *r1c_multi_to_pi) {
//  bool if_break;
//  double dif;
//  int sig;
//  int status = 0;
//  auto beg = high_resolution_clock::now();
//  auto end = beg;
//  double eps;
//  for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
//    int i = 1;
//    int cnt;
//    double sum_label;
//    STILL_EXIST1:
//    for (; i < Dim; ++i) {
//      auto &valid_num = IfExistExtraLabelsInBackwardSense[i][b].second;
//      if (!valid_num) continue;
//      auto &label_array = IfExistExtraLabelsInBackwardSense[i][b].first;
//      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
//        auto &ki = label_array[vec_index];
//        if (ki->if_extended) continue;
//        if_break = false;
//        double tmp_ki_rc_sub = ki->RC - RC_TOLERANCE;
//        for (int b4_b = b + 1; b4_b < NumBucketsPerVertex; ++b4_b) {
//          auto &b4_label_list = LabelArrayInBackwardSense[i][b4_b].first;
//          auto &b4_valid_num = LabelArrayInBackwardSense[i][b4_b].second;
//          if (RC2TillThisBinInBackwardSense[i][b4_b] > tmp_ki_rc_sub) break;
//          for (int vec_b4 = 0; vec_b4 < b4_valid_num; ++vec_b4) {
//            auto &b4_ki = b4_label_list[vec_b4];
//            if (b4_ki->RC > tmp_ki_rc_sub) break;
//            if (((ki->PI & b4_ki->PI) ^ (b4_ki->PI)).none()) {
//              dif = tmp_ki_rc_sub;
//              for (int l = 0; l < b4_ki->numValidRank1Cut; ++l) {
//                if (!ki->Rank1CutMem[b4_ki->validRank1Cut[l]]) dif += r1c_to_pi[b4_ki->validRank1Cut[l]];
//              }
//              for (int l = 0; l < b4_ki->numValidRank1Cut_multi; ++l) {
//                int tmp_cut = b4_ki->validRank1Cut_multi[l];
//                if (b4_ki->Rank1CutMem_multi[tmp_cut]
//                    > ki->Rank1CutMem_multi[tmp_cut]) {
//                  dif += r1c_multi_to_pi[tmp_cut];
//                }
//              }
//              if (dif > b4_ki->RC) {
//                if_break = true;
//                goto b4_break1;
//              }
//            }
//          }
//        }
//        b4_break1:
//        ki->if_extended = true;
//        if (if_break) continue;//don't have to ++vec_index
//        sig = extendKernel4ArcElimination_last_half_Backward(ki,
//                                                             i,
//                                                             ki->Sum_MainResource,
//                                                             node->AllBackwardBuckets[i][b].BucketArcs,
//                                                             r1c_to_pi,
//                                                             r1c_multi_to_pi);
//        if (sig == 2) {
//          status = 2;
//          goto QUIT;
//        }
//        sig = extendKernel4ArcElimination_last_half_Backward(ki,
//                                                             i,
//                                                             0,
//                                                             node->AllBackwardBuckets[i][b].JumpArcs,
//                                                             r1c_to_pi,
//                                                             r1c_multi_to_pi);
//        if (sig == 2) {
//          status = 2;
//          goto QUIT;
//        }
//      }
//      valid_num = 0;
//    }
//    //test if all labels are extended
//    end = high_resolution_clock::now();
//    eps = duration<double>(end - beg).count();
//    if (eps > CONFIG::HardTimeThresholdInArcElimination_last_half) {
//      status = 1;
//      Rollback = 1;
//      goto QUIT;
//    }
//    for (i = 1; i < Dim; ++i) {
//      if (IfExistExtraLabelsInBackwardSense[i][b].second)
//        goto STILL_EXIST1;
//    }
//    //write data
//    for (i = 1; i < Dim; ++i) {
//      std::stable_sort(LabelArrayInBackwardSense[i][b].first.begin(),
//                       LabelArrayInBackwardSense[i][b].first.begin() + LabelArrayInBackwardSense[i][b].second,
//                       CmpLabelRCLess);
//    }
//    if (b != NumBucketsPerVertex - 1) {
//      for (i = 1; i < Dim; ++i) {
//        if (LabelArrayInBackwardSense[i][b].second) {
//          RC2TillThisBinInBackwardSense[i][b] = min(RC2TillThisBinInBackwardSense[i][b + 1],
//                                                    LabelArrayInBackwardSense[i][b].first[0]->RC);
//        } else {
//          RC2TillThisBinInBackwardSense[i][b] = RC2TillThisBinInBackwardSense[i][b + 1];
//        }
//      }
//    } else {
//      for (i = 1; i < Dim; ++i) {
//        if (LabelArrayInBackwardSense[i][b].second) {
//          RC2TillThisBinInBackwardSense[i][b] = LabelArrayInBackwardSense[i][b].first[0]->RC;
//        } else {
//          RC2TillThisBinInBackwardSense[i][b] = LARGEFLOAT;
//        }
//      }
//    }
//  }
//  populateRC2TillThisBinNRC2Bin(node, 2);
//  QUIT:
//  return status;
//}
//
//int CVRP::backwardConcatenateInArcElimination(const double *r1c_to_pi, const double *r1c_multi_to_pi) {
//  int status = 0;
//  bool if_find;
//  double path_rc;
//  int arr_bj;
//  for (auto &label_list : concatenateLabelsInBackwardCG) {
//    int i = label_list.first.first;
//    int j = label_list.first.second;
//    auto &label_vec = label_list.second;
//    for (auto &pr : label_vec) {
//      auto &ki = pr.first;
//      auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
//      tmp_mainResource = pr.second;
//      auto &tmp_rc = AllLabel[IdxGlo].RC;
//      tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//      if_find = false;
//      arr_bj = int(tmp_mainResource / StepSize);
//      if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_till_this_bin
//        continue;
//
//      if (RC2BinInForwardSense[j][arr_bj] + tmp_rc < OptGap) {//most_negative_rc_in_this_bin
//        //add one more condition for testing capacity
//        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
//        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//          auto &kj = label_arr[vec_index];
//          path_rc = kj->RC + tmp_rc;
//          if (path_rc > OptGap) break;
//
//          if ((ki->PI & kj->PI).any()) continue;
//          if (tmp_mainResource < kj->Sum_MainResource) continue;
//
//          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//              if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//            }
//          }
//
//          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = ki->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut] +
//                  ki->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (ki->Rank1CutMem_multi[tmp_cut] +
//                  kj->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          }
//
//          if (path_rc < OptGap) {
//            if_find = true;
//            goto outside;
//          }
//        }
//      }
//      //bj-1
//      for (--arr_bj; arr_bj >= 0; --arr_bj) {
//        if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_till_this_bin
//          break;
//
//        if (RC2BinInForwardSense[j][arr_bj] + tmp_rc > OptGap)//most_negative_rc_in_this_bin
//          continue;
//
//        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
//
//        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//          auto &kj = label_arr[vec_index];
//          path_rc = kj->RC + tmp_rc;
//          if (path_rc > OptGap) break;
//
//          if ((ki->PI & kj->PI).any()) continue;
//
//          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//              if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//            }
//          }
//          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = ki->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut] +
//                  ki->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (ki->Rank1CutMem_multi[tmp_cut] +
//                  kj->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          }
//          if (path_rc < OptGap) {
//            if_find = true;
//            goto outside;
//          }
//        }
//      }
//
//      outside:
//      if (!if_find) continue;
//      //begin to extend by this direction
//      auto sig = extendKernel4ArcElimination_inner_Backward(ki, i, j, tmp_mainResource, r1c_to_pi, r1c_multi_to_pi);
//      if (sig == 1) continue;
//      if (sig == 2) {
//        status = 2;
//        goto QUIT;
//      }
//    }
//  }
//  QUIT:
//  return status;
//}
//
//void CVRP::eliminateBackwardArcs(BBNODE *node, const double *r1c_to_pi, const double *r1c_multi_to_pi,
//                                 int dim_sq,
//                                 bool *stateBetween2Buckets,
//                                 int *latest_bucket) {
//  memset(stateBetween2Buckets, 0, dim_sq * NumBucketsPerVertex * sizeof(bool));
//  fill_n(latest_bucket, dim_sq, NumBucketsPerVertex);
//
//  double tmp_rc, path_rc, tmp_mainResource;
//  int num_bucket_arcs = 0;
//  for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
//    for (int i = 1; i < Dim; ++i) {
//      for (auto pair : node->AllBackwardBuckets[i][b].JumpArcs) {
//        int j = pair.second;
//        int map = i * Dim + j;
//        int &latest = latest_bucket[map];
//        int old_latest = latest;
//        auto &label_array = LabelArrayInBackwardSense[i][b].first;
//        auto &label_num = LabelArrayInBackwardSense[i][b].second;
//        if (!decreaseMainResourceConsumption(pair.first, tmp_mainResource, i, j)) continue;
//        for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//          auto &ki = label_array[vec_index_i];
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//          //keep updating the latest_bucket
//          int arr_bj = min(int((tmp_mainResource) / StepSize), latest - 1);
//
//          for (int bj = 0; bj <= arr_bj; ++bj) {
//            //first test
//            if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//            for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//              auto &kj = label_arr[vec_index];
//              if ((kj->PI & ki->PI).any()) continue;
//              if (tmp_mainResource < kj->Sum_MainResource) continue;
//              path_rc = tmp_rc + kj->RC;
//              if (path_rc > OptGap) break;
//
//              if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//                  if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//                }
//              } else {
//                for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//                  if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//                }
//              }
//              if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = ki->validRank1Cut_multi[l];
//                  if (kj->Rank1CutMem_multi[tmp_cut] +
//                      ki->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              } else {
//                for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = kj->validRank1Cut_multi[l];
//                  if (ki->Rank1CutMem_multi[tmp_cut] +
//                      kj->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              }
//
//              if (path_rc < OptGap) {
//                latest = bj;
//                goto quit1;
//              }
//            }
//          }
//          quit1:
//          continue;
//        }
//        //update the bins
//        if (latest != old_latest) {
//          int bi = TellWhichBin4ArcEliminationInBackwardSense[map + latest * dim_sq];
//          int map2 = map + b * dim_sq;
//          for (int k = b; k >= bi; --k, map2 -= dim_sq) stateBetween2Buckets[map2] = true;
//        }
//      }
//
//      for (int j : node->AllBackwardBuckets[i][b].BucketArcs) {
//        int map = i * Dim + j;
//        int &latest = latest_bucket[map];
//        int old_latest = latest;
//        auto &label_array = LabelArrayInBackwardSense[i][b].first;
//        auto &label_num = LabelArrayInBackwardSense[i][b].second;
//        for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//          auto *ki = label_array[vec_index_i];
//          if (!decreaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//          //keep updating the latest_bucket
//
//          int arr_bj = min(int((tmp_mainResource) / StepSize), latest - 1);
//
//          for (int bj = 0; bj <= arr_bj; ++bj) {
//            //first test
//            if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//            for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//              auto &kj = label_arr[vec_index];
//              if ((kj->PI & ki->PI).any()) continue;
//
//              if (tmp_mainResource < kj->Sum_MainResource) continue;
//
//              path_rc = tmp_rc + kj->RC;
//              if (path_rc > OptGap) break;
//
//              if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//                  if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//                }
//              } else {
//                for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//                  if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//                }
//              }
//
//              if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = ki->validRank1Cut_multi[l];
//                  if (kj->Rank1CutMem_multi[tmp_cut] +
//                      ki->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              } else {
//                for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = kj->validRank1Cut_multi[l];
//                  if (ki->Rank1CutMem_multi[tmp_cut] +
//                      kj->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              }
//
//              if (path_rc < OptGap) {
//                latest = bj;
//                goto quit2;
//              }
//            }
//          }
//          quit2:
//          continue;
//        }
//        //update the bins
//        if (latest != old_latest) {
//          int bi = TellWhichBin4ArcEliminationInBackwardSense[map + latest * dim_sq];
//          int map2 = map + b * dim_sq;
//          for (int k = b; k >= bi; --k, map2 -= dim_sq) stateBetween2Buckets[map2] = true;
//        }
//      }
//    }
//  }
//
//  vector<int> tmp_vec;
//  tmp_vec.reserve(Dim);
//  int map1 = 0;
//  for (int b = 0; b < NumBucketsPerVertex; ++b, map1 += dim_sq) {
//    int map2 = map1 + Dim;
//    for (int i = 1; i < Dim; ++i, map2 += Dim) {
//      tmp_vec.clear();
//      for (int j : node->AllBackwardBuckets[i][b].BucketArcs) {
//        if (stateBetween2Buckets[map2 + j]) {
//          tmp_vec.emplace_back(j);
//          ++num_bucket_arcs;
//        }
//      }
//      node->AllBackwardBuckets[i][b].BucketArcs = tmp_vec;
//    }
//  }
//
//  cout << "Num of Backward BucketArcs= " << num_bucket_arcs << " prev.= "
//       << double(num_bucket_arcs) / node->NumBackwardBucketArcs * 100 << "%"
//       << " max.= " << double(num_bucket_arcs) / MaxNumBackwardGraphArc * 100 << "%"
//       << endl;
//  node->NumBackwardBucketArcs = num_bucket_arcs;
//}
//
//void CVRP::obtainBackwardJumpArcs(BBNODE *const node, bitset<2> **bitMap) const {
//  int num_jump_arcs = 0;
//  bool if_used;
//
//  for (int i = 1; i < Dim; ++i) {
//    //start to collect data
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      node->AllBackwardBuckets[i][b].JumpArcs.clear();
//      //2 represents the bin is not checked
//      //1 represents the bin is checked
//      //0 represents the bin has arc itself
//      for (int j = 1; j < Dim; ++j) bitMap[j][b] = 2;
//      bitMap[i][b] = 1;
//      for (int j : node->AllBackwardBuckets[i][b].BucketArcs) bitMap[j][b] = 0;
//    }
//    for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
//      for (int j = 1; j < Dim; ++j) {
//        if (bitMap[j][b] == 2) {//需要有jump arc
//          //开始找最近的step个有去j这个地方的bin
//          if_used = false;
//          for (int b4_i = b - 1; b4_i >= 0; --b4_i) {
//            if (bitMap[j][b4_i] == 0) {//发现有去的arc了
//              pair<int, int> map = make_pair((b4_i + 1) * StepSize - TOLERANCE, j);
//              for (int tmp_b = b; tmp_b > b4_i; --tmp_b) {//不取等号< qi
//                bitMap[j][tmp_b] = 1;
//                node->AllBackwardBuckets[i][tmp_b].JumpArcs.emplace_back(map);
//                ++num_jump_arcs;
//              }
//              if_used = true;
//              break;
//            }
//          }
//          if (!if_used) {
//            for (int tmp_b = b; tmp_b >= 0; --tmp_b) bitMap[j][tmp_b] = 1;
//          }
//        }
//      }
//    }
//  }
//
//  node->NumBackwardJumpArcs = num_jump_arcs;
//
//  cout << "Obtain Backward Jump Arcs= " << num_jump_arcs << endl;
//}
//#endif
//
//void CVRP::runLabeling4ArcElimination(BBNODE *const node, const PtrAllR1Cs &ptrAllR1Cs) {
//  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
//  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
//  double dif;
//  bool if_break, if_find, if_keep;
//  double path_rc;
//  int edgemap;
//  OptGap = calculateOptGap(node);
//
//  auto time_beg = chrono::high_resolution_clock::now();
//  auto time_end = time_beg;
//  double time_eps;
//
//  time_beg = chrono::high_resolution_clock::now();
////  if (forwardConcatenateInArcElimination(r1c_to_pi, r1c_multi_to_pi) == 2) goto QUIT;
//  concatenatePhaseInArcElimination<true, true>(r1c_to_pi, r1c_multi_to_pi);
//  if (Rollback == 1 || Rollback == 2) goto QUIT;
//  time_end = chrono::high_resolution_clock::now();
//  time_eps = chrono::duration<double>(time_end - time_beg).count();
//  cout << "Mid forward concatenate= " << time_eps << endl;
//
//  time_beg = chrono::high_resolution_clock::now();
////  if (lastHalfForwardInArcElimination(node, r1c_to_pi, r1c_multi_to_pi)) goto QUIT;
//  runLabeling<true, true, true>(node, ptrAllR1Cs);
//  if (Rollback == 1 || Rollback == 2) goto QUIT;
//  time_end = chrono::high_resolution_clock::now();
//  time_eps = chrono::duration<double>(time_end - time_beg).count();
//  cout << "Last half forward labeling= " << time_eps << endl;
//
//#ifdef SYMMETRY_PROHIBIT
//  time_beg = chrono::high_resolution_clock::now();
//  concatenatePhaseInArcElimination<false, false>(r1c_to_pi, r1c_multi_to_pi);
//  if (Rollback == 1 || Rollback == 2) goto QUIT;
//
////  if (backwardConcatenateInArcElimination(r1c_to_pi, r1c_multi_to_pi) == 2) goto QUIT;
//
//  time_end = chrono::high_resolution_clock::now();
//  time_eps = chrono::duration<double>(time_end - time_beg).count();
//  cout << "Mid Backward concatenate= " << time_eps << endl;
//
//  time_beg = chrono::high_resolution_clock::now();
////  if (lastHalfBackwardInArcElimination(node, r1c_to_pi, r1c_multi_to_pi)) goto QUIT;
//  runLabeling<false, true, true>(node, ptrAllR1Cs);
//  if (Rollback == 1 || Rollback == 2) goto QUIT;
//  time_end = chrono::high_resolution_clock::now();
//  time_eps = chrono::duration<double>(time_end - time_beg).count();
//  cout << "Last half Backward labeling= " << time_eps << endl;
//#endif
//  QUIT:
//  return;
//}
//
