//
// Created by Zhengzhong You on 7/22/22.
//

#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::eliminateArcs(BBNODE *const node) {
  //used only when cg convergence
  //first we run runLabeling4ArcElimination
  tellIfArcElimination(node);
  if (!final_decision_4_arc_elimination) {
    if_ArcEliminationSucceed = false;
    return;
  }
  final_decision_4_arc_elimination = false;
  if_ArcEliminationTriedButFailed = false;

  double time_labeling, time_elimination, time_obtainjump;

  cout << BIG_PHASE_SEPARATION;
  cout << "run ArcElimination..." << endl;

  populateTellWhichBin4ArcElimination<true>();
#ifdef SYMMETRY_PROHIBIT
  populateTellWhichBin4ArcElimination<false>();
#endif

  PtrAllR1Cs ptrAllR1Cs(node, this);

  auto beg = chrono::high_resolution_clock::now();
  auto end = chrono::high_resolution_clock::now();
  auto eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  Rollback = 0;

  runLabeling4ArcElimination(node, ptrAllR1Cs);

  if (Rollback == 2) {
    cout << "Rollback: " << Rollback << endl;
    if_ArcEliminationSucceed = false;
    if_ArcEliminationTriedButFailed = true;
    goto QUIT;
  } else if (Rollback == 1) {
    cout << "Rollback: " << Rollback << endl;
    if_ArcEliminationSucceed = false;
    if_ArcEliminationTriedButFailed = true;
    //disabled
//    Gap4LastFailDue2DifficultyInArcElimination =
//        min(Gap4LastFailDue2DifficultyInArcElimination, (UB - node->Val) / UB);
    cout << "Reach the hard limit!" << endl;
    goto QUIT;
  }
  safe_Hyperparameter(Rollback)

  end = chrono::high_resolution_clock::now();
  eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  time_labeling = (double) eps.count() * 1e-3;

  beg = chrono::high_resolution_clock::now();
  //eliminateBucketArcs
  eliminateBucketArcs(node, ptrAllR1Cs);

  end = chrono::high_resolution_clock::now();
  eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  time_elimination = (double) eps.count() * 1e-3;

  beg = chrono::high_resolution_clock::now();
  //obtainJumpArcs
  obtainJumpArcs(node);
  end = chrono::high_resolution_clock::now();
  eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  time_obtainjump = (double) eps.count() * 1e-3;

  cout << "time summary= " << endl;
  cout << "runLabeling4ArcElimination= " << time_labeling << endl;
  cout << "eliminateBucketArcs= " << time_elimination << endl;
  cout << "obtainJumpArcs= " << time_obtainjump << endl;
  if_ArcEliminationSucceed = true;
  QUIT:
//  cout << "这里是test！" << endl;
//  for (int i = 1; i < Dim; ++i) {
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      vector<bool> test(Dim, false);
//      for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//        test[j] = true;
//      }
//
//      for (auto j : node->AllForwardBuckets[i][b].JumpArcs) {
//        if (test[j.second] == true) {
//          cout << "error: jump arcs and bucket arcs overlap! forward!" << endl;
//          exit(1);
//        }
//      }
//      std::fill(test.begin(), test.end(), false);
//      for (int j : node->AllBackwardBuckets[i][b].BucketArcs) {
//        test[j] = true;
//      }
//      for (auto j : node->AllBackwardBuckets[i][b].JumpArcs) {
//        if (test[j.second] == true) {
//          cout << "error: jump arcs and bucket arcs overlap! backward!" << endl;
//          exit(1);
//        }
//      }
//    }
//  }

#ifdef find_missed_solution
  bool if_opt;
  vector<int> data;
  findWhySolutionDisappear(node, this, data, if_opt, true);
  if (data.empty()) {

  } else if (true || !if_use_heur_enumeration) {
    cout << "这里有检查！记得删除这个多余的检查！" << endl;
    vector<vector<int>> sol_ = {

        {0, 161, 71, 141, 182, 0},
        {0, 123, 77, 121, 73, 7, 128, 193, 0},
        {0, 167, 79, 100, 145, 0},
        {0, 68, 90, 105, 53, 181, 0},
        {0, 56, 151, 164, 80, 44, 143, 140, 0},
        {0, 34, 129, 78, 84, 126, 67, 0},
        {0, 106, 24, 183, 195, 8, 0},
        {0, 43, 13, 70, 2, 39, 0},
        {0, 49, 98, 5, 115, 52, 0},
        {0, 12, 180, 185, 134, 163, 87, 188, 0},
        {0, 21, 118, 60, 66, 104, 35, 64, 168, 0},
        {0, 107, 137, 32, 194, 36, 0},
        {0, 22, 184, 171, 114, 142, 166, 0},
        {0, 127, 65, 92, 97, 91, 103, 148, 31, 191, 0},
        {0, 42, 197, 122, 99, 28, 0},
        {0, 47, 138, 116, 95, 86, 146, 0},
        {0, 18, 169, 63, 26, 173, 196, 0},
        {0, 11, 179, 3, 62, 9, 0},
        {0, 6, 50, 10, 29, 102, 132, 0},
        {0, 23, 165, 74, 110, 33, 192, 136, 0},
        {0, 45, 108, 156, 153, 120, 0},
        {0, 61, 170, 119, 55, 0},
        {0, 93, 1, 178, 94, 69, 147, 41, 198, 0},
        {0, 112, 4, 113, 111, 48, 51, 17, 0},
        {0, 38, 46, 150, 162, 37, 85, 0},
        {0, 133, 177, 81, 72, 190, 176, 0},
        {0, 59, 57, 131, 109, 96, 0},
        {0, 155, 172, 130, 83, 25, 0},
        {0, 175, 159, 89, 144, 117, 54, 0},
        {0, 20, 139, 199, 14, 160, 186, 0},
        {0, 154, 189, 135, 19, 30, 75, 15, 0},
        {0, 101, 174, 152, 27, 187, 0},
        {0, 125, 124, 76, 157, 40, 0},
        {0, 16, 149, 88, 58, 82, 158, 0}

    };
    for (auto &s_ : sol_) {
      double res = 0;
      int b4 = s_[1];
      increaseMainResourceConsumption(res, res, 0, b4);
      for (auto i : s_) {
        if (i == 0 || i == b4) continue;
        int bin = int(res / StepSize);
        auto if_find = std::find(node->AllForwardBuckets[b4][bin].BucketArcs.begin(),
                                 node->AllForwardBuckets[b4][bin].BucketArcs.end(), i);
        if (if_find == node->AllForwardBuckets[b4][bin].BucketArcs.end()) {
          cout << "error! in forward" << endl;
          cout << "i= " << b4 << " bin= " << bin << " j= " << i << endl;
          cout << "res/StepSize= " << res / StepSize << endl;
          cout << "res= " << res << endl;
          increaseMainResourceConsumption(res, res, i, 0);
          cout << "res= " << res << endl;
          exit(1);
        }
        increaseMainResourceConsumption(res, res, b4, i);
        b4 = i;
      }
    }
#ifdef SYMMETRY_PROHIBIT
    for (auto &s_ : sol_) {
      double res = MaxMainResource;
      int b4 = s_[s_.size() - 2];
      decreaseMainResourceConsumption(res, res, 0, b4);
      for (int it = s_.size() - 3; it >= 1; --it) {
        int i = s_[it];
        int bin = int(res / StepSize);
        auto if_find = std::find(node->AllBackwardBuckets[b4][bin].BucketArcs.begin(),
                                 node->AllBackwardBuckets[b4][bin].BucketArcs.end(), i);
        if (if_find == node->AllBackwardBuckets[b4][bin].BucketArcs.end()) {
          cout << "error! in backward" << endl;
          cout << "i= " << b4 << " bin= " << bin << " j= " << i << endl;
          exit(1);
        }
        decreaseMainResourceConsumption(res, res, b4, i);
        b4 = i;
      }
    }
#endif

  }
#endif
//  TellWhichBin4ArcEliminationInForwardSense.clear();
//#ifdef SYMMETRY_PROHIBIT
//  TellWhichBin4ArcEliminationInBackwardSense.clear();
//#endif
#ifdef draw_residual_bucket_graph
  for (int i = 1; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      check_ij_bucket(node, i, j);
      check_ij_bucket(node, j, i);
    }
  }
  cout << "check bucket done" << endl;
  cout << "remain arc: " << endl;
  int num_bucket = 0, num_jump = 0;
  for (int i = 1; i < Dim; ++i) {
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      num_bucket += (int) node->AllForwardBuckets[i][b].BucketArcs.size();
      num_jump += (int) node->AllForwardBuckets[i][b].JumpArcs.size();
    }
  }
  cout << "num_bucket/max= " << double(num_bucket) / double(MaxNumForwardGraphArc) << endl;
  cout << "num_bucket/pri= " << double(num_bucket) / double(node->NumForwardBucketArcs) << endl;
  cout << "num_jump/max= " << double(num_jump) << endl;
#endif
  cout << BIG_PHASE_SEPARATION;
}

void CVRP::runLabeling4ArcElimination(BBNODE *const node, const PtrAllR1Cs &ptrAllR1Cs) {
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  double dif;
  bool if_break, if_find, if_keep;
  double path_rc;
  int edgemap;
  OptGap = calculateOptGap(node);

  auto time_beg = chrono::high_resolution_clock::now();
  auto time_end = time_beg;
  double time_eps;

#ifdef SYMMETRY_PROHIBIT
  concatenatePhaseInArcElimination<true, false>(r1c_to_pi, r1c_multi_to_pi);
  if (Rollback == 1 || Rollback == 2) goto QUIT;
  runLabeling<true, true, false>(node, ptrAllR1Cs);
  if (Rollback == 1 || Rollback == 2) goto QUIT;
  concatenatePhaseInArcElimination<false, false>(r1c_to_pi, r1c_multi_to_pi);
  if (Rollback == 1 || Rollback == 2) goto QUIT;
  runLabeling<false, true, false>(node, ptrAllR1Cs);
  if (Rollback == 1 || Rollback == 2) goto QUIT;
#else
  concatenatePhaseInArcElimination<true, true>(r1c_to_pi, r1c_multi_to_pi);
  if (Rollback == 1 || Rollback == 2) goto QUIT;
  runLabeling<true, true, true>(node, ptrAllR1Cs);
  if (Rollback == 1 || Rollback == 2) goto QUIT;
#endif

  time_end = chrono::high_resolution_clock::now();
  time_eps = chrono::duration<double>(time_end - time_beg).count();
  cout << "Last half forward labeling= " << time_eps << endl;

  QUIT:
  return;
}

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
//  if (forwardConcatenateInArcElimination(r1c_to_pi, r1c_multi_to_pi) == 2) goto QUIT;
//  time_end = chrono::high_resolution_clock::now();
//  time_eps = chrono::duration<double>(time_end - time_beg).count();
//  cout << "Mid forward concatenate= " << time_eps << endl;
//
//  time_beg = chrono::high_resolution_clock::now();
//  if (lastHalfForwardInArcElimination(node, r1c_to_pi, r1c_multi_to_pi)) goto QUIT;
//  time_end = chrono::high_resolution_clock::now();
//  time_eps = chrono::duration<double>(time_end - time_beg).count();
//  cout << "Last half forward labeling= " << time_eps << endl;
//
//  QUIT:
//  return;
//}

void CVRP::eliminateBucketArcs(BBNODE *const node, const PtrAllR1Cs &ptrAllR1Cs) {
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  int dim_sq = Dim * Dim;
  //false no，true yes
  auto stateBetween2Buckets = new bool[dim_sq * NumBucketsPerVertex];
  auto latest_bucket = new int[dim_sq];
  OptGap = calculateOptGap(node);

//  eliminatebuketArcs(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);

#ifdef SYMMETRY_PROHIBIT
  eliminatebuketArcs<true, false>(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);
  eliminatebuketArcs<false, false>(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);
#else
  eliminatebuketArcs<true, true>(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);
#endif


  /**
   * write out the residual graph
   */

//  cout << "here is writing the residual graph!" << endl;
//  ofstream residual_graph("residual_graph.txt", ios::out);
//  for (int i = 0; i < Dim; ++i) {
//    unordered_set<int> arcs;
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//        if (i == j)continue;
//        arcs.emplace(j);
//      }
//    }
//    residual_graph << i << "= ";
//    for (int j : arcs) {
//      residual_graph << j << " ";
//    }
//    residual_graph << endl;
//  }

/**
 * read in the residual graph
 */

//  ifstream residual_graph("bc_residual_graph.txt", ios::in);
//  string line;
//  while (getline(residual_graph, line)) {
//    if (line.empty()) continue;
//    stringstream ss(line);
//    int i;
//    ss >> i;
//    if (i == 0) continue;
//    int j;
//    unordered_set<int> in_cols;
//    while (!ss.eof()) {
//      ss >> j;
//      if (ss.fail()) {
//        ss.clear();
//        string temp;
//        ss >> temp;
//        continue;
//      }
//      in_cols.emplace(j);
//    }
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      for (auto it = node->AllForwardBuckets[i][b].BucketArcs.begin();
//           it != node->AllForwardBuckets[i][b].BucketArcs.end();) {
//        if (in_cols.find(*it) == in_cols.end()) {
//          it = node->AllForwardBuckets[i][b].BucketArcs.erase(it);
//        } else {
//          ++it;
//        }
//      }
//      for (auto it = node->AllForwardBuckets[i][b].JumpArcs.begin();
//           it != node->AllForwardBuckets[i][b].JumpArcs.end();) {
//        if (in_cols.find(it->second) == in_cols.end()) {
//          it = node->AllForwardBuckets[i][b].JumpArcs.erase(it);
//        } else {
//          ++it;
//        }
//      }
//    }
//  }

  delete[]stateBetween2Buckets;
  delete[]latest_bucket;
}

//void CVRP::eliminatebuketArcs(BBNODE *node, const double *r1c_to_pi,
//                              const double *r1c_multi_to_pi,
//                              int dim_sq,
//                              bool *stateBetween2Buckets,
//                              int *latest_bucket) {
//  memset(stateBetween2Buckets, 0, dim_sq * NumBucketsPerVertex * sizeof(bool));
//  fill_n(latest_bucket, dim_sq, NumBucketsPerVertex);
////  if constexpr (if_symmetry) {
////
////  } else {
////    memset(latest_bucket, -1, dim_sq * sizeof(int));
////  }
//  double tmp_rc, path_rc, tmp_mainResource;
//  int num_bucket_arcs = 0;
//  for (int b = 0; b < NumBucketsPerVertex; ++b) {
//    for (int i = 1; i < Dim; ++i) {
//      for (auto pair : node->AllForwardBuckets[i][b].JumpArcs) {
//        int j = pair.second;
//        int map = i * Dim + j;
//        int &latest = latest_bucket[map];
//        int old_latest = latest;
//        auto &label_array = LabelArrayInForwardSense[i][b].first;
//        auto &label_num = LabelArrayInForwardSense[i][b].second;
//        if (!increaseMainResourceConsumption(pair.first, tmp_mainResource, i, j)) continue;
//        for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//          auto &ki = label_array[vec_index_i];
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//          //keep updating the latest_bucket
//#ifdef SYMMETRY_PROHIBIT
//          int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
//          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
//            //first test
//            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
//#else
//          int arr_bj = min(int((MaxMainResource - tmp_mainResource) / StepSize), latest - 1);
//          for (int bj = 0; bj <= arr_bj; ++bj) {
//            //first test
//            if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//#endif
//            for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//              auto &kj = label_arr[vec_index];
//              if ((kj->PI & ki->PI).any()) continue;
//#ifdef SYMMETRY_PROHIBIT
//              if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//              if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
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
//                goto quit1;
//              }
//            }
//          }
//          quit1:
//          continue;
//        }
//        //update the bins
//        if (latest != old_latest) {
//#ifdef SYMMETRY_PROHIBIT
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + latest * dim_sq];
//#else
//          int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
//          if (chg_latest < 0) chg_latest = 0;
//          else if (chg_latest >= NumBucketsPerVertex) chg_latest = NumBucketsPerVertex - 1;
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq];
//#endif
//          int map2 = map + b * dim_sq;
//          for (int k = b; k <= bi; ++k, map2 += dim_sq) stateBetween2Buckets[map2] = true;
//        }
//      }
//      for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//        int map = i * Dim + j;
//        int &latest = latest_bucket[map];
//        int old_latest = latest;
//        auto &label_array = LabelArrayInForwardSense[i][b].first;
//        auto &label_num = LabelArrayInForwardSense[i][b].second;
//        for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//          auto *ki = label_array[vec_index_i];
//          if (!increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//          //keep updating the latest_bucket
//#ifdef SYMMETRY_PROHIBIT
//          int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
//          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
//            //first test
//            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
//#else
//          int arr_bj = min(int((MaxMainResource - tmp_mainResource) / StepSize), latest - 1);
//          for (int bj = 0; bj <= arr_bj; ++bj) {
//            //first test
//            if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//#endif
//            for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//              auto &kj = label_arr[vec_index];
//              if ((kj->PI & ki->PI).any()) continue;
//#ifdef SYMMETRY_PROHIBIT
//              if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//              if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
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
//#ifdef SYMMETRY_PROHIBIT
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + latest * dim_sq];
//#else
//          int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
//          if (chg_latest < 0) chg_latest = 0;
//          else if (chg_latest >= NumBucketsPerVertex) chg_latest = NumBucketsPerVertex - 1;
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq];
//#endif
//          int map2 = map + b * dim_sq;
//          for (int k = b; k <= bi; ++k, map2 += dim_sq) stateBetween2Buckets[map2] = true;
//        }
//      }
//    }
//  }
//#ifdef find_lost_arcs
//  int i = 37, b = 10, j = 32, map = i * Dim + j;
//  int latest = -1;
//  auto &label_array = LabelArrayInForwardSense[i][b].first;
//  auto &label_num = LabelArrayInForwardSense[i][b].second;
//  for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//    auto *ki = label_array[vec_index_i];
//    if (!increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
//    tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//    //keep updating the latest_bucket
//#ifdef SYMMETRY_PROHIBIT
//    int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
//          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
//            //first test
//            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
//#else
//    int arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//    for (int bj = 0; bj <= arr_bj; ++bj) {
//      //first test
//      if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//      //real test
//      auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//      auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//#endif
//      for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//        auto &kj = label_arr[vec_index];
//        if ((kj->PI & ki->PI).any()) continue;
//#ifdef SYMMETRY_PROHIBIT
//        if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//        if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
//        path_rc = tmp_rc + kj->RC;
//        if (path_rc > OptGap) break;
//
//        if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//          for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//            if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//          }
//        } else {
//          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//            if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//          }
//        }
//
//        if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//          for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = ki->validRank1Cut_multi[l];
//            if (kj->Rank1CutMem_multi[tmp_cut] +
//                ki->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        } else {
//          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = kj->validRank1Cut_multi[l];
//            if (ki->Rank1CutMem_multi[tmp_cut] +
//                kj->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        }
//
//        if (path_rc < OptGap) {
//          cout << "arc should exists!" << endl;
//          cout << "bin= " << bj << endl;
//          latest = bj;
//          break;
//        }
//      }
//    }
//  }
//  if (latest != -1) {
//    int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
//    cout << "tell= " << TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq] << endl;
//  }
//#endif
//
//  vector<int> tmp_vec;
//  tmp_vec.reserve(Dim);
//  int map1 = 0;
//  for (int b = 0; b < NumBucketsPerVertex; ++b, map1 += dim_sq) {
//    int map2 = map1 + Dim;
//    for (int i = 1; i < Dim; ++i, map2 += Dim) {
//      tmp_vec.clear();
//      for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//        if (stateBetween2Buckets[map2 + j]) {
//          tmp_vec.emplace_back(j);
//          ++num_bucket_arcs;
//        }
//      }
//      node->AllForwardBuckets[i][b].BucketArcs = tmp_vec;
//    }
//  }
//
//  cout << "Num of Forward BucketArcs= " << num_bucket_arcs << " prev.= "
//       << double(num_bucket_arcs) / node->NumForwardBucketArcs * 100 << "%"
//       << " max.= " << double(num_bucket_arcs) / MaxNumForwardGraphArc * 100 << "%"
//       << endl;
//  node->NumForwardBucketArcs = num_bucket_arcs;
//}

//void CVRP::eliminatebuketArcs(BBNODE *node, const double *r1c_to_pi,
//                              const double *r1c_multi_to_pi,
//                              int dim_sq,
//                              bool *stateBetween2Buckets,
//                              int *latest_bucket) {
//  memset(stateBetween2Buckets, 0, dim_sq * NumBucketsPerVertex * sizeof(bool));
//#ifdef SYMMETRY_PROHIBIT
//  memset(latest_bucket, -1, dim_sq * sizeof(int));
//#else
//  fill_n(latest_bucket, dim_sq, NumBucketsPerVertex);
//#endif
//  double tmp_rc, path_rc, tmp_mainResource;
//  int num_bucket_arcs = 0;
//  for (int b = 0; b < NumBucketsPerVertex; ++b) {
//    for (int i = 1; i < Dim; ++i) {
//      for (auto pair : node->AllForwardBuckets[i][b].JumpArcs) {
//        int j = pair.second;
//        int map = i * Dim + j;
//        int &latest = latest_bucket[map];
//        int old_latest = latest;
//        auto &label_array = LabelArrayInForwardSense[i][b].first;
//        auto &label_num = LabelArrayInForwardSense[i][b].second;
//        if (!increaseMainResourceConsumption(pair.first, tmp_mainResource, i, j)) continue;
//        for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//          auto &ki = label_array[vec_index_i];
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//          //keep updating the latest_bucket
//#ifdef SYMMETRY_PROHIBIT
//          int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
//          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
//            //first test
//            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
//#else
//            int arr_bj = min(int((MaxMainResource - tmp_mainResource) / StepSize), latest - 1);
//            for (int bj = 0; bj <= arr_bj; ++bj) {
//              //first test
//              if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//              //real test
//              auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//              auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//#endif
//            for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//              auto &kj = label_arr[vec_index];
//              if ((kj->PI & ki->PI).any()) continue;
//#ifdef SYMMETRY_PROHIBIT
//              if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//              if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
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
//                goto quit1;
//              }
//            }
//          }
//          quit1:
//          continue;
//        }
//        //update the bins
//        if (latest != old_latest) {
//#ifdef SYMMETRY_PROHIBIT
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + latest * dim_sq];
//#else
//          int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
//          if (chg_latest < 0) chg_latest = 0;
//          else if (chg_latest >= NumBucketsPerVertex) chg_latest = NumBucketsPerVertex - 1;
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq];
//#endif
//          int map2 = map + b * dim_sq;
//          for (int k = b; k <= bi; ++k, map2 += dim_sq) stateBetween2Buckets[map2] = true;
//        }
//      }
//      for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//        int map = i * Dim + j;
//        int &latest = latest_bucket[map];
//        int old_latest = latest;
//        auto &label_array = LabelArrayInForwardSense[i][b].first;
//        auto &label_num = LabelArrayInForwardSense[i][b].second;
//        for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//          auto *ki = label_array[vec_index_i];
//          if (!increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//          //keep updating the latest_bucket
//#ifdef SYMMETRY_PROHIBIT
//          int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
//          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
//            //first test
//            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
//#else
//            int arr_bj = min(int((MaxMainResource - tmp_mainResource) / StepSize), latest - 1);
//            for (int bj = 0; bj <= arr_bj; ++bj) {
//              //first test
//              if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//              //real test
//              auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//              auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//#endif
//            for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//              auto &kj = label_arr[vec_index];
//              if ((kj->PI & ki->PI).any()) continue;
//#ifdef SYMMETRY_PROHIBIT
//              if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//              if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
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
//#ifdef SYMMETRY_PROHIBIT
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + latest * dim_sq];
//#else
//          int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
//          if (chg_latest < 0) chg_latest = 0;
//          else if (chg_latest >= NumBucketsPerVertex) chg_latest = NumBucketsPerVertex - 1;
//          int bi = TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq];
//#endif
//          int map2 = map + b * dim_sq;
//          for (int k = b; k <= bi; ++k, map2 += dim_sq) stateBetween2Buckets[map2] = true;
//        }
//      }
//    }
//  }
//#ifdef find_lost_arcs
//  int i = 37, b = 10, j = 32, map = i * Dim + j;
//  int latest = -1;
//  auto &label_array = LabelArrayInForwardSense[i][b].first;
//  auto &label_num = LabelArrayInForwardSense[i][b].second;
//  for (int vec_index_i = 0; vec_index_i < label_num; ++vec_index_i) {
//    auto *ki = label_array[vec_index_i];
//    if (!increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
//    tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//    //keep updating the latest_bucket
//#ifdef SYMMETRY_PROHIBIT
//    int arr_bj = max(int((tmp_mainResource) / StepSize), latest + 1);
//          for (int bj = NumBucketsPerVertex - 1; bj >= arr_bj; --bj) {
//            //first test
//            if (tmp_rc + RC2BinInBackwardSense[j][bj] > OptGap)continue;
//            //real test
//            auto &label_arr = LabelArrayInBackwardSense[j][bj].first;
//            auto &label_valid_num = LabelArrayInBackwardSense[j][bj].second;
//#else
//    int arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//    for (int bj = 0; bj <= arr_bj; ++bj) {
//      //first test
//      if (tmp_rc + RC2BinInForwardSense[j][bj] > OptGap)continue;
//      //real test
//      auto &label_arr = LabelArrayInForwardSense[j][bj].first;
//      auto &label_valid_num = LabelArrayInForwardSense[j][bj].second;
//#endif
//      for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//        auto &kj = label_arr[vec_index];
//        if ((kj->PI & ki->PI).any()) continue;
//#ifdef SYMMETRY_PROHIBIT
//        if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//        if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
//        path_rc = tmp_rc + kj->RC;
//        if (path_rc > OptGap) break;
//
//        if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//          for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//            if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//          }
//        } else {
//          for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//            if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//          }
//        }
//
//        if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//          for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = ki->validRank1Cut_multi[l];
//            if (kj->Rank1CutMem_multi[tmp_cut] +
//                ki->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        } else {
//          for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//            int tmp_cut = kj->validRank1Cut_multi[l];
//            if (ki->Rank1CutMem_multi[tmp_cut] +
//                kj->Rank1CutMem_multi[tmp_cut]
//                >= R1C_multi_denominator_InCG[tmp_cut]
//                )
//              path_rc -= r1c_multi_to_pi[tmp_cut];
//          }
//        }
//
//        if (path_rc < OptGap) {
//          cout << "arc should exists!" << endl;
//          cout << "bin= " << bj << endl;
//          latest = bj;
//          break;
//        }
//      }
//    }
//  }
//  if (latest != -1) {
//    int chg_latest = int((MaxMainResource - latest * StepSize) / StepSize);
//    cout << "tell= " << TellWhichBin4ArcEliminationInForwardSense[map + chg_latest * dim_sq] << endl;
//  }
//#endif
//
//  vector<int> tmp_vec;
//  tmp_vec.reserve(Dim);
//  int map1 = 0;
//  for (int b = 0; b < NumBucketsPerVertex; ++b, map1 += dim_sq) {
//    int map2 = map1 + Dim;
//    for (int i = 1; i < Dim; ++i, map2 += Dim) {
//      tmp_vec.clear();
//      for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//        if (stateBetween2Buckets[map2 + j]) {
//          tmp_vec.emplace_back(j);
//          ++num_bucket_arcs;
//        }
//      }
//      node->AllForwardBuckets[i][b].BucketArcs = tmp_vec;
//    }
//  }
//
//  cout << "Num of Forward BucketArcs= " << num_bucket_arcs << " prev.= "
//       << double(num_bucket_arcs) / node->NumForwardBucketArcs * 100 << "%"
//       << " max.= " << double(num_bucket_arcs) / MaxNumForwardGraphArc * 100 << "%"
//       << endl;
//  node->NumForwardBucketArcs = num_bucket_arcs;
//}

void CVRP::obtainJumpArcs(BBNODE *const node) const {
  auto tmp = new bitset<2> *[Dim];
  for (int i = 0; i < Dim; ++i) tmp[i] = new bitset<2>[NumBucketsPerVertex];
  obtainjumpArcs<true>(node, tmp);
#ifdef SYMMETRY_PROHIBIT
  obtainjumpArcs<false>(node, tmp);
#endif
  for (int i = 0; i < Dim; ++i) {
    delete[]tmp[i];
  }
  delete[]tmp;
}

void CVRP::check_ij_bucket(BBNODE *node, int i, int j) const {
//  return;
  cout << "from i = " << i << " to j = " << j << endl;
  int num_buckets = 0, num_jump = 0;
  for (int bin = 0; bin < NumBucketsPerVertex; ++bin) {
    auto if_find = std::find(node->AllForwardBuckets[i][bin].BucketArcs.begin(),
                             node->AllForwardBuckets[i][bin].BucketArcs.end(), j);
    if (if_find == node->AllForwardBuckets[i][bin].BucketArcs.end()) {
      //jump arc is a pair, compare the second
      auto iff = std::find_if(node->AllForwardBuckets[i][bin].JumpArcs.begin(),
                              node->AllForwardBuckets[i][bin].JumpArcs.end(),
                              [&](const pair<double, int> &p) { return p.second == j; });
      if (iff == node->AllForwardBuckets[i][bin].JumpArcs.end()) {
        continue;
        cout << "(" << bin;
        cout << " n) ";
      } else {
        cout << "(" << bin;
        cout << " j) ";
        ++num_jump;
      }
    } else {
      cout << "(" << bin;
      cout << " B) ";
      ++num_buckets;
    }
  }
  cout << endl;
//  {
//    cout << "this is the test!" << endl;
//    if ((double) num_buckets / num_jump < 0.5) {
//      cout << " num_buckets/num_jump  < 0.2! " << endl;
//      for (int bin = 0; bin < NumBucketsPerVertex; ++bin) {
//        auto if_find = std::find(node->AllForwardBuckets[i][bin].BucketArcs.begin(),
//                                 node->AllForwardBuckets[i][bin].BucketArcs.end(), j);
//        if (if_find == node->AllForwardBuckets[i][bin].BucketArcs.end()) {
//          //jump arc is a pair, compare the second
//          auto iff = std::find_if(node->AllForwardBuckets[i][bin].JumpArcs.begin(),
//                                  node->AllForwardBuckets[i][bin].JumpArcs.end(),
//                                  [&](const pair<double, int> &p) { return p.second == j; });
//          if (iff == node->AllForwardBuckets[i][bin].JumpArcs.end()) {
//            continue;
//            cout << "(" << bin;
//            cout << " n)";
//          } else {
//            node->AllForwardBuckets[i][bin].JumpArcs.erase(iff);
//          }
//        } else {
//          node->AllForwardBuckets[i][bin].BucketArcs.erase(if_find);
//        }
//      }
//      for (int bin = 0; bin < NumBucketsPerVertex; ++bin) {
//        auto if_find = std::find(node->AllForwardBuckets[j][bin].BucketArcs.begin(),
//                                 node->AllForwardBuckets[j][bin].BucketArcs.end(), i);
//        if (if_find == node->AllForwardBuckets[j][bin].BucketArcs.end()) {
//          //jump arc is a pair, compare the second
//          auto iff = std::find_if(node->AllForwardBuckets[j][bin].JumpArcs.begin(),
//                                  node->AllForwardBuckets[j][bin].JumpArcs.end(),
//                                  [&](const pair<double, int> &p) { return p.second == i; });
//          if (iff == node->AllForwardBuckets[j][bin].JumpArcs.end()) {
//            continue;
//            cout << "(" << bin;
//            cout << " n)";
//          } else {
//            node->AllForwardBuckets[j][bin].JumpArcs.erase(iff);
//          }
//        } else {
//          node->AllForwardBuckets[j][bin].BucketArcs.erase(if_find);
//        }
//      }
//    }
//  }
}
