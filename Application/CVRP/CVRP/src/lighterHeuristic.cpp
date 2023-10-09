//
// Created by Zhengzhong You on 7/23/22.
//

#include "CVRP.hpp"

using namespace std;

template<typename T>
void CVRP::extendKernel4LightHeur(LABEL *&ki,
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
    if (ki->PI[j]) continue;
    auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;

    if (!increaseMainResourceConsumption(res, tmp_mainResource, i, j)) continue;

    int bj = int(tmp_mainResource / StepSize);
    auto &labelList_j = LabelArrayInForwardSense[j][bj].first;
    auto &valid_num = LabelArrayInForwardSense[j][bj].second;
    auto &tmp_rc = AllLabel[IdxGlo].RC;
    auto &tmp_PI = AllLabel[IdxGlo].PI;
    auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
    tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
    tmp_Rank1CutMem = ki->Rank1CutMem;
    auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
    copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);

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
      if (tmp_rc < kj->RC) {
        kj->if_extended = true;
      } else {
        continue;
      }
    }

    tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);

    for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;

    tmp_PI = (ki->PI) & (NGMem4Vertex[j]);
    tmp_PI.set(j);

    labelList_j[0] = AllLabel + IdxGlo;

    valid_num = 1;

    //Seq
    AllLabel[IdxGlo].PLabel = ki;
    //EndVertex
    AllLabel[IdxGlo].EndVertex = j;
    //if_extended
    AllLabel[IdxGlo].if_extended = false;
    //bucket
    auto &bucket = IfExistExtraLabelsInForwardSense[j][bj];
    bucket.first[0] = AllLabel + IdxGlo;
    bucket.second = 1;
    //state_b4
    double path_rc = tmp_rc + ChgCostMat4Vertex[j][0];
    if (if_return_to_depot) {
      addPathByRC(path_rc, AllLabel + IdxGlo, nullptr, CONFIG::MaxNumRoutesInLighterHeur);
    }
    ++IdxGlo;
  }
}

int CVRP::runLighterHeurLabeling(BBNODE *const node) {

  int index = NumCol - 1;

  PtrAllR1Cs ptrAllR1Cs(node, this);

  PriorPoolBeg4Pricing = PoolBeg4Pricing;

  initializeLabels(node, 1, true, {true, 1, true});

  // 开始遍历
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  LABEL *ki, *b4_ki;
  bool if_break;

  unordered_set<int> depot_set;
  depot_set.reserve(Dim);
#ifdef SYMMETRY_PROHIBIT
  for (auto j : node->AllBackwardBuckets[0][0].BucketArcs) depot_set.emplace(j);
#else
  for (auto j : node->AllForwardBuckets[0][0].BucketArcs) depot_set.emplace(j);
#endif

  //only one label are kept in one bucket
  for (int b = 0; b < NumBucketsPerVertex; ++b) {
    int i = 1;
    STILL_EXIST:
    for (; i < Dim; ++i) {
      auto &un_extend_num = IfExistExtraLabelsInForwardSense[i][b].second;
      if (!un_extend_num) continue;
      auto &label_array = IfExistExtraLabelsInForwardSense[i][b].first;
      ki = label_array[0];
      if (ki->if_extended) {
        un_extend_num = 0;
        continue;
      }
      if_break = false;
      for (int b4_b = b - 1; b4_b >= 0; --b4_b) {
        if (RC2TillThisBinInForwardSense[i][b4_b] > ki->RC) break;
        auto &b4_valid_num = LabelArrayInForwardSense[i][b4_b].second;
        if (!b4_valid_num) continue;
        auto &b4_label_list = LabelArrayInForwardSense[i][b4_b].first;
        b4_ki = b4_label_list[0];
        if (ki->RC > b4_ki->RC) {
          if_break = true;
          goto b4_break;
        }
      }
      b4_break:
      ki->if_extended = true;
      if (if_break) continue;
      bool if_return_to_depot = depot_set.find(i) != depot_set.end();
      extendKernel4LightHeur(ki,
                             i,
                             ki->Sum_MainResource,
                             node->AllForwardBuckets[i][b].BucketArcs,
                             r1c_to_pi,
                             r1c_multi_to_pi, if_return_to_depot);
      extendKernel4LightHeur(ki,
                             i,
                             0,
                             node->AllForwardBuckets[i][b].JumpArcs,
                             r1c_to_pi,
                             r1c_multi_to_pi,
                             if_return_to_depot);
      un_extend_num = 0;
    }
    //test if all labels are extended
    for (i = 1; i < Dim; ++i) {
      if (IfExistExtraLabelsInForwardSense[i][b].second)
        goto STILL_EXIST;
    }
    //write data
    for (i = 1; i < Dim; ++i) {
      std::stable_sort(LabelArrayInForwardSense[i][b].first.begin(),
                       LabelArrayInForwardSense[i][b].first.begin() + LabelArrayInForwardSense[i][b].second,
                       CmpLabelRCLess);
    }
    if (b) {
      for (i = 1; i < Dim; ++i) {
        if (LabelArrayInForwardSense[i][b].second) {
          RC2TillThisBinInForwardSense[i][b] = min(RC2TillThisBinInForwardSense[i][b - 1],
                                                   LabelArrayInForwardSense[i][b].first[0]->RC);
        } else {
          RC2TillThisBinInForwardSense[i][b] = RC2TillThisBinInForwardSense[i][b - 1];
        }
      }
    } else {
      for (i = 1; i < Dim; ++i) {
        if (LabelArrayInForwardSense[i][b].second) {
          RC2TillThisBinInForwardSense[i][b] = LabelArrayInForwardSense[i][b].first[0]->RC;
        } else {
          RC2TillThisBinInForwardSense[i][b] = LARGEFLOAT;
        }
      }
    }
  }

  writeColsInPricingPool(node, index);

  return index - NumCol + 1;
}

int CVRP::generateColsByLighterHeur(BBNODE *const node) {

  int ccnt = runLighterHeurLabeling(node);

  if (!ccnt) return 0;

  addCols(node, ccnt);

  return ccnt;
}

