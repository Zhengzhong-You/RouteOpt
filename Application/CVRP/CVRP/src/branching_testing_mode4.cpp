//
// Created by You, Zhengzhong on 5/26/23.
//

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

#ifdef MASTER_VALVE_ML

#if do_SB_mode == 4
void CVRP::do_SB(BBNODE *node, pair<int, int> &info) {
  if (node->if_terminated) return;
  Branch_pair.clear();

  ml.readEnuState(If_in_Enu_State);
  generateModelInPhase1_combinedFashion(node);
  ml.EdgeTmpInfo.clear();

  info = Branch_pair[0];
  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++BranchChoice[info];
  ++BranchTimes;
}
#endif

void CVRP::generateModelInPhase1_combinedFashion(BBNODE *node) {
  InitialScreening(node, true, true, CONFIG::ML_BranchPhase0, CONFIG::Frac4sudoCostBranPhase0);
  getTrainingDataInPhase1(node);
  simulateWriteLPPseudocost(node);

  CGTesting(node, true, true, false);

  for (auto &edge : Branch_Pair_Val) {
    ml.EdgeTmpInfo[edge.first].SB_scores = edge.second;
  }

  ml.writeTrainingLPFile_combinedFashion(false);
}

#define PseudoMark std::string("pseudo-")

void ML3::writeTrainingLPFile_combinedFashion(bool if_debug_features) {
  auto &path = if_in_enu ? enum_lp_path : lp_Output_path;
  ofstream trainingData;
  trainingData.open(path, ios::app);
  vector<pair<int, int>> record;
  record.reserve(EdgeTmpInfo.size());
  for (auto &tmp_info : EdgeTmpInfo) {
    trainingData << tmp_info.second.SB_scores;
    trainingData << " qid:" << QID;
    int cnt = 0;
    for (auto &feature : tmp_info.second.BasicFeatures) {
      trainingData << " " << cnt << ":" << (float) feature.second;
      ++cnt;
    }
    trainingData << endl;
  }

  ++QID;

  trainingData.close();
  if (if_debug_features) {
    printFeatures();
  }
}


#endif

