//
// Created by You, Zhengzhong on 5/26/23.
//

#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

#ifdef MASTER_VALVE_ML
#if do_SB_mode == 3

void CVRP::do_SB(BBNODE *node, pair<int, int> &info) {
  /**
   * test if the initial screening is sufficient
   */
  if (node->if_terminated) return;
  Branch_pair.clear();
#ifdef testAcc
  ml.readEnuState(If_in_Enu_State);
  testuseML(node);
  ml.EdgeTmpInfo.clear();
#else
  testInitialScreening(node);
#endif
  info = Branch_pair[0];
  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++BranchChoice[info];
  ++BranchTimes;
}

#endif

void CVRP::testuseML(BBNODE *node) {
  InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);

  auto tmp_br = Branch_pair;

  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
  useModelInPhase1_sep(node, num);
  If_in_Enu_State ? useModelInPhase2(node, CONFIG::BranPhase2InEnu) : useModelInPhase2(node, CONFIG::BranPhase2);
  CGTesting(node, false, false, true);

  auto final_tmp_br = Branch_pair;

  Branch_pair = tmp_br;
  CGTesting(node, true, true, false);
  writeUseModelAcc(final_tmp_br, 10);

  Branch_pair = final_tmp_br;
}

void CVRP::writeUseModelAcc(const std::vector<std::pair<int, int>> &tmp_br, int top_n) {
  //tackle the acc
  if (!is_sorted(Branch_Pair_Val.begin(), Branch_Pair_Val.end(),
                 [](const std::pair<std::pair<int, int>, double> &a,
                    const std::pair<std::pair<int, int>, double> &b) {
                   return a.second > b.second;
                 })) {
    std::sort(Branch_Pair_Val.begin(), Branch_Pair_Val.end(),
              [](const std::pair<std::pair<int, int>, double> &a,
                 const std::pair<std::pair<int, int>, double> &b) {
                return a.second > b.second;
              });
  }
  unordered_map<pair<int, int>, double, PairHasher> tmp_map;
  tmp_map.reserve(Branch_Pair_Val.size());
  for (auto &i : Branch_Pair_Val) {
    tmp_map[i.first] = i.second;
  }
  self_mkdir("ISAcc");
  string filename = "ISAcc/ISAcc_" + FileName + ".txt";
  ofstream out;
  out.open(filename, ios::app);
  out << "----------------------------------------" << endl;
  out << "#All= " << Branch_Pair_Val.size() << " Max= " << Branch_Pair_Val[0].second << endl;
  out << "top_n= " << top_n << ": ";
  for (int i = 0; i < top_n; ++i) {
    out << Branch_Pair_Val[i].second << " ";
  }
  out << endl;
  for (auto &pr : tmp_br) {
    out << tmp_map[pr] << " ";
  }
  out << endl;
  out.close();
}

void CVRP::testInitialScreening(BBNODE *node) {
  constructBranchingSet(node);
  auto tmp_br = Branch_pair;
  Branch_pair.clear();
  InitialScreening(node, true, true, CONFIG::ML_BranchPhase0, CONFIG::Frac4sudoCostBranPhase0);
  auto tmp_br2 = Branch_pair;
  simulateWriteLPPseudocost(node);
  Branch_pair = tmp_br;
  CGTesting(node, true, true, false);
  writeISAcc(tmp_br2);
}

void CVRP::writeISAcc(const std::vector<std::pair<int, int>> &tmp_br) {
  //tackle the acc
  if (!is_sorted(Branch_Pair_Val.begin(), Branch_Pair_Val.end(),
                 [](const std::pair<std::pair<int, int>, double> &a,
                    const std::pair<std::pair<int, int>, double> &b) {
                   return a.second > b.second;
                 })) {
    std::sort(Branch_Pair_Val.begin(), Branch_Pair_Val.end(),
              [](const std::pair<std::pair<int, int>, double> &a,
                 const std::pair<std::pair<int, int>, double> &b) {
                return a.second > b.second;
              });
  }
  unordered_map<pair<int, int>, double, PairHasher> tmp_map;
  tmp_map.reserve(Branch_Pair_Val.size());
  for (auto &i : Branch_Pair_Val) {
    tmp_map[i.first] = i.second;
  }
  self_mkdir("ISAcc");
  string filename = "ISAcc/ISAcc_" + FileName + ".txt";
  ofstream out;
  out.open(filename, ios::app);
  out << "----------------------------------------" << endl;
  out << "#All= " << Branch_Pair_Val.size() << " Max= " << Branch_Pair_Val[0].second << endl;
  out << "top_n= " << 10 << ": ";
  for (int i = 0; i < 10; ++i) {
    out << Branch_Pair_Val[i].second << " ";
  }
  for (auto &pr : tmp_br) {
    out << tmp_map[pr] << " ";
  }
  out << endl;
  out.close();
}

#endif