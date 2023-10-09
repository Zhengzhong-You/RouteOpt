//
// Created by You, Zhengzhong on 5/25/23.
//
#include "CVRP.hpp"

using namespace std;
using namespace std::chrono;

#ifdef MASTER_VALVE_ML
#if do_SB_mode == 2
void CVRP::do_SB(BBNODE *node, pair<int, int> &info) {
  if (node->if_terminated) return;
  Branch_pair.clear();
#ifdef MASTER_VALVE_ML
  /**
   * ranking 1>3>4>2(nodes) (367.17, 382.19, 480.22, 533.32), must need sep!
   */
  ml.readEnuState(If_in_Enu_State);
  if (ML_state != 3) throw runtime_error("do_SB: not in ML state");
  bool if_use_initial;
  bool if_use_sep;
  int plan = 4;
  switch (plan) {
    case 1:if_use_initial = true;
      if_use_sep = true;
      break;
    case 2:if_use_initial = true;
      if_use_sep = false;
      break;
    case 3:if_use_initial = false;
      if_use_sep = true;
      break;
    case 4:if_use_initial = false;
      if_use_sep = false;
      break;
    default:throw runtime_error("do_SB: plan not found");
  }
  ExperimentalMode(node, if_use_initial, if_use_sep);
  ml.EdgeTmpInfo.clear();
#endif
  info = Branch_pair[0];
  cout << SMALL_PHASE_SEPARATION;
  cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
  ++BranchChoice[info];
  ++BranchTimes;
}
#endif

void CVRP::ExperimentalMode(BBNODE *node, bool if_use_initial, bool if_use_sep) {
  /**
   * this mode works only under b4 enumeration
   * considering this is the first phase
   */
  if (If_in_Enu_State) throw runtime_error("ExperimentalMode: not in enu state");

  cout << "ExperimentalMode: if_use_initial = " << if_use_initial << ", if_use_sep = " << if_use_sep << endl;

  if (if_use_initial) {
    InitialScreening(node, true, true, ml.giveInitialScreeningNumCandidates(), CONFIG::Frac4sudoCostBranPhase0);
  } else {
    InitialScreening(node, true, true, numeric_limits<int>::max(), CONFIG::Frac4sudoCostBranPhase0);
  }
  int num = min(ml.giveTestingNumCandidates(node->TreeLevel), int(Branch_pair.size()));
  if (if_use_sep) {
    useModelInPhase1_sep(node, num);
  } else {
    useModelInPhase1(node, num);
  }
  LPTesting(node, 0, false);//write record!
  CGTesting(node, true, false, false);//do not write record!
  //use exact cg to tell who is the best!
}

void CVRP::useModelInPhase1_sep(BBNODE *node, int num) {
  if (Branch_pair.size() == 1) {
    cout << "useModelInPhase1_sep: Branch_pair.size() == 1, return!" << endl;
    return;
  }
  getTrainingDataInPhase1(node);
  auto ratio_pseudo =
      ((double) Branch_pair_from_pseudo.size()
          / (double) (Branch_pair_from_pseudo.size() + Branch_pair_from_fractional.size()));
  int pseudo_size, frac_size;
  if (num == 1) {
    pseudo_size = Branch_pair_from_pseudo.empty() ? 0 : 1;
    frac_size = 1 - pseudo_size;
  } else {
    pseudo_size = min((int) (num * ratio_pseudo), (int) Branch_pair_from_pseudo.size());
    frac_size = min(num - pseudo_size, (int) Branch_pair_from_fractional.size());
  }
//evaluate the model from pseudo cost
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_pseudo(Branch_pair_from_pseudo.size());
  transform(Branch_pair_from_pseudo.begin(),
            Branch_pair_from_pseudo.end(),
            Branch_Val_model_pseudo.begin(),
            [](const auto &edge) {
              return make_pair(edge, 0.0);
            });
  ml.predict(Branch_Val_model_pseudo, 1);
  transform(Branch_Val_model_pseudo.begin(),
            Branch_Val_model_pseudo.begin() + pseudo_size,
            Branch_pair.begin(),
            [](const auto &a) {
              return a.first;
            });
  std::vector<std::pair<std::pair<int, int>, double >> Branch_Val_model_fractional(
      Branch_pair_from_fractional.size());
  transform(Branch_pair_from_fractional.begin(),
            Branch_pair_from_fractional.end(),
            Branch_Val_model_fractional.begin(),
            [](const auto &edge) {
              return make_pair(edge, 0.0);
            });
  ml.predict(Branch_Val_model_fractional, 1);
  transform(Branch_Val_model_fractional.begin(),
            Branch_Val_model_fractional.begin() + frac_size,
            Branch_pair.begin() + pseudo_size,
            [](const auto &a) {
              return a.first;
            });
  Branch_pair.resize(num);
  cout << "useModelInPhase1_sep Branch_pair size= " << Branch_pair.size() << endl;
}

#endif