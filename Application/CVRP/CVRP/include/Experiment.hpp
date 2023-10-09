//
// Created by Zhengzhong You on 2/13/23.
//

#ifndef VRPTW_INCLUDE_EXPERIMENT_HPP_
#define VRPTW_INCLUDE_EXPERIMENT_HPP_

//#define specialBrInEnu

//#define OpenTest

#ifdef OpenTest
//#define TestIfPhase0IsWorthy
#define TestIfMLIsWorthy
class Experiment {
 public:
#ifdef TestIfPhase0IsWorthy
  static std::unordered_set<std::pair<int, int>, PairHasher> Phase0Candidates;
  static std::unordered_map<int, std::tuple<double, double, int>> Phase0TopNAcc; // top n & acc (sum_acc & sum_take_up & all)
  static void printPhase0TopNAcc();
#endif
  ///NOTE: ML is different with TestIfPhase0IsWorthy, the acc is defined differently.
  ///we use selected mapping rule to label the SB scores, and count if there exist 1 in the selected candidates.
  ///if there exist 1, then we take the acc as one! else it is zero!
#ifdef TestIfMLIsWorthy
  static std::pair<int, int> ML_1TopNAcc; // top n & acc (sum_acc & sum_take_up & all)
  static std::pair<int, int> ML_2TopNAcc; // top n & acc (sum_acc & sum_take_up & all)
  static void printMLTopNAcc();
  static void mappingfunctor(const std::vector<std::pair<std::pair<int, int>, double> > &Branch_Val,
                             std::unordered_map<std::pair<int, int>, int, PairHasher> &rank);
#endif
};
#endif

#endif //VRPTW_INCLUDE_EXPERIMENT_HPP_
