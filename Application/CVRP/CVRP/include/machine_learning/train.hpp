//
// Created by You, Zhengzhong on 5/5/24.
//

#ifndef INCLUDE_MACHINE_LEARNING_TRAIN_HPP_
#define INCLUDE_MACHINE_LEARNING_TRAIN_HPP_

/**
 * get training data
 */

#include "machine_learning.hpp"

#define MAX_TREE_LEVEL 5
class GetTrainingData : public MachineLearning {
 public:
  static int qid;
  static std::string lp_output_path;
  static std::string exact_output_path;
  static std::unordered_map<std::pair<int, int>, std::pair<double, double>, PairHasher> edge_cg_change;

  static void writeOneEdgeInfo(const std::pair<int, int> &edge, std::ofstream &trainingData, bool if_left);
  static void findCGScore4OneEdge(const std::pair<int, int> &edge, double dif, bool if_left);
  static void initOutputPath();
  static void generateModelPhase1();
  static void generateModelPhase2();
  static void simulateWriteLPPseudoCost();
  static void writeTrainingLPFile();
  static void writeTrainingExactFile();
};

#endif //INCLUDE_MACHINE_LEARNING_TRAIN_HPP_
