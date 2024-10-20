//
// Created by You, Zhengzhong on 5/5/24.
//

#ifndef INCLUDE_MACHINE_LEARNING_PREDICT_HPP_
#define INCLUDE_MACHINE_LEARNING_PREDICT_HPP_

/**
 * load model and predict
 */

#define MAX_R_SECOND_STAGE 2//scale the prediction value to 0-2!
#define MAX_RANK_DEVIATION (MAX_R_SECOND_STAGE +2)
#define BEST_PRE (-1)
#define TRUST_PRE_BAR (MAX_R_SECOND_STAGE)//equal is not trust!

#include "machine_learning.hpp"
struct LPNPreInfo {
  double lp;
  int pre;
};

class Predict : public MachineLearning {
 public:
  static int num_deep_testing;
  static BoosterHandle booster_1;
  static BoosterHandle booster_2;
  static std::unordered_map<std::pair<int, int>,
							std::pair<LPNPreInfo, LPNPreInfo>,
							PairHasher> edge_lp_pre;//pair<double, double>, first is lp, second is pre
  static std::vector<std::pair<std::pair<LPNPreInfo, LPNPreInfo>, std::pair<int, int>>> his_recording;
  static std::vector<std::pair<std::pair<LPNPreInfo, LPNPreInfo>, std::pair<int, int>>> his_recording_after_testing;
  static std::vector<std::pair<double, int>> aver_decrease_estimator;
  static double getDecreaseEstimator(int pre, bool if_over_estimate);
  static void testCGOneSide(std::pair<std::pair<int, int>, bool> &info);
  static bool compareTwoCandidates(const std::pair<int, int> &edge1, const std::pair<int, int> &edge2);
  static bool checkIfComparisonNecessary(const std::pair<int, int> &edge2);
  static void compareOnePair(std::pair<int, int> &current_bst, const std::pair<int, int> &edge2);
  static void deeperTesting(int num);
  static bool isDominated(const std::pair<int, int> &edge2,
						  const std::vector<std::pair<std::pair<LPNPreInfo, LPNPreInfo>,
													  std::pair<int, int>>> &recordings);
  static void loadModel(int phase);
  static void predict_model_1(std::vector<std::pair<std::pair<int, int>, double >> &branch_val);
  static void predict_model_2();
  static void useModelPhase1(int num);
  static void useModelPhase2(int num);
  static void useML1_N_LPTesting();
  static void useMLInGeneralFramework();
  static void freeModel();
};

#endif //INCLUDE_MACHINE_LEARNING_PREDICT_HPP_
