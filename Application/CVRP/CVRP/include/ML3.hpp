//
// Created by Zhengzhong You on 3/27/23.
//

#ifndef VRPTW_INCLUDE_ML3_HPP_
#define VRPTW_INCLUDE_ML3_HPP_

#include "MACRO.hpp"
#include "BBNODE.hpp"
#include <string>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <bitset>

#ifdef MASTER_VALVE_ML
#include <xgboost/c_api.h>

class CVRP;

class BBNODE;

struct tmp_edgeRelatedData {
  std::vector<std::pair<std::string, double>> BasicFeatures;
  std::vector<std::pair<std::string, double>>
      ResolvingLPFeatures;// 1 is the left branch, 3 is the right branch, and will be used!
  double SB_scores{};
  tmp_edgeRelatedData() = default;
};

struct tmp_deleteCols {
  std::vector<int> Idx4LeftBrCol;
  std::vector<int> Idx4RightBrCol;
  tmp_deleteCols() = default;
};

struct long_edgeRelatedData {
  std::pair<double, int> AverEdgeLP;
  long_edgeRelatedData() = default;
};

class ML3 {
 public:
  //support variable
  double MaxEdgeCost{};
  double MaxMidPointEdgeCord_2_depot{};
  std::vector<std::vector<std::pair<double, double>>> MidPointEdgeCord;
  std::vector<std::vector<double>> MidPointEdgeCord_2_depot;
  std::vector<std::vector<std::vector<int>>> NodeDensity_in_std_dis_vec_form;//
  std::vector<std::vector<double>> Edge_2_other_convert_dis;
  int QID{};
  std::string lp_Output_path{}, exact_Output_path{};
  std::unordered_map<std::pair<int, int>, tmp_edgeRelatedData, PairHasher> EdgeTmpInfo;
  std::unordered_map<std::pair<int, int>, long_edgeRelatedData, PairHasher> EdgeLongInfo;

  BoosterHandle booster_1{};
  BoosterHandle booster_2{};
  void freeModel();
  ~ML3() {
    freeModel();
  }
  void loadModel(const std::string &model_path, int phase);
  void predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx);
  bool if_in_enu{};
  [[nodiscard]] int giveTestingNumCandidates(int node_dep) const;
  [[nodiscard]] int giveInitialScreeningNumCandidates() const {
    return CONFIG::ML_BranchPhase0;
//    int out;
//    if (if_in_enu) {
//      out = CONFIG::ML_NumBrCandiEnuInPhase1;
//    } else {
//      out = CONFIG::ML_NumBrCandiInPhase1;
//    }
//    return out;
  }
  void writeTrainingLPFile(bool if_debug_features);
  void writeTrainingLPFile_combinedFashion(bool if_debug_features);
  void writeTrainingExactFile(bool if_debug_features);
  void printFeatures();

  void readEnuState(bool If_in_enu) {
    if_in_enu = If_in_enu;
  }
  void collectEdgeRelatedFeatures(CVRP *cvrp, BBNODE *node, double org_val);
  void collectVariableRelatedFeatures(CVRP *cvrp,
                                      BBNODE *node,
                                      std::pair<int, int> edge,
                                      const int *solver_ind,
                                      int BeforeNumRow,
                                      int numnz,
                                      double org_val);
  void collectResolvingFeatures(CVRP *cvrp,
                                BBNODE *node,
                                std::pair<int, int> edge,
                                int BeforeNumRow,
                                double tmp_val,
                                double org_val,
                                int numnz,
                                bool dir);

  void collectResolvingFeatures_runCG(CVRP *cvrp,
                                      BBNODE *node,
                                      std::pair<int, int> edge,
                                      double min_rc,
                                      double mean_rc,
                                      int add,
                                      int BeforeNumRow,
                                      double tmp_val,
                                      double org_val);

  void calculatePrerequisites(CVRP *cvrp);

  void getInfo(CVRP *cvrp);

  std::vector<double> if_in_solution;
  std::unordered_map<std::pair<int, int>, double, PairHasher> edge_val;

  /**
   * now this is the enumeration learning section, now I know a machine learning model is necessary!
   */
  std::string enum_lp_path{}, enum_exact_path{};

  std::unordered_map<std::pair<int, int>, tmp_deleteCols, PairHasher> EnumEdgeDeleteCols;
  BoosterHandle enum_booster1{}, enum_booster2{};
//  RowVectorXd rc;
  static void debugInputData(const std::pair<std::string, double> &fs);
  void collectExtraEdgeFeatures(CVRP *cvrp, BBNODE *node);

#ifdef testAcc
  std::string lp_Output_path2{}, enum_lp_path2{};
  std::string exact_Output_path2{}, enum_exact_path2{};
#endif

};

#endif
#endif //VRPTW_INCLUDE_ML3_HPP_
