////
//// Created by Zhengzhong You on 2/11/23.
////
//
//#ifndef VRPTW_INCLUDE_ML2_HPP_
//#define VRPTW_INCLUDE_ML2_HPP_
//
//#include "BBNODE.hpp"
//#include <string>
//#include <unordered_set>
//#include <iostream>
//#include <fstream>
//#include <bitset>
//#include <xgboost/c_api.h>
//
////#define MachineLearning
////#define GenerateTrainingData_1
////#define GenerateTrainingData_2
////#define UseModel
////#define DebugFeatures
////#define AccuracyTest //stage 1
//#define PseudoMark std::string("pseudo-")
////#define IfFindPseudoMark
//
////open zero phase
//#define openZeroPhase
//
////control if open cuts and enumeration
//#define openCutsAndEnumerationAtEachEndNode
//
////if this instance needs branching after enumeration, we abandon it!
//#define BranchInEnuNotAllowed
//
////if this instance needs rollback, we abandon it!
//#define CutRollBackNotAllowed
//
////not allow dynamic adjustment for all parameters
//#define DynamicAdjustmentNotAllowed
//
//#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2) || defined(Stage1)
////training details
//
////control training data tree level
//#define TrainingDataTreeLevel 5 //if over this number then skip this node
//#endif
//
//#ifdef UseModel
////#define OpenHeurCGInUseModel
////#define OnlyUseModel1 ///model1 -> heuristic -> res
////#define Model1_LP_Heuristic ///model1 -> LP -> heuristic -> res
//#endif
//
//#ifdef AccuracyTest
//#define AccuracyTest_alpha 0.8
//#define AccuracyTest_outputDir std::string("AccuracyTest")
//#define Stage1
//#ifdef Stage1
//#define NumCandidates 40
//#define NumCandidates_2 8
//#endif
////#define Stage2
////#define OverallStage
//#endif
//
//class CVRP;
//
//class BBNODE;
//
//#ifdef MachineLearning
//struct tmp_edgeRelatedData {
//  double mean_c{};
//  double mean_nonzeroAij{};
//  double mean_x_i_LP{};
//
//  std::vector<double> sum_vval;
//  std::vector<int> pass_ij;
//  double sum_1{};
//  double sum_2{};
//
//  std::vector<std::pair<std::string, double>> BasicFeatures;
//  std::vector<std::pair<std::string, double>> ResolvingLPFeatures;
//
//#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2)
//  double SB_scores{};
//#endif
//
//  tmp_edgeRelatedData() = default;
//};
//
//struct long_edgeRelatedData {
//  double NumEdgePositiveWhenCGConvergent{};
//  std::pair<double, int> AverEdgeLP;
//  std::pair<double, int> AverEdgePseudoCost_up;
//  std::pair<double, int> AverEdgePseudoCost_down;
//  std::pair<double, int> AverEdgePseudoCost_geomean;
//  std::pair<double, int> AverEdgePseudoCostInLPTesting_up;
//  std::pair<double, int> AverEdgePseudoCostInLPTesting_down;
//  std::pair<double, int> AverEdgePseudoCostInLPTesting_geomean;
//  std::pair<double, int> AverEdgeConvertedRC;
//  std::pair<double, int> AverEdgeOptColDensity;
//  long_edgeRelatedData() = default;
//};
//
//class ML2 {
// public:
//  ///TODO: This is the most naive version for selecting features
//  ///TODO: If Feature Selection is completed, this part should be changed
//  ///TODO: For every useful feature, we need to highly optimize it
//  //support variable
//  double RootValue{};
//  double NumCGConvergent{};
//  double MaxEdgeCost{};
//  double MaxMidPointEdgeCord_2_depot{};
//  double std_geo_dis{};
//  std::vector<std::vector<std::pair<double, double>>> MidPointEdgeCord;
//  std::vector<std::vector<double>> MidPointEdgeCord_2_depot;
//  std::vector<std::vector<yzzLong>> NodeDensity_in_std_dis;//
//  std::vector<std::vector<std::vector<int>>> NodeDensity_in_std_dis_vec_form;//
//  std::vector<std::vector<double>> Edge_2_other_convert_dis;
//  std::vector<std::vector<double>> approxi_EdgeRC;
//  std::vector<std::vector<double>> optColRatio;
//#if defined(GenerateTrainingData_1) || defined(GenerateTrainingData_2)
//  int QID{};
//#endif
//#ifdef GenerateTrainingData_1
//  std::string lp_Output_path{};
//  void writeTrainingLPFile();
//#endif
//
//#ifdef GenerateTrainingData_2
//  std::string exact_Output_path{};
//  void writeTrainingExactFile();
//#endif
//  //tmp info
//  std::unordered_map<std::pair<int, int>, tmp_edgeRelatedData, PairHasher> EdgeTmpInfo;//empty
//  //long info
//  std::unordered_map<std::pair<int, int>, long_edgeRelatedData, PairHasher> EdgeLongInfo;
//  std::vector<std::tuple<int, int, double>> all_fractional_edges;
//  void approximateEdgeRC(CVRP *cvrp, BBNODE *node);
//#if defined( UseModel) || defined (GenerateTrainingData_2) || defined (Stage1)
//  BoosterHandle booster_1{};
//  BoosterHandle booster_2{};
//  std::unordered_set<int> SelectedFeatures_model1{3, 6, 8, 52, 53, 55, 64, 66, 67, 71};
//  std::unordered_set<int> SelectedFeatures_model2{55, 3, 71, 89, 53, 52, 6, 8, 54, 9, 35, 38, 39, 44, 43};
//  void loadModel(const std::string &model_path);
//  void predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx);
//  void freeModel() const;
//  static int giveTestingNumCandidates(int node_dep);
//#endif
//
//#ifdef UseModel
//  BoosterHandle booster_3{};
//  std::unordered_set<int>
//      SelectedFeatures_model3{62, 103, 66, 120, 87, 64, 75, 71, 52, 83, 85, 76, 53, 3, 60, 100, 113};
//#endif
//
//#ifdef DebugFeatures
//  void printFeatures();
//#endif
//
//  ~ML2() {
//#if defined( UseModel) || defined (GenerateTrainingData_2) || defined (Stage1)
//    freeModel();
//#endif
//  }
//  void collectEdgeRelatedFeatures(CVRP *cvrp, BBNODE *node, double org_val);
//  void collectVariableRelatedFeatures(CVRP *cvrp,
//                                      BBNODE *node,
//                                      std::pair<int, int> edge,
//                                      int BeforeNumRow,
//                                      int numnz,
//                                      double org_val);
//  void collectResolvingFeatures(CVRP *cvrp,
//                                BBNODE *node,
//                                std::pair<int, int> edge,
//                                int BeforeNumRow,
//                                double tmp_val,
//                                double org_val,
//                                int numnz,
//                                bool dir);
//  void collectInteractFeatures(CVRP *cvrp);
//
//  void supp_findNonele(const std::vector<std::vector<int>> &frac_routes, int Dim);
//
//  void supp_calculateNoneleRatio(std::pair<int, int> &edge, bool if_basic_fea);
//
//  std::vector<size_t> vbeg4getX;
//  std::vector<int> vind4getX;
//  std::vector<double> vval4getX;
//  std::vector<double> Obj_coefficient;
//  std::vector<double> x_rc;
//  std::vector<double> if_in_solution;
//  std::vector<bool> if_activate;
//  std::vector<double> allAverageVval;
//  std::vector<double> cstrRhs;
//  std::vector<size_t> vbeg4getXconstrs;
//  std::vector<int> vind4getXconstrs;
//  std::vector<double> vval4getXconstrs;
//  std::unordered_map<std::pair<int, int>, double, PairHasher> edge_val;
//  std::vector<double> noneleRatio;
//  double allActivateCut{};
//  double old_num_opt_cols{};
//};
//#endif
//
//#endif //VRPTW_INCLUDE_ML2_HPP_
