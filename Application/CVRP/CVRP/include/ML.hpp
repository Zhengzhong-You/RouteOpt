//
// Created by Zhengzhong You on 10/6/22.
//

#ifndef CVRP_TRAINING_HPP
#define CVRP_TRAINING_HPP

#include "BBNODE.hpp"
#include <string>
#include <xgboost/c_api.h>
using data_T = float;

class ML_base {
 public:
  //feature section
  data_T f0{};
  data_T f1{};
  data_T f2{};
  data_T f3{};
  data_T f4{};
  data_T f5{};
  data_T f6{};
  data_T f7{};
  data_T f8{};
  data_T f9{};
  data_T f10{};
  data_T f11{};
  data_T f12{};
  data_T f13{};
  data_T f14{};
  data_T f15{};
  data_T f16{};
  data_T f17{};
  data_T f18{};
  data_T f19{};
  data_T f20{};

  data_T label{};
  //support section
  double MaxQ4Vertex{};
  double MaxCost4Edge{};
  double MinCost4Edge{};
  int Dim{};
  double RootLPVal{};

  std::string DataFileName;
  std::vector<std::vector<double>> EdgeQ;
  std::vector<std::vector<double>> CostMatrix;

  //dynamic support section
  std::vector<std::vector<int>> NumEdgeBranch;
  int AllNumEdgeBranch{};
  std::vector<std::vector<std::pair<double, int>>> AverEdgeRC;
  std::vector<std::vector<std::pair<double, int>>> AverEdgeXij;
  int QID{};

 public:

  void writeDynamicData(CVRP *cvrp, const std::pair<int, int> &edge, int updateMode);//

  void initializeSupportSection(CVRP *cvrp);//

  void giveFileName(const std::string &file_name);//

  void addQID();//

  bool tellIfStopGenerateBr(int UB) const;//

  bool tellIfEmptyData(int LB) const;//

  void readRootLPVal(double rootLPVal);//

  void writeBasicFeature(CVRP *cvrp, BBNODE *node, const std::pair<int, int> &edge);
};

class PhaseI : public ML_base {
 public:
  void calculateFeatureNLabel(CVRP *cvrp, BBNODE *node, const std::pair<int, int> &edge, double score);
  virtual void printRankData() const;
};

class PhaseII : public PhaseI {
 public:
  BoosterHandle booster{};
  bst_ulong output_length{};
  const float *output_result{};//no need to free
  int NumFeatures{};
  data_T f21{};
  data_T f22{};
  data_T f23{};
  data_T f24{};
  data_T f25{};
  data_T f26{};
  data_T f27{};
  data_T f28{};
  data_T f29{};
  data_T f30{};
  data_T f31{};
  data_T f32{};
  data_T f33{};
  data_T f34{};
  data_T f35{};
  data_T f36{};
  data_T f37{};
 public:
  void calculateFeatureNLabel(CVRP *cvrp,
                              BBNODE *node,
                              const std::pair<int, int> &edge,
                              int NumOldOptCols,
                              double OldValue,
                              double score
#ifdef ONLY_USE_MODEL
      , int type,
                              std::unordered_map<int, std::vector<int>> &map_node_lp
#endif
  );
  void printRankData() const override;
  void predictScore(float *data, int numData, std::vector<std::pair<std::pair<int, int>, double>> &rank);
  void writeIntoFloatData(float *data, int &beg);
  void readModel(const char *model_path, int NumFea);
  virtual void freeModel() const;
};

class PhaseIII : public PhaseII {
 public:
  BoosterHandle booster_phase3{};
  bst_ulong output_length_phase3{};
  const float *output_result_phase3{};//no need to free
  int NumFeatures_phase3{};
  void predictScore_phase3(float *data, int numData, std::vector<std::pair<std::pair<int, int>, double>> &rank);
  void writeIntoFloatData_phase3(float *data, int &beg);
  void readModel(const char *model_path, int NumFea, const char *model_path_phase3, int NumFea_phase3);
  void freeModel() const override;
};

#endif //CVRP_TRAINING_HPP
