//
// Created by Zhengzhong You on 12/11/22.
//

#include "CVRP.hpp"
#include "VRPTW.hpp"

using namespace std;

void CVRP::setResourceInBucketGraph() {
  //in CVRP the main resource is the capacity
  MaxMainResource = Cap;
  MeetPointResourceInBiDir = MaxMainResource / 2;
  MainResourceAcrossArcsInForwardSense.resize(Dim);
  for (auto &vertex : MainResourceAcrossArcsInForwardSense) vertex.resize(Dim);
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
#ifdef SYMMETRY_PROHIBIT
      MainResourceAcrossArcsInForwardSense[i][j] = InfoVertex[i][3];
#else
      MainResourceAcrossArcsInForwardSense[i][j] = (InfoVertex[i][3] + InfoVertex[j][3]) / 2;
#endif
    }
  }
#ifdef SYMMETRY_PROHIBIT
  MainResourceAcrossArcsInBackwardSense.resize(Dim);
  for (auto &vertex : MainResourceAcrossArcsInBackwardSense) vertex.resize(Dim);
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
      MainResourceAcrossArcsInBackwardSense[i][j] = InfoVertex[j][3];
    }
  }
#endif
  lb4Vertex.resize(Dim, 0);
  ub4Vertex.resize(Dim, MaxMainResource);
}

void CVRP::initialBucketGraph() {
  //map the q accumulation into arcs
  LabelArrayInForwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    LabelArrayInForwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      LabelArrayInForwardSense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
    }
  }

  IfExistExtraLabelsInForwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    IfExistExtraLabelsInForwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      IfExistExtraLabelsInForwardSense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
    }
  }
  RC2TillThisBinInForwardSense = new double *[Dim];
  for (int i = 0; i < Dim; ++i) {
    RC2TillThisBinInForwardSense[i] = new double[NumBucketsPerVertex];
  }
  RC2BinInForwardSense = new double *[Dim];
  for (int i = 0; i < Dim; ++i) {
    RC2BinInForwardSense[i] = new double[NumBucketsPerVertex];
  }
#ifdef SYMMETRY_PROHIBIT
  LabelArrayInBackwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    LabelArrayInBackwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      LabelArrayInBackwardSense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
    }
  }

  IfExistExtraLabelsInBackwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    IfExistExtraLabelsInBackwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      IfExistExtraLabelsInBackwardSense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
    }
  }
  RC2TillThisBinInBackwardSense = new double *[Dim];
  for (int i = 0; i < Dim; ++i) {
    RC2TillThisBinInBackwardSense[i] = new double[NumBucketsPerVertex];
  }
  RC2BinInBackwardSense = new double *[Dim];
  for (int i = 0; i < Dim; ++i) {
    RC2BinInBackwardSense[i] = new double[NumBucketsPerVertex];
  }
#endif
}
