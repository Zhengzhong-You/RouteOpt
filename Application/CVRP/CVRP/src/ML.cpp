//
// Created by Zhengzhong You on 10/6/22.
//

#include "ML.hpp"
#include "CVRP.hpp"
#include <fstream>

#define safe_xgboost(call) {  int Xerr = (call);\
if (Xerr != 0) { \
   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call + ":" + XGBGetLastError()); \
}\
}

using namespace std;

void ML_base::writeBasicFeature(CVRP *cvrp, BBNODE *node, const std::pair<int, int> &edge) {
  int ai = edge.first, aj = edge.second;
  auto ColPool4Mem = cvrp->ColPool4Mem;
  f0 = data_T(double(EdgeQ[ai][aj]) / MaxQ4Vertex);
  f1 = data_T(double(CostMatrix[0][ai] + CostMatrix[0][aj]) / 2 / MaxCost4Edge);
  f2 = data_T(sqrt_self(float(CostMatrix[0][ai] * CostMatrix[0][aj])) / MaxCost4Edge);
  f3 = data_T(double(CostMatrix[ai][aj] - MinCost4Edge) / (MaxCost4Edge - MinCost4Edge));
  f4 = data_T(node->TreeLevel);
  f5 = data_T(abs(node->Val - RootLPVal) / RootLPVal);
  if (f5 > 1000 || f5 < 0) f5 = 0;
  f6 = data_T(cvrp->ArcGraph[ai][aj]);
  f7 = data_T((f6 - 0.5) / f0);
  f8 = data_T(f6 * f3);
  f9 = data_T(0);
  for (auto &i : node->Idx4LPSolsInColPool) {
    for (auto j = i.first + 1;; ++j) {
      int curr_node = ColPool4Mem[j];
      if (!curr_node) break;
      if (curr_node == ai) {
        if (ColPool4Mem[j + 1] == aj || ColPool4Mem[j - 1] == aj) {
          ++f9;
          break;
        }
      }
    }
  }
  f9 /= data_T(double(node->Idx4LPSolsInColPool.size()));
  if (!AllNumEdgeBranch) {
    f10 = data_T(0);
  } else f10 = data_T(double(NumEdgeBranch[ai][aj]) / AllNumEdgeBranch);

  f11 = data_T(AverEdgeRC[ai][aj].first / AverEdgeRC[ai][aj].second / node->Val);
  f12 = data_T(AverEdgeXij[ai][aj].first / AverEdgeXij[ai][aj].second);
  size_t numnzP;
  double aver_ai, aver_aj;
  if (ai) {
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, nullptr, nullptr, nullptr, ai - 1, 1))
    aver_ai = double(numnzP) / cvrp->NumCol;
  } else {
    aver_ai = 2;
  }
  safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, nullptr, nullptr, nullptr, aj - 1, 1))
  aver_aj = data_T(double(numnzP) / cvrp->NumCol);
  f13 = data_T((aver_ai + aver_aj) / 2);
  f14 = data_T((aver_ai - aver_aj) / 2);
  f15 = data_T(sqrt_self(float(aver_ai * aver_aj)));
  f16 = data_T(aver_ai / aver_aj);
  int code = edge.first * Dim + edge.second;
  double aver_br_up = 0;
  int record;
  auto p = node->Ptr;
  while (p) {
    if (!p->Edge2Cols[code].empty()) {
      record = p->Edge2Cols[code][0];
      for (int i : p->Edge2Cols[code]) {
        if (record != i) {
          ++aver_br_up;
          record = i;
        }
      }
      ++aver_br_up;
    }
    p = p->PNode;
  }
  double aver_br_down = cvrp->NumCol - aver_br_up;
  aver_br_up /= cvrp->NumCol;
  aver_br_down /= cvrp->NumCol;
  f17 = data_T((aver_br_up + aver_br_down) / 2);
  f18 = data_T((aver_br_up - aver_br_down) / 2);
  auto test = float(aver_br_up * aver_br_down);
  if (test < 0) test = 0;
  f19 = data_T(sqrt_self(test));
  f20 = data_T(aver_br_up / aver_br_down);
}

void PhaseI::calculateFeatureNLabel(CVRP *const cvrp, BBNODE *const node, const std::pair<int, int> &edge,
                                    double score) {
//  cout << "After br choice has been made, please remember to call writeDynamicData before using it!" << endl;
  writeBasicFeature(cvrp, node, edge);
  label = data_T(score);
}

void ML_base::writeDynamicData(CVRP *const cvrp, const std::pair<int, int> &edge, int updateMode) {
  //used when br choice has been made
  int ai = edge.first, aj = edge.second;
  if (updateMode == 1) {
    ++NumEdgeBranch[ai][aj];
    ++AllNumEdgeBranch;
  } else if (updateMode == 2) {
    AverEdgeRC[ai][aj].first += cvrp->ChgCostMat4Vertex[ai][aj];
    ++AverEdgeRC[ai][aj].second;
    AverEdgeXij[ai][aj].first += cvrp->ArcGraph[ai][aj];
    ++AverEdgeXij[ai][aj].second;
  }
}

void PhaseI::printRankData() const {
  ofstream trainingData;
  trainingData.open(DataFileName, ios::app);
  trainingData << label;
  trainingData << " qid:" << QID;
  trainingData << " 0:" << f0;
  trainingData << " 1:" << f1;
  trainingData << " 2:" << f2;
  trainingData << " 3:" << f3;
  trainingData << " 4:" << f4;
  trainingData << " 5:" << f5;
  trainingData << " 6:" << f6;
  trainingData << " 7:" << f7;
  trainingData << " 8:" << f8;
  trainingData << " 9:" << f9;
  trainingData << " 10:" << f10;
  trainingData << " 11:" << f11;
  trainingData << " 12:" << f12;
  trainingData << " 13:" << f13;
  trainingData << " 14:" << f14;
  trainingData << " 15:" << f15;
  trainingData << " 16:" << f16;
  trainingData << " 17:" << f17;
  trainingData << " 18:" << f18;
  trainingData << " 19:" << f19;
  trainingData << " 20:" << f20;
  trainingData << endl;
  trainingData.close();
}

void ML_base::giveFileName(const std::string &file_name) {
  DataFileName = file_name;
}

void ML_base::initializeSupportSection(CVRP *const cvrp) {
  //use once at the very first
  MaxQ4Vertex = 0;
  MaxCost4Edge = 0;
  MinCost4Edge = MaxInt;
  Dim = cvrp->Dim;
  RootLPVal = 0;
  AllNumEdgeBranch = 0;
  QID = 0;
  EdgeQ = cvrp->MainResourceAcrossArcsInForwardSense;
  CostMatrix = cvrp->CostMat4Vertex;
  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      if (MinCost4Edge > CostMatrix[i][j]) MinCost4Edge = CostMatrix[i][j];
      if (MaxCost4Edge < CostMatrix[i][j]) MaxCost4Edge = CostMatrix[i][j];
    }
  }
  NumEdgeBranch.resize(Dim);
  for (auto &it : NumEdgeBranch) it.resize(Dim, 0);
  AverEdgeRC.resize(Dim);
  for (auto &it : AverEdgeRC) it.resize(Dim, make_pair(0, 0));
  AverEdgeXij.resize(Dim);
  for (auto &it : AverEdgeXij) it.resize(Dim, make_pair(0, 0));
  for (int i = 1; i < Dim; ++i) {
    if (MaxQ4Vertex < cvrp->InfoVertex[i][3]) {
      MaxQ4Vertex = cvrp->InfoVertex[i][3];
    }
  }
}

void ML_base::addQID() {
  ++QID;
}

void ML_base::readRootLPVal(double rootLPVal) {
  RootLPVal = rootLPVal;
}

bool ML_base::tellIfStopGenerateBr(int UB) const {
  if (QID >= UB) return true;
  return false;
}

bool ML_base::tellIfEmptyData(int LB) const {
  if (QID <= LB) {
    ofstream trainingData;
    trainingData.open(DataFileName, ios::out);//clear all content
    trainingData.close();
    return true;
  }
  return false;
}

void PhaseII::readModel(const char *model_path, int NumFea) {
  NumFeatures = NumFea;
  cout << "Load Model succeed! & NumFeatures= " << NumFeatures << endl;
  safe_xgboost(XGBoosterCreate(nullptr, 0, &booster))
  safe_xgboost(XGBoosterLoadModel(booster, model_path))
}

void PhaseII::predictScore(float *data, int numData, std::vector<std::pair<std::pair<int, int>, double>> &rank) {
  DMatrixHandle test;

  safe_xgboost(XGDMatrixCreateFromMat(data, numData, NumFeatures, 0, &test))

  safe_xgboost(XGBoosterPredict(booster, test, 0, 0, 0, &output_length, &output_result))

  rank.resize(output_length);

  for (unsigned int i = 0; i < output_length; i++) {
    rank[i].second = output_result[i];
  }
//  cout << "这里是测试！" << endl;
//  vector<double> a(output_result, output_result + output_length);
//  sort(a.begin(), a.end());
//  for (unsigned int i = 0; i < output_length; i++) {
//    cout << a[i] << " ";
//  }
  safe_xgboost(XGDMatrixFree(test))
}

void PhaseII::freeModel() const {
  XGBoosterFree(booster);
}

void PhaseII::writeIntoFloatData(float *data, int &beg) {
  data[beg++] = f0;
  data[beg++] = f1;
  data[beg++] = f2;
  data[beg++] = f3;
  data[beg++] = f4;
  data[beg++] = f5;
  data[beg++] = f6;
  data[beg++] = f7;
  data[beg++] = f8;
  data[beg++] = f9;
  data[beg++] = f10;
  data[beg++] = f11;
  data[beg++] = f12;
  data[beg++] = f13;
  data[beg++] = f14;
  data[beg++] = f15;
  data[beg++] = f16;
  data[beg++] = f17;
  data[beg++] = f18;
  data[beg++] = f19;
  data[beg++] = f20;
//  cout << "f0=" << f0 << endl;
//  cout << "f1=" << f1 << endl;
//  cout << "f2=" << f2 << endl;
//  cout << "f3=" << f3 << endl;
//  cout << "f4=" << f4 << endl;
//  cout << "f5=" << f5 << endl;
//  cout << "f6=" << f6 << endl;
//  cout << "f7=" << f7 << endl;
//  cout << "f8=" << f8 << endl;
//  cout << "f9=" << f9 << endl;
//  cout << "f10=" << f10 << endl;
//  cout << "f11=" << f11 << endl;
//  cout << "f12=" << f12 << endl;
//  cout << "f13=" << f13 << endl;
//  cout << "f14=" << f14 << endl;
//  cout << "f15=" << f15 << endl;
//  cout << "f16=" << f16 << endl;
//  cout << "f17=" << f17 << endl;
//  cout << "f18=" << f18 << endl;
//  cout << "f19=" << f19 << endl;
//  cout << "f20=" << f20 << endl;
}

#ifdef INCORPORATE_ML
void PhaseII::calculateFeatureNLabel(CVRP *const cvrp,
                                     BBNODE *const node,
                                     const std::pair<int, int> &edge,
                                     int NumOldOptCols,
                                     double OldValue,
                                     double score
#ifdef ONLY_USE_MODEL
    , int type,
                                     unordered_map<int, vector<int>> &map_node_lp
#endif
) {
  cvrp->Decision_Val.clear();
  double dif1, dif2;
  writeBasicFeature(cvrp, node, edge);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int ai = edge.first, aj = edge.second;
  auto ColPool4Mem = cvrp->ColPool4Mem;
  auto ColPool4Pricing = cvrp->ColPool4Pricing;
  size_t numnzP;
  //overall information about cols used
  auto cbeg = cvrp->solver_beg;
  auto cind = cvrp->solver_ind;
  auto cval = cvrp->solver_val;
  int BeforeNumRow = cvrp->NumRow;
  auto cols_i_j = new int[cvrp->NumCol]();//no pass 0, pass i only 1, pass j only 2, pass i and j 3
#ifdef GENERATE_ML_EXACT_DATA
  int NumPass_i = 0, NumPass_j = 0, NumPass_ij = (int) node->Ptr->Edge2Cols[ai * Dim + aj].size();
  for (int i : node->Ptr->Edge2Cols[ai * Dim + aj]) {
    cols_i_j[i] = 3;
  }
  if (ai) {
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, cbeg, cind, cval, ai - 1, 1))
    for (size_t i = 0; i < numnzP; ++i)
      if (!cols_i_j[cind[i]]) {
        cols_i_j[cind[i]] = 1;
        ++NumPass_i;
      }
  } else {
    for (int i = 0; i < cvrp->NumCol; ++i) {
      if (!cols_i_j[i]) {
        cols_i_j[i] = 1;
        ++NumPass_i;
      }
    }
  }
  safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, cbeg, cind, cval, aj - 1, 1))
  for (size_t i = 0; i < numnzP; ++i)
    if (!cols_i_j[cind[i]]) {
      cols_i_j[cind[i]] = 2;
      ++NumPass_j;
    }
#endif

#ifdef ONLY_USE_MODEL
  int NumPass_i = 0, NumPass_j = 0, NumPass_ij;
  if (type == 1) {
    NumPass_ij = (int) node->Ptr->Edge2Cols[ai * Dim + aj].size();
    for (int i : node->Ptr->Edge2Cols[ai * Dim + aj]) {
      cols_i_j[i] = 3;
    }
  } else if (type == 2) {
    NumPass_ij = (int) map_node_lp[ai * Dim + aj].size();
    NumPass_ij += (int) map_node_lp[aj * Dim + ai].size();
    for (int i : map_node_lp[ai * Dim + aj]) {
      cols_i_j[i] = 3;
    }
    for (int i : map_node_lp[aj * Dim + ai]) {
      cols_i_j[i] = 3;
    }
  }
  if (ai) {
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, cbeg, cind, cval, ai - 1, 1))
    for (size_t i = 0; i < numnzP; ++i)
      if (!cols_i_j[cind[i]]) {
        cols_i_j[cind[i]] = 1;
        ++NumPass_i;
      }
  } else {
    for (int i = 0; i < cvrp->NumCol; ++i) {
      if (!cols_i_j[i]) {
        cols_i_j[i] = 1;
        ++NumPass_i;
      }
    }
  }
  safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, cbeg, cind, cval, aj - 1, 1))
  for (size_t i = 0; i < numnzP; ++i)
    if (!cols_i_j[cind[i]]) {
      cols_i_j[cind[i]] = 2;
      ++NumPass_j;
    }

#endif
  //now we solve a LP to get the updated information
  double product = 1;
  int numnz;
#ifdef ONLY_USE_MODEL
  if (type == 1) {
    cvrp->getNewCstrCoeffByEdge(node, edge, cind, cval, numnz);
  } else if (type == 2) {
    vector<int> cols(cvrp->NumCol, 0);
    for (auto it : map_node_lp[ai * Dim + aj]) ++cols[it];
    for (auto it : map_node_lp[aj * Dim + ai]) ++cols[it];
    int ccnt = 0;
    numnz = 0;
    for (auto it : cols) {
      if (it) {
        cind[numnz] = ccnt;
        cval[numnz++] = it;
      }
      ++ccnt;
    }
  }
#endif

#ifdef GENERATE_ML_EXACT_DATA
  cvrp->getNewCstrCoeffByEdge(node, edge, cind, cval, numnz);
#endif

  safe_solver(cvrp->addBranchconstr(numnz, cind, cval, SOLVER_LESS_EQUAL, 0, nullptr, node->solver))
  ++cvrp->NumRow;
  safe_Hyperparameter(cvrp->checkCST_LIMIT())

  double dual;
  double obj;
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetDual(BeforeNumRow, 1, &dual))
  safe_solver(node->solver.SOLVERgetObjVal(&obj))
  safe_solver(node->solver.SOLVERgetX(0, cvrp->NumCol, cvrp->X))
  dif1 = obj - OldValue;
  auto X = cvrp->X;
  bool if_update_ub = true;
  for (int i = 0; i < cvrp->NumCol; ++i) {
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      if_update_ub = false;
      break;
    }
  }
  if (if_update_ub) {
    if (cvrp->ceil_transformed_number_related(obj - TOLERANCE) + TOLERANCE < cvrp->UB) {
      cvrp->UB = cvrp->ceil_transformed_number_related(obj - TOLERANCE);

      auto IPOptSol = cvrp->IPOptSol;
      IPOptSol[0] = 0;
      int *mem;
#ifdef GENERATE_ML_EXACT_DATA
      mem = ColPool4Mem;
#endif
#ifdef ONLY_USE_MODEL
      if (type == 1) {
        mem = ColPool4Mem;
      } else if (type == 2) {
        mem = ColPool4Pricing;
      }
#endif
      for (int i = 0; i < cvrp->NumCol; ++i) {
        if (X[i] > TOLERANCE) {
          IPOptSol[++IPOptSol[0]] = 0;
          for (size_t j = node->IdxCols[i] + 1;; ++j) {
            if (!mem[j])
              break;
            IPOptSol[++IPOptSol[0]] = mem[j];
          }
          IPOptSol[++IPOptSol[0]] = 0;
        }
      }
    }
  }

  f21 = data_T(dual / obj);

  int NumNewOptCols = 0;
  int NumOptCol_pass_i_j = 0;
  for (int i = 0; i < cvrp->NumCol; ++i) {
    if (cvrp->X[i] > TOLERANCE) {
      ++NumNewOptCols;
      if (cols_i_j[i] == 3) ++NumOptCol_pass_i_j;
    }
  }
  f22 = data_T(double(NumNewOptCols - NumOldOptCols) / NumOldOptCols);
  f23 = data_T((obj - OldValue) / OldValue);
  product *= (obj - OldValue);
  f24 = data_T(double(NumPass_ij) / cvrp->NumCol);
  f25 = data_T(double(NumPass_i) / cvrp->NumCol);
  f26 = data_T(double(NumPass_j) / cvrp->NumCol);
  f27 = data_T(double(NumOptCol_pass_i_j) / NumNewOptCols);
  safe_solver(node->solver.SOLVERgetRC(0, cvrp->NumCol, cvrp->RC))
  double sum_rc_i_j = 0, sum_rc_i = 0, sum_rc_j = 0;
  for (int i = 0; i < cvrp->NumCol; ++i) {
    if (cols_i_j[i] == 3) {
      sum_rc_i_j += cvrp->RC[i];
    } else if (cols_i_j[i] == 1) {
      sum_rc_i += cvrp->RC[i];
    } else if (cols_i_j[i] == 2) {
      sum_rc_j += cvrp->RC[i];
    }
  }
  f28 = data_T(sum_rc_i_j / obj);
  f29 = data_T(sum_rc_i / obj);
  f30 = data_T(sum_rc_j / obj);

  safe_solver(cvrp->inverseLastBranchconstr(SOLVER_GREATER_EQUAL, 1, node->solver))
  safe_solver(node->solver.SOLVERreoptimize())

  safe_solver(node->solver.SOLVERgetDual(BeforeNumRow, 1, &dual))
  safe_solver(node->solver.SOLVERgetObjVal(&obj))
  safe_solver(node->solver.SOLVERgetX(0, cvrp->NumCol, cvrp->X))

  dif2 = obj - OldValue;
  cvrp->Decision_Val[edge] = {dif1 + OldValue, dif2 + OldValue};
  cout << "dif1: " << dif1 << " dif2: " << dif2 << endl;
  X = cvrp->X;
  if_update_ub = true;
  for (int i = 0; i < cvrp->NumCol; ++i) {
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      if_update_ub = false;
      break;
    }
  }
  if (if_update_ub) {
    if (cvrp->ceil_transformed_number_related(obj - TOLERANCE) + TOLERANCE < cvrp->UB) {
      cvrp->UB = cvrp->ceil_transformed_number_related(obj - TOLERANCE);

      auto IPOptSol = cvrp->IPOptSol;
      IPOptSol[0] = 0;
      int *mem;
#ifdef GENERATE_ML_EXACT_DATA
      mem = ColPool4Mem;
#endif
#ifdef ONLY_USE_MODEL
      if (type == 1) {
        mem = ColPool4Mem;
      } else if (type == 2) {
        mem = ColPool4Pricing;
      }
#endif
      for (int i = 0; i < cvrp->NumCol; ++i) {
        if (X[i] > TOLERANCE) {
          IPOptSol[++IPOptSol[0]] = 0;
          for (size_t j = node->IdxCols[i] + 1;; ++j) {
            if (!mem[j])
              break;
            IPOptSol[++IPOptSol[0]] = mem[j];
          }
          IPOptSol[++IPOptSol[0]] = 0;
        }
      }
    }
  }

  f31 = data_T(dual / obj);

  NumNewOptCols = 0;
  NumOptCol_pass_i_j = 0;
  for (int i = 0; i < cvrp->NumCol; ++i) {
    if (cvrp->X[i] > TOLERANCE) {
      ++NumNewOptCols;
      if (cols_i_j[i] == 1 || cols_i_j[i] == 2) ++NumOptCol_pass_i_j;//meaning "or"
    }
  }
  f32 = data_T(double(NumNewOptCols - NumOldOptCols) / NumOldOptCols);
  f33 = data_T((obj - OldValue) / OldValue);

  product *= (obj - OldValue);
  int tmp = ai * Dim + aj;
  auto edge_iter = cvrp->LPTestingBranch.find(tmp);
  if (edge_iter == cvrp->LPTestingBranch.end()) cvrp->LPTestingBranch[tmp] = product;
  else cvrp->LPTestingBranch[tmp] = (cvrp->LPTestingBranch[tmp] + product) / 2;

  f34 = data_T(double(NumOptCol_pass_i_j) / NumNewOptCols);
  safe_solver(node->solver.SOLVERgetRC(0, cvrp->NumCol, cvrp->RC))
  sum_rc_i_j = 0;
  sum_rc_i = 0;
  sum_rc_j = 0;
  for (int i = 0; i < cvrp->NumCol; ++i) {
    if (cols_i_j[i] == 3) {
      sum_rc_i_j += cvrp->RC[i];
    } else if (cols_i_j[i] == 1) {
      sum_rc_i += cvrp->RC[i];
    } else if (cols_i_j[i] == 2) {
      sum_rc_j += cvrp->RC[i];
    }
  }
  f35 = data_T(sum_rc_i_j / obj);
  f36 = data_T(sum_rc_i / obj);
  f37 = data_T(sum_rc_j / obj);

  label = data_T(score);

  safe_solver(node->solver.SOLVERdelconstrs(1, &BeforeNumRow))
  safe_solver(node->solver.SOLVERreoptimize())
  --cvrp->NumRow;

  delete[]cols_i_j;
}

void PhaseII::printRankData() const {
  ofstream trainingData;
  trainingData.open(DataFileName, ios::app);
  trainingData << label;
  trainingData << " qid:" << QID;
  trainingData << " 0:" << f0;
  trainingData << " 1:" << f1;
  trainingData << " 2:" << f2;
  trainingData << " 3:" << f3;
  trainingData << " 4:" << f4;
  trainingData << " 5:" << f5;
  trainingData << " 6:" << f6;
  trainingData << " 7:" << f7;
  trainingData << " 8:" << f8;
  trainingData << " 9:" << f9;
  trainingData << " 10:" << f10;
  trainingData << " 11:" << f11;
  trainingData << " 12:" << f12;
  trainingData << " 13:" << f13;
  trainingData << " 14:" << f14;
  trainingData << " 15:" << f15;
  trainingData << " 16:" << f16;
  trainingData << " 17:" << f17;
  trainingData << " 18:" << f18;
  trainingData << " 19:" << f19;
  trainingData << " 20:" << f20;
  trainingData << " 21:" << f21;
  trainingData << " 22:" << f22;
  trainingData << " 23:" << f23;
  trainingData << " 24:" << f24;
  trainingData << " 25:" << f25;
  trainingData << " 26:" << f26;
  trainingData << " 27:" << f27;
  trainingData << " 28:" << f28;
  trainingData << " 29:" << f29;
  trainingData << " 30:" << f30;
  trainingData << " 31:" << f31;
  trainingData << " 32:" << f32;
  trainingData << " 33:" << f33;
  trainingData << " 34:" << f34;
  trainingData << " 35:" << f35;
  trainingData << " 36:" << f36;
  trainingData << " 37:" << f37;
  trainingData << endl;
  trainingData.close();
}

void PhaseIII::predictScore_phase3(float *data,
                                   int numData,
                                   std::vector<std::pair<std::pair<int, int>, double>> &rank) {
  DMatrixHandle test;

  safe_xgboost(XGDMatrixCreateFromMat(data, numData, NumFeatures_phase3, 0, &test))

  safe_xgboost(XGBoosterPredict(booster_phase3, test, 0, 0, 0, &output_length_phase3, &output_result_phase3))

  rank.resize(output_length_phase3);

  for (unsigned int i = 0; i < output_length_phase3; i++) {
    rank[i].second = output_result_phase3[i];
  }
  safe_xgboost(XGDMatrixFree(test))
}

void PhaseIII::writeIntoFloatData_phase3(float *data, int &beg) {
  data[beg++] = f0;
  data[beg++] = f1;
  data[beg++] = f2;
  data[beg++] = f3;
  data[beg++] = f4;
  data[beg++] = f5;
  data[beg++] = f6;
  data[beg++] = f7;
  data[beg++] = f8;
  data[beg++] = f9;
  data[beg++] = f10;
  data[beg++] = f11;
  data[beg++] = f12;
  data[beg++] = f13;
  data[beg++] = f14;
  data[beg++] = f15;
  data[beg++] = f16;
  data[beg++] = f17;
  data[beg++] = f18;
  data[beg++] = f19;
  data[beg++] = f20;
  data[beg++] = f21;
  data[beg++] = f22;
  data[beg++] = f23;
  data[beg++] = f24;
  data[beg++] = f25;
  data[beg++] = f26;
  data[beg++] = f27;
  data[beg++] = f28;
  data[beg++] = f29;
  data[beg++] = f30;
  data[beg++] = f31;
  data[beg++] = f32;
  data[beg++] = f33;
  data[beg++] = f34;
  data[beg++] = f35;
  data[beg++] = f36;
  data[beg++] = f37;
}

void PhaseIII::readModel(const char *model_path, int NumFea, const char *model_path_phase3, int NumFea_phase3) {
  NumFeatures = NumFea;
  NumFeatures_phase3 = NumFea_phase3;
  safe_xgboost(XGBoosterCreate(nullptr, 0, &booster))
  safe_xgboost(XGBoosterLoadModel(booster, model_path))
  safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_phase3))
  safe_xgboost(XGBoosterLoadModel(booster_phase3, model_path_phase3))
  cout << "Load Model succeed! & NumFeatures= " << NumFeatures << endl;
  cout << "Load Model succeed! & NumFeatures_phase3= " << NumFeatures_phase3 << endl;
}

void PhaseIII::freeModel() const {
  XGBoosterFree(booster);
  XGBoosterFree(booster_phase3);
}
#endif

