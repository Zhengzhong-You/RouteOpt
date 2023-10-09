////
//// Created by Zhengzhong You on 2/11/23.
////
//
//#include "ML2.hpp"
//#include "CVRP.hpp"
//using namespace std;
//
//#define safe_xgboost(call) {  int Xerr = (call);\
//if (Xerr != 0) { \
//   throw  std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) +": error in " + #call + ":" + XGBGetLastError()); \
//}\
//}
//
//#ifdef GenerateTrainingData_1
//void ML2::writeTrainingLPFile() {
//  ofstream trainingData;
//  trainingData.open(lp_Output_path, ios::app);
//  for (auto &tmp_info : EdgeTmpInfo) {
//    trainingData << tmp_info.second.SB_scores;
//    trainingData << " qid:" << QID;
//    int cnt = 0;
//    for (auto &feature : tmp_info.second.BasicFeatures) {
//      trainingData << " " << cnt << ":" << feature.second;
//      ++cnt;
//    }
//    trainingData << endl;
//  }
//  trainingData.close();
//#ifdef DebugFeatures
//  printFeatures();
//#endif
//}
//#endif
//#ifdef GenerateTrainingData_2
//void ML2::writeTrainingExactFile() {
//  ofstream trainingData;
//  trainingData.open(exact_Output_path, ios::app);
//  for (auto &tmp_info : EdgeTmpInfo) {
//    if (tmp_info.second.ResolvingLPFeatures.empty())
//      continue;
//    trainingData << tmp_info.second.SB_scores;
//    trainingData << " qid:" << QID;
//    int cnt = 0;
//    for (auto &feature : tmp_info.second.BasicFeatures) {
//      trainingData << " " << cnt << ":" << feature.second;
//      ++cnt;
//    }
//    for (auto &feature : tmp_info.second.ResolvingLPFeatures) {
//      trainingData << " " << cnt << ":" << feature.second;
//      ++cnt;
//    }
//    trainingData << endl;
//  }
//#ifdef DebugFeatures
//  printFeatures();
//#endif
//}
//#endif
//
//#if defined(UseModel) || defined(GenerateTrainingData_2) || defined(Stage1)
//
//void ML2::loadModel(const std::string &model_path) {
//  //load model
//  auto path1 = model_path + "/model_1";
//  safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_1))
//  safe_xgboost(XGBoosterLoadModel(booster_1, path1.c_str()))
//  cout << "model_phase1 loaded successfully" << endl;
//  auto path2 = model_path + "/model_2";
//  safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_2))
//  safe_xgboost(XGBoosterLoadModel(booster_2, path2.c_str()))
//  cout << "model_phase2 loaded successfully" << endl;
//#ifdef UseModel
//  auto path3 = model_path + "/model_3";
//  safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_3))
//  safe_xgboost(XGBoosterLoadModel(booster_3, path3.c_str()))
//  cout << "model_phase3 loaded successfully" << endl;
//#endif
//}
//void ML2::predict(std::vector<std::pair<std::pair<int, int>, double >> &Branch_Val, int model_idx) {
//  if (Branch_Val.empty())
//    return;
//  DMatrixHandle test;
//
//  int numFeatures;
//  BoosterHandle booster;
//  bst_ulong output_length;
//  const float *output_result;
//
//  auto &edge = EdgeTmpInfo[Branch_Val[0].first];
//  if (model_idx == 1) {
//    numFeatures = (int) (edge.BasicFeatures.size());
//    booster = booster_1;
//  } else if (model_idx == 2) {
//    numFeatures = (int) (edge.BasicFeatures.size());
//    booster = booster_2;
//  }
//#ifdef UseModel
//  else if (model_idx == 3) {
//    numFeatures = (int) (edge.BasicFeatures.size() + edge.ResolvingLPFeatures.size());
//    booster = booster_3;
//  }
//#endif
//  else {
//    throw std::runtime_error("model_idx is not valid");
//  }
//
//  auto data = new float[Branch_Val.size() * numFeatures];
//
//  if (model_idx == 1) {
//    for (int i = 0; i < Branch_Val.size(); i++) {
//      auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
//      int j = 0;
//      for (auto &fs : tmp_edge.BasicFeatures) {
//        if (SelectedFeatures_model1.find(j) == SelectedFeatures_model1.end())
//          data[i * numFeatures + j] = 0.;
//        else
//          data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//    }
//  } else if (model_idx == 2) {
//    for (int i = 0; i < Branch_Val.size(); i++) {
//      auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
//      int j = 0;
//      for (auto &fs : tmp_edge.BasicFeatures) {
//        if (SelectedFeatures_model2.find(j) == SelectedFeatures_model2.end())
//          data[i * numFeatures + j] = 0.;
//        else
//          data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//    }
//  }
//#ifdef UseModel
//  else if (model_idx == 3) {
//    for (int i = 0; i < Branch_Val.size(); i++) {
//      auto &tmp_edge = EdgeTmpInfo[Branch_Val[i].first];
//      int j = 0;
//      for (auto &fs : tmp_edge.BasicFeatures) {
//        if (SelectedFeatures_model3.find(j) == SelectedFeatures_model3.end())
//          data[i * numFeatures + j] = 0.;
//        else
//          data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//      for (auto &fs : tmp_edge.ResolvingLPFeatures) {
//        if (SelectedFeatures_model3.find(j) == SelectedFeatures_model3.end())
//          data[i * numFeatures + j] = 0.;
//        else
//          data[i * numFeatures + j] = (float) fs.second;
//        ++j;
//      }
//    }
//  }
//#endif
//  else {
//    throw std::runtime_error("model_idx is not valid");
//  }
//
//  safe_xgboost(XGDMatrixCreateFromMat(data, (int) Branch_Val.size(), numFeatures, 0, &test))
//
//  safe_xgboost(XGBoosterPredict(booster, test, 0, 0, 0, &output_length, &output_result))
//
//  for (unsigned int i = 0; i < output_length; i++) {
//    Branch_Val[i].second = output_result[i];
//  }
//  safe_xgboost(XGDMatrixFree(test))
//  delete[] data;
//}
//
//void ML2::freeModel() const {
//  safe_xgboost(XGBoosterFree(booster_1))
//  safe_xgboost(XGBoosterFree(booster_2))
//#ifdef UseModel
//  safe_xgboost(XGBoosterFree(booster_3))
//#endif
//}
//
//int ML2::giveTestingNumCandidates(int node_dep) {
//  for (auto &it : CONFIG::TreeLevel_Num_Vec) {
//    if (node_dep <= it.first)
//      return it.second;
//  }
//  return 0;
//}
//#endif
//
//#ifdef MachineLearning
//void ML2::approximateEdgeRC(CVRP *cvrp, BBNODE *node) {
//  auto NumRow = cvrp->NumRow;
//  auto Pi = cvrp->Pi;
//  auto Dim = cvrp->Dim;
//  auto RealDim = cvrp->RealDim;
//  auto &CostMat4Vertex = cvrp->CostMat4Vertex;
//  auto &map_rank1_multiplier = cvrp->map_rank1_multiplier;
//
//  safe_solver(node->solver.SOLVERupdatemodel())
//  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//  for (int i = 1; i < Dim; ++i) {
//    for (int j = i + 1; j < Dim; ++j) {
//      approxi_EdgeRC[i][j] = CostMat4Vertex[i][j] - 0.5 * (Pi[i - 1] + Pi[j - 1]);
//    }
//  }
//  for (int i = 1; i < Dim; ++i) {
//    approxi_EdgeRC[0][i] = CostMat4Vertex[0][i] - 0.5 * (Pi[i - 1] + Pi[RealDim]);
//  }
//  for (int i = 1; i < Dim; ++i) {
//    for (int j = i + 1; j < Dim; ++j) {
//      approxi_EdgeRC[j][i] = approxi_EdgeRC[i][j];
//    }
//  }
//  for (int i = 1; i < Dim; ++i) {
//    approxi_EdgeRC[i][0] = approxi_EdgeRC[0][i];
//  }
//
//  //deal with RCC
//  double rc;
//  for (auto &rcc : node->RCCs) {
//    if (rcc.FormRCC) {
//      auto &info = rcc.InfoRCCCustomer;
//      rc = Pi[rcc.IdxRCC];
//      for (auto i = info.begin(); i != info.end(); ++i) {
//        auto j = i;
//        ++j;
//        for (; j != info.end(); ++j) {
//          approxi_EdgeRC[*i][*j] -= rc;
//          approxi_EdgeRC[*j][*i] -= rc;
//        }
//      }
//    } else {
//      auto &outside_customer_info = rcc.InfoRCCOutsideCustomer;
//      auto &customer_info = rcc.InfoRCCCustomer;
//      rc = Pi[rcc.IdxRCC];
//      for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
//        auto j = i;
//        ++j;
//        for (; j != outside_customer_info.end(); ++j) {
//          approxi_EdgeRC[*i][*j] -= rc;
//          approxi_EdgeRC[*j][*i] -= rc;
//        }
//      }
//      double half_rc = 0.5 * rc;
//      for (auto it : outside_customer_info) {
//        approxi_EdgeRC[0][it] -= half_rc;
//        approxi_EdgeRC[it][0] -= half_rc;
//      }
//      for (auto it : customer_info) {
//        approxi_EdgeRC[0][it] += half_rc;
//        approxi_EdgeRC[it][0] += half_rc;
//      }
//    }
//  }
//  //deal with Branch
//  for (auto &brc : node->BrCs) {
//    //use real data, not transformed data
//    approxi_EdgeRC[brc.Edge.first][brc.Edge.second] -= Pi[brc.IdxBrC];
//    approxi_EdgeRC[brc.Edge.second][brc.Edge.first] -= Pi[brc.IdxBrC];
//  }
//  //deal with r1cs
//  for (auto &r1c : node->R1Cs) {
//    auto trans_dual = Pi[r1c.IdxR1C] * 0.5;
//    for (auto i : r1c.InfoR1C) {
//      for (int j = 0; j < Dim; ++j) {
//        approxi_EdgeRC[i][j] -= trans_dual;
//      }
//    }
//  }
//  //deal with r1cs_multi
//  for (auto &r1c : node->R1Cs_multi) {
//    auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
//    auto &multi = get<0>(plan);
//    int denominator = get<1>(plan);
//    int cnt = 0;
//    for (auto i : r1c.InfoR1C.first) {
//      auto trans_dual = Pi[r1c.IdxR1C] * multi[cnt] / (double) denominator;
//      for (int j = 0; j < Dim; ++j) {
//        approxi_EdgeRC[i][j] -= trans_dual;
//      }
//      ++cnt;
//    }
//  }
//}
//void ML2::collectEdgeRelatedFeatures(CVRP *cvrp, BBNODE *node, double org_val) {
//  auto NumCol = cvrp->NumCol;
//  auto X = cvrp->X;
//  auto NumRow = cvrp->NumRow;
//  auto Pi = cvrp->Pi;
//  auto Dim = cvrp->Dim;
//  auto RealDim = cvrp->RealDim;
//  auto &Branch_pair = cvrp->Branch_pair;
//  auto &CostMat4Vertex = cvrp->CostMat4Vertex;
//  auto MaxMainResource = cvrp->MaxMainResource;
//  auto Slack = cvrp->Slack;
//  auto ColPool4Mem = cvrp->ColPool4Mem;
//  auto &MainResourceAcrossArcsInForwardSense = cvrp->MainResourceAcrossArcsInForwardSense;
//#ifdef SYMMETRY_PROHIBIT
//  auto &MainResourceAcrossArcsInBackwardSense= cvrp->MainResourceAcrossArcsInBackwardSense;
//#endif
//
//  safe_solver(node->solver.SOLVERoptimize());
//  safe_solver(node->solver.SOLVERgetX(0, NumCol, X));
//  EdgeTmpInfo.clear();
//  old_num_opt_cols = 0;
//  for (int i = 0; i < NumCol; ++i) {
//    if (X[i] > TOLERANCE) {
//      old_num_opt_cols++;
//    }
//  }
//  edge_val.clear();
//  for (int i = 1; i <= node->NumEdges; ++i) {
//    edge_val[make_pair(node->EdgeTail[i], node->EdgeHead[i])] = node->EdgeVal[i];
//  }
//
//  vector<vector<int>> frac_routes(node->Idx4LPSolsInColPool.size());
//  int col_cnt = 0;
//  for (auto &col : node->Idx4LPSolsInColPool) {
//    for (auto j = col.first + 1;; ++j) {
//      if (!ColPool4Mem[j])break;
//      frac_routes[col_cnt].emplace_back(ColPool4Mem[j]);
//    }
//    ++col_cnt;
//  }
//
//  supp_findNonele(frac_routes, Dim);
//  //depth of the node
//  int edge_cnt = 0;
//  for (auto &edge : Branch_pair) {
//    auto &e = EdgeTmpInfo[edge];
//    supp_calculateNoneleRatio(edge, true);
//    e.BasicFeatures.emplace_back("TreeLevel", node->TreeLevel);
//    e.BasicFeatures.emplace_back("NodeImproved", abs(node->Val - RootValue) / RootValue);
//    double sum_frac = 0;
//    for (auto &i : node->Idx4LPSolsInColPool) {
//      sum_frac += min(i.second, 1 - i.second);
//    }
//    e.BasicFeatures.emplace_back("Aver_frac", sum_frac / (int) node->Idx4LPSolsInColPool.size());
//    e.BasicFeatures.emplace_back("EdgeCost_ratio", CostMat4Vertex[edge.first][edge.second] / MaxEdgeCost);
//    e.BasicFeatures.emplace_back("EdgeRC_ratio", approxi_EdgeRC[edge.first][edge.second] / org_val);
//    //update the rc
//    EdgeLongInfo[edge].AverEdgeConvertedRC.first += approxi_EdgeRC[edge.first][edge.second];
//    ++EdgeLongInfo[edge].AverEdgeConvertedRC.second;
//#ifdef SYMMETRY_PROHIBIT
//    auto aver_res=(MainResourceAcrossArcsInForwardSense[edge.first][edge.second]+MainResourceAcrossArcsInBackwardSense[edge.first][edge.second])/2;
//    e.BasicFeatures.emplace_back("EdgeRes_ratio", aver_res/ MaxMainResource);
//#else
//    e.BasicFeatures.emplace_back("EdgeRes_ratio",
//                                 MainResourceAcrossArcsInForwardSense[edge.first][edge.second] / MaxMainResource);
//#endif
//    e.BasicFeatures.emplace_back("EdgeDis2Depot_ratio",
//                                 MidPointEdgeCord_2_depot[edge.first][edge.second] / MaxMidPointEdgeCord_2_depot);
//    double num_frac_edge_related = 0;
//    double aver_dis = 0;
//    double aver_include_frac = 0;
//    auto &mid = MidPointEdgeCord;
//    for (auto &e2 : all_fractional_edges) {
//      pair<int, int> new_edge = {get<0>(e2), get<1>(e2)};
//      int ai = edge.first, aj = edge.second, ii = get<0>(e2), jj = get<1>(e2);
//      if (NodeDensity_in_std_dis[ai][aj].test(ii) && NodeDensity_in_std_dis[ai][aj].test(jj)) {
//        ++num_frac_edge_related;
//        auto dis = sqrt_self(float((mid[ii][jj].first - mid[ai][aj].first) * (mid[ii][jj].first - mid[ai][aj].first) +
//            (mid[ii][jj].second - mid[ai][aj].second) * (mid[ii][jj].second - mid[ai][aj].second)));
//        aver_dis += dis / std_geo_dis;
//        aver_include_frac += get<2>(e2);
//      }
//    }
//    pair<double, int> num_lp_count_down;
//    pair<double, int> num_lp_count_up;
//    pair<double, int> num_exact_count_down;
//    pair<double, int> num_exact_count_up;
//    auto &LPTestingImprovement_down = cvrp->LPTestingImprovement_down;
//    auto &LPTestingImprovement_up = cvrp->LPTestingImprovement_up;
//    auto &RealImprovement_down = cvrp->RealImprovement_down;
//    auto &RealImprovement_up = cvrp->RealImprovement_up;
//
//    for (int k = 0; k < NodeDensity_in_std_dis_vec_form[edge.first][edge.second].size(); ++k) {
//      for (int n = k + 1; n < NodeDensity_in_std_dis_vec_form[edge.first][edge.second].size(); ++n) {
//        auto e2 = make_pair(NodeDensity_in_std_dis_vec_form[edge.first][edge.second][k],
//                            NodeDensity_in_std_dis_vec_form[edge.first][edge.second][n]);
//        if (LPTestingImprovement_down.find(e2) != cvrp->LPTestingImprovement_down.end()) {
//          num_lp_count_down.first += LPTestingImprovement_down[e2].first / LPTestingImprovement_down[e2].second;
//          ++num_lp_count_down.second;
//        }
//        if (LPTestingImprovement_up.find(e2) != cvrp->LPTestingImprovement_up.end()) {
//          num_lp_count_up.first += LPTestingImprovement_up[e2].first / LPTestingImprovement_up[e2].second;
//          ++num_lp_count_up.second;
//        }
//        if (RealImprovement_down.find(e2) != cvrp->RealImprovement_down.end()) {
//          num_exact_count_down.first += RealImprovement_down[e2].first / RealImprovement_down[e2].second;
//          ++num_exact_count_down.second;
//        }
//        if (RealImprovement_up.find(e2) != cvrp->RealImprovement_up.end()) {
//          num_exact_count_up.first += RealImprovement_up[e2].first / RealImprovement_up[e2].second;
//          ++num_exact_count_up.second;
//        }
//      }
//    }
//
//    if (num_lp_count_down.second) {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_lp_down_improved",
//                                   num_lp_count_down.first / num_lp_count_down.second);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_lp_down_improved", 0);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 0);
//    }
//    if (num_lp_count_up.second) {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_lp_up_improved", num_lp_count_up.first / num_lp_count_up.second);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_lp_up_improved", 0);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 0);
//    }
//    if (num_exact_count_down.second) {
//      e.BasicFeatures.emplace_back("Aver_exact_down_improved",
//                                   num_exact_count_down.first / num_exact_count_down.second);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back("Aver_exact_down_improved", 0);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 0);
//    }
//    if (num_exact_count_up.second) {
//      e.BasicFeatures.emplace_back("Aver_exact_up_improved", num_exact_count_up.first / num_exact_count_up.second);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back("Aver_exact_up_improved", 0);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 0);
//    }
//    pair<double, int> adjacent_edge_count_LP_down;
//    pair<double, int> adjacent_edge_count_LP_up;
//    pair<double, int> adjacent_edge_count_exact_down;
//    pair<double, int> adjacent_edge_count_exact_up;
//    for (int i = 1; i < Dim; ++i) {
//      if (i == edge.first || i == edge.second) continue;
//      pair<int, int> tmp_edge;
//      if (edge.first < i) tmp_edge = make_pair(edge.first, i);
//      else tmp_edge = make_pair(i, edge.first);
//      if (LPTestingImprovement_down.find(tmp_edge) != LPTestingImprovement_down.end()) {
//        adjacent_edge_count_LP_down.first +=
//            LPTestingImprovement_down[tmp_edge].first / LPTestingImprovement_down[tmp_edge].second;
//        ++adjacent_edge_count_LP_down.second;
//      }
//      if (LPTestingImprovement_up.find(tmp_edge) != LPTestingImprovement_up.end()) {
//        adjacent_edge_count_LP_up.first +=
//            LPTestingImprovement_up[tmp_edge].first / LPTestingImprovement_up[tmp_edge].second;
//        ++adjacent_edge_count_LP_up.second;
//      }
//      if (RealImprovement_down.find(tmp_edge) != RealImprovement_down.end()) {
//        adjacent_edge_count_exact_down.first +=
//            RealImprovement_down[tmp_edge].first / RealImprovement_down[tmp_edge].second;
//        ++adjacent_edge_count_exact_down.second;
//      }
//      if (RealImprovement_up.find(tmp_edge) != RealImprovement_up.end()) {
//        adjacent_edge_count_exact_up.first +=
//            RealImprovement_up[tmp_edge].first / RealImprovement_up[tmp_edge].second;
//        ++adjacent_edge_count_exact_up.second;
//      }
//    }
//    if (adjacent_edge_count_LP_down.second) {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_adjacent_lp_down_improved",
//                                   adjacent_edge_count_LP_down.first / adjacent_edge_count_LP_down.second);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_adjacent_lp_down_improved", 0);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 0);
//    }
//    if (adjacent_edge_count_LP_up.second) {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_adjacent_lp_up_improved",
//                                   adjacent_edge_count_LP_up.first / adjacent_edge_count_LP_up.second);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back(PseudoMark + "Aver_adjacent_lp_up_improved", 0);
//      e.BasicFeatures.emplace_back(PseudoMark + "If_ever_evaluated", 0);
//    }
//    if (adjacent_edge_count_exact_down.second) {
//      e.BasicFeatures.emplace_back("Aver_adjacent_exact_down_improved",
//                                   adjacent_edge_count_exact_down.first / adjacent_edge_count_exact_down.second);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back("Aver_adjacent_exact_down_improved", 0);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 0);
//    }
//    if (adjacent_edge_count_exact_up.second) {
//      e.BasicFeatures.emplace_back("Aver_adjacent_exact_up_improved",
//                                   adjacent_edge_count_exact_up.first / adjacent_edge_count_exact_up.second);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 1);
//    } else {
//      e.BasicFeatures.emplace_back("Aver_adjacent_exact_up_improved", 0);
//      e.BasicFeatures.emplace_back("If_ever_evaluated", 0);
//    }
//    vector<double> dis_already_br_down;
//    vector<double> dis_already_br_up;
//    int ai = edge.first, aj = edge.second;
//    for (auto &brc : node->BrCs) {
//      int ii = brc.Edge.first, jj = brc.Edge.second;
//      double x_dif = MidPointEdgeCord[ai][aj].first - MidPointEdgeCord[ii][jj].first;
//      double y_dif = MidPointEdgeCord[ai][aj].second - MidPointEdgeCord[ii][jj].second;
//      double dis = sqrt_self(float(x_dif * x_dif + y_dif * y_dif));
//      if (brc.BrDir) {
//        dis_already_br_down.emplace_back(dis);
//      } else {
//        dis_already_br_up.emplace_back(dis);
//      }
//    }
//    //min, max, aver
//    double min_dis_already_br_down = 0, max_dis_already_br_down = 0, aver_dis_already_br_down = 0;
//    bool if_already_br_down = false;
//    if (!dis_already_br_down.empty()) {
//      min_dis_already_br_down = *min_element(dis_already_br_down.begin(), dis_already_br_down.end()) / std_geo_dis;
//      max_dis_already_br_down = *max_element(dis_already_br_down.begin(), dis_already_br_down.end()) / std_geo_dis;
//      aver_dis_already_br_down =
//          accumulate(dis_already_br_down.begin(), dis_already_br_down.end(), 0.0) / (double) dis_already_br_down.size()
//              / std_geo_dis;
//      if_already_br_down = true;
//    }
//    e.BasicFeatures.emplace_back("min_dis_already_br_down", min_dis_already_br_down);
//    e.BasicFeatures.emplace_back("max_dis_already_br_down", max_dis_already_br_down);
//    e.BasicFeatures.emplace_back("aver_dis_already_br_down", aver_dis_already_br_down);
//    e.BasicFeatures.emplace_back("valid", if_already_br_down);
//    double min_dis_already_br_up = 0, max_dis_already_br_up = 0, aver_dis_already_br_up = 0;
//    bool if_already_br_up = false;
//    if (!dis_already_br_up.empty()) {
//      min_dis_already_br_up = *min_element(dis_already_br_up.begin(), dis_already_br_up.end()) / std_geo_dis;
//      max_dis_already_br_up = *max_element(dis_already_br_up.begin(), dis_already_br_up.end()) / std_geo_dis;
//      aver_dis_already_br_up =
//          accumulate(dis_already_br_up.begin(), dis_already_br_up.end(), 0.0) / (double) dis_already_br_up.size()
//              / std_geo_dis;
//    }
//    e.BasicFeatures.emplace_back("min_dis_already_br_up", min_dis_already_br_up);
//    e.BasicFeatures.emplace_back("max_dis_already_br_up", max_dis_already_br_up);
//    e.BasicFeatures.emplace_back("aver_dis_already_br_up", aver_dis_already_br_up);
//    e.BasicFeatures.emplace_back("valid", if_already_br_up);
//    bool is_included = false;
//    if (num_frac_edge_related > TOLERANCE) {
//      aver_dis /= num_frac_edge_related;
//      aver_include_frac /= num_frac_edge_related;
//      is_included = true;
//    } else {
//      aver_dis = 0;
//      aver_include_frac = 0;
//    }
//    e.BasicFeatures.emplace_back("num_frac_edge_related_ratio",
//                                 num_frac_edge_related / (double) all_fractional_edges.size());
//    e.BasicFeatures.emplace_back("aver_dis", aver_dis);
//    e.BasicFeatures.emplace_back("aver_include_frac", aver_include_frac);
//    e.BasicFeatures.emplace_back("is_included", is_included);
//    e.BasicFeatures.emplace_back("NodeDensity_in_std_dis_vec_form_ratio",
//                                 (double) pow(NodeDensity_in_std_dis_vec_form[edge.first][edge.second].size(), 2) / Dim
//                                     / Dim);
//    e.BasicFeatures.emplace_back("Edge_2_other_convert_dis", Edge_2_other_convert_dis[edge.first][edge.second]);
//    e.BasicFeatures.emplace_back("optColPercent", optColRatio[edge.first][edge.second]);
//    //update
//    EdgeLongInfo[edge].AverEdgeOptColDensity.first += optColRatio[edge.first][edge.second];
//    ++EdgeLongInfo[edge].AverEdgeOptColDensity.second;
//    e.sum_vval.resize(NumRow);
//    e.pass_ij.resize(NumCol);
//    ++edge_cnt;
//  }
//  //find the index of fractional variable
//  if_in_solution.resize(NumCol);
//fill(if_in_solution.begin(), if_in_solution.end(), 0);
//  unordered_map<size_t, int> idx_in_solution;
//  idx_in_solution.reserve(NumCol);
//  for (int i = 0; i < NumCol; ++i) {
//    idx_in_solution[node->IdxCols[i]] = i;
//  }
//  for (auto &i : node->Idx4LPSolsInColPool) {
//    if_in_solution[idx_in_solution[i.first]] = i.second;
//  }
//  size_t numnz4getX;
//  vbeg4getX.resize(NumCol + 1);
//  vind4getX.resize(NumRow);
//  vval4getX.resize(NumRow);
//  Obj_coefficient.resize(NumCol);
//  x_rc.resize(NumCol);
//  allAverageVval.resize(NumRow, 0);
//fill(allAverageVval.begin(), allAverageVval.end(), 0);
//  cstrRhs.resize(NumRow);
//  if_activate.resize(NumRow, false);
//fill(if_activate.begin(), if_activate.end(), false);
//  safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack))
//  safe_solver(node->solver.SOLVERgetObj(0, NumCol, Obj_coefficient.data()))
//  safe_solver(node->solver.SOLVERgetRHS(0, NumRow, cstrRhs.data()))
//  safe_solver(node->solver.SOLVERgetRC(0, NumCol, x_rc.data()))
//  transform(cstrRhs.begin(), cstrRhs.end(), cstrRhs.begin(), [](double a) {
//    if (abs(a) > TOLERANCE) return a;
//    else
//      return TOLERANCE;
//  });
//  allActivateCut = 0;
//  for (int i = 0; i < NumRow; ++i) {
//    if (abs(Slack[i]) < TOLERANCE) {
//      if_activate[i] = true;
//      ++allActivateCut;
//    }
//  }
//  size_t numnz4getXconstrs;
//  safe_solver(node->solver.SOLVERXgetconstrs(&numnz4getXconstrs,
//                                             nullptr,
//                                             nullptr,
//                                             nullptr,
//                                             0,
//                                             NumRow))
//  vbeg4getXconstrs.resize(NumRow + 1);
//  vind4getXconstrs.resize(numnz4getXconstrs);
//  vval4getXconstrs.resize(numnz4getXconstrs);
//  safe_solver(node->solver.SOLVERXgetconstrs(&numnz4getXconstrs,
//                                             vbeg4getXconstrs.data(),
//                                             vind4getXconstrs.data(),
//                                             vval4getXconstrs.data(),
//                                             0,
//                                             NumRow))
//  vbeg4getXconstrs[NumRow] = numnz4getXconstrs;
//  for (int i = 0; i < NumRow; ++i) {
//    for (auto j = vbeg4getXconstrs[i]; j < vbeg4getXconstrs[i + 1]; ++j) {
//      allAverageVval[i] += vval4getXconstrs[j];
//    }
//    allAverageVval[i] /= NumCol;
//  }
//}
//
//void ML2::collectVariableRelatedFeatures(CVRP *cvrp,
//                                         BBNODE *node,
//                                         pair<int, int> edge,
//                                         int BeforeNumRow,
//                                         int numnz,
//                                         double org_val) {
//  auto NumCol = cvrp->NumCol;
//  auto NumRow = cvrp->NumRow;
//  auto Dim = cvrp->Dim;
//  auto &solver_ind = cvrp->solver_ind;
//  double sum_p = 0;
//  double sum_c = 0;
//  double sum_abs_Aij = 0;
//  double sum_ij_rc = 0;
//  auto &e = EdgeTmpInfo[edge];
//  int ai = edge.first, aj = edge.second;
//  size_t numnz4getX;
//  auto &pass_ij = e.pass_ij;
//  auto &sum_vval = e.sum_vval;
//  auto &sum_1 = e.sum_1;
//  auto &sum_2 = e.sum_2;
//  //if 0, no related to this edge, if 1, only related to i, if 2, only related to j, if 3, related to both i and j
//  memset(pass_ij.data(), 0, sizeof(int) * NumCol);
//  memset(sum_vval.data(), 0, sizeof(double) * BeforeNumRow);
//  vector<double> q;
//  vector<double> frac_down;
//  vector<double> frac_up;
//  q.reserve(numnz);
//  frac_down.reserve(numnz);
//  frac_up.reserve(numnz);
//  for (int i = 0; i < numnz; ++i) {
//    double one_Aij = 0;
//    double activate = 0;
//    int col_idx = solver_ind[i];
//    pass_ij[col_idx] = 3;
//    sum_ij_rc += x_rc[col_idx];
//    safe_solver(node->solver.SOLVERXgetvars(&numnz4getX,
//                                            vbeg4getX.data(),
//                                            vind4getX.data(),
//                                            vval4getX.data(),
//                                            col_idx,
//                                            1))
//    for (int j = 0; j < numnz4getX; ++j) {
//      one_Aij += abs(vval4getX[j]);
//      sum_vval[vind4getX[j]] += vval4getX[j];
//    }
//    sum_p += (double) numnz4getX / BeforeNumRow;
//    sum_c += Obj_coefficient[col_idx];
//    sum_abs_Aij += (double) one_Aij / (double) numnz4getX;
//    if (abs(if_in_solution[col_idx]) > TOLERANCE) {
//      frac_down.emplace_back(if_in_solution[col_idx]);
//      frac_up.emplace_back(1 - if_in_solution[col_idx]);
//      for (int j = 0; j < numnz4getX; ++j) {
//        if (if_activate[vind4getX[j]]) {
//          ++activate;
//        }
//      }
//      q.emplace_back(activate / allActivateCut);
//    }
//  }
//  sum_1 = 0;
//  sum_2 = 0;
//  //first we look j
//  for (auto j = vbeg4getXconstrs[aj - 1]; j < vbeg4getXconstrs[aj]; ++j) {
//    if (pass_ij[vind4getXconstrs[j]] != 3) {
//      pass_ij[vind4getXconstrs[j]] = 2;
//      ++sum_2;
//    }
//  }
//  //then we look i
//  if (ai) {
//    for (auto j = vbeg4getXconstrs[ai - 1]; j < vbeg4getXconstrs[ai]; ++j) {
//      if (pass_ij[vind4getXconstrs[j]] != 3) {
//        pass_ij[vind4getXconstrs[j]] = 1;
//        ++sum_1;
//      }
//    }
//  } else {
//    for (auto i : pass_ij) {
//      if (!i) {
//        i = 1;
//      }
//    }
//    sum_1 = NumCol - sum_2 - numnz;
//  }
//  //average rc
//  e.BasicFeatures.emplace_back("Aver_rc_pass_ij_ratio", sum_ij_rc / numnz / org_val);
//  e.BasicFeatures.emplace_back("Frac_pass_ij", double(numnz) / NumCol);
//  e.BasicFeatures.emplace_back("Frac_pass_i", double(sum_1) / NumCol);
//  e.BasicFeatures.emplace_back("Frac_pass_j", double(sum_2) / NumCol);
//  //mean p
//  e.BasicFeatures.emplace_back("Mean_p", sum_p / numnz);
//  //mean row_density
//  transform(sum_vval.begin(), sum_vval.end(), sum_vval.begin(), [numnz](const auto &a) {
//    return a / (numnz);
//  });
//  transform(sum_vval.begin(),
//            sum_vval.end(),
//            allAverageVval.begin(),
//            sum_vval.begin(),
//            [](const auto &a, const auto &b) {
//              if (abs(a) < TOLERANCE && abs(b) < TOLERANCE) {
//                return 0.0;
//              } else {
//                return (a - b) / (abs(a) + abs(b));
//              }
//            });
//  e.BasicFeatures.emplace_back("Mean_row_density", accumulate(sum_vval.begin(), sum_vval.end(), 0.0) / BeforeNumRow);
//  e.BasicFeatures.emplace_back("Min_row_density", *min_element(sum_vval.begin(), sum_vval.end()));
//  e.BasicFeatures.emplace_back("Max_row_density", *max_element(sum_vval.begin(), sum_vval.end()));
//  double mean_frac_down = accumulate(frac_down.begin(), frac_down.end(), 0.0) / (int) frac_down.size();
//  e.BasicFeatures.emplace_back("Mean_frac_down", mean_frac_down);
//  e.BasicFeatures.emplace_back("Min_frac_down", *min_element(frac_down.begin(), frac_down.end()));
//  e.BasicFeatures.emplace_back("Max_frac_down", *max_element(frac_down.begin(), frac_down.end()));
//  e.BasicFeatures.emplace_back("Mean_frac_up", accumulate(frac_up.begin(), frac_up.end(), 0.0) / (int) frac_up.size());
//  e.BasicFeatures.emplace_back("Min_frac_up", *min_element(frac_up.begin(), frac_up.end()));
//  e.BasicFeatures.emplace_back("Max_frac_up", *max_element(frac_up.begin(), frac_up.end()));
//  e.mean_c = sum_c / numnz;
//  e.mean_nonzeroAij = sum_abs_Aij / numnz;
//  e.mean_x_i_LP = mean_frac_down;
//  auto frac_edge_down = edge_val[edge];
//  e.BasicFeatures.emplace_back("Frac_edge", frac_edge_down);
//  auto &RealImprovement_up = cvrp->RealImprovement_up;
//  auto &RealImprovement_down = cvrp->RealImprovement_down;
//  auto frac_edge_up = 1 - frac_edge_down;
//  auto it_up = RealImprovement_up.find(edge);
//  auto it_down = RealImprovement_down.find(edge);
//  bool if_up = it_up == RealImprovement_up.end() ? false : true;
//  bool if_down = it_down == RealImprovement_down.end() ? false : true;
//  double improvement_up = if_up ? (it_up->second.first / it_up->second.second) : 0;
//  double improvement_down = if_down ? (it_down->second.first / it_down->second.second) : 0;
//  double pseudo_cost_up = improvement_up * frac_edge_up;
//  double pseudo_cost_down = improvement_down * frac_edge_down;
//  double pseudo_cost_mean = sqrt(pseudo_cost_up * pseudo_cost_down);
//  if (if_up) {
//    EdgeLongInfo[edge].AverEdgePseudoCost_up.first += pseudo_cost_up;
//    ++EdgeLongInfo[edge].AverEdgePseudoCost_up.second;
//  }
//  if (if_down) {
//    EdgeLongInfo[edge].AverEdgePseudoCost_down.first += pseudo_cost_down;
//    ++EdgeLongInfo[edge].AverEdgePseudoCost_down.second;
//  }
//  if (if_up && if_down) {
//    EdgeLongInfo[edge].AverEdgePseudoCost_geomean.first += pseudo_cost_mean;
//    ++EdgeLongInfo[edge].AverEdgePseudoCost_geomean.second;
//  }
//  e.BasicFeatures.emplace_back("pseudo_cost_up_ratio", pseudo_cost_up / org_val);
//  e.BasicFeatures.emplace_back("ever_up", if_up);
//  e.BasicFeatures.emplace_back("pseudo_cost_down_ratio", pseudo_cost_down / org_val);
//  e.BasicFeatures.emplace_back("ever_down", if_down);
//  e.BasicFeatures.emplace_back("pseudo_cost_geomean_ratio", pseudo_cost_mean / org_val);
//  e.BasicFeatures.emplace_back("ever_geomean", if_up && if_down);
//  auto &LPTestingImprovement_up = cvrp->LPTestingImprovement_up;
//  auto &LPTestingImprovement_down = cvrp->LPTestingImprovement_down;
//  auto it_lp_up = LPTestingImprovement_up.find(edge);
//  auto if_lp_up = it_lp_up == LPTestingImprovement_up.end() ? false : true;
//  auto it_lp_down = LPTestingImprovement_down.find(edge);
//  auto if_lp_down = it_lp_down == LPTestingImprovement_down.end() ? false : true;
//  double improvement_lp_up = if_lp_up ? (it_lp_up->second.first / it_lp_up->second.second) : 0;
//  double improvement_lp_down = if_lp_down ? (it_lp_down->second.first / it_lp_down->second.second) : 0;
//  double pseudo_cost_lp_up = improvement_lp_up * frac_edge_up;
//  double pseudo_cost_lp_down = improvement_lp_down * frac_edge_down;
//  double pseudo_cost_lp_mean = sqrt(pseudo_cost_lp_up * pseudo_cost_lp_down);
//  if (if_lp_up) {
//    EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_up.first += pseudo_cost_lp_up;
//    ++EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_up.second;
//  }
//  if (if_lp_down) {
//    EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_down.first += pseudo_cost_lp_down;
//    ++EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_down.second;
//  }
//  if (if_lp_up && if_lp_down) {
//    EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_geomean.first += pseudo_cost_lp_mean;
//    ++EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_geomean.second;
//  }
//  e.BasicFeatures.emplace_back(PseudoMark + "pseudo_cost_lp_up_ratio", pseudo_cost_lp_up / org_val);
//  e.BasicFeatures.emplace_back(PseudoMark + "ever_lp_up", if_lp_up);
//  e.BasicFeatures.emplace_back(PseudoMark + "pseudo_cost_lp_down_ratio", pseudo_cost_lp_down / org_val);
//  e.BasicFeatures.emplace_back(PseudoMark + "ever_lp_down", if_lp_down);
//  e.BasicFeatures.emplace_back(PseudoMark + "pseudo_cost_lp_geomean_ratio", pseudo_cost_lp_mean / org_val);
//  e.BasicFeatures.emplace_back(PseudoMark + "ever_lp_geomean", if_lp_up && if_lp_down);
//  e.BasicFeatures.emplace_back("mean_q", accumulate(q.begin(), q.end(), 0.0) / (int) q.size());
//  e.BasicFeatures.emplace_back("min_q", *min_element(q.begin(), q.end()));
//  e.BasicFeatures.emplace_back("max_q", *max_element(q.begin(), q.end()));
//  e.BasicFeatures.emplace_back("branch_times", cvrp->BranchChoice[edge]);
//  e.BasicFeatures.emplace_back("num_edge_positive_when_cg_convergent_ratio",
//                               EdgeLongInfo[edge].NumEdgePositiveWhenCGConvergent / NumCGConvergent);
//  e.BasicFeatures.emplace_back("improvement_up_ratio", improvement_up / org_val);
//  e.BasicFeatures.emplace_back("improvement_down_ratio", improvement_down / org_val);
//  e.BasicFeatures.emplace_back(PseudoMark + "improvement_lp_up_ratio", improvement_lp_up / org_val);
//  e.BasicFeatures.emplace_back(PseudoMark + "improvement_lp_down_ratio", improvement_lp_down / org_val);
//  auto &up = EdgeLongInfo[edge].AverEdgePseudoCost_up;
//  if (!up.second) {
//    e.BasicFeatures.emplace_back("aver_edge_pseudo_cost_up_ratio", 0);
//    e.BasicFeatures.emplace_back("ever_edge_pseudo_cost_up", 0);
//  } else {
//    e.BasicFeatures.emplace_back("aver_edge_pseudo_cost_up_ratio", up.first / up.second / org_val);
//    e.BasicFeatures.emplace_back("ever_edge_pseudo_cost_up", 1);
//  }
//  //average improvement of the pseudo cost on the down branch
//  auto &down = EdgeLongInfo[edge].AverEdgePseudoCost_down;
//  if (!down.second) {
//    e.BasicFeatures.emplace_back("aver_edge_pseudo_cost_down_ratio", 0);
//    e.BasicFeatures.emplace_back("ever_edge_pseudo_cost_down", 0);
//  } else {
//    e.BasicFeatures.emplace_back("aver_edge_pseudo_cost_down_ratio", down.first / down.second / org_val);
//    e.BasicFeatures.emplace_back("ever_edge_pseudo_cost_down", 1);
//  }
//  //average improvement of the pseudo cost on the geometric mean branch
//  auto &geomean = EdgeLongInfo[edge].AverEdgePseudoCost_geomean;
//  if (!geomean.second) {
//    e.BasicFeatures.emplace_back("aver_edge_pseudo_cost_geomean_ratio", 0);
//    e.BasicFeatures.emplace_back("ever_edge_pseudo_cost_geomean", 0);
//  } else {
//    e.BasicFeatures.emplace_back("aver_edge_pseudo_cost_geomean_ratio",
//                                 geomean.first / geomean.second / org_val);
//    e.BasicFeatures.emplace_back("ever_edge_pseudo_cost_geomean", 1);
//  }
//  //average improvement of the pseudo cost on the up branch in LP testing
//  auto &up_lp = EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_up;
//  if (!up_lp.second) {
//    e.BasicFeatures.emplace_back(PseudoMark + "aver_edge_pseudo_cost_lp_up_ratio", 0);
//    e.BasicFeatures.emplace_back(PseudoMark + "ever_edge_pseudo_cost_lp_up", 0);
//  } else {
//    e.BasicFeatures.emplace_back(PseudoMark + "aver_edge_pseudo_cost_lp_up_ratio",
//                                 up_lp.first / up_lp.second / org_val);
//    e.BasicFeatures.emplace_back(PseudoMark + "ever_edge_pseudo_cost_lp_up", 1);
//  }
//  //average improvement of the pseudo cost on the down_lp branch in LP testing
//  auto &down_lp = EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_down;
//  if (!down_lp.second) {
//    e.BasicFeatures.emplace_back(PseudoMark + "aver_edge_pseudo_cost_lp_down_ratio", 0);
//    e.BasicFeatures.emplace_back(PseudoMark + "ever_edge_pseudo_cost_lp_down", 0);
//  } else {
//    e.BasicFeatures.emplace_back(PseudoMark + "aver_edge_pseudo_cost_lp_down_ratio",
//                                 down_lp.first / down_lp.second / org_val);
//    e.BasicFeatures.emplace_back(PseudoMark + "ever_edge_pseudo_cost_lp_down", 1);
//  }
//  //average improvement of the pseudo cost on the geometric mean branch in LP testing
//  auto &geomean_lp = EdgeLongInfo[edge].AverEdgePseudoCostInLPTesting_geomean;
//  if (!geomean_lp.second) {
//    e.BasicFeatures.emplace_back(PseudoMark + "aver_edge_pseudo_cost_lp_geomean_ratio", 0);
//    e.BasicFeatures.emplace_back(PseudoMark + "ever_edge_pseudo_cost_lp_geomean", 0);
//  } else {
//    e.BasicFeatures.emplace_back(PseudoMark + "aver_edge_pseudo_cost_lp_geomean_ratio",
//                                 geomean_lp.first / geomean_lp.second / org_val);
//    e.BasicFeatures.emplace_back(PseudoMark + "ever_edge_pseudo_cost_lp_geomean", 1);
//  }
//  auto &lp = EdgeLongInfo[edge].AverEdgeLP;
//  e.BasicFeatures.emplace_back("aver_edge_lp", lp.first / lp.second);
//}
//void ML2::collectResolvingFeatures(CVRP *cvrp,
//                                   BBNODE *node,
//                                   pair<int, int> edge,
//                                   int BeforeNumRow,
//                                   double tmp_val,
//                                   double org_val,
//                                   int numnz,
//                                   bool dir) {
//  auto &e = EdgeTmpInfo[edge];
//  auto NumCol = cvrp->NumCol;
//  auto X = cvrp->X;
//  auto RC = cvrp->RC;
//  //warning for enumeration!
//  auto col_pool = cvrp->ColPool4Mem;
//  double dual;
//  auto &pass_ij = e.pass_ij;
//  auto &sum_1 = e.sum_1;
//  auto &sum_2 = e.sum_2;
//  auto dif = cvrp->calculateDif(tmp_val, org_val);
//  safe_solver(node->solver.SOLVERgetDual(BeforeNumRow, 1, &dual))
//  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
//  vector<vector<int>> frac_routes;
//  vector<double> old_rc_of_new_x;
//  pair<double, int> sum_frac_x;
//  for (int i = 0; i < NumCol; ++i) {
//    if (X[i] > TOLERANCE) {
//      frac_routes.emplace_back(vector<int>());
//      for (auto j = node->IdxCols[i] + 1;; ++j) {
//        if (!col_pool[j]) break;
//        frac_routes.back().emplace_back(col_pool[j]);
//      }
//      sum_frac_x.first += min(X[i], 1 - X[i]);
//      ++sum_frac_x.second;
//      old_rc_of_new_x.emplace_back(x_rc[i]);
//    }
//  }
//  e.ResolvingLPFeatures.emplace_back("aver_frac_x", sum_frac_x.first / sum_frac_x.second);
//  e.ResolvingLPFeatures.emplace_back("min_frac_x",
//                                     *min_element(old_rc_of_new_x.begin(), old_rc_of_new_x.end()) / org_val);
//  e.ResolvingLPFeatures.emplace_back("max_frac_x",
//                                     *max_element(old_rc_of_new_x.begin(), old_rc_of_new_x.end()) / org_val);
//  e.ResolvingLPFeatures.emplace_back("aver_frac_x",
//                                     accumulate(old_rc_of_new_x.begin(), old_rc_of_new_x.end(), 0.0) /
//                                         (double) old_rc_of_new_x.size() / org_val);
//  supp_findNonele(frac_routes, cvrp->Dim);
//  supp_calculateNoneleRatio(edge, false);
//  e.ResolvingLPFeatures.emplace_back("dual/tmp_val", dual / tmp_val);
//  int num_opt_col = 0;
//  int num_pass_ij_or_i_or_j = 0;
//  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
//  if (dir) {
//    for (int i = 0; i < NumCol; ++i) {
//      if (X[i] > TOLERANCE) {
//        ++num_opt_col;
//        if (pass_ij[i] == 3) ++num_pass_ij_or_i_or_j;
//      }
//    }
//  } else {
//    for (int i = 0; i < NumCol; ++i) {
//      if (X[i] > TOLERANCE) {
//        ++num_opt_col;
//        if (pass_ij[i] == 1 || pass_ij[i] == 2) ++num_pass_ij_or_i_or_j;
//      }
//    }
//  }
//  e.ResolvingLPFeatures.emplace_back("num_opt_col - old_num_opt_cols_ratio",
//                                     double(num_opt_col - old_num_opt_cols) / (num_opt_col + old_num_opt_cols));
//  e.ResolvingLPFeatures.emplace_back("num_pass_ij_or_i_or_j_ratio", double(num_pass_ij_or_i_or_j) / num_opt_col);
//  e.ResolvingLPFeatures.emplace_back("dif / (tmp_val + org_val)", dif / (tmp_val + org_val));
//  safe_solver(node->solver.SOLVERgetRC(0, NumCol, RC))
//  double sum_rc_ij = 0;
//  double sum_rc_i = 0;
//  double sum_rc_j = 0;
//  for (int i = 1; i < NumCol; ++i) {
//    if (pass_ij[i] == 1) sum_rc_i += RC[i];
//    else if (pass_ij[i] == 2) sum_rc_j += RC[i];
//    else if (pass_ij[i] == 3) sum_rc_ij += RC[i];
//  }
//  bool if_numnz_nonzero = numnz > 0;
//  bool if_sum_1_nonzero = sum_1 > 0;
//  bool if_sum_2_nonzero = sum_2 > 0;
//  sum_rc_ij = if_numnz_nonzero ? sum_rc_ij / (numnz * tmp_val) : 0;
//  sum_rc_i = if_sum_1_nonzero ? sum_rc_i / (sum_1 * tmp_val) : 0;
//  sum_rc_j = if_sum_2_nonzero ? sum_rc_j / (sum_2 * tmp_val) : 0;
//  e.ResolvingLPFeatures.emplace_back("sum_rc_ij", sum_rc_ij);
//  e.ResolvingLPFeatures.emplace_back("sum_rc_i", sum_rc_i);
//  e.ResolvingLPFeatures.emplace_back("sum_rc_j", sum_rc_j);
//  e.ResolvingLPFeatures.emplace_back("if_numnz_nonzero", if_numnz_nonzero);
//  e.ResolvingLPFeatures.emplace_back("if_sum_1_nonzero", if_sum_1_nonzero);
//  e.ResolvingLPFeatures.emplace_back("if_sum_2_nonzero", if_sum_2_nonzero);
//}
//
//void ML2::collectInteractFeatures(CVRP *cvrp) {
//  //get interactive features
//  auto &Branch_pair = cvrp->Branch_pair;
//  vector<double> mean_c(Branch_pair.size());
//  for (auto &edge : Branch_pair) {
//    mean_c.emplace_back(EdgeTmpInfo[edge].mean_c);
//  }
//  auto max_mean_c = *max_element(mean_c.begin(), mean_c.end());
//  auto min_mean_c = *min_element(mean_c.begin(), mean_c.end());
//  for (auto &edge : Branch_pair) {
//    auto &e = EdgeTmpInfo[edge];
//    double discount = (e.mean_c - min_mean_c) / (max_mean_c - min_mean_c);
//    e.BasicFeatures.emplace_back("interact", discount);
//    e.BasicFeatures.emplace_back("interact/e.mean_nonzeroAij", discount / e.mean_nonzeroAij);
//    e.BasicFeatures.emplace_back("interact*e.mean_x_i_LP", discount * e.mean_x_i_LP);
//  }
//}
//
//void ML2::supp_findNonele(const vector<vector<int>> &frac_routes, int Dim) {
//  //frac_routes cannot contain 0
//  if (noneleRatio.empty()) {
//    noneleRatio.resize(Dim, 0);
//  } else {
//    memset(noneleRatio.data(), 0, sizeof(double) * Dim);
//  }
//  std::vector<std::unordered_map<int, pair<int, bool>>> nonele_info(Dim);
//  int cnt = 0;
//  for (auto &route : frac_routes) {
//    yzzLong PI = 0;
//    for (auto n : route) {
//      if (PI[n]) {
//        ++nonele_info[n][cnt].first;
//      } else {
//        PI.set(n);
//      }
//      nonele_info[n][cnt].second = true;
//    }
//    ++cnt;
//  }
//  for (int i = 1; i < Dim; ++i) {
//    double sum = accumulate(nonele_info[i].begin(), nonele_info[i].end(), 0.0,
//                            [](double a, const std::pair<int, std::pair<int, bool>> &b) {
//                              return a + b.second.first;
//                            });
//    //only depot could be 0 in the denominator
//    noneleRatio[i] = nonele_info[i].empty() ? 0 : sum / (double) nonele_info[i].size();
//  }
//}
//
//void ML2::supp_calculateNoneleRatio(pair<int, int> &edge, bool if_basic_fea) {
//  auto &e = EdgeTmpInfo[edge];
//  if (if_basic_fea) {
//    e.BasicFeatures.emplace_back("noneleRatio[edge.first] for basic", noneleRatio[edge.first]);
//    e.BasicFeatures.emplace_back("noneleRatio[edge.second] for basic", noneleRatio[edge.second]);
//    e.BasicFeatures.emplace_back("edge.first != 0 for basic", edge.first != 0);
//  } else {
//    e.ResolvingLPFeatures.emplace_back("noneleRatio[edge.first] for resolving", noneleRatio[edge.first]);
//    e.ResolvingLPFeatures.emplace_back("noneleRatio[edge.second] for resolving", noneleRatio[edge.second]);
//    e.ResolvingLPFeatures.emplace_back("edge.first != 0 for resolving", edge.first != 0);
//  }
//}
//#endif
//
//#ifdef DebugFeatures
//void ML2::printFeatures() {
//  for (auto &tmp_info : EdgeTmpInfo) {
//    int cnt = 0;
//    for (auto &feature : tmp_info.second.BasicFeatures) {
//      cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
//      ++cnt;
//    }
//    for (auto &feature : tmp_info.second.ResolvingLPFeatures) {
//      cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
//      ++cnt;
//    }
//    cout << BIG_PHASE_SEPARATION;
//  }
//#ifdef IfFindPseudoMark
//  auto &edge = *EdgeTmpInfo.begin();
//  int cnt = 0;
//  vector<int> f_set;
//  for (auto &feature : edge.second.BasicFeatures) {
//    if (feature.first.find(PseudoMark) != string::npos) {
//      cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
//      f_set.emplace_back(cnt);
//    }
//    ++cnt;
//  }
//  cout << "f_set: " << endl;
//  cout << "[";
//  for (auto &tmp : f_set) {
//    cout << "," << tmp;
//  }
//  cout << "]";
//#endif
//}
//#endif
//
//
