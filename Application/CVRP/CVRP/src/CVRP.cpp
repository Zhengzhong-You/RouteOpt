//
// Created by Zhengzhong You on 1/20/22.


#include "CVRP.hpp"
#include "templateFunctors.hpp"
//#include <Python.h>

using namespace std;
using namespace chrono;

//---------------------------------------------------------------------------------------------------------------------

CVRP::CVRP(const InstanceData &instanceData) {
  Dim = instanceData.Dim;
  RealDim = Dim - 1;
  K = instanceData.K;
  Cap = instanceData.Cap;
  InfoVertex = instanceData.InfoVertex;
  FileName = instanceData.name;
  safe_Hyperparameter(checkMaxNum_Customers())
}

void CVRP::lateProcessing() {
  Count4Tolerance4tryEnumerationWhenArcEliminationFails =
      CONFIG::InitialTolerance4tryEnumerationWhenArcEliminationFails;
  MaxNumRoute4MIP = CONFIG::MaxNumRoute4MIP;
  Rank1MemSizeLimit = int(Dim * CONFIG::MaxCutMemFactor);
  CutGenTimeThresholdInPricing = CONFIG::CutGenTimeThresholdInPricingInitial;
  SizeNGMem = CONFIG::InitialNGSize;
  GapTolerance4ArcEliminationNEnumeration = CONFIG::InitGapTolerance4ArcEliminationNEnumeration;
  NumBucketsPerVertex = CONFIG::InitialNumBuckets;
  CostMat4Vertex.resize(Dim, vector<double>(Dim, 0));
  Vertex2AllInOneLPR1Cs.resize(Dim);
  Vertex2ActiveInOnePricingR1Cs.resize(Dim);
  NGMem4Vertex.resize(Dim, 0);
  Rank1SepHeurMem4Vertex.resize(Dim, 0);
  SizeNGMem4Vertex.resize(Dim, CONFIG::InitialNGSize);
  PriorMemSets.resize(Dim, CONFIG::InitialNGSize);
  ChgCostMat4Vertex.resize(Dim, vector<double>(Dim));
  ColSetUsedInRootMIP.reserve(InitialNumCol4MapInRootMIP);
  roundUpTolerance = -1. / pow_self(10, transformed_number) + MIP_TOLERANCE;

#ifdef readEnumerationTrees
  if (if_only_read_enumerationTree) {
    tree_path = tree_folder + "/" + CONFIG::tree_path;
    colPool_path = col_pool_folder + "/" + CONFIG::colPool_path;
//	tree_path = CONFIG::tree_path;
//	colPool_path = CONFIG::colPool_path;
	auto pos = tree_path.find_last_of('/')+1;
	auto end=tree_path.find_last_of('.');
	string new_path = tree_path.substr(pos, end);
	FileName= new_path;
  }
#endif

  //initial cost matrix
  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      CostMat4Vertex[i][j] = transformCost(
          sqrt_self(float((InfoVertex[i][1] - InfoVertex[j][1]) * (InfoVertex[i][1] - InfoVertex[j][1]) +
              (InfoVertex[i][2] - InfoVertex[j][2]) * (InfoVertex[i][2] - InfoVertex[j][2]))));
    }
  }
  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      CostMat4Vertex[j][i] = CostMat4Vertex[i][j];
    }
  }

  vector<pair<int, double>> cost(Dim);
  vector<vector<int>> tex(Dim, vector<int>(SizeNGMem, 0));
  for (int i = 1; i < Dim; ++i) {
    for (int j = 1; j < Dim; ++j) {
      cost[j].first = j;
      cost[j].second = CostMat4Vertex[i][j];
    }
    cost[0].second = INFINITY;
    std::stable_sort(cost.begin(), cost.end(),
                     [](const pair<int, double> &a, const pair<int, double> &b) {
                       return a.second < b.second;
                     });
    //ng
//    cout << "i= " << i << " ";
    yzzLong &vst = NGMem4Vertex[i];
    for (int k = 0; k < SizeNGMem; ++k) {
      vst.set(cost[k].first);
//      cout << cost[k].first << " ";
    }//    cout << endl;
    //rank1
    yzzLong &vst2 = Rank1SepHeurMem4Vertex[i];
    for (int k = 0; k < CONFIG::MaxHeurSepMem4RowRank1; ++k) {
      vst2.set(cost[k].first);
    }
//    for (int k = 1; k < Dim; ++k) {
//      vst2.set(1);
//    }
  }
  //initial demand
  Demand = new double[Dim];
  for (int i = 1; i < Dim; ++i) {
    Demand[i] = InfoVertex[i][3];
  }

  //rank1 multiplier
  generateOptimalMultiplier4R1C_multi();
  setResourceInBucketGraph();

  specialize_MaxLengthEleRoute();
  cout << "MaxLengthRoute= " << MaxLengthEleRoute << endl;
  //set the main resource
  StepSize = 2 * int(ceil(MaxMainResource / NumBucketsPerVertex / 2));
  NumBucketsPerVertex = (int) floor(MaxMainResource / StepSize) + 1;

  assignMem();

  if (!if_only_read_enumerationTree) initialBucketGraph();

  //need to change this in vrptw
  if (!if_only_read_enumerationTree)initialLabels();

  int candi_size = Dim * Dim / 2;
  RealImprovement_up.reserve(candi_size);
  RealImprovement_down.reserve(candi_size);
  HeuristicImprovement_up.reserve(candi_size);
  HeuristicImprovement_down.reserve(candi_size);
  LPTestingImprovement_up.reserve(candi_size);
  LPTestingImprovement_down.reserve(candi_size);

  getLowerBoundofMinimumNumCars();

#ifdef MASTER_VALVE_ML
  ml.calculatePrerequisites(this);
#endif
}

void CVRP::initialLabels() {
  SeqBeg = 0;
  AllLabel[0].EndVertex = 0;
  AllLabel[0].PLabel = nullptr;
  //column initialization
  for (int i = 1; i < Dim; ++i) {
    AllLabel[i].PI.set(i);
    AllLabel[i].Cost = CostMat4Vertex[0][i];
    AllLabel[i].Seq = AllSeq + SeqBeg;
    *(AllLabel[i].Seq) = 0;
    *(AllLabel[i].Seq + 1) = i;
    AllLabel[i].IdxEndSeq = 1;
    SeqBeg += AllLabel[i].IdxEndSeq + 1;//
    AllLabel[i].EndVertex = i;
    AllLabel[i].PLabel = AllLabel;
    increaseMainResourceConsumption(0, AllLabel[i].Sum_MainResource, 0, i);
  }
#ifdef SYMMETRY_PROHIBIT
  int max_num = 2 * Dim - 1;
  for (int i = Dim; i < max_num; ++i) {
    int point = i - Dim + 1;
    AllLabel[i].PI.set(point);
    AllLabel[i].Cost = CostMat4Vertex[0][point];
    AllLabel[i].Seq = AllSeq + SeqBeg;
    *(AllLabel[i].Seq) = 0;
    *(AllLabel[i].Seq + 1) = point;
    AllLabel[i].IdxEndSeq = 1;
    SeqBeg += AllLabel[i].IdxEndSeq + 1;//
    AllLabel[i].EndVertex = point;
    AllLabel[i].PLabel = AllLabel;
    decreaseMainResourceConsumption(MaxMainResource, AllLabel[i].Sum_MainResource, 0, point);
  }
#endif
}

///---------------------------------------------------------------------------------------------------------------------

void CVRP::buildModel() {
  //initialize a binary tree p0
  auto *node = new BBNODE(MaxNumEdge, Dim, this);
  initialBucketGraph4Node(node);
  auto p0 = new BidirLinkList(nullptr);
  node->Ptr = new BidirLinkList(p0);
  p0->LNode = node->Ptr;

  for (int i = 1; i < Dim; ++i) {
    p0->Edge2Cols[i].emplace_back(i);
    p0->Edge2Cols[i].emplace_back(i);
  }

  int tmp_cnt = 0;

  PoolBeg4Mem = 0;
  node->IdxCols[tmp_cnt++] = PoolBeg4Mem;

  for (int i = 0; i < Dim; ++i) {
    ColPool4Mem[i] = i;
  }
  ColPool4Mem[Dim] = 0;
  PoolBeg4Mem = Dim + 1;
  for (int i = 1; i < Dim; ++i) {
    node->IdxCols[tmp_cnt++] = PoolBeg4Mem;
    ColPool4Mem[PoolBeg4Mem++] = 0;
    ColPool4Mem[PoolBeg4Mem++] = i;
    ColPool4Mem[PoolBeg4Mem++] = 0;
  }

  double obj_sum = 0;
  for (int i = 1; i < Dim; ++i) {
    solver_obj[i] = 2 * CostMat4Vertex[0][i];
    obj_sum += solver_obj[i];
  }

  UB = CONFIG::UB < obj_sum ? CONFIG::UB : obj_sum;
  cout << "UB= " << UB << endl;

  double increase_UB = int(UB * 1.5) + 1;
  Obj4FirstCol = min(increase_UB, obj_sum);
  solver_obj[0] = Obj4FirstCol;

  //build the initial model
  safe_solver(Solver.SOLVERloadenv(nullptr))
  safe_solver(Solver.SOLVERsetenvThreads(NUM_THREAD_LP, true))
  safe_solver(Solver.SOLVERsetenvOutputFlag(0, true))
  safe_solver(Solver.SOLVERsetenvInfUnbdInfo(1, true))
  safe_solver(Solver.SOLVERsetenvMIPGap(MIPGap, true))
  node->solver.SOLVERgetenv(&Solver);
  const char *model_name = "CVRP.lp";

  cout << SMALL_PHASE_SEPARATION;
  cout << "<Instance  " << FileName << "  Capacity  " << Cap << ">" << endl;

  //we add one more constraint that is the number of vehicles is no less than the LB
  auto *rhs = new double[Dim];
  char *sense = new char[Dim];
  for (int i = 0; i < RealDim; ++i) {
    rhs[i] = 1;
    sense[i] = SOLVER_EQUAL;
    solver_beg[i] = 2 * i;
    solver_ind[2 * i] = 0;
    solver_ind[2 * i + 1] = i + 1;
    solver_val[2 * i] = 1;
    solver_val[2 * i + 1] = 1;
  }
  rhs[RealDim] = K;
  sense[RealDim] = SOLVER_GREATER_EQUAL;
  solver_beg[RealDim] = 2 * RealDim;
  auto i_sup_it = solver_beg[RealDim];
  solver_ind[i_sup_it] = 0;
  solver_val[i_sup_it++] = rhs[RealDim];
  for (int i = 1; i < Dim; ++i, ++i_sup_it) {
    solver_ind[i_sup_it] = i;
    solver_val[i_sup_it] = 1;
  }

  safe_solver(node->solver.SOLVERnewmodel(model_name, Dim, solver_obj, nullptr, nullptr, nullptr, nullptr))
  safe_solver(node->solver.SOLVERXaddconstrs(Dim,
                                             2 * RealDim + Dim,
                                             solver_beg,
                                             solver_ind,
                                             solver_val,
                                             sense,
                                             rhs,
                                             nullptr))
  safe_solver(node->solver.SOLVERupdatemodel())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
  safe_Hyperparameter(checkCST_LIMIT())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  safe_solver(node->solver.SOLVERoptimize())

  IdxNode = 0;
  node->TreeLevel = 0;
  node->Val = 0;
  node->Idx = IdxNode;
  node->NumParentCols = Dim;
//  node->VBasis.resize(CONFIG::MaxNumCols);
//  node->CBasis.resize(CST_LIMIT);

  BBT.push(node);
  LB = 0;
  LB_transformed = 0;

#ifdef MASTER_VALVE_ML
#ifdef if_force_read_sol
  readSolFile(true);
#else
  readSolFile(false);
#endif
#endif

  NumStrongArtiVars = 1;
  NumArtiVars = Dim;
#ifndef SYMMETRY_PROHIBIT
  SeqSizeArtiVars = 2 * RealDim;
#else
  SeqSizeArtiVars = 4 * RealDim;
#endif

  yzzLong tmp_pi = 0;
  for (int i = 1; i < Dim; ++i) {
    tmp_pi.set(i);
    ColSetUsedInRootMIP[tmp_pi] = {vector<int>(1, i), 2 * CostMat4Vertex[0][i]};
  }

  delete[] rhs;
  delete[] sense;
}

///---------------------------------------------------------------------------------------------------------------------

void CVRP::solveLPInLabeling(BBNODE *node, bool if_open_heur, bool if_open_exact, bool if_record_sol) {
  node->if_Int = false;

  if (!node->Idx) {
    Ratio_DominanceChecks_NonDominant = {};
  }

  optimizeLP4OneIter(node, node->Val);

  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.SOLVERgetenvMethod(&env_method))
  if (env_method != SOLVER_PRIMAL_SIMPLEX) {
    safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_PRIMAL_SIMPLEX))
    if_changed = true;
  }

  if (if_open_heur) {
    //Lighter Heuristic Phase
    if (runColumnGenerationType(node, 1) || runColumnGenerationType(node, 2))
      goto QUIT;
  }

  if (if_open_exact) {
//    cout << "hrewrwe!" << endl;
//    int num_col;
//    vector<size_t> idx;
//    {
//      num_col = NumCol;
//      idx.assign(node->IdxCols, node->IdxCols + num_col);
//    }

    //find aver dual value as initial value
    if (CONFIG::InitialGapInPricing != 0) {
      vector<double> rank1_dual(node->R1Cs.size() + node->R1Cs_multi.size(), 0);
      int cnt = 0;
      for (auto &r1c : node->R1Cs) {
        if (-Pi4Labeling[r1c.IdxR1C] > TOLERANCE)
          rank1_dual[cnt++] = -Pi4Labeling[r1c.IdxR1C];
      }
      for (auto &r1c : node->R1Cs_multi) {
        if (-Pi4Labeling[r1c.IdxR1C] > TOLERANCE)
          rank1_dual[cnt++] = -Pi4Labeling[r1c.IdxR1C];
      }
      if (cnt) {
        //aver
        GapBetweenLastSmallestRCAndRCThreshold = accumulate(rank1_dual.begin(), rank1_dual.end(), 0.0) / cnt;
      } else GapBetweenLastSmallestRCAndRCThreshold = 0;
    } else {
      GapBetweenLastSmallestRCAndRCThreshold = 0;
    }

//    if_exact = true;

//    int num_col;
//    GRBmodel *model;
//    vector<size_t> idx;
//    {
//      num_col = NumCol;
//      idx.assign(node->IdxCols, node->IdxCols + num_col);
//      GapBetweenLastSmallestRCAndRCThreshold = 0;
//      if_ban_convertDual = true;
//      CONFIG::convert_dual_pricingTime_vs_LPTime = 1000;
//      model = GRBcopymodel(node->solver.model);
//    }
//
//    auto beg = high_resolution_clock::now();
    if (runColumnGenerationType(node, 3)) {
      goto QUIT;
    }
//    auto end = high_resolution_clock::now();
//    cout << "CG used " << duration<double>(end - beg).count() << " s" << endl;

/// jun 27 删除
//    cout << " success= " << succ_ratio.first / succ_ratio.second << endl;
//    cout << " fail= " << fail_ratio.first / fail_ratio.second << endl;
//    succ_ratio = {0, 0};
//    fail_ratio = {0, 0};

//    {
//      cout << "rerun..." << endl;
//      CONFIG::convert_dual_pricingTime_vs_LPTime = 0;
//      GapBetweenLastSmallestRCAndRCThreshold = 0;
//      if_ban_convertDual = false;
//      for (int i = 0; i < num_col; ++i) {
//        node->IdxCols[i] = idx[i];
//      }
//      NumCol = num_col;
//      GRBfreemodel(node->solver.model);
//      node->solver.model = model;
//      safe_solver(node->solver.SOLVERupdatemodel())
//      safe_solver(node->solver.SOLVERoptimize())
//      optimizeLP4OneIter(node, node->Val);
//      beg = high_resolution_clock::now();
//      if (runColumnGenerationType(node, 3)) {
//        goto QUIT;
//      }
//      end = high_resolution_clock::now();
//      cout << "Good CG used " << duration<double>(end - beg).count() << " s" << endl;
//    }

//    if (mode == 0) {
//      GapBetweenLastSmallestRCAndRCThreshold = 0;
//      for (int i = 0; i < num_col; ++i) {
//        node->IdxCols[i] = idx[i];
//      }
//      int len = NumCol - num_col;
//      vector<int> tmp(len);
//      iota(tmp.begin(), tmp.end(), num_col);
//      NumCol = num_col;
//      GRBdelvars(node->solver.model, len, tmp.data());
//      safe_solver(node->solver.SOLVERupdatemodel())
//      safe_solver(node->solver.SOLVERoptimize())
//      mode = 1;
////      Pi4Labeling = tmp_rc;
//      optimizeLP4OneIter(node, node->Val);
//      if (Pi4Labeling != tmp_rc) {
//        cout << "No error!" << endl;
//      }
//      mode = 0;
//      cout << "rerun..." << endl;
//      if (runColumnGenerationType(node, 3)) {
//        cout << "quit!" << endl;
//        exit(0);
//      }
//    }

    if (NumCol > LPCol_FINAL_LIMIT) cleanIdxCol4Node(node, node->NumParentCols);
    if_can_arc_eliminationByExactCG = true;
  }

  if (if_record_sol && Rollback != 1) {
    recordOptCol(node);
  }

  QUIT:
  if (if_changed) {
    safe_solver(node->solver.SOLVERsetenvMethod(env_method))
  }
  if (node->if_terminated) return;
  if (!node->Idx && if_open_exact) {
    double aver_ratio = Ratio_DominanceChecks_NonDominant.first / Ratio_DominanceChecks_NonDominant.second;
    cout << "aver_ratio: " << aver_ratio << endl;
    if (aver_ratio > CONFIG::BucketResizeFactor_Ratio_DominanceChecks_NonDominant) {
      auto numArcs = node->NumForwardBucketArcs + node->NumForwardJumpArcs;
#ifdef SYMMETRY_PROHIBIT
      numArcs += node->NumBackwardBucketArcs + node->NumBackwardJumpArcs;
#endif
      if (double(numArcs) / Dim < CONFIG::BucketResizeFactor_Num_BucketArcsPerVertex
          && !if_force_not_regenerate_bucket_graph) {
        regenerateBucketGraph(node);
      }
    }
  }
}

///---------------------------------------------------------------------------------------------------------------------

///---------------------------------------------------------------------------------------------------------------------

void CVRP::changeModel4BetterDual(BBNODE *node) {
//  safe_solver(node->solver.SOLVERsetenvOutputFlag(1, false))
//  safe_solver(node->solver.SOLVERreoptimize())
//  //change the coefficient of the objective function of the first column
//  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//  LPVal -= 1e-6;
//print LPVal for 6 digits
//  cout << "-------------------" << endl;
//  cout << fixed << setprecision(10);
//  cout << "LPVal: " << LPVal << endl;
//  safe_solver(node->solver.SOLVERreoptimize())
//  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
//  cout << fixed << setprecision(10);
//  cout << "LPVal: " << LPVal << endl;
  safe_solver(node->solver.SOLVERchgObj(0, 1, &LPVal))
//  safe_solver(node->solver.SOLVERupdatemodel())
//  double obj;
//  safe_solver(node->solver.SOLVERgetObj(0, 1, &obj))
//  cout << "obj: " << obj << endl;
//  cout << "-------------------" << endl;
  //change the rhs
  vector<double> old_rhs(NumRow);
  safe_solver(node->solver.SOLVERgetRhs(0, NumRow, old_rhs.data()))
  vector<double> rhs(NumRow, 0);
  for (auto &r1c : node->R1Cs) {
    rhs[r1c.IdxR1C] = int(r1c.InfoR1C.size() + r1c.Mem.size());
//    rhs[r1c.IdxR1C] = 1;
//    rhs[r1c.IdxR1C] = int(r1c.Mem.size());
  }
  for (auto &r1c : node->R1Cs_multi) {
    rhs[r1c.IdxR1C] = int(r1c.InfoR1C.first.size() + r1c.Mem.size());
//    rhs[r1c.IdxR1C] = 1;
//    rhs[r1c.IdxR1C] = int(r1c.Mem.size());
  }
//  vector<int> mem_size(node->R1Cs.size() + node->R1Cs_multi.size(), 0);
//  int cnt = 0;
//  for (auto &r1c : node->R1Cs) {
//    mem_size[cnt++] = int(r1c.InfoR1C.size() + r1c.Mem.size());
////    rhs[r1c.IdxR1C] = 1;
////    rhs[r1c.IdxR1C] = int(r1c.Mem.size());
//  }
//  for (auto &r1c : node->R1Cs_multi) {
//    mem_size[cnt++] = int(r1c.InfoR1C.first.size() + r1c.Mem.size());
////    rhs[r1c.IdxR1C] = 1;
////    rhs[r1c.IdxR1C] = int(r1c.Mem.size());
//  }
//  sort(mem_size.begin(), mem_size.end());
//  int mem_size_threshold;
//  if (mem_size.empty()) mem_size_threshold = 0;
//  else
//    mem_size_threshold = mem_size[int(mem_size.size() * 0.5)];
//  for (auto &r1c : node->R1Cs) {
//    rhs[r1c.IdxR1C] = int(r1c.InfoR1C.size() + r1c.Mem.size()) > mem_size_threshold ? 1 : 0;
////    rhs[r1c.IdxR1C] = 1;
////    rhs[r1c.IdxR1C] = int(r1c.Mem.size());
//  }
//  for (auto &r1c : node->R1Cs_multi) {
//    rhs[r1c.IdxR1C] = int(r1c.InfoR1C.first.size() + r1c.Mem.size()) > mem_size_threshold ? 1 : 0;
////    rhs[r1c.IdxR1C] = 1;
////    rhs[r1c.IdxR1C] = int(r1c.Mem.size());
//  }
//  vector<unordered_set<int>> vertex_2_cuts(Dim);
//  for (auto &r1c : node->R1Cs) {
//    for (auto &n : r1c.InfoR1C) {
//      vertex_2_cuts[n].insert(r1c.IdxR1C);
//    }
//    for (auto &mem : r1c.Mem) {
//      vertex_2_cuts[mem].insert(r1c.IdxR1C);
//    }
//  }
//  for (auto &r1c : node->R1Cs_multi) {
//    for (auto &n : r1c.InfoR1C.first) {
//      vertex_2_cuts[n].insert(r1c.IdxR1C);
//    }
//    for (auto &mem : r1c.Mem) {
//      vertex_2_cuts[mem].insert(r1c.IdxR1C);
//    }
//  }
//  //sort by the vector size
//  sort(vertex_2_cuts.begin(), vertex_2_cuts.end(), [](const unordered_set<int> &a, const unordered_set<int> &b) {
//    return a.size() > b.size();
//  });
//  int start_from = Dim + 1;
//  for (auto &v : vertex_2_cuts) {
//    for (auto &c : v) {
//      if (rhs[c] != 0.) rhs[c] = start_from;
//    }
//    --start_from;
//  }
  safe_solver(node->solver.SOLVERsetRhs(0, NumRow, rhs.data()))
  safe_solver(node->solver.SOLVERsetColUpper(0, 0.))
  safe_solver(node->solver.SOLVERremoveColLower(0))
  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.SOLVERgetenvMethod(&env_method))
  if (env_method != SOLVER_DUAL_SIMPLEX) {
    safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_DUAL_SIMPLEX))
    if_changed = true;
  }
  safe_solver(node->solver.SOLVERreoptimize())
//  double test_val;
//  safe_solver(node->solver.SOLVERgetObjVal(&test_val))
  //change back model
  bool if_new_dual = true;
  int status;
  safe_solver(node->solver.SOLVERgetStatus(&status))
  if (status == SOLVER_UNBOUNDED || status == SOLVER_INF_OR_UNBD) {
//    cout << "find dual Model is unbound after pre_solving! use the old dual!" << endl;
//    cout << "status: " << status << endl;
    if_new_dual = false;
  }
  if (if_new_dual) {
    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi4Labeling.data()))
    safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack4Labeling.data()))
  }
  //safe_solver(node->solver.SOLVERchgObj(0, 1, &UB))
  //safe_solver(node->solver.SOLVERsetRHS(0, NumRow, old_rhs.data()))
  //safe_solver(node->solver.SOLVERsetColLower(0, 0.))
  //safe_solver(node->solver.SOLVERremoveColUpper(0))
  safe_solver(node->solver.SOLVERchgObj(0, 1, &Obj4FirstCol))
  safe_solver(node->solver.SOLVERsetRhs(0, NumRow, old_rhs.data()))
  safe_solver(node->solver.SOLVERsetColLower(0, 0.))
  safe_solver(node->solver.SOLVERremoveColUpper(0))
//  safe_solver(node->solver.SOLVERupdatemodel())
  safe_solver(node->solver.SOLVERreoptimize())
  if (!if_new_dual) {
    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi4Labeling.data()))
    safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack4Labeling.data()))
  }
//  safe_solver(node->solver.SOLVERsetenvOutputFlag(0, false))
  if (if_changed) {
    safe_solver(node->solver.SOLVERsetenvMethod(env_method))
  }
//  double pre_dicted = 0;
//  for (int i = 0; i < NumRow; ++i) {
//    pre_dicted += old_rhs[i] * Pi4Labeling[i];
//  }
//
//  if (abs(LPVal - pre_dicted) > 2e-6) {
//    cout << "LPVal= " << LPVal << " pre_dicted= " << pre_dicted << " test_val= " << test_val << endl;
////    cout << "LPVal= " << LPVal << " pre_dicted= " << pre_dicted << endl;
//    exit(0);
//  }


  /**
   * compare!
   */

//  get rhs and coefficient of the first column and compare them!
//  vector<double> coef(NumRow);
//  for (int i = 0; i < NumRow; ++i) {
//    GRBgetcoeff(node->solver.model, i, 0, &coef[i]);
//    if (coef[i] != old_rhs[i]) cout << "i= " << i << " coef[i]= " << coef[i] << " old_rhs[i]= " << old_rhs[i] << endl;
//  }

//  for (auto &r1c : node->R1Cs) {
//    double dif = Pi[r1c.IdxR1C] - Pi4Labeling[r1c.IdxR1C];
//    if (abs(dif) > 1e-6)
//      cout << dif << " ";
//  }
//  cout << endl;
//  for (auto &r1c : node->R1Cs_multi) {
//    double dif = Pi[r1c.IdxR1C] - Pi4Labeling[r1c.IdxR1C];
//    if (abs(dif) > 1e-6)
//      cout << dif << " ";
//  }
//  cout << endl;
}

void CVRP::priceLabeling(BBNODE *node) {
//use LPVal
/**
 * in this function, we get slack and dual to find the Pi4Labeling and slack used in labeling and cutting
 */
  Pi4Labeling.resize(NumRow);
  Slack4Labeling.resize(NumRow);

  if (
//      false &&
      if_exact_labeling_CG && !if_ban_convertDual && GapBetweenLastSmallestRCAndRCThreshold < 1
          && (!node->R1Cs.empty() || !node->R1Cs_multi.empty())
      ) {
    /**
   * if at first iteration, the pricing time is larger than 10 times resolving lp time,
   * we do the convert
   *, and at next iteration, if the resolving time is even larger then pricing time,
   * we stop the convert
   */
    if (TimePricing4Iter > CONFIG::convert_dual_pricingTime_vs_LPTime * TimeResolveLP4Iter) {
      changeModel4BetterDual(node);
      cout << "change model" << endl;
    } else {
      if (TimeResolveLP4Iter > TimePricing4Iter) {
        if_ban_convertDual = true;
        cout << "ban convert" << endl;
      }
      goto here;
    }
  } else {
    here:
    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi4Labeling.data()))
    safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack4Labeling.data()))
  }

//  if (mode == 1) {
//    Pi4Labeling = tmp_rc;
//    cout << "change dual!" << endl;
//  }


//  cout << "Numrow= " << NumRow << " readdim= " << RealDim << endl;

//  safe_solver(node->solver.SOLVERupdatemodel())
//  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi4Labeling))
  for (int i = 1; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      ChgCostMat4Vertex[i][j] = CostMat4Vertex[i][j] - 0.5 * (Pi4Labeling[i - 1] + Pi4Labeling[j - 1]);
    }
  }
  for (int i = 1; i < Dim; ++i) {
    ChgCostMat4Vertex[0][i] = CostMat4Vertex[0][i] - 0.5 * (Pi4Labeling[i - 1] + Pi4Labeling[RealDim]);
  }
  for (int i = 1; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      ChgCostMat4Vertex[j][i] = ChgCostMat4Vertex[i][j];
    }
  }
  for (int i = 1; i < Dim; ++i) {
    ChgCostMat4Vertex[i][0] = ChgCostMat4Vertex[0][i];
  }

  //deal with RCC
  double rc;
  for (auto &rcc : node->RCCs) {
    if (rcc.FormRCC) {
      auto &info = rcc.InfoRCCCustomer;
      rc = Pi4Labeling[rcc.IdxRCC];
      for (auto i = info.begin(); i != info.end(); ++i) {
        auto j = i;
        ++j;
        for (; j != info.end(); ++j) {
          ChgCostMat4Vertex[*i][*j] -= rc;
          ChgCostMat4Vertex[*j][*i] -= rc;
        }
      }
    } else {
      auto &outside_customer_info = rcc.InfoRCCOutsideCustomer;
      auto &customer_info = rcc.InfoRCCCustomer;
      rc = Pi4Labeling[rcc.IdxRCC];
      for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
        auto j = i;
        ++j;
        for (; j != outside_customer_info.end(); ++j) {
          ChgCostMat4Vertex[*i][*j] -= rc;
          ChgCostMat4Vertex[*j][*i] -= rc;
        }
      }
      double half_rc = 0.5 * rc;
      for (auto it : outside_customer_info) {
        ChgCostMat4Vertex[0][it] -= half_rc;
        ChgCostMat4Vertex[it][0] -= half_rc;
      }
      for (auto it : customer_info) {
        ChgCostMat4Vertex[0][it] += half_rc;
        ChgCostMat4Vertex[it][0] += half_rc;
      }
    }
  }
  //deal with Branch
  for (auto &brc : node->BrCs) {
    if (!brc.BrDir) {
      ChgCostMat4Vertex[brc.Edge.first][brc.Edge.second] = numeric_limits<float>::max();
      ChgCostMat4Vertex[brc.Edge.second][brc.Edge.first] =
          numeric_limits<float>::max();//do not use double since the number will overflow
    } else {
      ChgCostMat4Vertex[brc.Edge.first][brc.Edge.second] -= Pi4Labeling[brc.IdxBrC];
      ChgCostMat4Vertex[brc.Edge.second][brc.Edge.first] -= Pi4Labeling[brc.IdxBrC];
    }
  }
}



///---------------------------------------------------------------------------------------------------------------------


///---------------------------------------------------------------------------------------------------------------------

void CVRP::assignMem() {

  RouteInPricingAssign = MAX_ROUTE_PRICING;
  AverRouteLength = int(ceil(Dim / K)) + 5;
  Mem4Pricing = AverRouteLength * (RouteInPricingAssign);
  ColPool4Pricing = new int[Mem4Pricing];

  if (!if_only_read_enumerationTree) {
    LabelAssign = LABEL_ASSIGN;
    RouteInMemAssign = MAX_ROUTE_MEM;
    AllLabel = new LABEL[LabelAssign];
    AllSeq = new int[(LabelAssign) * MaxLengthEleRoute];
    Mem4Mem = AverRouteLength * (RouteInMemAssign);
    ColPool4Mem = new int[Mem4Mem];
    AllValidR1Cs = new int[(LabelAssign) * MaxNum_R1Cs];
    AllValidR1C_multi = new int[(LabelAssign) * MaxNum_R1C_multi];
    AllR1C_multi_Mem = new int[(LabelAssign) * MaxNum_R1C_multi];

    for (size_t i = 0; i < LabelAssign; ++i) {
      AllLabel[i].validRank1Cut = AllValidR1Cs + i * MaxNum_R1Cs;
      AllLabel[i].validRank1Cut_multi = AllValidR1C_multi + i * MaxNum_R1C_multi;
      AllLabel[i].Rank1CutMem_multi = AllR1C_multi_Mem + i * MaxNum_R1C_multi;
    }
  }

  MaxNonZeroEntry = CONFIG::MaxNumCols * MaxLengthEleRoute;
  solver_ind = new int[MaxNonZeroEntry];
  solver_ind2 = new int[MaxNonZeroEntry];
  solver_val = new double[MaxNonZeroEntry];

  ///---------------------------------------------------------------------------------------------------------------------
  //will not be reallocated!

  solver_beg = new size_t[CONFIG::MaxNumCols];
  solver_obj = new double[CONFIG::MaxNumCols];
  X = new double[CONFIG::MaxNumCols];
  const_for_branching = new int[CONFIG::MaxNumCols];
  std::iota(const_for_branching, const_for_branching + CONFIG::MaxNumCols, 0);
  RC = new double[CONFIG::MaxNumCols];
  Pi = new double[CST_LIMIT];
  Slack = new double[CST_LIMIT];
  //#################################################################################################################
  //BBNode section
  ArcGraph = new double *[Dim];
  for (int i = 0; i < Dim; ++i) {
    ArcGraph[i] = new double[Dim];
  }
  ArcGraph_revised = new double *[Dim];//
  for (int i = 0; i < Dim; ++i) {
    ArcGraph_revised[i] = new double[Dim];
  }
}

///---------------------------------------------------------------------------------------------------------------------

CVRP::~CVRP() {
  delete[]AllLabel;
  delete[]AllSeq;
  delete[]X;
  delete[]RC;
  delete[]Pi;
  delete[]Slack;
  delete[]ColPool4Pricing;
  delete[]copyColPool4Pricing;
  delete[]ColPool4Mem;

  delete[]solver_beg;
  delete[]solver_ind;
  delete[]solver_ind2;
  delete[]solver_obj;
  delete[]solver_val;

  delete[]const_for_branching;
  delete[]Demand;
  delete[]AllValidR1Cs;
  delete[]AllValidR1C_multi;
  delete[]AllR1C_multi_Mem;

  if (!if_only_read_enumerationTree) {
    for (int i = 0; i < Dim; ++i) {
      delete[]RC2TillThisBinInForwardSense[i];
    }
    delete[]RC2TillThisBinInForwardSense;
    for (int i = 0; i < Dim; ++i) {
      delete[]RC2BinInForwardSense[i];
    }
    delete[]RC2BinInForwardSense;
    for (int i = 0; i < Dim; ++i) {
      delete[]LabelArrayInForwardSense[i];
    }
    delete[]LabelArrayInForwardSense;

    for (int i = 0; i < Dim; ++i) {
      delete[]IfExistExtraLabelsInForwardSense[i];
    }
    delete[]IfExistExtraLabelsInForwardSense;
  }

#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < Dim; ++i) {
    delete[]RC2TillThisBinInBackwardSense[i];
  }
  delete[]RC2TillThisBinInBackwardSense;
  for (int i = 0; i < Dim; ++i) {
    delete[]RC2BinInBackwardSense[i];
  }
  delete[]RC2BinInBackwardSense;
  for (int i = 0; i < Dim; ++i) {
    delete[]LabelArrayInBackwardSense[i];
  }
  delete[]LabelArrayInBackwardSense;
  for (int i = 0; i < Dim; ++i) {
    delete[]IfExistExtraLabelsInBackwardSense[i];
  }
  delete[]IfExistExtraLabelsInBackwardSense;
#endif

  for (int i = 0; i < Dim; ++i) {
    delete[]ArcGraph[i];
  }
  delete[]ArcGraph;
  for (int i = 0; i < Dim; ++i) {
    delete[]ArcGraph_revised[i];
  }
  delete[]ArcGraph_revised;
}

///---------------------------------------------------------------------------------------------------------------------
bool operator==(const RCC &lhs, const RCC &rhs) {
  if (lhs.RHS != rhs.RHS) return false;
  if (lhs.FormRCC != rhs.FormRCC) {
    return false;
  } else if (lhs.FormRCC) {
    if (lhs.InfoRCCCustomer != rhs.InfoRCCCustomer) {
      return false;
    }
  } else {
    if (lhs.InfoRCCOutsideCustomer != rhs.InfoRCCOutsideCustomer) {
      return false;
    }
  }
  return true;
}

bool CmpLabelRCLess(const LABEL *l1, const LABEL *l2) {
  return (l1->RC < l2->RC);
}

///#####################################################################################
void CVRP::recordOptCol(BBNODE *node, bool if_force_rewrite) {
  safe_solver(node->solver.SOLVERupdatemodel())
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  node->Val = LPVal;
  if (ceil_transformed_number_related(node->Val - TOLERANCE) + TOLERANCE >= UB) {
    node->if_terminated = true;
    return;
  }

  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))

  int parent_column = If_in_Enu_State ? Dim : node->NumParentCols;
  //record the sequence of columns
  vector<pair<size_t, double>> tmp_solindex;

  for (int i = 0; i < parent_column; ++i) {
    if (X[i] > TOLERANCE) {
      tmp_solindex.emplace_back(node->IdxCols[i], X[i]);
    }
  }

  int break_record = (int) tmp_solindex.size();
  node->NumParentColsInLPSols = break_record;

  for (int i = parent_column; i < NumCol; ++i) {
    if (X[i] > TOLERANCE) {
      tmp_solindex.emplace_back(node->IdxCols[i], X[i]);
    }
  }

  if (node->Idx4LPSolsInColPool == tmp_solindex) {
    cout << "sol remains the same" << endl;
    if (if_force_rewrite) {
      goto REWRITE;
    } else {
      return;
    }
  }

  node->Idx4LPSolsInColPool = tmp_solindex;

#ifdef MachineLearning
  auto &optColRatio = ml.optColRatio;
  for (int i = 0; i < Dim; ++i) {
    memset(optColRatio[i].data(), 0, sizeof(double) * Dim);
  }
  unordered_map<pair<int, int>, unordered_set<size_t>, PairHasher> edge2cols;
  for (int i = 0; i < break_record; ++i) {
    for (size_t j = tmp_solindex[i].first;;) {
      int from = ColPool4Mem[j];
      int to = ColPool4Mem[j + 1];
      if (from > to) swap(from, to);
      edge2cols[{from, to}].insert(i);
      if (!ColPool4Mem[++j]) break;
    }
  }
  for (int i = break_record; i < tmp_solindex.size(); ++i) {
    for (size_t j = tmp_solindex[i].first;;) {
      int from = ColPool4Pricing[j];
      int to = ColPool4Pricing[j + 1];
      if (from > to) swap(from, to);
      edge2cols[{from, to}].insert(i);
      if (!ColPool4Pricing[++j]) break;
    }
  }
  for (auto &e : edge2cols) {
    int from = e.first.first;
    int to = e.first.second;
    optColRatio[from][to] += (double) e.second.size() / (double) tmp_solindex.size();
  }
#endif

  REWRITE:
  //update the edge info
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
      ArcGraph[i][j] = 0;
      ArcGraph_revised[i][j] = 0;
    }
  }

  int beg = 0;
  int end = (int) tmp_solindex.size();
  int *colpool = ColPool4Pricing;
  if (!If_in_Enu_State) {
    end = break_record;
    colpool = ColPool4Mem;
  }

  AGAIN:
  for (int i = beg; i < end; ++i) {
    for (size_t j = tmp_solindex[i].first;;) {
      ArcGraph[colpool[j]][colpool[j + 1]] += tmp_solindex[i].second;
      if (!colpool[++j]) break;
    }
  }

  if (end != tmp_solindex.size()) {
    beg = end;
    end = (int) tmp_solindex.size();
    colpool = ColPool4Pricing;
    goto AGAIN;
  }

  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      ArcGraph[i][j] += ArcGraph[j][i];
      //construct revised ArcGraph
      ArcGraph_revised[i][j] = ArcGraph[i][j];
    }
  }


  //revise the arti_cols from 1 to fx_Dimension
  //the Dim cols in front will never be touched, so it's ok
  //but u should realize that this is not the case for enumeration
  //so there is a particular treatment for enumeration
  vector<pair<int, double>> tmp_first_dim;
  for (int i = 0; i < Dim; ++i) if (X[i] > TOLERANCE) tmp_first_dim.emplace_back(i, X[i]);
  if (If_in_Enu_State) {
    for (auto &i : tmp_first_dim) {
      if (i.first) {
        int valid_length = 0, curr_node;
        for (size_t j = node->IdxCols[i.first] + 1;; ++j) {
          curr_node = ColPool4Pricing[j];
          if (!curr_node) break;
          ++valid_length;
        }
        curr_node = ColPool4Pricing[node->IdxCols[i.first] + 1];
        if (valid_length == 1)ArcGraph_revised[0][curr_node] -= i.second;
        else break;
      }
    }
  } else {
    for (auto &i : tmp_first_dim)
      if (i.first) ArcGraph_revised[0][i.first] -= i.second;
  }

  bool if_reformulate = true;
  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      if (ArcGraph_revised[i][j] > TOLERANCE && abs(ArcGraph_revised[i][j] - 1) > TOLERANCE) {
        if_reformulate = false;
        goto outside;
      }
    }
  }

  outside:
  if (if_reformulate) {
    reformulateIPSol(node);
  }

  NumEdge = 0;

  if (if_force_rewrite) {
    for (int i = 0; i < Dim; ++i) {
      for (int j = i + 1; j < Dim; ++j) {
        if (ArcGraph_revised[i][j] > TOLERANCE) {
          ++NumEdge;
        }
      }
    }
  } else {
    for (int i = 0; i < Dim; ++i) {
      for (int j = i + 1; j < Dim; ++j) {
        if (ArcGraph_revised[i][j] > TOLERANCE) {
          ++NumEdge;
#ifdef MachineLearning
          ++ml.EdgeLongInfo[{i, j}].NumEdgePositiveWhenCGConvergent;
#endif
        }
#ifdef MachineLearning
        auto &tmp = ml.EdgeLongInfo[{i, j}].AverEdgeLP;
        tmp.first += ArcGraph_revised[i][j];
        ++tmp.second;
#endif
#ifdef MASTER_VALVE_ML
        auto &tmp = ml.EdgeLongInfo[{i, j}].AverEdgeLP;
        tmp.first += ArcGraph_revised[i][j];
        ++tmp.second;
#endif
      }
    }
#ifdef MachineLearning
    ++ml.NumCGConvergent;
    ml.approximateEdgeRC(this, node);
#endif
  }

  if (NumEdge > MaxNumEdge) {
    MaxNumEdge = 2 * NumEdge;
  }
}

int CVRP::optimizeLP4OneIter(BBNODE *node, double prior_value) {

//  if (NumCol == 5366)
//  cout << "-----------------------------------" << endl;
//  safe_solver(node->solver.SOLVERsetenvOutputFlag(1, false))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  if (NumCol - node->NumParentCols > LIMIT_NEW_ADDED_COLS && abs(prior_value - LPVal) > TOLERANCE)
    cleanIdxCol4Node(node, node->NumParentCols);
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))

  priceLabeling(node);

  node->if_Int = true;
  for (int i = 0; i < NumCol; ++i) {
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      node->if_Int = false;
      break;
    }
  }
  if (node->if_Int) {
    if (ceil_transformed_number_related(LPVal - TOLERANCE) + TOLERANCE < UB) {

      UB = ceil_transformed_number_related(LPVal - TOLERANCE);

      IPOptSol.clear();
      int tmp_p_col = node->NumParentCols;
      for (int i = 0; i < tmp_p_col; ++i) {
        if (X[i] > TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(MaxLengthEleRoute);
          tmp.emplace_back(0);
          for (size_t j = node->IdxCols[i] + 1;; ++j) {
            if (!ColPool4Mem[j])
              break;
            tmp.emplace_back(ColPool4Mem[j]);
          }
          tmp.emplace_back(0);
          IPOptSol.emplace_back(std::move(tmp));
        }
      }

      for (int i = tmp_p_col; i < NumCol; ++i) {
        if (X[i] > TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(MaxLengthEleRoute);
          tmp.emplace_back(0);
          for (size_t j = node->IdxCols[i] + 1;; ++j) {
            if (!ColPool4Pricing[j])
              break;
            tmp.emplace_back(ColPool4Pricing[j]);
          }
          tmp.emplace_back(0);
          IPOptSol.emplace_back(std::move(tmp));
        }
      }

#ifdef MASTER_VALVE_ML
      updateUB_EdgeSol();
#endif

      if (ceil_transformed_number_related(LB_transformed - TOLERANCE) + TOLERANCE >= UB) {
        node->Val = LPVal;
        node->if_terminated = true;
        cout << TERMINATED_MESSAGE_PROMISING_UPDATE_UB;
        return 1;
      }
    }
  }

  return 0;
}

void CVRP::cleanIdxCol4Node(BBNODE *node, int beg, bool if_only_rcfixing) {
  if (beg >= NumCol) return;
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetRC(0, NumCol, RC))
  vector<double> rc_copy(RC + beg, RC + NumCol);
  OptGap = calculateOptGap(node);
  double threshold;
  if (if_only_rcfixing) {
    threshold = OptGap;
  } else {
    int n = int((NumCol - beg) * COL_KEPT_FRAC);
    nth_element(rc_copy.begin(), rc_copy.begin() + n, rc_copy.end());
    threshold = max(min(rc_copy[n], OptGap), TOLERANCE);
  }
  int len = 0, keep = beg;
  for (int i = beg; i < NumCol; ++i) {
    if (RC[i] > threshold) {
      solver_ind[len++] = i;
    } else {
      node->IdxCols[keep++] = node->IdxCols[i];
    }
  }
  safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
}

void CVRP::getNewCstrCoeffByEdge(BBNODE *node, std::pair<int, int> edge, int *c_ind, double *c_val, int &num_nz) {
  //parents and their own maps can be read casually,
  //because this function is only used when selecting candidates in the SB stage after adding cuts
  if (If_in_Enu_State) {
    vector<int> cols(NumCol, 0);
    for (auto it : Map_Edge_ColIdx_in_Enu[edge]) ++cols[it];//will be automatically populated!
    int ccnt = 0;
    int numnz = 0;
    for (auto it : cols) {
      if (it) {
        c_ind[numnz] = ccnt;
        c_val[numnz++] = it;
      }
      ++ccnt;
    }
    num_nz = numnz;
  } else {
    auto p = node->Ptr;
    int code = edge.first * Dim + edge.second;
    int cnt = 0;
    int times;
    int record;
    while (p) {
      if (!p->Edge2Cols[code].empty()) {
        times = 0;
        record = p->Edge2Cols[code][0];
        for (int i : p->Edge2Cols[code]) {
          if (record == i) {
            ++times;
          } else {
            c_ind[cnt] = record;
            if (record < NumArtiVars) {
              c_val[cnt++] = times - 1;
              safe_Hyperparameter(times - 2)
            } else
              c_val[cnt++] = times;
            record = i;
            times = 1;
          }
        }
        c_ind[cnt] = record;
        if (record < NumArtiVars) {
          c_val[cnt++] = 1;
          safe_Hyperparameter(times - 2)
        } else
          c_val[cnt++] = times;
      }
      p = p->PNode;
    }
    num_nz = cnt;
  }
}

void CVRP::getCoeffRCC(BBNODE *node, RCC &rcc, int *c_ind, double *c_val, int &num_nz) const {
  //first column: artificial variable
  unordered_map<int, vector<int>> seq_map;
  seq_map.reserve(Dim * Dim);
  int curr_node, next_node;
  auto supp = new double[NumCol]();
  //also write 0 into the data set
  for (int n = node->NumParentCols; n < NumCol; ++n) {
    for (size_t m = node->IdxCols[n];;) {
      curr_node = ColPool4Pricing[m];
      next_node = ColPool4Pricing[++m];
      seq_map[curr_node * Dim + next_node].emplace_back(n);
      if (!next_node) break;
    }
  }

  if (rcc.FormRCC) {//normal case
    //customer_info needed to be std::stable_sorted
    auto &customer_info = rcc.InfoRCCCustomer;
    auto code = new int[customer_info.size() * customer_info.size() / 2 + 1]();
    //write seq into map
    int cnt = 0;
    for (auto i = customer_info.begin(); i != customer_info.end(); ++i) {
      auto j = i;
      ++j;
      for (; j != customer_info.end(); ++j) {
        int tmp1 = (*i) * Dim + (*j);
        int tmp2 = (*j) * Dim + (*i);
        for (auto it : seq_map[tmp1])++supp[it];
        for (auto it : seq_map[tmp2])++supp[it];
        code[cnt++] = min(tmp1, tmp2);
      }
    }

    //use support array to operate
    auto p = node->Ptr->PNode;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->Edge2Cols[code[i]]) {
          ++supp[j];
        }
      }
      p = p->PNode;
    }
    delete[]code;
  } else {//another case
    //x(\bar s, \bar s)+1/2(x(0,\bar s)-x(0,s))\leq RHS
    //no need to multiply two
    auto &customer_info = rcc.InfoRCCCustomer;
    auto &outside_customer_info = rcc.InfoRCCOutsideCustomer;
    auto code = new int[Dim * Dim]();
    //write seq into map
    int cnt = 0;
    for (auto i = outside_customer_info.begin(); i != outside_customer_info.end(); ++i) {
      auto j = i;
      ++j;
      for (; j != outside_customer_info.end(); ++j) {
        int tmp1 = (*i) * Dim + (*j);
        int tmp2 = (*j) * Dim + (*i);
        for (auto it : seq_map[tmp1])++supp[it];
        for (auto it : seq_map[tmp2])++supp[it];
        code[cnt++] = min(tmp1, tmp2);
      }
    }

    //use support array to operate
    auto p = node->Ptr->PNode;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->Edge2Cols[code[i]]) {
          ++supp[j];
        }
      }
      p = p->PNode;
    }

    //calculate x(0:\bar s) and x(0:s)
    cnt = 0;
    //for \bar s
    for (auto customer_it : outside_customer_info) {
      int tmp1 = customer_it;
      int tmp2 = customer_it * Dim;
      for (auto it : seq_map[tmp1])supp[it] += 0.5;
      for (auto it : seq_map[tmp2])supp[it] += 0.5;
      code[cnt++] = min(tmp1, tmp2);
    }
    p = node->Ptr->PNode;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->Edge2Cols[code[i]]) {
          supp[j] += 0.5;
        }
      }
      p = p->PNode;
    }

    //for s
    cnt = 0;
    for (auto customer_it : customer_info) {
      int tmp1 = customer_it;
      int tmp2 = customer_it * Dim;
      for (auto it : seq_map[tmp1])supp[it] -= 0.5;
      for (auto it : seq_map[tmp2])supp[it] -= 0.5;
      code[cnt++] = min(tmp1, tmp2);
    }
    p = node->Ptr->PNode;
    while (p) {
      for (int i = 0; i < cnt; ++i) {
        for (int j : p->Edge2Cols[code[i]]) {
          supp[j] -= 0.5;
        }
      }
      p = p->PNode;
    }
    delete[]code;
  }

  int cnt = 0;
  supp[0] = rcc.RHS;
  for (int i = 0; i < NumCol; ++i) {
    if (supp[i] != 0) {
      c_val[cnt] = supp[i];
      c_ind[cnt++] = i;
    }
  }
  num_nz = cnt;

  delete[]supp;
}

void CVRP::writeCol2Mem(BBNODE *node) {
  //start to replace all col index after parent column (including) with the value in in_memory
  if (checkMemPool()) reallocateMemPool();
  for (int i = node->NumParentCols; i < NumCol; ++i) {
    size_t start = PoolBeg4Mem;
    ColPool4Mem[PoolBeg4Mem++] = 0;
    for (size_t j = node->IdxCols[i] + 1;; ++j) {
      if (!(ColPool4Pricing[j]))break;
      ColPool4Mem[PoolBeg4Mem++] = ColPool4Pricing[j];
    }
    ColPool4Mem[PoolBeg4Mem++] = 0;
    node->IdxCols[i] = start;
  }
  safe_solver(node->solver.SOLVERupdatemodel())
//  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
//  vector<pair<size_t, double>> tmp_solindex;
//  for (int i = 0; i < NumCol; ++i)
//    if (X[i] > TOLERANCE) {
//      tmp_solindex.emplace_back(node->IdxCols[i], X[i]);
//    }
//  node->Idx4LPSolsInColPool = tmp_solindex;
//  node->NumParentColsInLPSols = (int) node->Idx4LPSolsInColPool.size();
}

void CVRP::constructMap(BBNODE *node, int beg) const {
  auto edge_map = &(node->Ptr->Edge2Cols);
  for (int i = beg; i < NumCol; ++i) {
    for (size_t j = node->IdxCols[i];; ++j) {
      int ai = ColPool4Mem[j], aj = ColPool4Mem[j + 1];
      if (ai > aj) {
        (*edge_map)[aj * Dim + ai].emplace_back(i);
      } else {
        (*edge_map)[ai * Dim + aj].emplace_back(i);
      }
      if (!aj) break;
    }
  }
  node->NumParentCols = NumCol;
}

///#####################################################################################
/// section
void CVRP::reformulateIPSol(BBNODE *node) {
  cout << "reformulateIPSol" << endl;
  node->if_Int = true;
  node->if_terminated = true;
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  if (ceil_transformed_number_related(LPVal - TOLERANCE) + TOLERANCE < UB) {
    UB = ceil_transformed_number_related(LPVal - TOLERANCE);
    IPOptSol.clear();
    for (int i = 0; i < node->NumParentColsInLPSols; ++i) {
      if ((node->Idx4LPSolsInColPool[i].second - 1) < TOLERANCE) {
        vector<int> tmp;
        tmp.reserve(MaxLengthEleRoute);
        tmp.emplace_back(0);
        for (auto j = node->Idx4LPSolsInColPool[i].first + 1;; ++j) {
          if (!ColPool4Mem[j]) break;
          tmp.emplace_back(ColPool4Mem[j]);
        }
        tmp.emplace_back(0);
        IPOptSol.emplace_back(std::move(tmp));
        for (auto j = node->Idx4LPSolsInColPool[i].first;;) {
          int p = ColPool4Mem[j], q = ColPool4Mem[j + 1];
          if (p > q)
            ArcGraph_revised[q][p] = 0;
          else
            ArcGraph_revised[p][q] = 0;
          if (!ColPool4Mem[++j]) break;
        }
      }
    }
    for (int i = node->NumParentColsInLPSols; i < node->Idx4LPSolsInColPool.size(); ++i) {
      if ((node->Idx4LPSolsInColPool[i].second - 1) < TOLERANCE) {
        vector<int> tmp;
        for (auto j = node->Idx4LPSolsInColPool[i].first + 1;; ++j) {
          if (!ColPool4Pricing[j]) break;
          tmp.emplace_back(ColPool4Pricing[j]);
        }
        tmp.emplace_back(0);
        IPOptSol.emplace_back(std::move(tmp));
        for (auto j = node->Idx4LPSolsInColPool[i].first;;) {
          int p = ColPool4Pricing[j], q = ColPool4Pricing[j + 1];
          if (p > q)
            ArcGraph_revised[q][p] = 0;
          else
            ArcGraph_revised[p][q] = 0;
          if (!ColPool4Pricing[++j]) break;
        }
      }
    }
    for (int j = 1; j < Dim; ++j) {
      if (ArcGraph_revised[0][j] == 0) continue;
      vector<int> tmp;
      tmp.reserve(MaxLengthEleRoute);
      tmp.emplace_back(0);
      tmp.emplace_back(j);
      ArcGraph_revised[0][j] = 0;
      int tmp_j = j;
      while (true) {
        bool ifbreak = false;
        for (int k = 0; k < tmp_j; ++k) {
          if (ArcGraph_revised[k][tmp_j] != 0) {
            ArcGraph_revised[k][tmp_j] = 0;
            tmp.emplace_back(k);
            tmp_j = k;
            ifbreak = true;
            break;
          }
        }
        if (!ifbreak) {
          for (int k = tmp_j + 1; k < Dim; ++k) {
            if (ArcGraph_revised[tmp_j][k] != 0) {
              ArcGraph_revised[tmp_j][k] = 0;
              tmp.emplace_back(k);
              tmp_j = k;
              break;
            }
          }
        }
        if (tmp_j == 0) {
          break;
        }
      }
      IPOptSol.emplace_back(std::move(tmp));
    }
#ifdef MASTER_VALVE_ML
    updateUB_EdgeSol();
#endif
  }
}

///#####################################################################################
/// bi_directional

void CVRP::convertVertex2R1CsInPricing(BBNODE *node) {
  Vertex2ActiveInOnePricingR1Cs.assign(Dim, tuple<vector<int>, R1CMem, vector<int >>());
  int cnt = 0;
  //rank the dual according to the dual min is the best!
  unordered_map<int, double> cnt_dual;
  for (auto &r1c : node->R1Cs) {
    if (abs(Pi4Labeling[r1c.IdxR1C]) < DUAL_TOLERANCE) continue;
    cnt_dual.emplace(cnt, Pi4Labeling[r1c.IdxR1C]);
    for (auto i : r1c.InfoR1C) {
      get<0>(Vertex2ActiveInOnePricingR1Cs[i]).emplace_back(cnt);
      get<1>(Vertex2ActiveInOnePricingR1Cs[i]).set(cnt);
      get<2>(Vertex2ActiveInOnePricingR1Cs[i]).emplace_back(cnt);
    }
    for (auto &v : r1c.Mem) {
      get<1>(Vertex2ActiveInOnePricingR1Cs[v]).set(cnt);
      get<2>(Vertex2ActiveInOnePricingR1Cs[v]).emplace_back(cnt);
    }
    ++cnt;
  }
  for (int i = 1; i < Dim; ++i) {
    // sort the get<2>(Vertex2ActiveInOnePricingR1Cs[i]) by cnt_dual
    vector<pair<int, double>> tmp(get<2>(Vertex2ActiveInOnePricingR1Cs[i]).size());
    transform(get<2>(Vertex2ActiveInOnePricingR1Cs[i]).begin(),
              get<2>(Vertex2ActiveInOnePricingR1Cs[i]).end(), tmp.begin(),
              [&cnt_dual](const auto &a) { return make_pair(a, cnt_dual[a]); });
    sort(tmp.begin(), tmp.end(), [](const auto &a, const auto &b) {
      return a.second < b.second;
    });
    transform(tmp.begin(), tmp.end(), get<2>(Vertex2ActiveInOnePricingR1Cs[i]).begin(),
              [](const auto &a) { return a.first; });
  }

//  cout << "-----------------------------------------------" << endl;
  cnt_dual.clear();
  cnt = 0;
  Vertex2ActiveInOnePricingR1C_multi.assign(Dim, tuple<vector<tuple<int, int, int>>, vector<int>, vector<int>>());
  R1C_multi_denominator_InCG.clear();
  for (auto &r1c : node->R1Cs_multi) {
    if (abs(Pi4Labeling[r1c.IdxR1C]) < DUAL_TOLERANCE) continue;
    const auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    R1C_multi_denominator_InCG.emplace_back(denominator);
    int count = 0;
    for (auto &i : r1c.InfoR1C.first) {
      get<0>(Vertex2ActiveInOnePricingR1C_multi[i]).emplace_back(cnt, multi[count], denominator);
      get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).emplace_back(cnt);
      ++count;
    }
    for (auto &v : r1c.Mem) {
      get<2>(Vertex2ActiveInOnePricingR1C_multi[v]).emplace_back(cnt);
    }
    cnt_dual.emplace(cnt, Pi4Labeling[r1c.IdxR1C]);
    ++cnt;
  }
  NumValidR1C_multi_InCG = cnt;
  auto if_use = new bool[NumValidR1C_multi_InCG];
  for (int i = 1; i < Dim; ++i) {
    memset(if_use, 0, sizeof(bool) * NumValidR1C_multi_InCG);
    for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[i])) {
      if_use[l] = true;
    }
    for (int j = 0; j < NumValidR1C_multi_InCG; ++j) {
      if (!if_use[j]) get<1>(Vertex2ActiveInOnePricingR1C_multi[i]).emplace_back(j);
    }
  }
  delete[] if_use;

  for (int i = 1; i < Dim; ++i) {
    // sort the get<2>(Vertex2ActiveInOnePricingR1C_multi[i]) by cnt_dual
    vector<pair<int, double>> tmp(get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).size());
    transform(get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).begin(),
              get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).end(), tmp.begin(),
              [&cnt_dual](const auto &a) { return make_pair(a, cnt_dual[a]); });
    sort(tmp.begin(), tmp.end(), [](const auto &a, const auto &b) {
      return a.second < b.second;
    });
    transform(tmp.begin(), tmp.end(), get<2>(Vertex2ActiveInOnePricingR1C_multi[i]).begin(),
              [](const auto &a) { return a.first; });
  }
}

void CVRP::convertVertex2R1CsInOneLP(BBNODE *node) {
  Vertex2AllInOneLPR1Cs.assign(Dim, pair<vector<int>, R1CMem>());
  int cnt = 0;
  for (auto &r1c : node->R1Cs) {
    for (auto i : r1c.InfoR1C) {
      Vertex2AllInOneLPR1Cs[i].first.emplace_back(cnt);
      Vertex2AllInOneLPR1Cs[i].second.set(cnt);
    }
    for (auto &v : r1c.Mem) {
      Vertex2AllInOneLPR1Cs[v].second.set(cnt);
    }
    ++cnt;
  }
  cnt = 0;
  Vertex2AllInOneLPR1C_multi.assign(Dim, pair<vector<pair<int, int>>, vector<int>>());
  R1C_multi_denominator_InLP.clear();
  vector<vector<int>> tmp(Dim, vector<int>());
  for (auto &r1c : node->R1Cs_multi) {
    const auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    R1C_multi_denominator_InLP.emplace_back(denominator);
    int count = 0;
    for (auto &i : r1c.InfoR1C.first) {
      Vertex2AllInOneLPR1C_multi[i].first.emplace_back(cnt, multi[count]);
      tmp[i].emplace_back(cnt);
      ++count;
    }
    for (auto &v : r1c.Mem) {
      tmp[v].emplace_back(cnt);
    }
    ++cnt;
  }
  auto if_use = new bool[cnt];
  for (int i = 1; i < Dim; ++i) {
    memset(if_use, 0, sizeof(bool) * cnt);
    for (auto l : tmp[i]) {
      if_use[l] = true;
    }
    for (int j = 0; j < cnt; ++j) {
      if (!if_use[j]) Vertex2AllInOneLPR1C_multi[i].second.emplace_back(j);
    }
  }
  delete[] if_use;
#ifdef checkR1CTotally
  cout << "check if cuts are not well defined!" << endl;
  for (auto &r1c : node->R1Cs) {
    yzzLong tmp2 = 0;
    for (auto &v : r1c.InfoR1C) {
      if (tmp2.test(v)) {
        throw std::runtime_error("duplicate error in cut self!");
      }
      tmp2.set(v);
    }
    for (auto &v : r1c.Mem) {
      if (tmp2.test(v)) {
        throw std::runtime_error("duplicate error in mem!");
      }
      tmp2.set(v);
    }
  }
  for (auto &r1c : node->R1Cs_multi) {
    yzzLong tmp2 = 0;
    for (auto &v : r1c.InfoR1C.first) {
      if (tmp2.test(v)) {
        throw std::runtime_error("duplicate error in mul_cut self!");
      }
      tmp2.set(v);
    }
    for (auto &v : r1c.Mem) {
      if (tmp2.test(v)) {
        throw std::runtime_error("duplicate error in mul_mem!");
      }
      tmp2.set(v);
    }
  }
#endif
}

void CVRP::cleanAllPtr(BBNODE *const node, int mode, bool if_clear_concatenate) {
  if (mode == 1) {
    for (int i = 0; i < Dim; ++i) {
      for (int b = 0; b < NumBucketsPerVertex; ++b) {
        LabelArrayInForwardSense[i][b].second = 0;
        IfExistExtraLabelsInForwardSense[i][b].second = 0;
      }
    }
  }
#ifdef SYMMETRY_PROHIBIT
  else if (mode == 2) {
    for (int i = 0; i < Dim; ++i) {
      for (int b = 0; b < NumBucketsPerVertex; ++b) {
        LabelArrayInBackwardSense[i][b].second = 0;
        IfExistExtraLabelsInBackwardSense[i][b].second = 0;
      }
    }
  }
#endif
  if (if_clear_concatenate) {
    concatenateLabelsInForwardCG.clear();
#ifdef SYMMETRY_PROHIBIT
    concatenateLabelsInBackwardCG.clear();
#endif
  }
}

double CVRP::calculateOptGap(BBNODE *const node) const {
  return (UB - node->Val + roundUpTolerance);
}

void CVRP::tellIfArcElimination(BBNODE *node) {
  double now_gap, required_gap;
  if (!if_can_arc_eliminationByExactCG) {
    final_decision_4_arc_elimination = false;
    cerr << "Error: if_can_arc_eliminationByExactCG is false, but tries to call arcElimination!" << endl;
    exit(-1);
  } else if_can_arc_eliminationByExactCG = false;
  if (if_stopArcElimination) {
    final_decision_4_arc_elimination = false;
    if_stopArcElimination = false;
    cout << "ArcEliminationNEnumeration is banned!" << endl;
    goto QUIT;
  }
  if (final_decision_4_arc_elimination) {
    cout << "ArcEliminationNEnumeration must be performed!" << endl;
    goto QUIT;
  }
  if (abs(OldUB - UB) < TOLERANCE) {
    cout << "OldUB is updated from " << OldUB << " to " << UB << endl;
    cout << "ArcEliminationNEnumeration must be performed!" << endl;
    OldUB = UB;
    final_decision_4_arc_elimination = true;
    goto QUIT;
  }
  now_gap = (UB - node->Val) / UB;
//  required_gap = Gap4LastFailDue2DifficultyInArcElimination * CONFIG::GapRequirementReduction4Arc;
//  if (required_gap < now_gap) {
//    cout << "now_gap= " << now_gap << " larger than Gap4LastFailDue2DifficultyInArcElimination= " << required_gap
//         << endl;
//    goto QUIT;
//  }
  if (node->LastGap / now_gap > CONFIG::GapImproved4ArcEliminationNEnumeration
      || now_gap < GapTolerance4ArcEliminationNEnumeration) {
    node->LastGap = now_gap;
    cout << "gap updated as " << now_gap << endl;
    final_decision_4_arc_elimination = true;
    final_decision_4_enumeration = true;
    goto QUIT;
  }
  QUIT:
  return;
}

void CVRP::tellIfEnumeration(BBNODE *node) {
  double gap;
  if (final_decision_4_enumeration) {
    cout << "Enumeration must be performed!" << endl;
    enumerationMode = !if_ArcEliminationSucceed;
    if (!if_ArcEliminationSucceed) {
      cout << "But ArcElimination is failed! So skip enumeration!" << endl;
      final_decision_4_enumeration = false;
    }
    goto QUIT;
  }
  if (!if_ArcEliminationSucceed && !if_ArcEliminationTriedButFailed) {
    cout << "Arc elimination is even not tried!" << endl;
    goto QUIT;
  }

  if (abs(OldUB - UB) < TOLERANCE) {
    cout << "OldUB is updated from " << OldUB << " to " << UB << endl;
    cout << "ArcEliminationNEnumeration must be performed!" << endl;
    OldUB = UB;
    final_decision_4_enumeration = true;
    enumerationMode = !if_ArcEliminationSucceed;
    goto QUIT;
  }

  gap = (UB - node->Val) / UB;
//  if (node->LastGap / gap > CONFIG::GapImproved4ArcEliminationNEnumeration
//      || gap < GapTolerance4ArcEliminationNEnumeration) {
//    node->LastGap = gap;
//    cout << "gap updated as " << gap << endl;
//    final_decision_4_enumeration = true;
//    enumerationMode = !if_ArcEliminationSucceed;
//    goto QUIT;
//  }
  if (gap > CONFIG::MaxGap2TryEnumeration || gap > LastEnumerationFailGap) {
    goto QUIT;
  }
  final_decision_4_enumeration = true;
  enumerationMode = !if_ArcEliminationSucceed;
  if (if_ArcEliminationTriedButFailed) {
    if (Count4Tolerance4tryEnumerationWhenArcEliminationFails
        <= CONFIG::InitialTolerance4tryEnumerationWhenArcEliminationFails) {
      ++Count4Tolerance4tryEnumerationWhenArcEliminationFails;//if succeed in enumeration, this will be reset!
    } else {
      final_decision_4_enumeration = false;
    }
  }
  QUIT:
  return;
}

//void CVRP::rollback2PreState(BBNODE *const node, const std::vector<RCC> &old_rcc,
//                             const std::vector<R1C> &old_r1c,
//                             const std::vector<R1C_multi> &old_r1c_multi,
//                             const std::vector<BrC> &old_brc, const std::vector<size_t> &old_col_idx) {
//  throw runtime_error("rollback2PreState has been banned! Please use rollbackEasyWay!");
//  HardRollBackFactor = CONFIG::HardTimeThresholdInPricing / LastMaxTimeLabeling * 1.1;
//  cout << "HardRollBackFactor: " << HardRollBackFactor << endl;
//#ifdef CutRollBackNotAllowed
//  cout << "Error: CutRollBackNotAllowed is defined, but rollback2PreState is called!" << endl;
//  exit(0);
//#endif
//  cout << "Very bad behavior! because the lp is rollback2preState!" << endl;
//  cout << "here we construct a new lp from scratch, and we solve it using column generation to get the optimal dual!"
//       << endl;
//  cout << "and if we even fail here, we conclude fail to this instance!" << endl;
//  cout << "here we only not move the first col" << endl;
//
//  node->RCCs = old_rcc;
//  node->R1Cs = old_r1c;
//  node->R1Cs_multi = old_r1c_multi;
//  node->BrCs = old_brc;
//  copy(old_col_idx.begin(), old_col_idx.end(), node->IdxCols);
//  NumRow = Dim + int(node->RCCs.size() + node->R1Cs.size() + node->R1Cs_multi.size() + node->BrCs.size());
//  NumCol = (int) old_col_idx.size();
//  //re-generate this (necessary)
//  convertVertex2R1CsInOneLP(node);
//  ////////////////////////////////
//  RowMatrixXd mat = RowMatrixXd::Zero(NumRow, NumCol);
//  vector<double> obj(NumCol);
//  unordered_map<int, vector<int>> vertex_map;
//  vertex_map.reserve((int) pow_self(Dim, 2));
//  vector<int> r1c_eff(node->R1Cs.size(), 0);
//  vector<int> r1c_multi_eff(node->R1Cs_multi.size(), 0);
//  vector<int> r1c_multi_state(node->R1Cs_multi.size(), 0);
//  R1CMem r1cMem = 0;
//  //basic constraints & r1c
//  int ccnt = 0;
//  for (int i = 0; i < node->NumParentCols; ++i) {
//    int past_node = 0;
//    double cost_sum = 0;
//    for (auto j = node->IdxCols[i] + 1;; ++j) {
//      int curr_node = ColPool4Mem[j];
//      vertex_map[past_node * Dim + curr_node].emplace_back(ccnt);
//      cost_sum += CostMat4Vertex[past_node][curr_node];
//      if (!curr_node) break;
//      ++mat(curr_node - 1, ccnt);
//      for (auto l : Vertex2AllInOneLPR1Cs[curr_node].first) {
//        if (r1cMem[l]) {
//          r1cMem[l] = false;
//          ++r1c_eff[l];
//        } else r1cMem[l] = true;
//      }
//      r1cMem &= Vertex2AllInOneLPR1Cs[curr_node].second;
//
//      //r1c_multi
//      for (auto &l : Vertex2AllInOneLPR1C_multi[curr_node].first) {
//        int tmp_cut = l.first;
//        r1c_multi_state[tmp_cut] += l.second;
//        if (r1c_multi_state[tmp_cut] >= R1C_multi_denominator_InLP[tmp_cut]) {
//          r1c_multi_state[tmp_cut] -= R1C_multi_denominator_InLP[tmp_cut];
//          ++r1c_multi_eff[tmp_cut];
//        }
//      }
//
//      for (auto l : Vertex2AllInOneLPR1C_multi[curr_node].second) r1c_multi_state[l] = 0;
//
//      past_node = curr_node;
//    }
//    for (int k = 0; k < r1c_eff.size(); ++k) {
//      mat(node->R1Cs[k].IdxR1C, ccnt) = r1c_eff[k];
//    }
//    for (int j = 0; j < r1c_multi_eff.size(); ++j) {
//      if (r1c_multi_eff[j]) {
//        mat(node->R1Cs_multi[j].IdxR1C, ccnt) = r1c_multi_eff[j];
//      }
//    }
//    memset(r1c_eff.data(), 0, sizeof(int) * r1c_eff.size());
//    memset(r1c_multi_eff.data(), 0, sizeof(int) * r1c_multi_eff.size());
//    memset(r1c_multi_state.data(), 0, sizeof(int) * r1c_multi_state.size());
//    r1cMem = 0;
//    obj[ccnt] = cost_sum;
//    ++ccnt;
//  }
//  for (int i = node->NumParentCols; i < NumCol; ++i) {
//    int past_node = 0;
//    double cost_sum = 0;
//    for (auto j = node->IdxCols[i] + 1;; ++j) {
//      int curr_node = ColPool4Pricing[j];
//      vertex_map[past_node * Dim + curr_node].emplace_back(ccnt);
//      cost_sum += CostMat4Vertex[past_node][curr_node];
//      if (!curr_node) break;
//      ++mat(curr_node - 1, ccnt);
//      for (auto l : Vertex2AllInOneLPR1Cs[curr_node].first) {
//        if (r1cMem[l]) {
//          r1cMem[l] = false;
//          ++r1c_eff[l];
//        } else r1cMem[l] = true;
//      }
//      r1cMem &= Vertex2AllInOneLPR1Cs[curr_node].second;
//
//      //r1c_multi
//      for (auto &l : Vertex2AllInOneLPR1C_multi[curr_node].first) {
//        int tmp_cut = l.first;
//        r1c_multi_state[tmp_cut] += l.second;
//        if (r1c_multi_state[tmp_cut] >= R1C_multi_denominator_InLP[tmp_cut]) {
//          r1c_multi_state[tmp_cut] -= R1C_multi_denominator_InLP[tmp_cut];
//          ++r1c_multi_eff[tmp_cut];
//        }
//      }
//
//      for (auto l : Vertex2AllInOneLPR1C_multi[curr_node].second) r1c_multi_state[l] = 0;
//
//      past_node = curr_node;
//    }
//    for (int k = 0; k < r1c_eff.size(); ++k) {
//      mat(node->R1Cs[k].IdxR1C, ccnt) = r1c_eff[k];
//    }
//    for (int j = 0; j < r1c_multi_eff.size(); ++j) {
//      mat(node->R1Cs_multi[j].IdxR1C, ccnt) = r1c_multi_eff[j];
//    }
//    memset(r1c_eff.data(), 0, sizeof(int) * r1c_eff.size());
//    memset(r1c_multi_eff.data(), 0, sizeof(int) * r1c_multi_eff.size());
//    memset(r1c_multi_state.data(), 0, sizeof(int) * r1c_multi_state.size());
//    r1cMem = 0;
//    obj[ccnt] = cost_sum;
//    ++ccnt;
//  }
//  //vehicle constraintscnt
//  mat.row(RealDim) = RowVectorXd::Ones(NumCol);
//  //rcc constraints
//  addRCC2LP(node, vertex_map, mat);
//  //brc constraints
//  for (auto &br : node->BrCs) {
//    int ai = br.Edge.first;
//    int aj = br.Edge.second;
//    int idx = br.IdxBrC;
//    for (auto it_map : vertex_map[ai * Dim + aj]) ++mat(idx, it_map);
//    for (auto it_map : vertex_map[aj * Dim + ai]) ++mat(idx, it_map);
//  }
//  ///////////////////////
//  //basic constraints
//  vector<char> sense(NumRow, SOLVER_EQUAL);
//  vector<double> rhs(NumRow, 1);
//  //vehicle
//  sense[RealDim] = SOLVER_GREATER_EQUAL;
//  rhs[RealDim] = K;
//  //cuts
//  for (auto &r1c : node->R1Cs) {
//    sense[r1c.IdxR1C] = SOLVER_LESS_EQUAL;
//    rhs[r1c.IdxR1C] = r1c.RHS;
//  }
//  for (auto &r1c_multi : node->R1Cs_multi) {
//    sense[r1c_multi.IdxR1C] = SOLVER_LESS_EQUAL;
//    rhs[r1c_multi.IdxR1C] = r1c_multi.RHS;
//  }
//  for (auto &rcc : node->RCCs) {
//    sense[rcc.IdxRCC] = SOLVER_LESS_EQUAL;
//    rhs[rcc.IdxRCC] = rcc.RHS;
//  }
//  for (auto &brc : node->BrCs) {
//    if (brc.BrDir) {
//      sense[brc.IdxBrC] = SOLVER_GREATER_EQUAL;
//      rhs[brc.IdxBrC] = 1;
//    } else {
//      sense[brc.IdxBrC] = SOLVER_LESS_EQUAL;
//      rhs[brc.IdxBrC] = 0;
//    }
//  }
//  //first we free the model.
//  auto name = "RollBack_" + to_string(CGRollBackTimes++) + ".lp";
//  cout << "new built model name is " << name << endl;
//  node->solver.SOLVERgetenv(&Solver);//need load environment
//  safe_solver(node->solver.SOLVERfreemodel())
//  //then we create a new model.
//  safe_solver(node->solver.SOLVERnewmodel(name.c_str(), 0, nullptr, nullptr, nullptr, nullptr, nullptr))
//  safe_solver(node->solver.SOLVERaddconstrs(NumRow, 0, nullptr, nullptr, nullptr, sense.data(), rhs.data(), nullptr))
//  //add cols in MIP
//  size_t nzcnt = 0;
//  //set the first col as the same with the rhs
//  for (int i = 0; i < NumRow; ++i) mat(i, 0) = rhs[i];
//  for (int i = 0; i < NumCol; ++i) {
//    solver_beg[i] = nzcnt;
//    solver_obj[i] = obj[i];
//    for (int j = 0; j < NumRow; ++j) {
//      if (mat(j, i) != 0) {
//        solver_ind[nzcnt] = j;
//        solver_val[nzcnt] = mat(j, i);
//        ++nzcnt;
//      }
//    }
//  }
//  solver_beg[NumCol] = nzcnt;
//  safe_solver(node->solver.SOLVERXaddvars(ccnt,
//                                          nzcnt,
//                                          solver_beg,
//                                          solver_ind,
//                                          solver_val,
//                                          solver_obj,
//                                          nullptr,
//                                          nullptr,
//                                          nullptr,
//                                          nullptr))
//  safe_solver(node->solver.SOLVERupdatemodel())
////  safe_solver(node->solver.SOlVERsetVBasis(0, NumCol, node->VBasis.data()))
////  safe_solver(node->solver.SOLVERsetCBasis(0, NumRow, node->CBasis.data()))
//  safe_solver(node->solver.SOLVERoptimize())
//
//  Rollback = 0;
//  solveLPInLabeling(node);
//  cout << "Rollback: " << Rollback << endl;
//  cout << "node->Val= " << node->Val << endl;
//  safe_Hyperparameter(Rollback == 1)
//
//  CutGenTimeThresholdInPricing *= CONFIG::CutGenTimeThresholdInPricing_reduce_factor;
//}

void CVRP::rollbackEaseWay(BBNODE *const node, int old_num) {
  if (If_in_Enu_State) throw runtime_error("Error: rollbackEaseWay is called in Enu_State!");
  cout << "rollbackEaseWay" << endl;
  HardRollBackFactor = CONFIG::HardTimeThresholdInPricing / LastMaxTimeLabeling * 1.1;
  CutGenTimeThresholdInPricing *= CONFIG::CutGenTimeThresholdInPricing_reduce_factor;
  cout << "HardRollBackFactor: " << HardRollBackFactor << " CutGenTimeThresholdInPricing: "
       << CutGenTimeThresholdInPricing << endl;
#ifdef CutRollBackNotAllowed
  throw runtime_error("Error: CutRollBackNotAllowed is defined, but rollback2PreState is called!");
#endif
  vector<int> delete_cuts(NumRow - old_num);
  iota(delete_cuts.begin(), delete_cuts.end(), old_num);
  deleteNonactiveCuts(node, delete_cuts);

  // reset the memory
  for (auto &i : ResetCutMem) {
    auto &d = get<2>(i);
    std::vector<int> tmp(d.count());
    int idx = 0;
    for (int j = 1; j < Dim; ++j) {
      if (d.test(j)) tmp[idx++] = j;
    }
    if (get<0>(i)) {
      node->R1Cs[get<1>(i)].Mem = std::move(tmp);
    } else {
      node->R1Cs_multi[get<1>(i)].Mem = std::move(tmp);
    }
  }

  ResetCutMem.clear();

  convertVertex2R1CsInOneLP(node);
  safe_solver(node->solver.SOLVERoptimize())
  Rollback = 0;
  // need resolve since we have changed the mem of old cuts, we do not know if rollback happens again
  solveLPInLabeling(node);
  safe_Hyperparameter(Rollback == 1)
}

void CVRP::initialBucketGraph4Node(BBNODE *node) {
  vector<int> bucketArc(RealDim - 1);
  for (int i = 1; i < Dim; ++i) {
    int tmp = 0;
    for (int j = 1; j < i; ++j) bucketArc[tmp++] = j;
    for (int j = i + 1; j < Dim; ++j) bucketArc[tmp++] = j;
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      auto &bucket = node->AllForwardBuckets[i][b];
      bucket.BucketArcs = bucketArc;
      bucket.i = i;
    }
  }
  auto &bucket = node->AllForwardBuckets[0][0];
  bucket.BucketArcs.resize(RealDim);
  iota(bucket.BucketArcs.begin(), bucket.BucketArcs.end(), 1);

  MaxNumForwardGraphArc = NumBucketsPerVertex * (RealDim - 1) * RealDim;
  node->NumForwardBucketArcs = MaxNumForwardGraphArc;
#ifdef SYMMETRY_PROHIBIT
  for (int i = 1; i < Dim; ++i) {
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      node->AllBackwardBuckets[i][b] = node->AllForwardBuckets[i][0];
    }
  }
  auto &bucket_back = node->AllBackwardBuckets[0][0];
    bucket_back=bucket;
  MaxNumBackwardGraphArc = NumBucketsPerVertex * (RealDim - 1) * RealDim;
  node->NumBackwardBucketArcs = MaxNumBackwardGraphArc;
#endif
}

double CVRP::ceil_transformed_number_related(double x) const {
  for (int i = 0; i < transformed_number; ++i) {
    x *= 10;
  }
  x = std::ceil(x);
  for (int i = 0; i < transformed_number; ++i) {
    x /= 10;
  }
  return x;
}

//split string
void splitString(const string &str, const string &splits, vector<string> &res) {
  if (str.empty()) return;
  auto strs = str + splits;
  size_t pos = strs.find(splits);
  int step = (int) splits.size();
  while (pos != std::string::npos) {
    string tmp = strs.substr(0, pos);
    res.push_back(tmp);
    strs = strs.substr(pos + step, strs.size());
    pos = strs.find(splits);
  }
}

string readInstanceFile(const std::string &file_name, int line) {
  safe_Hyperparameter(line == ReadNoLine)
  fstream file(file_name);
  if (!file.is_open()) {
    cout << "File not found!" << endl;
    exit(EXIT_FAILURE);
  }
  file.seekg(std::ios::beg);
  for (int i = 0; i < line; ++i) {
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  string line_str;
  getline(file, line_str);
  vector<string> strList;
  splitString(line_str, " ", strList);
  if (strList.size() == 2) {
    CONFIG::UB = strtod(strList[1].c_str(), nullptr);
  }
#ifdef readEnumerationTrees
  else if (strList.size() == 4) {
    CONFIG::UB = strtod(strList[1].c_str(), nullptr);
    CONFIG::colPool_path = strList[2].c_str();
    CONFIG::tree_path = strList[3].c_str();
  }
#endif
  file.close();
  return strList[0];
}

string generateInstancePath(int argc, char *argv[]) {
  int n = ReadNoLine;
  string d;
  bool if_type1 = false;
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-d" && i + 1 < argc) {
      d = argv[++i];
      if_type1 = true;
    } else if (arg == "-n" && i + 1 < argc) {
      stringstream convert(argv[++i]);
      if (!(convert >> n)) {
        // error handling: not a number
        printf("Invalid number: %s\n", argv[i]);
        continue;
      }
      if_type1 = true;
    } else if (arg == "-u" && i + 1 < argc) {
      stringstream convert(argv[++i]);
      if (!(convert >> CONFIG::UB)) {
        // error handling: not a number
        printf("Invalid number: %s\n", argv[i]);
        continue;
      }
    }
#ifdef readEnumerationTrees
      else if (arg == "-c" && i + 1 < argc) {
        stringstream convert(argv[++i]);
        if (!(convert >> CONFIG::colPool_path)) {
          // error handling: not a number
          printf("Invalid string: %s\n", argv[i]);
          continue;
        }
      } else if (arg == "-t" && i + 1 < argc) {
        stringstream convert(argv[++i]);
        if (!(convert >> CONFIG::tree_path)) {
          // error handling: not a number
          printf("Invalid string: %s\n", argv[i]);
          continue;
        }
      }
#endif
    else {
      printf("Unknown option: %s\n", argv[i]);
    }
  }

  if (if_type1) {
    return readInstanceFile(d, n);
  } else {
    if (argc > 1) {
      return argv[1];
    } else {
      // Handle case where no extra arguments were passed
      throw std::invalid_argument("Not enough arguments");
    }
  }
}

int CVRP::inverseLastBranchconstr(char sense, double rhs, SOLVER &solver) {
  int error = solver.SOLVERupdatemodel();
  if (error) return error;
  int idxconstrs;
  error += solver.SOLVERgetNumRow(&idxconstrs);
  if (error) return error;
  --idxconstrs;
  error += solver.SOLVERsetRhs(idxconstrs, 1, &sense, &rhs);
  if (error) return error;
  int vind = 0;
  error += solver.SOLVERXchgcoeffs(1, &idxconstrs, &vind, &rhs);
  return error;
}

int CVRP::addBranchconstr(int numnz,
                          int *cind,
                          double *cval,
                          char sense,
                          double rhs,
                          const char *constrname,
                          SOLVER &solver) {
  int error = solver.SOLVERaddconstr(numnz, cind, cval, sense, rhs, constrname);
  if (error) return error;
  error += solver.SOLVERupdatemodel();
  if (error) return error;
  int idxconstrs;
  error += solver.SOLVERgetNumRow(&idxconstrs);
  if (error) return error;
  --idxconstrs;
  //only consider the first column
  int vind = 0;
  error += solver.SOLVERXchgcoeffs(1, &idxconstrs, &vind, &rhs);
  return error;
}

void CVRP::chgBranchconstr(double *val,
                           int *cind,
                           int *vind,
                           int numnz,
                           int *old_ind,
                           double *old_val,
                           char sense,
                           double rhs,
                           SOLVER &solver) {
  memset(val, 0, sizeof(double) * NumCol);
  for (int i = 0; i < numnz; ++i) {
    val[old_ind[i]] = old_val[i];
  }
  safe_solver(solver.SOLVERchgcoeffs(NumCol, cind, vind, val))
  int BeforeNumRow = cind[0];
  safe_solver(solver.SOLVERsetRhs(BeforeNumRow, 1, &sense, &rhs))
  int vind2 = 0;
  safe_solver(solver.SOLVERXchgcoeffs(1, &BeforeNumRow, &vind2, &rhs))
}

float sqrt_self(float x) {
  float xhalf = 0.5f * x;
  int i = *(int *) &x; // get bits for floating VALUE
  i = 0x5f375a86 - (i >> 1); // gives initial guess y0
  x = *(float *) &i; // convert bits BACK to float
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
  return 1 / x;
}

double pow_self(double x, int n) {
  //x^n and n needs to be positive integer >= 0
  if (n == 0) return 1;
  double y = x;
  for (int i = 1; i < n; ++i) {
    y *= x;
  }
  return y;
}

///this is version that select the same number of cols in each vertex
//void CVRP::giveInitialUBByMIP(BBNODE *&node) {
//
//  cout << "ColSetUsedInRootMIP  Size: " << ColSetUsedInRootMIP.size() << endl;
//  int numcol = (int) ColSetUsedInRootMIP.size(), numrow = NumRow;
//  vector<pair<vector<int>, double >> vec_ColSetUsedInRootMIP(numcol);
//  int ccnt = 0;
//  for (auto &col : ColSetUsedInRootMIP) {
//    vec_ColSetUsedInRootMIP[ccnt++] = col.second;
//  }
//  if (ColSetUsedInRootMIP.size() > CONFIG::MaxNumRoute4RootMIP) {
//    //now we only consider use the dual from the essential routes
//    RowMatrixXd mat = RowMatrixXd::Zero(numrow, numcol);
//    RowVectorXd obj(numcol);
//    unordered_map<int, vector<int>> vertex_map;
//    vertex_map.reserve((int) pow_self(Dim, 2));
//    vector<int> r1c_eff(node->R1Cs.size(), 0);
//    R1CMem r1cMem = 0;
//    //basic constraints & r1c3
//    ccnt = 0;
//    size_t sum_length = 0;
//    for (auto &col : vec_ColSetUsedInRootMIP) {
//      auto &vec = col.first;
//      sum_length += vec.size();
//      obj[ccnt] = col.second;
//      int past_node = 0;
//      for (auto i : vec) {
//        mat(i - 1, ccnt) = 1;
//        vertex_map[past_node * Dim + i].emplace_back(ccnt);
//        for (auto l : Vertex2AllInOneLPR1Cs[i].first) {
//          if (r1cMem[l]) {
//            r1cMem[l] = false;
//            ++r1c_eff[l];
//          } else r1cMem[l] = true;
//        }
//        r1cMem &= Vertex2AllInOneLPR1Cs[i].second;
//        past_node = i;
//      }
//      for (int i = 0; i < node->R1Cs.size(); ++i) {
//        mat(node->R1Cs[i].IdxR1C, ccnt) = r1c_eff[i];
//      }
//      fill(r1c_eff.begin(), r1c_eff.end(), 0);
//      r1cMem = 0;
//      vertex_map[past_node * Dim].emplace_back(ccnt);
//      ++ccnt;
//    }
//    sum_length /= ccnt;
//    size_t base_line = CONFIG::MaxNumRoute4RootMIP * sum_length / Dim;
//    //vehicle constraints
//    mat.row(RealDim) = RowVectorXd::Ones(numcol);
//    //rcc constraints
//    addRCC2LP(node, vertex_map, mat);
//    //r1c1 constraints (all zero)
//    RowVectorXd rc = obj - DualInRootMIP * mat;
//    vector<vector<pair<int, double>>> vertex_statistic(Dim);
//    ccnt = 0;
//    for (auto &col : vec_ColSetUsedInRootMIP) {
//      auto &vec = col.first;
//      for (int i : vec) {
//        vertex_statistic[i].emplace_back(ccnt, rc[ccnt]);
//      }
//      ++ccnt;
//    }
//    set<int> cols_2_kept;
//    for (auto &vec : vertex_statistic) {
//      sort(vec.begin(), vec.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
//        return a.second < b.second;
//      });
//      auto limit = min(vec.size(), base_line);
//      for (size_t i = 0; i < limit; ++i) {
//        cols_2_kept.emplace(vec[i].first);
//      }
//    }
//    vector<int> keep_cols(cols_2_kept.begin(), cols_2_kept.end());
//    std::sort(keep_cols.begin(), keep_cols.end());
//    int keep = 0;
//    for (int i = keep; i < keep_cols.size(); ++i) {
//      vec_ColSetUsedInRootMIP[keep++] = vec_ColSetUsedInRootMIP[keep_cols[i]];
//    }
//    vec_ColSetUsedInRootMIP.resize(keep);
//    cout << "ColSetUsedInRootMIP  Size: " << vec_ColSetUsedInRootMIP.size() << endl;
//    numcol = keep;
//  }
//  numrow = Dim;
//  SOLVER solver{};
//  solver.SOLVERgetenv(&Solver);
//  auto xtype = new char[numcol];
//  auto sense = new char[numrow];
//  auto rhs = new double[numrow];
//  fill_n(xtype, numcol, SOLVER_BINARY);
//  fill_n(sense, numrow, SOLVER_EQUAL);
//  sense[RealDim] = SOLVER_GREATER_EQUAL;
//  fill_n(rhs, numrow, 1);
//  rhs[RealDim] = K;
//  safe_solver(solver.SOLVERnewmodel("RootMIP", 0, nullptr, nullptr, nullptr, nullptr, nullptr))
//  safe_solver(solver.SOLVERaddconstrs(numrow, 0, nullptr, nullptr, nullptr, sense, rhs, nullptr))
//  safe_solver(solver.SOLVERsetenvCutoff(UB + roundUpTolerance))
//  //add cols in MIP
//  ccnt = 0;
//  size_t nzcnt = 0;
//  for (auto &col : vec_ColSetUsedInRootMIP) {
//    auto &vec = col.first;
//    solver_beg[ccnt] = nzcnt;
//    solver_obj[ccnt++] = col.second;
//    for (auto i : vec) {
//      solver_ind[nzcnt] = i - 1;
//      solver_val[nzcnt++] = 1;
//    }
//    solver_ind[nzcnt] = RealDim;
//    solver_val[nzcnt++] = 1;
//  }
//  solver_beg[ccnt] = nzcnt;
//  safe_solver(solver.SOLVERXaddvars(ccnt,
//                                    nzcnt,
//                                    solver_beg,
//                                    solver_ind,
//                                    solver_val,
//                                    solver_obj,
//                                    nullptr,
//                                    nullptr,
//                                    xtype,
//                                    nullptr))
//  safe_solver(solver.SOLVERsetenvTimeLimit(100))
//  auto beg = high_resolution_clock::now();
//  safe_solver(solver.SOLVERoptimize())
//  auto end = high_resolution_clock::now();
//  auto eps = duration_cast<milliseconds>(end - beg);
//  cout << "MIP solved in " << (double) eps.count() * 1e-3 << " s!" << endl;
//  //    writeMIP(solver);
//  int status;
//  safe_solver(solver.SOLVERgetStatus(&status))
//  double obj;
//  if (status == SOLVER_INFEASIBLE || status == SOLVER_INF_OR_UNBD) {
//    cout << "Model is infeasible after pre_solving!" << endl;
//    cout << "status: " << status << endl;
//    obj = 1e100;
//  } else {
//    safe_solver(solver.SOLVERgetObjVal(&obj))
//  }
//  bool if_updated = false;
//  if (UB > ceil_transformed_number_related(obj - TOLERANCE) + TOLERANCE) {
//    if_updated = true;
//    safe_solver(solver.SOLVERgetX(0, numcol, X))
//    UB = ceil_transformed_number_related(obj - TOLERANCE);
//    cout << "solve MIP get " << UB << endl;
//    IPOptSol.clear();
//    for (int i = 0; i < numcol; ++i) {
//      if (abs(X[i] - 1) < MIP_TOLERANCE) {
//        vector<int> tmp;
//        tmp.reserve(MaxLengthEleRoute);
//        tmp.emplace_back(0);
//        for (auto vertex : vec_ColSetUsedInRootMIP[i].first) {
//          tmp.emplace_back(vertex);
//        }
//        tmp.emplace_back(0);
//        IPOptSol.emplace_back(std::move(tmp));
//      }
//    }
//
//#ifdef MASTER_VALVE_ML
//    updateUB_EdgeSol();
//#endif
//
//  }
//  delete[]xtype;
//  delete[]sense;
//  delete[]rhs;
//  //recover environment
//  safe_solver(solver.SOLVERsetenvCutoff(MAXCUTOFF))
//  safe_solver(solver.SOLVERsetenvTimeLimit(MAXTIMELIMT4MIP))
//  //free the memory
//  safe_solver(solver.SOLVERfreemodel())
//  ColSetUsedInRootMIP.clear();
//  if (if_updated) {
//    cout << "scheme 1 update UB= " << UB << endl;
//    cout << "resolve LP by CG!" << endl;
//    solveLPInLabeling(node);
//    eliminateArcs(node);
//    enumerateMIP(node);
//  }
//}

void CVRP::copyColGeneratedInRoot4MIP() {
  int curr_node, past_node = 0;
  double cost_sum = 0;
  yzzLong tmp_pi = 0;
  auto seq = new int[Dim];
  int size_seq = 0;
  for (size_t start = PriorPoolBeg4Pricing + 1; start < PoolBeg4Pricing; ++start) {
    curr_node = ColPool4Pricing[start];
    cost_sum += CostMat4Vertex[past_node][curr_node];
    if (!curr_node) {
      //add col into map
      if (ColSetUsedInRootMIP.find(tmp_pi) == ColSetUsedInRootMIP.end()) {
        ColSetUsedInRootMIP[tmp_pi] = {vector<int>(seq, seq + size_seq), cost_sum};
      } else {
        if (ColSetUsedInRootMIP[tmp_pi].second > cost_sum + TOLERANCE) {
          ColSetUsedInRootMIP[tmp_pi] = {vector<int>(seq, seq + size_seq), cost_sum};
        }
      }
      cost_sum = 0;
      past_node = 0;
      tmp_pi = 0;
      size_seq = 0;
      ++start;
      continue;
    }
    if (tmp_pi[curr_node]) {
      //give up this col
      for (; start < PoolBeg4Pricing; ++start) {
        if (!ColPool4Pricing[start]) {
          cost_sum = 0;
          past_node = 0;
          tmp_pi = 0;
          size_seq = 0;
          ++start;
          break;
        }
      }
      continue;
    } else {
      tmp_pi.set(curr_node);
    }
    seq[size_seq++] = curr_node;
    past_node = curr_node;
  }
  delete[]seq;
}

void CVRP::collectDualInRoot4MIP(BBNODE *const node) {
  auto dual = new double[NumRow];
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetDual(0, NumRow, dual))
  DualInRootMIP.resize(NumRow);
  for (int i = 0; i < NumRow; ++i) {
    DualInRootMIP[i] = dual[i];
  }
  delete[]dual;
}

#ifdef DEBUG_MIP_FILE
void CVRP::writeMIP(SOLVER &solver, const string &file_type) const {
  cout << "WARNING: we generate IP file here!" << endl;
  string name;
  if (file_type == "LP") {
    name =
        "LPs/LP_CutOff_" + to_string(UB + roundUpTolerance) + "_col_" + to_string(NumCol)
            + ".rlp";
  } else if (file_type == "IP") {
    name =
        "IPs/CutOff_" + to_string(UB + roundUpTolerance) + "_col_" + to_string(NumCol)
            + ".rlp";
  } else if (file_type == "MPS") {
    name =
        "MPSs/CutOff_" + to_string(UB + roundUpTolerance) + "_col_" + to_string(NumCol)
            + ".mps";
  }
  if (NumCol > 5000) safe_solver(solver.SOLVERwrite(name.c_str()))
  cout << "write file " << name << endl;
}

#endif

#ifdef DEBUG_LP_FILE
void CVRP::writeLP_N_tell_if_LP_corrected(SOLVER &solver) const {
  string name = "tmp.rlp";
  cout << "write file " << name << endl;
  safe_solver(solver.SOLVERwrite(name.c_str()))
  cout << "write file " << name << " done!" << endl;
  ///
  cout << "check if LP is corrected" << endl;
  cout << "we invoke the function in tell_if_lp_file_contains_same_variables.py" << endl;
  //test in python
  Py_Initialize();
  // add .py path
  PyRun_SimpleString("import os");
  PyRun_SimpleString("import re");
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('./')");
  PyObject *pModule = PyImport_ImportModule("tell_if_lp_file_contains_same_variables");
  // 导入要运行的函数
  PyObject *pFunc = PyObject_GetAttrString(pModule, "search4log");
  // 构造传入参数
  PyObject *args = PyTuple_New(1);
  PyTuple_SetItem(args, 0, Py_BuildValue("s", name.c_str()));
  // 运行函数，并获取返回值
  PyObject *pRet = PyObject_CallObject(pFunc, args);
  bool result = false;
  if (pRet) {
    result = PyLong_AsLong(pRet);
    cout << "if_wrong:" << result << endl;
  }
  Py_Finalize();
  if (result) {
    cout << "WARNING: the LP file is not corrected!" << endl;
    exit(0);
  }
}
#endif

void CVRP::specialize_MaxLengthEleRoute() {
  //need to be irrelevant to the number of vehicles
  //MaxLengthRoute
  vector<double> demand(RealDim);
  for (int i = 1; i < Dim; ++i) {
    demand[i - 1] = InfoVertex[i][3];
  }
  std::stable_sort(demand.begin(), demand.end());
  double acc_cap = 0;
  MaxLengthEleRoute = RealDim + 3;
  for (int i = 0; i < RealDim; ++i) {
    acc_cap += demand[i];
    if (acc_cap > Cap - TOLERANCE) {
      MaxLengthEleRoute = i + 3;
      break;
    }
  }
}

//void CVRP::writeAllColsInIP4Check(BBNODE *node) const {
//  //write all cols in IP
//  SOLVER solver{};
//  solver.SOLVERgetsolver(&node->solver);
//  vector<int> Col_added(node->SizeEnuColPool);
//  iota(Col_added.begin(), Col_added.end(), 0);
//  auto &ptr = node->IdxColsInEnuColPool;
//  auto &ptr_cost = node->Cost4ColsInEnuColPool;
//  auto &mat = node->MatInEnu;
//  size_t nzcnt = 0;
//  int ccnt = 0;
//  int index = NumCol;
//
//  for (auto col : Col_added) {
//    solver_obj[ccnt] = ptr_cost[col];
//    solver_beg[ccnt++] = nzcnt;
//    for (int row = 0; row < NumRow; ++row) {
//      if (mat(row, col) != 0) {
//        solver_ind[nzcnt] = row;
//        solver_val[nzcnt++] = mat(row, col);
//      }
//    }
//  }
//
//  safe_Hyperparameter(ccnt >= CONFIG::MaxNumCols)
//  safe_Hyperparameter(nzcnt >= size_t(CONFIG::MaxNumCols) * MaxLengthEleRoute)
//
//  safe_solver(solver.SOLVERXaddvars(ccnt,
//                                    nzcnt,
//                                    solver_beg,
//                                    solver_ind,
//                                    solver_val,
//                                    solver_obj,
//                                    nullptr,
//                                    nullptr,
//                                    nullptr,
//                                    nullptr))
//  safe_solver(solver.SOLVERupdatemodel())
//  int numcol;
//  safe_solver(solver.SOLVERgetNumCol(&numcol))
//  vector<char> vtype(numcol, SOLVER_BINARY);
//  safe_solver(solver.SOLVERsetVTypearray(0, numcol, vtype.data()))
//  string name = "MIP_Rows_" + to_string(NumRow) + "_Cols_" + to_string(numcol) + ".mps";
//  cout << "write file " << name << endl;
//  safe_solver(solver.SOLVERwrite(name.c_str()))
//  cout << "write file " << name << " done!" << endl;
//  //free
//  safe_solver(solver.SOLVERfreemodel())
//}

//void CVRP::addRCC2LP(BBNODE *const node, unordered_map<int, vector<int>> &vertex_map, sparseRowMatrixXd &mat) const {
//  for (auto &rcc : node->RCCs) {
//    if (rcc.FormRCC) {
//      int idx = rcc.IdxRCC;
//      auto &info = rcc.InfoRCCCustomer;
//      for (auto it = info.begin(); it != info.end(); ++it) {
//        int ai = *it;
//        auto it_inner = it;
//        ++it_inner;
//        for (; it_inner != info.end(); ++it_inner) {
//          int aj = *it_inner;
//          for (auto it_map : vertex_map[ai * Dim + aj]) ++mat(idx, it_map);
//          for (auto it_map : vertex_map[aj * Dim + ai]) ++mat(idx, it_map);
//        }
//      }
//    } else {
//      int idx = rcc.IdxRCC;
//      auto &customer_info = rcc.InfoRCCCustomer;
//      auto &outside_customer_info = rcc.InfoRCCOutsideCustomer;
//      for (auto it = outside_customer_info.begin(); it != outside_customer_info.end(); ++it) {
//        int ai = *it;
//        auto it_inner = it;
//        ++it_inner;
//        for (; it_inner != outside_customer_info.end(); ++it_inner) {
//          int aj = *it_inner;
//          for (auto it_map : vertex_map[ai * Dim + aj]) ++mat(idx, it_map);
//          for (auto it_map : vertex_map[aj * Dim + ai]) ++mat(idx, it_map);
//        }
//      }
//      for (int aj : outside_customer_info) {
//        for (auto it_map : vertex_map[aj]) mat(idx, it_map) += 0.5;
//        for (auto it_map : vertex_map[aj * Dim]) mat(idx, it_map) += 0.5;
//      }
//      for (int aj : customer_info) {
//        for (auto it_map : vertex_map[aj]) mat(idx, it_map) -= 0.5;
//        for (auto it_map : vertex_map[aj * Dim]) mat(idx, it_map) -= 0.5;
//      }
//    }
//  }
//}

bool CVRP::increaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end) {
  newMainResource = nowMainResource + MainResourceAcrossArcsInForwardSense[start][end];
#ifdef CAPACITY_AS_MAIN_RESOURCE
#else
  if (newMainResource > ub4Vertex[end]) return false;
  newMainResource = newMainResource > lb4Vertex[end] ? newMainResource : lb4Vertex[end];
#endif
  if (newMainResource + MainResourceAcrossArcsInForwardSense[end][0] > MaxMainResource) return false;
  return true;
}

bool CVRP::decreaseMainResourceConsumption(double nowMainResource, double &newMainResource, int start, int end) {
  newMainResource = nowMainResource - MainResourceAcrossArcsInBackwardSense[start][end];
#ifdef CAPACITY_AS_MAIN_RESOURCE
#else
  if (newMainResource < lb4Vertex[end]) return false;
  newMainResource = newMainResource < ub4Vertex[end] ? newMainResource : ub4Vertex[end];
#endif
  if (newMainResource - MainResourceAcrossArcsInBackwardSense[end][0] < 0) return false;
  return true;
}

bool CVRP::runColumnGenerationType(BBNODE *node, int mode) {
  if (node->if_terminated) return true;
  int status = 0;
  bool switch_is_on = true;
  switch (mode) {
    case 1:cout << "LighterHeur phase begin...\n";
      break;
    case 2:cout << "HeavierHeur phase begin...\n";
      break;
    case 3:cout << "Exact phase begin...\n";
      break;
    default:cerr << "None of these modes are used!\n";
      exit(-1);
  }
  int iter = 0, ccnt = 0, old_ncol = NumCol, tag;
  double eps_CG = 0, eps_LP = 0, b4_node_val = node->Val;
  auto beg = chrono::high_resolution_clock::now();
  auto end = beg;
  bool if_tail_off = false;
  bool old_force;
  bool if_gap_0 = false;
//  GapBetweenLastSmallestRCAndRCThreshold = 0;
  if (mode == 3) {
    if_exact_labeling_CG = true;
    old_force = ForceNotRollback;
    if_ban_convertDual = false;
  }

  while (true) {

    ++iter;

    beg = chrono::high_resolution_clock::now();

    switch (mode) {
      case 1:ccnt = generateColsByLighterHeur(node);
        break;
      case 2:ccnt = generateColsByHeavierHeur(node);
        break;
      case 3:
        if (abs(GapBetweenLastSmallestRCAndRCThreshold) < TOLERANCE) {
          if_gap_0 = true;
        }
#ifdef DETAILED_EXACT_PRINT_INFO
        cout << "---------------------------" << endl;
        cout << "gap= " << GapBetweenLastSmallestRCAndRCThreshold << endl;
#endif
        ccnt = generateColsByBidir<
#ifdef SYMMETRY_PROHIBIT
            false
#else
            true
#endif
        >(node);
        if (if_gap_0) {
          ForceNotRollback = true;
        }
        if (Rollback == 3) {
          if_tail_off = true;
        }

        break;
      default:cerr << "None of these modes are used!\n";
        exit(-1);
    }

    end = chrono::high_resolution_clock::now();
    TimePricing4Iter = chrono::duration<double>(end - beg).count();
    eps_CG += TimePricing4Iter;

    if (!ccnt) {
      if (mode != 3 || Rollback == 1 || if_gap_0) {
        --iter;
        break;
      }
    }
    if (node->if_terminated) {
      --iter;
      break;
    }

    beg = chrono::high_resolution_clock::now();

    //tag mark if the CG is converged
    tag = optimizeLP4OneIter(node, b4_node_val);
    b4_node_val = LPVal;

    end = chrono::high_resolution_clock::now();
    TimeResolveLP4Iter = chrono::duration<double>(end - beg).count();
    eps_LP += TimeResolveLP4Iter;

    if (tag) {
      status = 1;//QUIT
      goto QUIT;
    }

    if (!(iter % PRINT_LABELING_STEP_SIZE)) {
      REPORT:
      GloEnd = chrono::high_resolution_clock::now();
      GloEps = chrono::duration<double>(GloEnd - GloBeg).count();
      printInfoLabeling(iter, NumCol - old_ncol, NumCol, NumRow, eps_LP,
                        eps_CG, GloEps,
                        LPVal, LB, UB);
      if (!switch_is_on) {
        goto QUIT;
      }
      eps_CG = 0;
      eps_LP = 0;
      old_ncol = NumCol;
    }
  }
  if (!iter || (switch_is_on && iter % PRINT_LABELING_STEP_SIZE)) {
    switch_is_on = false;
    goto REPORT;
  }
  QUIT:
  if (mode == 3) {
    if_exact_labeling_CG = false;
    ForceNotRollback = old_force;
    if (if_tail_off && !Rollback)
      Rollback = 3;
  }
  return status;
}

double CVRP::transformCost(double x) {
  return std::floor(x + 0.5);
}

void CVRP::setTailOffStd_N_RollBackStd() const {
  for (int b = 0; b < NumBucketsPerVertex; ++b) {
    int min_num_labels = MaxInt;
    for (int i = 0; i < Dim; ++i) {
      if (LabelArrayInForwardSense[i][b].second && min_num_labels > LabelArrayInForwardSense[i][b].second) {
        min_num_labels = LabelArrayInForwardSense[i][b].second;
      }
    }
    if (min_num_labels == MaxInt) {
      min_num_labels = 1;
    }
    for (int i = 0; i < Dim; ++i) {
      int hard_max = max(LabelArrayInForwardSense[i][b].second, min_num_labels) * FACTOR_NUM_LABEL;
      //take the power of 2
      hard_max = int(pow_self(2.0, int(ceil(log(hard_max) / log(2)) + TOLERANCE)));
      LabelArrayInForwardSense[i][b].first.resize(hard_max);
      IfExistExtraLabelsInForwardSense[i][b].first.resize(hard_max);
    }
#ifdef SYMMETRY_PROHIBIT
    min_num_labels = MaxInt;
    for (int i = 0; i < Dim; ++i) {
      if (LabelArrayInBackwardSense[i][b].second && min_num_labels > LabelArrayInBackwardSense[i][b].second) {
        min_num_labels = LabelArrayInBackwardSense[i][b].second;
      }
    }
    if (min_num_labels == MaxInt) {
      min_num_labels = 1;
    }
    for (int i = 0; i < Dim; ++i) {
      int hard_max = max(LabelArrayInBackwardSense[i][b].second, min_num_labels) * FACTOR_NUM_LABEL;
      //take the power of 2
      hard_max = (int) (pow_self(2.0, (int) ceil(log(hard_max) / log(2))) + TOLERANCE);
      LabelArrayInBackwardSense[i][b].first.resize(hard_max);
      IfExistExtraLabelsInBackwardSense[i][b].first.resize(hard_max);
    }
#endif
  }
}

void CVRP::deleteNewAddedNonActiveCutsBySlack(BBNODE *node, int oldNum, bool if_keep_rcc) {
  int cnt = 0;

  safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack))
  int *cstr_index = new int[NumRow];
  for (int i = 0; i < NumRow; ++i) cstr_index[i] = i;
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);
  unordered_set<int> new_rcc_index;
  if (if_keep_rcc) {
    new_rcc_index.reserve(NumRow);
    //check rcc index reverse sequence
    for (auto i = node->RCCs.rbegin(); i < node->RCCs.rend(); ++i) {
      if (i->IdxRCC < oldNum) {
        break;
      }
      new_rcc_index.emplace(i->IdxRCC);
    }
  }
  for (int i = oldNum; i < NumRow; ++i) {
    if (abs(Slack[i]) > TOLERANCE) {
      if (new_rcc_index.find(i) != new_rcc_index.end()) continue;
      solver_ind[cnt++] = i;
      cstr_index[i] = -1;
      deleted_cstrs.emplace_back(i);
    }
  }

  if (deleted_cstrs.empty()) {
    delete[]cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))

  //rcc
  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
//      cout << "delete rcc" << endl;
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //r1c
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //r1c_multi
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //BrC finally, BrC does not need to delete, only need to reassign
  for (auto &i : node->BrCs) {
    i.IdxBrC = cstr_index[i.IdxBrC];
  }

  convertVertex2R1CsInOneLP(node);

  delete[]cstr_index;
}

void CVRP::deleteNonActiveCutsBySlack(BBNODE *node, bool if_rcc, bool if_r1c) {
  int cnt = 0;
  if (!if_rcc && !if_r1c) return;
  safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack))
  int *cstr_index = new int[NumRow];
  iota(cstr_index, cstr_index + NumRow, 0);
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);

  if (if_rcc) {
    for (auto &i : node->RCCs) {
      if (abs(Slack[i.IdxRCC]) > TOLERANCE) {
        solver_ind[cnt++] = i.IdxRCC;
        cstr_index[i.IdxRCC] = -1;
        deleted_cstrs.emplace_back(i.IdxRCC);
      }
    }
  }

  if (if_r1c) {
    for (auto &i : node->R1Cs) {
      if (abs(Slack[i.IdxR1C]) > TOLERANCE) {
        solver_ind[cnt++] = i.IdxR1C;
        cstr_index[i.IdxR1C] = -1;
        deleted_cstrs.emplace_back(i.IdxR1C);
      }
    }
    //check R1Cs_multi
    for (auto &i : node->R1Cs_multi) {
      if (abs(Slack[i.IdxR1C]) > TOLERANCE) {
        solver_ind[cnt++] = i.IdxR1C;
        cstr_index[i.IdxR1C] = -1;
        deleted_cstrs.emplace_back(i.IdxR1C);
      }
    }
  }

  if (deleted_cstrs.empty()) {
    delete[]cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))

  //start the delete operation and the assignment operation synchronously,
  //but BrC does not need it, BrC only needs to perform the reassignment operation.
  //RCCs are the first

  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //R1C3s next
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //R1C3s_multi next
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //BrC finally, BrC does not need to delete, only need to reassign
  for (auto &i : node->BrCs) i.IdxBrC = cstr_index[i.IdxBrC];

  convertVertex2R1CsInOneLP(node);
  delete[]cstr_index;
}

void CVRP::deleteNonActiveCutsByDual(BBNODE *node, bool if_rcc_by_slack) {
  int cnt = 0;
  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
  int *cstr_index = new int[NumRow];
  iota(cstr_index, cstr_index + NumRow, 0);
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);

  if (if_rcc_by_slack) {
    safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack))
    for (auto &i : node->RCCs) {
      if (abs(Slack[i.IdxRCC]) > TOLERANCE) {
        solver_ind[cnt++] = i.IdxRCC;
        cstr_index[i.IdxRCC] = -1;
        deleted_cstrs.emplace_back(i.IdxRCC);
      }
    }
  } else {
    for (auto &i : node->RCCs) {
      if (abs(Pi[i.IdxRCC]) < DUAL_TOLERANCE) {
        solver_ind[cnt++] = i.IdxRCC;
        cstr_index[i.IdxRCC] = -1;
        deleted_cstrs.emplace_back(i.IdxRCC);
      }
    }
  }
  //check R1Cs
  for (auto &i : node->R1Cs) {
    if (abs(Pi[i.IdxR1C]) < DUAL_TOLERANCE) {
      solver_ind[cnt++] = i.IdxR1C;
      cstr_index[i.IdxR1C] = -1;
      deleted_cstrs.emplace_back(i.IdxR1C);
    }
  }
//check r1c_multi
  for (auto &i : node->R1Cs_multi) {
    if (abs(Pi[i.IdxR1C]) < DUAL_TOLERANCE) {
      solver_ind[cnt++] = i.IdxR1C;
      cstr_index[i.IdxR1C] = -1;
      deleted_cstrs.emplace_back(i.IdxR1C);
    }
  }

  if (deleted_cstrs.empty()) {
    delete[]cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))

  //start the delete operation and the assignment operation synchronously,
  //but BrC does not need it, BrC only needs to perform the reassignment operation.
  //RCCs are the first
  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //R1C3s next
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //R1C_multi next
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //BrC finally, BrC does not need to delete, only need to reassign
  for (auto &i : node->BrCs) i.IdxBrC = cstr_index[i.IdxBrC];

  convertVertex2R1CsInOneLP(node);
  delete[]cstr_index;
}

//void CVRP::deleteNonActiveCutsByDual(BBNODE *node) {
//  int cnt = 0;
//  vector<int> CBasis(NumRow);
//  cout << "关于cbasis的测试！" << endl;
//  safe_solver(GRBgetintattrarray(node->solver.model, "CBasis", 0, NumRow, CBasis.data()))
//  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//  int *cstr_index = new int[NumRow];
//  iota(cstr_index, cstr_index + NumRow, 0);
//  vector<int> deleted_cstrs;
//  deleted_cstrs.reserve(NumRow);
//
//  for (auto &i : node->RCCs) {
//    if (abs(Pi[i.IdxRCC]) < TOLERANCE && CBasis[i.IdxRCC]) {
//      solver_ind[cnt++] = i.IdxRCC;
//      cstr_index[i.IdxRCC] = -1;
//      deleted_cstrs.emplace_back(i.IdxRCC);
//    }
//  }
//  //check R1Cs
//  for (auto &i : node->R1Cs) {
//    if (abs(Pi[i.IdxR1C]) < TOLERANCE && CBasis[i.IdxR1C]) {
//      solver_ind[cnt++] = i.IdxR1C;
//      cstr_index[i.IdxR1C] = -1;
//      deleted_cstrs.emplace_back(i.IdxR1C);
//    }
//  }
////check r1c_multi
//  for (auto &i : node->R1Cs_multi) {
//    if (abs(Pi[i.IdxR1C]) < TOLERANCE && CBasis[i.IdxR1C]) {
//      solver_ind[cnt++] = i.IdxR1C;
//      cstr_index[i.IdxR1C] = -1;
//      deleted_cstrs.emplace_back(i.IdxR1C);
//    }
//  }
//
//  if (deleted_cstrs.empty()) {
//    delete[]cstr_index;
//    return;
//  }
//  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
//  int delta = 0;
//  auto stop_sign = deleted_cstrs.end() - 1;
//  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
//    ++delta;
//    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
//  }
//  ++delta;
//  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;
//
//  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
//  safe_solver(node->solver.SOLVERreoptimize())
//  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
//
//  //start the delete operation and the assignment operation synchronously,
//  //but BrC does not need it, BrC only needs to perform the reassignment operation.
//  //RCCs are the first
//  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
//    if (cstr_index[i->IdxRCC] == -1) {
//      i = node->RCCs.erase(i);
//    } else {
//      i->IdxRCC = cstr_index[i->IdxRCC];
//      ++i;
//    }
//  }
//
//  //R1C3s next
//  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
//    if (cstr_index[i->IdxR1C] == -1) {
//      i = node->R1Cs.erase(i);
//    } else {
//      i->IdxR1C = cstr_index[i->IdxR1C];
//      ++i;
//    }
//  }
//
//  //R1C_multi next
//  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
//    if (cstr_index[i->IdxR1C] == -1) {
//      i = node->R1Cs_multi.erase(i);
//    } else {
//      i->IdxR1C = cstr_index[i->IdxR1C];
//      ++i;
//    }
//  }
//
//  //BrC finally, BrC does not need to delete, only need to reassign
//  for (auto &i : node->BrCs) i.IdxBrC = cstr_index[i.IdxBrC];
//
//  convertVertex2R1CsInOneLP(node);
//  delete[]cstr_index;
//}

void CVRP::getInfoEdge(BBNODE *node, bool if_br) const {
  //Need LPOpt information to activate
  if (node->if_Int) return;
  node->NumEdges = 0;
  double **arc_graph;
  if (if_br) arc_graph = ArcGraph_revised;
  else arc_graph = ArcGraph;
  if (node->SizeAllocatedMem < NumEdge) {
    node->allocateMem(MaxNumEdge);
    node->SizeAllocatedMem = MaxNumEdge;
  }

  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      if (arc_graph[i][j] > TOLERANCE) {
        ++node->NumEdges;
        node->EdgeTail[node->NumEdges] = i;
        node->EdgeHead[node->NumEdges] = j;
        node->EdgeVal[node->NumEdges] = arc_graph[i][j];
      }
    }
  }
#ifdef debugMLInputData
  if (if_br) {
    cout << "check MLInputData by solution!" << endl;
    vector<double> x(NumCol);
    safe_solver(node->solver.SOLVERgetX(0, NumCol, x.data()))
    vector<vector<double>> ArcGraph2(Dim, vector<double>(Dim, 0));
    auto colpool = If_in_Enu_State ? ColPool4Pricing : ColPool4Mem;
    for (int i = 0; i < NumCol; ++i) {
      if (x[i] > TOLERANCE) {
        int past = 0, cur;
        for (auto j = node->IdxCols[i] + 1;; ++j) {
          cur = colpool[j];
          ArcGraph2[past][cur] += x[i];
          ArcGraph2[cur][past] += x[i];
          past = cur;
          if (cur == 0) break;
        }
      }
    }
    for (int i = 0; i < Dim; ++i) {
      for (int j = i + 1; j < Dim; ++j) {
        if (abs(ArcGraph[i][j] - ArcGraph2[i][j]) > TOLERANCE && ArcGraph2[i][j] < 1) {
          cerr << "error!" << endl;
          cout << "error! in debugMLInputData" << endl;
          cout << ArcGraph[i][j] << " " << ArcGraph2[i][j] << endl;
        }
      }
    }
  }
#endif
}

#ifdef DEBUG
void CVRP::checksol_(BBNODE *node) const {
  vector<vector<int>> sol_ =
      {
          {0, 28, 26, 12, 68, 80, 24, 29, 78, 34, 35, 65, 71, 66, 20, 51, 9, 81, 33, 79, 3, 77, 76, 50, 1, 69, 27, 0},
          {0, 13, 87, 42, 43, 15, 57, 2, 41, 22, 74, 75, 56, 23, 67, 39, 25, 55, 54, 4, 72, 73, 21, 40, 58, 53, 0},
          {0, 94, 95, 97, 92, 98, 37, 100, 91, 44, 14, 38, 86, 16, 61, 85, 93, 59, 99, 96, 6, 89, 0},
          {0, 31, 88, 7, 82, 48, 19, 11, 62, 10, 70, 30, 32, 90, 63, 64, 49, 36, 47, 46, 8, 45, 17, 84, 5, 60, 83, 18,
           52, 0}
      };
  //first, branch cannot be contradicted with right solution
  set<int> edge;
  for (auto &sol : sol_) {
    int old_seq = 0;
    for (int i = 1; i < sol.size() - 1; ++i) {
      int new_seq = sol[i];
      if (new_seq == 0) break;
      if (old_seq < new_seq) edge.insert(old_seq * Dim + new_seq);
      else edge.insert(new_seq * Dim + old_seq);
      old_seq = new_seq;
    }
  }
  for (auto &brc : node->BrCs) {
    int ai = brc.Edge.first, aj = brc.Edge.second;
    int tmp = ai * Dim + aj;
    if (brc.BrDir) {
      if (edge.find(tmp) == edge.end()) {
        cout << "must use branch " << ai << " " << aj << " is wrong" << endl;
        cout << "we dont check this !" << endl;
        return;
      }
    } else {
      if (edge.find(tmp) != edge.end()) {
        cout << "dont use branch " << ai << " " << aj << " is wrong" << endl;
        cout << "we dont check this !" << endl;
        return;
      }
    }
  }
  bool if_quit = false;
  for (auto &sol : sol_) {
    double q = 0;
    int old_seq = 0;
    for (int i = 0; i < sol.size() - 1; ++i) {
      int new_seq = sol[i];
      //check if there is an arc for old_seq to new_seq
      int bin = int(q / StepSize);
      if (old_seq && new_seq) {
        auto if_find = std::find(node->AllForwardBuckets[old_seq][bin].BucketArcs.begin(),
                                 node->AllForwardBuckets[old_seq][bin].BucketArcs.end(),
                                 new_seq);
        if (if_find == node->AllForwardBuckets[old_seq][bin].BucketArcs.end()) {
          cout << "error" << endl;
          cout << "missing " << old_seq << "->" << new_seq << endl;
          if_quit = true;
        }
      }
      q += MainResourceAcrossArcsInForwardSense[old_seq][new_seq];
      old_seq = new_seq;
    }
  }
  if (if_quit) exit(0);
}
#endif

#ifdef writeIP
void CVRP::writeRCFIP(double ub, const SOLVER &solver) const {
  cout << "write lp ...." << endl;
  self_mkdir(writeIP);
  int seq = 0;
  string filename;
  do {
    filename = string(writeIP) + "/ins_" + FileName + "_dim_" + to_string(Dim) + "_ub_" + to_string(ub) + "_seq_"
        + to_string(seq++) + ".mps";
  } while (std::filesystem::exists(filename));
  safe_solver(solver.SOLVERwrite(filename.c_str()))
  cout << "write lp done!" << endl;
}
#endif

#ifdef DELUXING_APPLIED
void CVRP::applyRCF(BBNODE *node, int round, bool if_verbose) {
  cout << "WARNING: RCF can only be applied by GRB model!" << endl;
  SOLVER solver{};
  solver.SOLVERgetsolver(&node->solver);
  int numcol = NumCol;
  int numrow = NumRow;
  size_t numnz = 0;
  int ccnt = 0;
  auto &mat = node->MatInEnu;
  auto &Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
  auto &cost = node->Cost4ColsInEnuColPool;
  size_t num = mat.count();
  if (checkSolverPtr(num)) {
    reallocateSolverPtr(num);
  }
  vector<int> map_new_old;
  map_new_old.reserve(node->SizeEnuColPool);
  for (int i = 0; i < node->SizeEnuColPool; ++i) {
    if (Deleted_ColsInEnuPool[i]) continue;
    solver_beg[ccnt] = numnz;
    solver_obj[ccnt++] = cost[i];
    for (int j = 0; j < numrow; ++j) {
      if (mat(j, i) != 0) {
        solver_ind[numnz] = j;
        solver_val[numnz++] = mat(j, i);
      }
    }
    map_new_old.emplace_back(i);
  }
  safe_solver(solver.SOLVERXaddvars(ccnt,
                                    numnz,
                                    solver_beg,
                                    solver_ind,
                                    solver_val,
                                    solver_obj,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr))
  safe_solver(solver.SOLVERupdatemodel())
  vector<int> idx;
  idx.reserve(node->SizeEnuColPool);
  cout << "RCF applying..." << endl;
//  string name = to_string(UB + roundUpTolerance) + ".rcf.mps";
//  solver.SOLVERwrite(name.c_str());
//  exit(0);
#ifdef writeIP
  writeRCFIP(UB + roundUpTolerance, solver);
#endif
  auto beg = high_resolution_clock::now();
  int beta1 = 1000, beta2 = 100000;
  const double timelimit = 100000;
  deLuxing(solver.model,
           UB + roundUpTolerance,
           RealDim,
           round,
           beta1,
           beta2,
           idx,
           timelimit,
           MIP_TOLERANCE,
           if_verbose);
  cout << "beta1= " << beta1 << " beta2= " << beta2 << endl;
  auto end = high_resolution_clock::now();
  cout << "RCF time: " << duration<double>(end - beg).count() << " s" << endl;
  cout << "RCF over!" << endl;
  int new_col;
  safe_solver(solver.SOLVERgetNumCol(&new_col))
  if (!new_col) {
    idx.resize(numcol + ccnt);
    iota(idx.begin(), idx.end(), 0);
  }
  safe_solver(solver.SOLVERfreemodel())

  auto if_del = new bool[numcol]();
  int cnt = 0;
  for (auto i : idx) {
    if (i >= numcol) break;
    if_del[i] = true;
    ++cnt;
  }

  //delete columns in lp but always keep the first column!
  int len = 0, keep = 1;
  for (int i = keep; i < numcol; ++i) {
    if (if_del[i]) {
      solver_ind[len++] = i;
    } else {
      node->IdxCols[keep++] = node->IdxCols[i];
    }
  }
  delete[]if_del;
  safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  //delete columns in enumeration pool
  for (int i = cnt; i < idx.size(); ++i) {
    node->Deleted_ColsInEnuPool[map_new_old[idx[i] - numcol]] = true;
  }
  cleanColsInPool(node, nullptr);
}
#endif

void CVRP::regenerateBucketGraph(BBNODE *node) {
  if (abs(StepSize / 2) < 1 - TOLERANCE) return;
  if_stopArcElimination = true;
  StepSize /= 2;
  NumBucketsPerVertex *= 2;
  MaxNumForwardGraphArc *= 2;
  node->NumForwardBucketArcs *= 2;
  node->NumForwardJumpArcs *= 2;
#ifdef SYMMETRY_PROHIBIT
  node->NumBackwardBucketArcs *= 2;
  node->NumBackwardJumpArcs *= 2;
  MaxNumBackwardGraphArc *= 2;
#endif
  //read old bucket graph
  auto new_LabelArrayInForwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    new_LabelArrayInForwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      //keep the half space!
      new_LabelArrayInForwardSense[i][j].first.resize(LabelArrayInForwardSense[i][j / 2].first.size() / 2);
    }
  }

  auto new_IfExistExtraLabelsInForwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    new_IfExistExtraLabelsInForwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      new_IfExistExtraLabelsInForwardSense[i][j].first.resize(
          LabelArrayInForwardSense[i][j / 2].first.size() / 2);
    }
  }

  for (int i = 0; i < Dim; ++i) {
    delete[]RC2TillThisBinInForwardSense[i];
    RC2TillThisBinInForwardSense[i] = new double[NumBucketsPerVertex];
  }

  for (int i = 0; i < Dim; ++i) {
    delete[]RC2BinInForwardSense[i];
    RC2BinInForwardSense[i] = new double[NumBucketsPerVertex];
  }
#ifdef SYMMETRY_PROHIBIT
  auto new_LabelArrayInBackwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    new_LabelArrayInBackwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      new_LabelArrayInBackwardSense[i][j].first.resize(LabelArrayInBackwardSense[i][j / 2].first.size() / 2);
    }
  }

  auto new_IfExistExtraLabelsInBackwardSense = new VecLabel *[Dim];
  for (int i = 0; i < Dim; ++i) {
    new_IfExistExtraLabelsInBackwardSense[i] = new VecLabel[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      new_IfExistExtraLabelsInBackwardSense[i][j].first.resize(LabelArrayInBackwardSense[i][j / 2].first.size() / 2);
    }
  }
  for (int i = 0; i < Dim; ++i) {
    delete[]RC2TillThisBinInBackwardSense[i];
    RC2TillThisBinInBackwardSense[i] = new double[NumBucketsPerVertex];
  }
  for (int i = 0; i < Dim; ++i) {
    delete[]RC2BinInBackwardSense[i];
    RC2BinInBackwardSense[i] = new double[NumBucketsPerVertex];
  }
#endif
  //delete old mem and assign new mem
  for (int i = 0; i < Dim; ++i) {
    delete[]LabelArrayInForwardSense[i];
  }
  delete[]LabelArrayInForwardSense;
  LabelArrayInForwardSense = new_LabelArrayInForwardSense;

  for (int i = 0; i < Dim; ++i) {
    delete[]IfExistExtraLabelsInForwardSense[i];
  }
  delete[]IfExistExtraLabelsInForwardSense;
  IfExistExtraLabelsInForwardSense = new_IfExistExtraLabelsInForwardSense;

#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < Dim; ++i) {
    delete[]LabelArrayInBackwardSense[i];
  }
  delete[]LabelArrayInBackwardSense;
  LabelArrayInBackwardSense = new_LabelArrayInBackwardSense;
  for (int i = 0; i < Dim; ++i) {
    delete[]IfExistExtraLabelsInBackwardSense[i];
  }
  delete[]IfExistExtraLabelsInBackwardSense;
  IfExistExtraLabelsInBackwardSense = new_IfExistExtraLabelsInBackwardSense;
#endif

  //now change the bucket graph in the node!
  auto new_AllForwardBuckets = new Bucket *[Dim];
  for (int i = 0; i < Dim; ++i) {
    new_AllForwardBuckets[i] = new Bucket[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      new_AllForwardBuckets[i][j] = node->AllForwardBuckets[i][j / 2];
    }
  }
#ifdef SYMMETRY_PROHIBIT
  auto new_AllBackwardBuckets = new Bucket *[Dim];
  for (int i = 0; i < Dim; ++i) {
    new_AllBackwardBuckets[i] = new Bucket[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      new_AllBackwardBuckets[i][j] = node->AllBackwardBuckets[i][j / 2];
    }
  }
#endif
  for (int i = 0; i < Dim; ++i) {
    delete[]node->AllForwardBuckets[i];
  }
  delete[]node->AllForwardBuckets;
  node->AllForwardBuckets = new_AllForwardBuckets;
#ifdef SYMMETRY_PROHIBIT
  for (int i = 0; i < Dim; ++i) {
    delete[]node->AllBackwardBuckets[i];
  }
  delete[]node->AllBackwardBuckets;
  node->AllBackwardBuckets = new_AllBackwardBuckets;
#endif
  TellWhichBin4ArcEliminationInForwardSense.clear();
  populateTellWhichBin4ArcElimination<true>();

#ifdef SYMMETRY_PROHIBIT
  TellWhichBin4ArcEliminationInBackwardSense.clear();
  populateTellWhichBin4ArcElimination<false>();
#endif

  cout << "new generated bucket graph: NumBucketsPerVertex= " << NumBucketsPerVertex << " StepSize= " << StepSize
       << endl;
  cout << "we cannot use arc elimination next!" << endl;
}

void CVRP::getLowerBoundofMinimumNumCars() {
  //calculate the minimum number of K by Capacity
  double sum_demand = accumulate(Demand + 1, Demand + Dim, 0.0);
  int cap_k = ceil(sum_demand / Cap);
  K = cap_k;
  cout << "LBNumVehicle= " << K << endl;
}

void CVRP::deleteNewAddedActiveCutsByDual_N_Mem(BBNODE *node, int oldNum) {
  int cnt = 0;

  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
  int *cstr_index = new int[NumRow];
  iota(cstr_index, cstr_index + NumRow, 0);

  vector<pair<int, double>> record;
  record.reserve(NumRow);

//  for (int i = 0; i < node->R1Cs.size(); ++i) {
//    int row_idx = node->R1Cs[i].IdxR1C;
//    if (row_idx >= oldNum) {
//      record.emplace_back(i,
//                          abs(Pi[row_idx])
//                              * pow(node->R1Cs[i].Mem.size() + node->R1Cs[i].InfoR1C.size(), CONFIG::MemFactor));
//    }
//  }
//
//  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
//    int row_idx = node->R1Cs_multi[i].IdxR1C;
//    if (row_idx >= oldNum) {
//      record.emplace_back(i + node->R1Cs.size(),
//                          abs(Pi[row_idx])
//                              * pow(node->R1Cs_multi[i].Mem.size() + node->R1Cs_multi[i].InfoR1C.first.size(),
//                                    CONFIG::MemFactor));
//    }
//  }

  int MemControl = (int) (Dim * CONFIG::MaxCutMemFactor);
  for (int i = 0; i < node->R1Cs.size(); ++i) {
    int row_idx = node->R1Cs[i].IdxR1C;
    if (row_idx >= oldNum) {
      double tmp;
      if (node->R1Cs[i].Mem.size() + node->R1Cs[i].InfoR1C.size() > MemControl) {
        tmp = numeric_limits<double>::max();
      } else {
        tmp = abs(Pi[row_idx])
            * pow(node->R1Cs[i].Mem.size() + node->R1Cs[i].InfoR1C.size(), CONFIG::MemFactor);
      }
      record.emplace_back(i, tmp);
    }
  }

  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
    int row_idx = node->R1Cs_multi[i].IdxR1C;
    if (row_idx >= oldNum) {
      double tmp;
      if (node->R1Cs_multi[i].Mem.size() + node->R1Cs_multi[i].InfoR1C.first.size() > MemControl) {
        tmp = numeric_limits<double>::max();
      } else {
        tmp = abs(Pi[row_idx])
            *
                pow(node->R1Cs_multi[i].Mem.size() + node->R1Cs_multi[i].InfoR1C.first.size(), CONFIG::MemFactor);
      }
      record.emplace_back(i + node->R1Cs.size(), tmp);
    }
  }

  sort(record.begin(), record.end(), [](auto &a, auto &b) { return a.second < b.second; });

  int idx = int((double) record.size() * CONFIG::CutsKeptFactor);

  for (int i = idx; i >= 0; --i) {
    if (record[i].second < 1e100) {
      idx = i;
      break;
    }
  }

  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);
  for (int i = idx; i < record.size(); ++i) {
    int ind = record[i].first;
    int row_idx;
    if (ind < node->R1Cs.size()) {
      row_idx = node->R1Cs[ind].IdxR1C;
    } else {
      row_idx = node->R1Cs_multi[ind - node->R1Cs.size()].IdxR1C;
    }
    solver_ind[cnt++] = row_idx;
    cstr_index[row_idx] = -1;
    deleted_cstrs.emplace_back(row_idx);
  }

  if (deleted_cstrs.empty()) {
    delete[]cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))

  //rcc
  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //r1c
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //r1c_multi
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //BrC finally, BrC does not need to delete, only need to reassign
  for (auto &i : node->BrCs) {
    i.IdxBrC = cstr_index[i.IdxBrC];
  }

  convertVertex2R1CsInOneLP(node);

  delete[]cstr_index;
}

void CVRP::findBetterUsingHeurEnumeration(BBNODE *&old_node) {
  if (if_findBetterUsingHeurEnumeration_suc) return;
  if (abs(guessed_UB) < TOLERANCE) {
    guessed_UB = ceil_transformed_number_related(old_node->Val * (1 + CONFIG::EarlyHeurEnumerationGuessGap));
  } else {
    guessed_UB =
        ceil_transformed_number_related(
            (guessed_UB + old_node->Val * (1 + CONFIG::EarlyHeurEnumerationGuessGapChgFactor)) / 2);
    guessed_UB =
        max(guessed_UB,
            ceil_transformed_number_related(old_node->Val) +
                2 / pow_self(10, transformed_number));//at least two units!
  }
  if (UB <= guessed_UB) {
    if_findBetterUsingHeurEnumeration_suc = true;
    return;
  }
  cout << "guessed_UB: " << guessed_UB << endl;
  auto old_UB = UB;
  UB = guessed_UB;
  bool *is_heuristic_arc;
  BBNODE *node;
  double all;
  int max_route;
  double last_half;
#ifdef Resolve_Ins_with_Optimal_Val
  cout << "Resolve_Ins_with_Optimal_Val" << endl;
  cout << "we assume this is the found UB!" << endl;
  goto Here;
#endif
  ///Zero we find the heuristic arcs
  is_heuristic_arc = new bool[Dim * Dim]();
//  for (int i = 1; i < Dim; ++i) {
//    for (int j = i + 1; j < Dim; ++j) {
//      if (ArcGraph_revised[i][j] > TOLERANCE) {
//        is_heuristic_arc[i * Dim + j] = true;
//        is_heuristic_arc[j * Dim + i] = true;
//      }
//    }
//  }
  for (int i = 1; i < old_node->NumParentCols; ++i) {
    int past_node = 0, curr_node;
    for (size_t j = old_node->IdxCols[i] + 1;; ++j) {
      curr_node = ColPool4Mem[j];
      is_heuristic_arc[past_node * Dim + curr_node] = true;
      is_heuristic_arc[curr_node * Dim + past_node] = true;
      if (curr_node == 0) break;
      past_node = curr_node;
    }
  }
  for (int i = old_node->NumParentCols; i < NumCol; ++i) {
    int past_node = 0, curr_node;
    for (size_t j = old_node->IdxCols[i] + 1;; ++j) {
      curr_node = ColPool4Pricing[j];
      is_heuristic_arc[past_node * Dim + curr_node] = true;
      is_heuristic_arc[curr_node * Dim + past_node] = true;
      if (curr_node == 0) break;
      past_node = curr_node;
    }
  }
//  cout << "here!" << endl;
//  fill(is_heuristic_arc, is_heuristic_arc + Dim * Dim, true);

  //we copy the old_node
  node = new BBNODE(old_node, NumBucketsPerVertex, NumCol, is_heuristic_arc);
  delete[] is_heuristic_arc;

  //change hyperparameters!
  all = CONFIG::HardTimeThresholdInAllEnumeration;
  max_route = CONFIG::MaxNumRouteInEnumeration_half;
  last_half = CONFIG::HardTimeThresholdInArcElimination_last_half;
  CONFIG::HardTimeThresholdInAllEnumeration = CONFIG::HardTimeThresholdInHeurEnumeration;
  CONFIG::MaxNumRouteInEnumeration_half = CONFIG::MaxNumRouteInHeurEnumeration_half;
  CONFIG::HardTimeThresholdInArcElimination_last_half = CONFIG::HardTimeThresholdInHeurArcElimination_last_half;

  if_use_heur_enumeration = true;
  ///zero we run exact CG
  if_force_not_regenerate_bucket_graph = true;
  solveLPInLabeling(node, false, true, false);
  if_force_not_regenerate_bucket_graph = false;
  ///first we run heuristic arc elimination
  final_decision_4_arc_elimination = true;
  eliminateArcs(node);
  ///second we run heuristic enumeration
  final_decision_4_enumeration = true;
  enumerateMIP(node);

  delete node;

  //reset environment
  PoolBeg4Pricing = PoolBeg4copyPricing;
  memcpy(ColPool4Pricing, copyColPool4Pricing, sizeof(int) * PoolBeg4Pricing);
  delete[] copyColPool4Pricing;
  copyColPool4Pricing = nullptr;
  PoolBeg4copyPricing = 0;

  CONFIG::HardTimeThresholdInAllEnumeration = all;
  CONFIG::MaxNumRouteInEnumeration_half = max_route;
  CONFIG::HardTimeThresholdInArcElimination_last_half = last_half;
  safe_solver(old_node->solver.SOLVERgetNumCol(&NumCol))
  safe_solver(old_node->solver.SOLVERgetNumRow(&NumRow))
  if_use_heur_enumeration = false;

  if (abs(UB - guessed_UB) > TOLERANCE && UB < old_UB) {
    Here:
    if_findBetterUsingHeurEnumeration_suc = true;
    cout << "findBetterUsingHeurEnumeration: " << UB << endl;
    //change hyperparameters back!
    if_force_not_regenerate_bucket_graph = true;
    solveLPInLabeling(old_node);
    if_force_not_regenerate_bucket_graph = false;
    final_decision_4_arc_elimination = true;
    eliminateArcs(old_node);
    final_decision_4_enumeration = true;
    enumerateMIP(old_node);
  } else {
    double chg_val = guessed_UB * CONFIG::controlGuessGap;
    chg_val = max(chg_val, 1 / pow_self(CONFIG::MinimumGapChgUnit, transformed_number));
    chg_val = (if_enumeration_suc && if_ArcEliminationSucceed) ? chg_val : -chg_val;
    guessed_UB = ceil_transformed_number_related(guessed_UB + chg_val);
    UB = old_UB;
    cout << "findBetterUsingHeurEnumeration: no better solution found" << endl;
    cout << "UB: " << UB << endl;
    cout << "update guessed_UB: " << guessed_UB << endl;
  }
}

void CVRP::augmentNGRound(BBNODE *node) {
  if (node->Idx) return;
  //we clean the ng sets
  vector<yzzLong> old_mem = NGMem4Vertex;
  for (int i = 1; i < Dim; ++i) {
    set<int> tmp;
    set<int> tmp2;
    for (int j = 1; j < i; ++j) {
      if (ArcGraph_revised[j][i] > TOLERANCE) {
        tmp.emplace(j);
      }
    }
    for (int j = i + 1; j < Dim; ++j) {
      if (ArcGraph_revised[i][j] > TOLERANCE) {
        tmp.emplace(j);
      }
    }
    for (int j = 1; j < Dim; ++j) {
      if (NGMem4Vertex[i][j]) tmp2.emplace(j);
    }
    set<int> tmp3;
    set_intersection(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end(), inserter(tmp3, tmp3.begin()));
    NGMem4Vertex[i] = 0;
    NGMem4Vertex[i].set(i);
    for (auto n : tmp3) {
      NGMem4Vertex[i].set(n);
    }
  }
  auto last = LastMaxTimeLabeling;
  //we find the NGMemSets
  int round = 0;
  bool if_empty;
  double old_val = node->Val;
  double cost_time = 0;
  while (true) {
    ++round;
    cout << "DSSR: " << round << endl;
    findNGMemSets(node, if_empty);//do not test if empty
    //print the NGMemSets
    auto beg = chrono::high_resolution_clock::now();
    ForceNotRollback = true;
    solveLPInLabeling(node);
    ForceNotRollback = false;
    auto end = chrono::high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    if (eps > cost_time) cost_time = eps;
    if (LastMaxTimeLabeling > last and node->Val + TOLERANCE < old_val) {
      cout << "DSSR: Reach the hard limit! We use ng Rollback!" << endl;
      NGMem4Vertex = old_mem;
      deleteColByNGMem(node);
      ForceNotRollback = true;
      solveLPInLabeling(node);
      ForceNotRollback = false;
      return;
    }
    old_mem = NGMem4Vertex;
    if (node->Val + TOLERANCE > old_val) break;
  }
  round = 0;
  old_val = node->Val;
//  NGMem4Vertex = old_mem;
//  cout << "here is a test!" << endl;
//  a = accumulate(NGMem4Vertex.begin(), NGMem4Vertex.end(), 0, [](int a, yzzLong b) { return a + b.count(); }
//  );
//  cout << "aver= " << a / Dim << endl;
  double std = TOLERANCE;
  while (true) {
    ++round;
    cout << "NGAugmentation: " << round << endl;
    findNGMemSets(node, if_empty);
    if (if_empty) break;
    auto beg = chrono::high_resolution_clock::now();
    solveLPInLabeling(node);
    auto end = chrono::high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    if (eps > CONFIG::NGAugTimeHardThresholdFactor * cost_time) {
      Rollback = 1;
    } else if (eps > CONFIG::NGAugTimeSoftThresholdFactor * cost_time) {
      Rollback = 3;
    }
    cout << "eps= " << eps << endl;
    cout << "cost_time= " << cost_time << endl;
    if (Rollback == 3) break;
    else if (Rollback == 1) {
      cout << "NGAugmentation: Reach the hard limit! We use ng Rollback!" << endl;
      NGMem4Vertex = old_mem;
      deleteColByNGMem(node);
      ForceNotRollback = true;
      solveLPInLabeling(node);
      ForceNotRollback = false;
      return;
    }
    old_mem = NGMem4Vertex;
    double tmp = abs(node->Val - old_val) / node->Val;
    if (tmp > std) std = tmp;
    else if (tmp < std * CONFIG::NGAugTailOff) break;
    else {
      old_val = node->Val;
    }
  }
}

void CVRP::findNGMemSets(BBNODE *node, bool &if_empty) {
  if (node->Idx) {
    cerr << "findNGMemSets: node->Idx!=0" << endl;
    exit(1);
  }
  unordered_map<int, vector<pair<int, yzzLong>>> map_size_cycle;
  vector<size_t> tmp(Dim);
  vector<pair<size_t, double>>
      OptCols(node->Idx4LPSolsInColPool.begin() + node->NumParentColsInLPSols, node->Idx4LPSolsInColPool.end());
  sort(OptCols.begin(), OptCols.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
    return a.second > b.second;
  });
  vector<size_t> OptColIdx(OptCols.size());
  transform(OptCols.begin(),
            OptCols.end(),
            OptColIdx.begin(),
            [](const pair<int, double> &a) { return a.first; });

  int re_cols = 0;
  for (auto idx : OptColIdx) {
    memset(tmp.data(), 0, sizeof(size_t) * Dim);
    bool if_re_col = false;
    for (auto j = idx + 1;; ++j) {
      int curr_node = ColPool4Pricing[j];
      if (!curr_node) break;
      if (tmp[curr_node]) {
        if_re_col = true;
        auto length = int(j - tmp[curr_node] - 1);
        map_size_cycle[length].emplace_back(curr_node, 0);
        auto &mem = map_size_cycle[length].back().second;
        for (auto k = tmp[curr_node] + 1; k < j; ++k) {
          mem.set(ColPool4Pricing[k]);
        }
      }
      tmp[curr_node] = j;//not in else!
    }
    if (if_re_col) {
      ++re_cols;
      if (re_cols > CONFIG::MaxNumColsInNGAug) break;
    }
  }
  vector<pair<int, vector<pair<int, yzzLong>>>> size_cycle(map_size_cycle.begin(), map_size_cycle.end());
  sort(size_cycle.begin(), size_cycle.end(), [](const pair<int, vector<pair<int, yzzLong>>> &a,
                                                const pair<int, vector<pair<int, yzzLong>>> &b) {
    return a.first < b.first;
  });
  if (size_cycle.empty()) {
    cout << "size_cycle is empty!" << endl;
    if_empty = true;
    return;
  } else if_empty = false;
  for (auto &i : size_cycle[0].second) {
    int re = i.first;
    for (int j = 1; j < Dim; ++j) {
      if (i.second[j]) {
        auto &mem = NGMem4Vertex[j];
        mem.set(re);
        if (mem.count() == CONFIG::MaxNGSize) break;
      }
    }
  }
  for (int i = 1; i < size_cycle.size(); ++i) {
    if (size_cycle[i].first > CONFIG::CycleSize) break;
    for (auto &j : size_cycle[i].second) {
      int re = j.first;
      for (int k = 1; k < Dim; ++k) {
        if (j.second[k]) {
          auto &mem = NGMem4Vertex[k];
          mem.set(re);
          if (mem.count() == CONFIG::MaxNGSize) break;
        }
      }
    }
  }
  deleteColByNGMem(node);
}

void CVRP::deleteColByNGMem(BBNODE *node) {
  int len = 0, keep = node->NumParentCols, current_node;
  bool if_break;
  yzzLong PI;
  for (int i = keep; i < NumCol; ++i) {
    if_break = false;
    PI = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      current_node = ColPool4Pricing[j];
      if (!current_node) break;
      if (PI[current_node]) {
        solver_ind[len++] = i;
        if_break = true;
        break;
      }
      PI = PI & NGMem4Vertex[current_node];
      PI.set(current_node);
    }
    if (!if_break) {
      node->IdxCols[keep++] = node->IdxCols[i];
    }
  }

  if (len) {
    safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  }
}

double CVRP::calculateGapImprovement(double nowVal, double b4Val) const {
  double guess_UB = UB > 2 * b4Val ? b4Val * 1.008 : UB;
  return (nowVal - b4Val) / (guess_UB - b4Val);
}

void CVRP::deleteArcByBrC_false(BBNODE *node, std::pair<int, int> edge) const {
  int i = edge.first, j = edge.second;
  bool if_rep = false;
  AGAIN:
  for (int bin = 0; bin < NumBucketsPerVertex; ++bin) {
    auto if_find = std::find(node->AllForwardBuckets[i][bin].BucketArcs.begin(),
                             node->AllForwardBuckets[i][bin].BucketArcs.end(), j);
    if (if_find == node->AllForwardBuckets[i][bin].BucketArcs.end()) {
      //jump arc is a pair, compare the second
      auto iff = std::find_if(node->AllForwardBuckets[i][bin].JumpArcs.begin(),
                              node->AllForwardBuckets[i][bin].JumpArcs.end(),
                              [&](const pair<double, int> &p) { return p.second == j; });
      if (iff != node->AllForwardBuckets[i][bin].JumpArcs.end()) {
        node->AllForwardBuckets[i][bin].JumpArcs.erase(iff);
      }
    } else {
      node->AllForwardBuckets[i][bin].BucketArcs.erase(if_find);
    }
  }
  if (!if_rep) {
    if_rep = true;
    std::swap(i, j);
    goto AGAIN;
  }
}

void CVRP::checkIfCutsLegal(BBNODE *node) const {
  vector<bool> if_takeup(NumRow, false);
  for (int i = 0; i < Dim; ++i)if_takeup[i] = true;
  for (auto &rcc : node->RCCs) {
    if (if_takeup[rcc.IdxRCC]) throw std::runtime_error("Error in RCCs");
    if_takeup[rcc.IdxRCC] = true;
  }
  for (auto &r1c : node->R1Cs) {
    if (if_takeup[r1c.IdxR1C]) throw std::runtime_error("Error in R1Cs");
    if_takeup[r1c.IdxR1C] = true;
  }
  for (auto &r1c_multi : node->R1Cs_multi) {
    if (if_takeup[r1c_multi.IdxR1C]) throw std::runtime_error("Error in R1CMulti");
    if_takeup[r1c_multi.IdxR1C] = true;
  }
  for (auto &brc : node->BrCs) {
    if (brc.IdxBrC != -1) {
      if (if_takeup[brc.IdxBrC])throw std::runtime_error("Error in BrCs");
      if_takeup[brc.IdxBrC] = true;
    }
  }
  for (int i = 0; i < NumRow; ++i) {
    if (!if_takeup[i]) {
      cout << "i: " << i << endl;
      throw std::runtime_error("Error in check_if_cuts_are_legal");
    }
  }

  safe_solver(node->solver.SOLVERupdatemodel())
  for (int i = 0; i < NumRow; i++) {
    double cof, rhs;
    safe_solver(node->solver.SOLVERgetcoeff(i, 0, cof))
    safe_solver(node->solver.SOLVERgetRhs(i, 1, &rhs))
    if (cof != rhs) {
      cout << "i: " << i << endl;
      for (auto &rcc : node->RCCs) {
        if (i == rcc.IdxRCC) {
          cout << "rcc: " << rcc.RHS << endl;
          break;
        }
      }
      for (auto &r1c : node->R1Cs) {
        if (i == r1c.IdxR1C) {
          cout << "r1c: " << r1c.RHS << endl;
          break;
        }
      }
      for (auto &r2c : node->R1Cs_multi) {
        if (i == r2c.IdxR1C) {
          cout << "r2c: " << r2c.RHS << endl;
          break;
        }
      }
      cout << "cof: " << cof << " rhs: " << rhs << endl;
      throw std::runtime_error("Error in not equal");
    }
  }
  /**
   * more test about the coefficients
   */

  cout << "check_if_cuts_are_legal more!" << endl;
  cout << "tmp not check rcc..." << endl;
  //check r1c

}

void CVRP::updateUB_EdgeSol() {
  EdgeSol.clear();
  EdgeSol.reserve(2 * Dim);
  for (auto &r : IPOptSol) {
    int past_node = 0;
    for (int i = 1; i < r.size(); ++i) {
      int current_node = r[i];
      auto pr =
          past_node < current_node ? make_pair(past_node, current_node) : make_pair(current_node, past_node);
      EdgeSol.emplace(pr);
      past_node = current_node;
    }
  }
  //print_edge_sol();
//  for (auto &e : EdgeSol) {
//    if (e.first == 0 && e.second == 0) throw std::runtime_error("Error in updateUB_EdgeSol");
//    else cout << "(" << e.first << "," << e.second << ")" << endl;
//  }
}

void CVRP::readSolFile(bool if_force) {
  string fileName = "sol_used/" + FileName + ".sol";
  ifstream fin(fileName);
  if (!fin && if_force) {
    throw std::runtime_error("Error in readSolFile");
  }
  if (!fin) return;
  IPOptSol.clear();
  double optval = 0;
  string line, item;
  while (std::getline(fin, line)) {
    if (line.find("<Optval") != std::string::npos) {
      stringstream ss(line);
      string optval_str;
      ss >> optval_str >> optval;
    } else if (line[0] == '0' && line[1] == '-') {
      stringstream ss(line);
      vector<int> tmp;
      tmp.reserve(MaxLengthEleRoute);
      while (getline(ss, item, '-')) {
        tmp.emplace_back(std::stoi(item));
      }
      IPOptSol.emplace_back(std::move(tmp));
    }
  }
  fin.close();
  if (IPOptSol.empty() && if_force) {
    throw std::runtime_error("Error in No solution!");
  }
  cout << "optval: " << optval << " UB: " << UB << endl;
  if (abs(optval - UB) < TOLERANCE || optval < UB) {
    cout << "UB is updated to: " << optval << endl;
  } else {
    cout << "we abandoned the old UB (though better, we want the solution) and update to " << optval << endl;
  }
  UB = optval;
  updateUB_EdgeSol();
}

void CVRP::reviseBranch_pair() {
  if (EdgeSol.empty()) {
    cout << "EdgeSol is empty" << endl;
    return;
  }
  vector<pair<int, int>> new_branch_pair;
  for (auto &edge : Branch_pair) {
    if (EdgeSol.find(edge) != EdgeSol.end()) {
      new_branch_pair.emplace_back(edge);
    }
  }
  if (new_branch_pair.empty()) {
    cout << "new_branch_pair is empty" << endl;
    return;
  } else {
    Branch_pair = new_branch_pair;
    cout << "Branch_pair is revised and its size is " << Branch_pair.size() << endl;
  }
}

void CVRP::updateLB(double val) {
  LB = val;
  LB_transformed = ceil_transformed_number_related(LB - TOLERANCE);
}