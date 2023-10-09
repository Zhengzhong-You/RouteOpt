#include "CVRP.hpp"
#include "VRPTW.hpp"
#include "InstanceData.hpp"
#include "solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
#ifdef SOLVER_VRPTW
  cout << "SOLVER_VRPTW" << endl;
  cout << "MaxCustomerNum= " << MaxNum_Customers << endl;
  cout << "Initial ng= " << CONFIG::InitialNGSize << endl;
#else
  cout << "SOLVER_CVRP" << endl;
#endif
  //to-do need print all parameters
  cout << "InitialNumBuckets= " << CONFIG::InitialNumBuckets << endl;
  CONFIG::checkParameters();

  string p = generateInstancePath(argc, argv);
  cout << "Instance Path: " << p << endl;

  InstanceData data;
  takeDataFromFile(p, data);
#ifdef Resolve_Ins_with_Optimal_Val
  CONFIG::MaxNumLabelInEnumeration = 500000;
  CONFIG::MIPInEnumerationTimeLimit = 20;
  CONFIG::MaxNumRouteInEnumeration_half = 5000;
  CONFIG::MaxNumRouteInEnumeration = 1000000;
  bool if_rep = false;
#endif

#ifdef SOLVER_VRPTW
  auto *solver = new VRPTW(data);

//  GRBenv *env = NULL;
//  int error = GRBloadenv(&env, NULL);
//  if (error) {
//    printf("无法创建Gurobi环境: %d\n", error);
//    exit(1);
//  }
//
//  error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
//  if (error) {
//    printf("无法设置输出标志: %d\n", error);
//    exit(1);
//  }
//
//  GRBmodel *model = nullptr;
//  error = GRBnewmodel(env, &model, "model.lp", 0, NULL, NULL, NULL, NULL, NULL);
//  if (error) {
//    printf("无法创建模型: %d\n", error);
//    exit(1);
//  }
//
//  // 从MPS文件加载模型
//  error = GRBreadmodel(env, "42200.800100.rcf.mps", &model);
//  if (error) {
//    printf("无法从MPS文件加载模型: %d\n", error);
//    exit(1);
//  }
//
//  vector<int> idx;
//
//  auto beg = chrono::high_resolution_clock::now();
////  heavyRCF(solver.model, UB + roundUpTolerance, round, idx, MIP_TOLERANCE, if_verbose);
//  heavyRCF(model, CONFIG::UB + solver->roundUpTolerance, 200, idx, MIP_TOLERANCE, true);
//  auto end = chrono::high_resolution_clock::now();
//  cout << "RCF time: " << chrono::duration<double>(end - beg).count() << " s" << endl;
//  cout << "RCF over!" << endl;
//
//  exit(0);

#else
  auto *solver = new CVRP(data);
#endif

#ifdef HGS_APPLIED
  cout << "HGS_APPLIED" << endl;
#ifndef SOLVER_VRPTW
//  CONFIG::MaxNumLabelInEnumeration = 500000;
//  CONFIG::MIPInEnumerationTimeLimit = 20;
//  CONFIG::MaxNumRouteInEnumeration_half = 5000;
//  CONFIG::MaxNumRouteInEnumeration = 1000000;
  vector<vector<int>> sol;
  double heur_UB;
  HGS(p, true, 1, 10 + 2 * pow((solver->Dim / 20.), 2), sol, heur_UB);
  if (CONFIG::UB > heur_UB + TOLERANCE) {
    CONFIG::UB = heur_UB;
    cout << "UB updated by HGS: " << CONFIG::UB << endl;
    solver->IPOptSol.reserve(sol.size());
    for (auto &route : sol) {
      if (route.empty()) continue;
      vector<int> tmp(route.size() + 2);
      tmp[0] = 0;
      tmp[route.size() + 1] = 0;
      for (int i = 0; i < route.size(); ++i) {
        tmp[i + 1] = route[i];
      }
      solver->IPOptSol.emplace_back(std::move(tmp));
    }
  }
#else
  cout << "HGS_APPLIED cannot run in VRPTW" << endl;
#endif
#endif

  here:
#ifdef MASTER_VALVE_ML
  solver->ml.getInfo(solver);
#endif

#ifdef readEnumerationTrees
  solver->solveEnuTree();
#else
#ifdef BranchFashion_MemSaving
  solver->solveModel();
#else
  solver->solveSBModel();
#endif
#endif
  delete solver;

#ifdef Resolve_Ins_with_Optimal_Val
  if (!if_rep) {
    if_rep = true;
    cout << "changed here!" << endl;
    CONFIG::MaxNumLabelInEnumeration = 50000;
    CONFIG::MIPInEnumerationTimeLimit = 100;
    CONFIG::MaxNumRouteInEnumeration_half = 10000;
    CONFIG::MaxNumRouteInEnumeration = 10000;
    ++CONFIG::Marker4FindLargeGapIns;
    if (!CONFIG::Marker4FindLargeGapIns) return 0;// if the old one is -1, then it is the bad instance!
#ifdef SOLVER_VRPTW
    solver = new VRPTW(data);
#else
    solver = new CVRP(data);
#endif
    goto here;
  }
#endif

#ifdef SOLVER_STATISTICS
  STATISTICS::giveStatistics(data.name);
#endif

  return 0;
}

