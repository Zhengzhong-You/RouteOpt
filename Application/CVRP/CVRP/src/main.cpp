#include "cvrp.hpp"
#ifdef SOLVER_VRPTW
#include "vrptw.hpp"
#endif
#ifdef SOLVER_RVRPSTW
#include "rvrpstw.hpp"
#endif
#include "instance_data.hpp"
#include "solver_interface/solver.hpp"
#if SOLUTION_TYPE == 1
#include "best_bound_first_branching.hpp"
#elif SOLUTION_TYPE == 2
#include "depth_first_branching.hpp"
#endif

#ifdef READ_ENUMERATION_TREES
#include "read_enumeration_tree.hpp"
#endif

// #include <sys/resource.h>

using namespace std;

void solveRVRPSTW(int argc, char *argv[]);

void solveCVRP_N_VRPTW(int argc, char *argv[]);

int main(int argc, char *argv[]) {
#ifdef SOLVER_RVRPSTW
  solveRVRPSTW(argc, argv);
#else
    solveCVRP_N_VRPTW(argc, argv);
#endif
    return 0;
}

#ifdef SOLVER_RVRPSTW
void solveRVRPSTW(int argc, char *argv[]) {
  cout << "SOLVER_RVRPSTW" << endl;
  string p = generateInstancePath(argc, argv);
  cout << "Instance Path: " << p << endl;
  InstanceData data;
  takeDataFromFile(p, data);
  auto *solver = new RVRPSTW(data);
  solver->seed = std::atoi(argv[argc - 1]);
  cout << "seed: " << solver->seed << endl;
  BaseBranching::init(solver);
  BaseBranching::prepare();
  while (true) {
#if SOLUTION_TYPE == 1
	BestBoundFirstBranching::solve(BestBoundFirstBranching::bbt);
#elif SOLUTION_TYPE == 2
	DepthFirstBranching::solve(DepthFirstBranching::bbt);
#endif
	if (!solver->separateRobustifyingCut()) break;
  }
  delete solver;
}
#endif

void solveCVRP_N_VRPTW(int argc, char *argv[]) {
#ifdef SOLVER_VRPTW
    cout << "SOLVER_VRPTW | solution mode= " << SOLVER_VRPTW << endl;
#else
    cout << "SOLVER_CVRP" << endl;
#endif

    // rusage usage{};
    // getrusage(RUSAGE_SELF, &usage);
    // cout << "Memory usage: " << usage.ru_maxrss << endl;

    string p = generateInstancePath(argc, argv);
    cout << "Instance Path: " << p << endl;
    InstanceData data;
    takeDataFromFile(p, data);

    Config::ub += 1;
    cout << "add one unit to ub: " << Config::ub << endl;

#ifdef SOLVER_VRPTW
    auto *solver = new VRPTW(data);
    solver->tryGetTravelTimeMatrix(p);
#else
    auto *solver = new CVRP(data);
#endif


#ifdef HGS_APPLIED
  cout << "HGS_APPLIED" << endl;
#ifndef SOLVER_VRPTW
#ifdef FIND_ALL_SOLUTIONS
  vector<pair<vector<vector<int>>, double>> sol;
  double heur_UB;
  double time = max(1200., 10. + 2 * pow((solver->dim / 20.), 2));
  HGSWithAllSolution(p, true, 1, time, sol);
  for (auto &route : sol) {
	for (auto &r : route.first) {
	  if (r.empty()) continue;
	  solver->ip_opt_sol.emplace_back(r);
	}
  }
#else
  vector<vector<int>> sol;
  double heur_UB;
  double time = max(1200., 10. + 2 * pow((solver->dim / 20.), 2));
  HGS(p, true, 1, time, sol, heur_UB);
  if (Config::ub > heur_UB + TOLERANCE) {
	Config::ub = heur_UB;
	cout << "ub updated by HGS: " << Config::ub << endl;
	solver->ip_opt_sol.reserve(sol.size());
	for (auto &route : sol) {
	  if (route.empty()) continue;
	  solver->ip_opt_sol.emplace_back(route);
	}
  }
#endif
#else
  cout << "HGS_APPLIED cannot run in VRPTW" << endl;
#endif
#endif

    BaseBranching::init(solver);
    BaseBranching::prepare();

#ifdef READ_ENUMERATION_TREES
#if SOLUTION_TYPE == 1
  BestBoundFirstBranching::solve(BestBoundFirstBranching::sub_bbt);
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::solve(DepthFirstBranching::sub_bbt);
#endif
#else
#if SOLUTION_TYPE == 1
    BestBoundFirstBranching::solve(BestBoundFirstBranching::bbt);
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::solve(DepthFirstBranching::bbt);
#endif
#endif

    delete solver;

#ifdef SOLVER_STATISTICS
  STATISTICS::giveStatistics(data.name);
#endif
}
