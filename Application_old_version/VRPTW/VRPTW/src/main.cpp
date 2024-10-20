#include "CVRP.hpp"
#include "VRPTW.hpp"
#include "InstanceData.hpp"
#include "Solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
#ifdef SOLVER_VRPTW
  cout << "SOLVER_VRPTW" << endl;
#else
  cout << "SOLVER_CVRP" << endl;
#endif
  Config::checkParameters();
  string p = generateInstancePath(argc, argv);
  cout << "Instance Path: " << p << endl;
  InstanceData data;
  takeDataFromFile(p, data);

#ifdef SOLVER_VRPTW
  auto *solver = new VRPTW(data);
#else
  auto *solver = new CVRP(data);
#endif

#ifdef HGS_APPLIED
  cout << "HGS_APPLIED" << endl;
#ifndef SOLVER_VRPTW
  vector<vector<int>> sol;
  double heur_UB;
  HGS(p, true, 1, 10 + 2 * pow((solver->dim / 20.), 2), sol, heur_UB);
  if (Config::ub > heur_UB + TOLERANCE) {
	Config::ub = heur_UB;
	cout << "ub updated by HGS: " << Config::ub << endl;
	solver->ip_opt_sol.reserve(sol.size());
	for (auto &route : sol) {
	  if (route.empty()) continue;
	  vector<int> tmp(route.size() + 2);
	  tmp[0] = 0;
	  tmp[route.size() + 1] = 0;
	  for (int i = 0; i < route.size(); ++i) {
		tmp[i + 1] = route[i];
	  }
	  solver->ip_opt_sol.emplace_back(std::move(tmp));
	}
  }
#else
  cout << "HGS_APPLIED cannot run in VRPTW" << endl;
#endif
#endif

#ifdef MASTER_VALVE_ML
  solver->ml.get_info(solver);
#endif

#ifdef READ_ENUMERATION_TREES
  solver->solveEnumTree();
#else
#ifdef BRANCH_FASHION_MEM_SAVING
  solver->solveModel();
#else
  solver->solveSBModel();
#endif
#endif
  delete solver;

#ifdef SOLVER_STATISTICS
  STATISTICS::giveStatistics(data.name);
#endif

  return 0;
}

