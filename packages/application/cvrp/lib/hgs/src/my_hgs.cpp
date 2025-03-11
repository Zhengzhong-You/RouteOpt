/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "Genetic.h"
#include "InstanceCVRPLIB.h"
#include "AlgorithmParameters.h"
#include "Params.h"
#include <numeric>
#include "my_hgs.hpp"

/**
* https://github.com/vidalt/HGS-CVRP
 * thanks Vidal T. N. for sharing this code
*/

/**
 * Now this code only uses very little of the whole function
 */

void HGS(const std::string &pathInstance,
		 bool verbose,
		 int seed,
		 double timeLimit,
		 std::vector<std::vector<int>> &sol,
		 double &ub) {
  InstanceCVRPLIB cvrp(pathInstance, true);

  AlgorithmParameters ap = default_algorithm_parameters();
  ap.seed = seed;
  ap.timeLimit = timeLimit;
  Params params(cvrp.x_coords,
				cvrp.y_coords,
				cvrp.dist_mtx,
				cvrp.service_time,
				cvrp.demands,
				cvrp.vehicleCapacity,
				cvrp.durationLimit,
				std::numeric_limits<int>::max(),
				cvrp.isDurationConstraint,
				verbose,
				ap);

  // Running HGS
  Genetic solver(params);
  solver.run();

  if (solver.population.getBestFound() != NULL) {
	auto &indiv = *solver.population.getBestFound();
	sol = indiv.chromR;
	ub = indiv.eval.penalizedCost;
  } else {
	sol = {};
	ub = std::numeric_limits<double>::max();
  }
}

void HGSWithAllSolution(const std::string &pathInstance,
						bool verbose,
						int seed,
						double timeLimit,
						std::vector<SolutionGroup> &sol) {
  InstanceCVRPLIB cvrp(pathInstance, true);

  AlgorithmParameters ap = default_algorithm_parameters();
  ap.seed = seed;
  ap.timeLimit = timeLimit;
  Params params(cvrp.x_coords,
				cvrp.y_coords,
				cvrp.dist_mtx,
				cvrp.service_time,
				cvrp.demands,
				cvrp.vehicleCapacity,
				cvrp.durationLimit,
				std::numeric_limits<int>::max(),
				cvrp.isDurationConstraint,
				verbose,
				ap);

  // Running HGS
  Genetic solver(params);
  solver.run();

  sol = {};

  if (!solver.population.feasibleSubpop.empty()) {
	for (auto &ind : solver.population.feasibleSubpop) {
	  auto &indiv = *ind;
	  SolutionGroup sg;
	  sg.first = indiv.chromR;
	  sg.second = indiv.eval.penalizedCost;
	  sol.emplace_back(sg);
	}
  }
}
