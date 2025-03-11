/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <vector>
#include <string>

/**
* https://github.com/vidalt/HGS-CVRP
 * thanks Vidal T. N. for sharing this code
*/

/**
 * Now this code only uses very little of the whole function
 */

typedef std::pair<std::vector<std::vector<int>>, double> SolutionGroup;


void HGS(const std::string &pathInstance,
		 bool verbose,
		 int seed,
		 double timeLimit,
		 std::vector<std::vector<int>> &sol,
		 double &ub);

void HGSWithAllSolution(const std::string &pathInstance,
						bool verbose,
						int seed,
						double timeLimit,
						std::vector<SolutionGroup> &sol);