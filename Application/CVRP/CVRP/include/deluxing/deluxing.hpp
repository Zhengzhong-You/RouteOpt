#ifndef DELUXING_H
#define DELUXING_H

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "gurobi_c.h"


double getWallTime();


void deLuxing(GRBmodel *orig, double UB, int NClust, int beta1, int beta2, std::vector<int> &Idxdel, double Timelimit,
              double Tolerance, bool Verbose);

#endif
