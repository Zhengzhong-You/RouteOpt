//
// Created by Zhengzhong You on 7/7/22.
//

#ifndef CVRP_CUTS_HPP
#define CVRP_CUTS_HPP

#include "CONFIG.hpp"
#include <vector>
#include <set>

struct RCC {
  bool FormRCC;
  std::vector<int> InfoRCCCustomer;
  std::vector<int> InfoRCCOutsideCustomer;
  double RHS;
  int IdxRCC;
};

struct R1C {
  //we know what kind of cuts by size of InfoR1C
  std::vector<int> InfoR1C;
  std::vector<int> Mem;
  int IdxR1C;
  int RHS;// = int(InfoR1C.size()/2)
};

struct R1C_multi {
  //treat these cuts different due to performance
  std::pair<std::vector<int>, int> InfoR1C;// cut and plan
  std::vector<int> Mem;
  int IdxR1C;
  int RHS;// get<2>map[cut.size()][plan_idx]
};

struct BrC {
  std::pair<int, int> Edge;
  int IdxBrC;
  bool BrDir;
};

struct NBrC {
  //except for the cols that already tagged,
  //only index will be needed!
  int IdxNBrC;
};

#endif //CVRP_CUTS_HPP
