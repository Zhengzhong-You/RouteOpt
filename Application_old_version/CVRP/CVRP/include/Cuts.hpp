
#ifndef CVRP_CUTS_HPP
#define CVRP_CUTS_HPP

#include "Config.hpp"
#include <vector>
#include <set>

struct Rcc {
  bool form_rcc;
  std::vector<int> info_rcc_customer;
  std::vector<int> info_rcc_outside_customer;
  double rhs;
  int idx_rcc;
#ifdef SOLVER_VRPTW
  bool if_keep{};// default: false;
#endif
};

struct R1c {
  std::vector<int> info_r1c;
  std::vector<int> mem;
  int idx_r1c;
  int rhs;// = int(info_r1c.size()/2)
};

struct R1cMulti {
  std::pair<std::vector<int>, int> info_r1c;// cut and plan
  std::vector<int> mem;
  int idx_r1c;
  int rhs;// get<2>map[cut.size()][plan_idx]
};

struct Brc {
  std::pair<int, int> edge;
  int idx_br_c;
  bool br_dir;
};
#endif //CVRP_CUTS_HPP
