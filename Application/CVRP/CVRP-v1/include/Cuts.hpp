
#ifndef CVRP_CUTS_HPP
#define CVRP_CUTS_HPP

#include "Config.hpp"
#include <vector>
#include <set>
#include <unordered_set>

struct Rcc {
  /**
   * x(S: S) \leq|S|-k(S)
   * x(\bar{S}: \bar{S})+\frac{1}{2} x(\{0\}: \bar{S})-\frac{1}{2} x(\{0\}: S) \leq|\bar{S}|-k(S)
   *
   * in enumeration, the rcc can be strengthened by
   * \sum_{p \in P} I_{(|p \cap S| \geq 1)} x_p \geq\left\lceil\frac{1}{Q} \sum_{i \in S} q_i\right\rceil
   */
  bool form_rcc;
  std::vector<int> info_rcc_customer;
  std::vector<int> info_rcc_outside_customer;
  double rhs;
  int idx_rcc;
#if  SOLVER_VRPTW == 1
  bool if_keep{};// default: false;
#endif
};

struct R1c {
  std::pair<std::vector<int>, int> info_r1c;// cut and plan
  int idx_r1c;
  int rhs;// get<2>map[cut.size()][plan_idx]
  /**
   * arc memory: any vertex to info_r1c.first + some vertexes to mem
   */
  std::vector<int> mem;
  std::vector<std::pair<std::vector<int>, int>>
	  arc_mem;// some vertexes to mem, never include some to info_r1c, like in-flow
  double convert_dual_coeff{1};
};

struct Brc {
  std::pair<int, int> edge;
  int idx_br_c;
  bool br_dir;
};
#endif //CVRP_CUTS_HPP
