#ifndef CVRP_CUTS_HPP
#define CVRP_CUTS_HPP

#include "config.hpp"
#include "rank1_cuts_separator.hpp"
#include <vector>
#include <set>
#include <unordered_set>

struct R1c;

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

struct Brc {
    std::pair<int, int> edge;
    int idx_br_c;
    bool br_dir;
};
#endif //CVRP_CUTS_HPP
