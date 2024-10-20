
#ifndef CVRP_HEAVIERHEURISTIC_HPP
#define CVRP_HEAVIERHEURISTIC_HPP
#include "cvrp.hpp"
#include "template_functors.hpp"

template<bool if_symmetry>
int CVRP::generateColumnsByHeavierHeuristic(BbNode *node) {

  if_exact_labeling_cg = false;
  getRank1DualsInCG(node, pi4_labeling);

  runLabeling<true, false, false, if_symmetry, 1>(node);
  if (if_short_memory) {
	throw std::runtime_error("Lack of memory even for heuristic labeling");
  }

  if (!if_symmetry) {
	runLabeling<false, false, false, if_symmetry, 1>(node);
	if (if_short_memory) {
	  throw std::runtime_error("Lack of memory even for heuristic labeling");
	}
  }

  auto tmp = Config::MaxNumRoutesInExact;
  Config::MaxNumRoutesInExact = Config::MaxNumRoutesInHeavierHeur;

  int ccnt = concatenateCols_prior_forward<if_symmetry>(node);

  Config::MaxNumRoutesInExact = tmp;

  if (!ccnt) return 0;

  addColumns(node, ccnt);

  return ccnt;
}

#endif //CVRP_HEAVIERHEURISTIC_HPP

