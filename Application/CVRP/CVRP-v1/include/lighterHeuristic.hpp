
#ifndef CVRP_LIGHTERHEURISTIC_HPP
#define CVRP_LIGHTERHEURISTIC_HPP
#include "CVRP.hpp"
#include "templateFunctors.hpp"

template<bool if_symmetry>
int CVRP::generateColumnsByLighterHeuristic(BbNode *node) {

  if_exact_labeling_cg = false;
  priceLabeling(node, pi4_labeling);

  runLabeling<true, false, if_symmetry, 2>(node);

  if (!if_symmetry) {
	runLabeling<false, false, if_symmetry, 2>(node);
  }

  auto tmp = Config::MaxNumRoutesInExact;
  Config::MaxNumRoutesInExact = Config::MaxNumRoutesInLighterHeur;

  int ccnt = concatenateCols_prior_forward<if_symmetry>(node);

  Config::MaxNumRoutesInExact = tmp;

  if (!ccnt) return 0;

  addColumns(node, ccnt);

  return ccnt;
}

#endif //CVRP_LIGHTERHEURISTIC_HPP

