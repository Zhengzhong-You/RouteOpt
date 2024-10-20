//
// Created by You, Zhengzhong on 5/5/24.
//

#ifndef INCLUDE_BRANCHING_SETTING2_HPP_
#define INCLUDE_BRANCHING_SETTING2_HPP_

/**
 * setting 2 is used for use less rank1 cuts and no enumeration to focus on the performance of the branching section
 */


#include <iostream>

class Setting2 {
 public:
  static void easeRank1Cuts(int &max_num_r1c3_per_round,
							int &max_num_r1c_per_round,
							int &cuts_tolerance,
							double &cuts_tail_off) {
	max_num_r1c3_per_round = 60;
	max_num_r1c_per_round = 40;
	cuts_tolerance = 1;
	cuts_tail_off = 0.15;
  }

  static void banCutsAtNonRoot(bool root_node, bool &if_ban_cuts) {
	if_ban_cuts = !root_node;
	if (if_ban_cuts)
	  std::cout << "In setting II, we do not do cuts for non root node!" << std::endl;
  }

  static void banEnumeration(bool &if_ban_enu) {
	if_ban_enu = true;
	std::cout << "In setting II, we do not do enumeration!" << std::endl;
  }

  static void banArcElimination(bool root_node, bool &if_ban_arc_elimination) {
	if_ban_arc_elimination = !root_node;
	if (if_ban_arc_elimination)
	  std::cout << "In setting II, we do not do arc elimination for non root node!" << std::endl;
  }
};

#endif //INCLUDE_BRANCHING_SETTING2_HPP_
