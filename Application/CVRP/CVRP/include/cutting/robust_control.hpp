//
// Created by You, Zhengzhong on 7/23/24.
//

#ifndef CUTTING_ROBUST_CONTROL_HPP_
#define CUTTING_ROBUST_CONTROL_HPP_

#include "cvrp.hpp"

#define MINIMUM_CUTS_ADDED 50

class RobustControl {
 public:
  static CVRP *cvrp;
  static BbNode *node;
  static int initial_max_num_r1c_per_round;
  static int initial_max_num_r1c3_per_round;
  static bool rank1_cuts_mode;//0-node based, 1-arc based
  static std::vector<R1c> new_generated_cuts;
  static bool if_fix_resource_point;
  static int pricing_hard_level;//0:extremely easy, 1:be careful; 2. hard;
  static bool if_ever_roll_back;
  static double val_b4_rank1;
  static std::vector<yzzLong> solution_arcs;

 public:
  static void init(CVRP *pr_cvrp);
  static void getInitialMaxNumR1CPerRound();
  static void estimatedNumberOfNonRobustCutsCanBeAdded();

  static void updateNode(BbNode *pr_node);
  static void recordNewGeneratedCuts(int oldNum);
  static bool restoreCuts();
  static void recordSolutionArcs();
  static void convert2ArcBasedMemory();
  static void robustControl();
  static void robustControl(int oldNum, double &prior_nodeVal, int &goto_state);
};

#endif //CUTTING_ROBUST_CONTROL_HPP_
