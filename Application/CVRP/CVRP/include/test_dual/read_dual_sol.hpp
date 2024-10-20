//
// Created by You, Zhengzhong on 7/10/24.
//

#ifndef INCLUDE_TEST_DUAL_READ_DUAL_SOL_HPP_
#define INCLUDE_TEST_DUAL_READ_DUAL_SOL_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"

class ReadDualSol {
 public:
  static std::vector<double> dual_sol;
  static double coeff;
  static void resetCoeff();
  static void readDualSol(const std::string &file_name);
  static void changeDualSol(std::vector<double> &dual);
  static void checkIfMisPrice(int &cnt);
  static void changeRCStd(BbNode *node, double &rc_std, const std::vector<double> &dual);
};

#endif //INCLUDE_TEST_DUAL_READ_DUAL_SOL_HPP_
