//
// Created by You, Zhengzhong on 6/5/24.
//

#include "cvrp.hpp"

using namespace std;

void getNewRhsByMem(BbNode *node, vector<double> &new_model_rhs);

void CVRP::getNewRhs4ChangingDual(BbNode *node, vector<double> &new_model_rhs) {
  getNewRhsByMem(node, new_model_rhs);
}

bool CVRP::stopChangingDualCondition(BbNode *node) {
  if (node->r1cs.empty() || !if_exact_labeling_cg) return true;
  return false;
}

void getNewRhsByMem(BbNode *node, vector<double> &new_model_rhs) {
  for (auto &r1c : node->r1cs) {
	new_model_rhs[r1c.idx_r1c] = pow((double)r1c.arc_mem.size(), 2);
  }
}
