//
// Created by You, Zhengzhong on 5/2/24.
//

#ifndef INCLUDE_PRICING_DUAL_DUAL_SMOOTHING_HPP_
#define INCLUDE_PRICING_DUAL_DUAL_SMOOTHING_HPP_

#include "cvrp.hpp"
#include "branching_node.hpp"

#define MODE_IN_BEST_LAGRANGE_BOUND 1
#define MODE_IN_LAST_ITERATION 2
#define NORM_TOLERANCE 1e-4

class DualSmoothing {
 public:
  static CVRP *cvrp;
  static BbNode *node;
  static int num_k;
  static double alpha;
  static double lp_value;
  static double best_lagrange_bound;
  static double most_negative_reduced_cost_in_lp;
  static int most_negative_reduced_cost_col_idx;
  static std::vector<double> most_negative_reduced_cost_col;
  static std::vector<double> rhs_lp_matrix;
  static denseColMatrixXd lp_matrix;
  static RowVectorXd lp_obj;
  static std::vector<double> best_lagrange_multipliers;
  static std::vector<double> pi_in;
  static std::vector<double> pi_out;
  static std::vector<double> pi_price;
  static bool if_mis_pricing;
  static int counter_mis_pricing;

  static void reset();
  static void init(CVRP *pr_cvrp);
  static void initNumK();
  static void updateNode(BbNode *pr_node);
  static void updateBestLagrangeBound();
  static void updatePiIn(int mode);
  static void updatePiOut();
  static void updatePiPrice();
  static void tellIfMisPrice(int &ccnt);
  static void buildLPMatrix();
  static void calculateMostNegativeReducedCostInLP();
  static void addCol2LPMatrix(int ccnt_cnt, std::vector<size_t> &solver_beg,
							  std::vector<int> &solver_ind, std::vector<double> &solver_val,
							  std::vector<double> &solver_obj);
  static void deleteCol4LPMatrix(std::vector<int> &col_idx);
  static void adjustAlpha();
  static void recordLPValue();
  static void freeLPMatrix();
  static void changeRCStd();
  static void printInfo();
};

#endif //INCLUDE_PRICING_DUAL_DUAL_SMOOTHING_HPP_
