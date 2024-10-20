//
// Created by You, Zhengzhong on 5/7/24.
//


#include "combine_faster_smoothing.hpp"

using namespace std;

CVRP *CombineFasterSmoothing::cvrp{};
BbNode *CombineFasterSmoothing::node{};
double CombineFasterSmoothing::beta{1};
double CombineFasterSmoothing::gamma{0.002};
double CombineFasterSmoothing::lp2_value{};
double CombineFasterSmoothing::lp2_norm{};
std::vector<double> CombineFasterSmoothing::lp2_obj{};
std::vector<double> CombineFasterSmoothing::pi_lp1{};
std::vector<double> CombineFasterSmoothing::pi_lp2{};
void calculateNormDistance(const std::vector<double> &vec1, const std::vector<double> &vec2, double &norm);

void CombineFasterSmoothing::init(CVRP *pr_cvrp) {
  cvrp = pr_cvrp;
  DualSmoothing::init(pr_cvrp);
  DualSmoothing::initNumK();
  FasterDual::init(pr_cvrp);
}

void CombineFasterSmoothing::updateNode(BbNode *pr_node) {
  node = pr_node;
  DualSmoothing::updateNode(pr_node);
  FasterDual::updateNode(pr_node);
}

void CombineFasterSmoothing::getLP2Obj(const std::vector<double> &pr_lp2_obj) {
  lp2_obj = pr_lp2_obj;
  lp2_norm = std::sqrt(std::inner_product(lp2_obj.begin(), lp2_obj.end(), lp2_obj.begin(), 0.0,
										  std::plus<>(),
										  [](double a, double b) {
											return pow(b, 2);
										  }));
}

void CombineFasterSmoothing::getLP2Value(double &pr_lp2_value) {
  lp2_value = pr_lp2_value;
}

void CombineFasterSmoothing::updatePiLP1() {
  pi_lp1 = cvrp->pi4_labeling;
}

void CombineFasterSmoothing::updatePiLP2() {
  pi_lp2 = cvrp->pi4_labeling;
}

void CombineFasterSmoothing::updateSmoothingState() {
  if (FasterDual::sub_pricing_solver.model) {
	double norm1, norm2;
	calculateNormDistance(DualSmoothing::pi_in, pi_lp1, norm1);
	calculateNormDistance(DualSmoothing::pi_in, pi_lp2, norm2);
	if (norm1 / norm2 < beta) {
	  vector<double> smooth_dual(DualSmoothing::pi_in.size());
	  transform(DualSmoothing::pi_in.begin(),
				DualSmoothing::pi_in.end(),
				pi_lp1.begin(),
				smooth_dual.begin(),
				[](double a, double b) { return DualSmoothing::alpha * a + (1 - DualSmoothing::alpha) * b; });
	  double smooth_point_val = std::inner_product(lp2_obj.begin(),
												   lp2_obj.end(),
												   smooth_dual.begin(),
												   0.0,
												   plus<>(),
												   [](double a, double b) { return a * b; });
	  DualSmoothing::pi_out = pi_lp1;
	  if (smooth_point_val > lp2_value) {
		cout << "smooth_point_val = " << smooth_point_val << " lp2_value = " << lp2_value << endl;
	  } else {
		double norm_smooth_point = std::sqrt(std::accumulate(smooth_dual.begin(),
															 smooth_dual.end(),
															 0.0,
															 [](double a, double b) { return a + pow(b, 2); }));
		double gamma_ratio = (lp2_value - smooth_point_val) / (norm_smooth_point * lp2_norm);
		cout << "gamma_ratio= " << gamma_ratio << endl;
		if (gamma_ratio > gamma) {
		  if (DualSmoothing::alpha == 0 || DualSmoothing::if_mis_pricing) {
			DualSmoothing::pi_out = pi_lp2;
			cout << "use lp2" << endl;
			goto OUT;
		  }
		  double diff_lp_in = std::inner_product(lp2_obj.begin(),
												 lp2_obj.end(),
												 DualSmoothing::pi_in.begin(),
												 0.0,
												 plus<>(),
												 [](double a, double b) { return a * b; });
		  double diff_lp1 = std::inner_product(lp2_obj.begin(),
											   lp2_obj.end(),
											   pi_lp1.begin(),
											   0.0,
											   plus<>(),
											   [](double a, double b) { return a * b; });

		  double lambda = (lp2_value - diff_lp_in) / (diff_lp1 - diff_lp_in);
		  std::transform(DualSmoothing::pi_in.begin(),
						 DualSmoothing::pi_in.end(),
						 pi_lp1.begin(),
						 DualSmoothing::pi_price.begin(),
						 [lambda](double a, double b) { return a + lambda * (b - a); });
		  cout << "use combine with lp2" << endl;
		  goto HERE;
		} else {
		  cout << "use combine with lp1" << endl;
		}
	  }
	} else {
	  cout << "use lp2 in combine fashion" << endl;
	  DualSmoothing::pi_out = pi_lp2;
	}
  } else {
	DualSmoothing::updatePiOut();
  }
OUT:
  DualSmoothing::updatePiPrice();
HERE:
  DualSmoothing::calculateMostNegativeReducedCostInLP();
  DualSmoothing::recordLPValue();
}

void calculateNormDistance(const std::vector<double> &vec1, const std::vector<double> &vec2, double &norm) {
  norm = std::sqrt(std::inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.0,
									  std::plus<>(),
									  [](double a, double b) {
										return pow(a - b, 2);
									  }));
}