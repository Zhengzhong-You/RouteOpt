//
// Created by You, Zhengzhong on 7/10/24.
//

#include "read_dual_sol.hpp"

double ReadDualSol::coeff = 0.8;//initialize the coeff to 0.8
std::vector<double> ReadDualSol::dual_sol;

using namespace std;

void ReadDualSol::resetCoeff() {
  coeff = 0.8;
}

void ReadDualSol::changeRCStd(BbNode *node, double &rc_std, const std::vector<double> &dual) {
  auto &model = node->getSolver().model;
/**
 * step 1, get the number of non-zeros of the model (gurobi)
 * step 2, resize eigen triplets by the number of non-zeros
 * step 3, the lp matrix by gurobi model
 * step 4, get the objective value of the model
 * step 5, calculate rc by cost - dual * matrix (sparse matrix)
 * step 6, find the one with minimum rc to be rc_std
 */

  int num_constrs = (int)dual.size();
  int num_non_zeros;
  node->getSolver().getConstraints(&num_non_zeros, nullptr, nullptr, nullptr, 0, num_constrs);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(num_non_zeros);

  int num_vars;
  node->getSolver().getNumCol(&num_vars);

  Eigen::SparseMatrix<double> A(num_constrs, num_vars);
  vector<int> cbeg(num_constrs + 1);
  vector<int> cind(num_non_zeros);
  vector<double> cval(num_non_zeros);
  int cnt = 0;
  node->getSolver().getConstraints(&num_non_zeros, cbeg.data(), cind.data(), cval.data(), 0, num_constrs);
  cbeg[num_constrs] = num_non_zeros;
  for (int i = 0; i < num_constrs; ++i) {
	for (int j = cbeg[i]; j < cbeg[i + 1]; ++j) {
	  triplets[cnt++] = Eigen::Triplet<double>(i, cind[j], cval[j]);
	}
  }
  A.setFromTriplets(triplets.begin(), triplets.end());

  vector<double> obj_coeffs(num_vars);
  node->getSolver().getObj(0, num_vars, obj_coeffs.data());

  Eigen::RowVectorXd dual_eigen = Eigen::Map<const Eigen::RowVectorXd>(dual.data(), (int)dual.size());
  Eigen::RowVectorXd
	  rc = Eigen::Map<const Eigen::RowVectorXd>(obj_coeffs.data(), num_vars) - dual_eigen * A;
  rc_std = min(rc.minCoeff(), -0.1);
}

void ReadDualSol::checkIfMisPrice(int &cnt) {
  if (cnt == 0) {
	coeff /= 2;
	cout << "adjust the coeff to " << coeff << endl;
	if (coeff < 0.1) {
	  coeff = 0;
	  cout << "The coeff is less than 0.1, stop the dual combination!" << endl;
	} else {
	  cnt = 1;//tag 1 represents the misprice;
	}
  }
}

void ReadDualSol::changeDualSol(std::vector<double> &dual) {
/**
 * change the dual solution by multiplying the coeff;
 */
  int size = (int)dual.size();
  if (size != dual_sol.size()) {
	throw runtime_error("The size of dual solution is not equal to the size of dual vector!");
  }
  for (int i = 0; i < size; ++i) {
	dual[i] = dual[i] * (1 - coeff) + dual_sol[i] * coeff;
  }
}

void ReadDualSol::readDualSol(const std::string &file_name) {
/**
 * read the file, named file_name_dual.txt;
 * and find the the last 3 lines,
 * there is a line starts with "pi:" like,
 * pi: 221 458.314 292.103,
 * and read 221 458.314 292.103 and so on, till the next line;
 * this is the dual solution;
 */

  std::ifstream file("dual/" + file_name + "_dual.txt");
  if (!file.is_open()) {
	throw runtime_error("Cannot open file " + file_name + "_dual.txt");
  }

  std::string line;
  std::vector<std::string> last_three_lines;

  while (std::getline(file, line)) {
	last_three_lines.push_back(line);
	if (last_three_lines.size() > 3) {
	  last_three_lines.erase(last_three_lines.begin());
	}
  }

  file.close();

  for (const auto &l : last_three_lines) {
	if (l.rfind("pi:", 0) == 0) {  // line starts with "pi:"
	  std::istringstream iss(l.substr(3));  // skip "pi:"
	  double value;
	  while (iss >> value) {
		dual_sol.emplace_back(value);
	  }
	  break;  // exit after finding the first "pi:" line
	}
  }
  for (const auto &d : dual_sol) {
	std::cout << d << " ";
  }
  std::cout << std::endl;
}
