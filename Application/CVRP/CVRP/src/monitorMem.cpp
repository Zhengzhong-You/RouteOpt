//
// Created by Zhengzhong You on 7/14/22.
//

#include "CVRP.hpp"

using namespace std;

int CVRP::checkMaxNum_Customers() const {
  if (Dim >= MaxNum_Customers) {
    cerr << "The Dim is " << Dim << " the solver cannot handle Customers over " << MaxNum_Customers << endl;
    return -1;
  }
  return 0;
}

int CVRP::checkMaxNum_R1Cs(int num_r1cs) {
  if (num_r1cs >= MaxNum_R1Cs) {
    cerr << "num_r1cs= " << num_r1cs << " but MaxNum_SR3= " << MaxNum_R1Cs << endl;
    return -2;
  }
  return 0;
}

int CVRP::checkCST_LIMIT() const {
  if (NumRow >= CST_LIMIT) {
    cerr << "NumRow= " << NumRow << " but CST_LIMIT= " << CST_LIMIT << endl;
    return -3;
  }
  return 0;
}

int CVRP::checkAssignGreaterMAXINT(const size_t &value, const string &str) {
  if (value >= size_t(MaxInt)) {
    cerr << str << " value= " << value << " but MaxInt= " << MaxInt << endl;
    return -4;
  }
  return 0;
}

void CVRP::reallocateLabel() {
  cout << "bad behavior! because labels are reallocated! and labeling will be rolled back!" << endl;
  delete[]AllLabel;
  delete[]AllSeq;
  delete[]AllValidR1Cs;

  LabelAssign *= 2;
  safe_Hyperparameter(checkAssignGreaterMAXINT(LabelAssign, "LabelAssign"))
  AllLabel = new LABEL[LabelAssign];
  AllSeq = new int[LabelAssign * MaxLengthEleRoute];
  AllValidR1Cs = new int[LabelAssign * MaxNum_R1Cs];
  for (size_t i = 0; i < LabelAssign; ++i) {
    AllLabel[i].validRank1Cut = AllValidR1Cs + i * MaxNum_R1Cs;
  }
}

int CVRP::checkMemPool() const {
  if (PoolBeg4Mem >= size_t(CONFIG::FracMemTolerance * double(Mem4Mem))) {
    cout << "reallocate!" << endl;
    return true;
  }
  return false;
}

void CVRP::reallocateMemPool() {
  cout << "bad behavior! because MemPool are reallocated!" << endl;
  RouteInMemAssign *= 2;
  safe_Hyperparameter(checkAssignGreaterMAXINT(RouteInMemAssign, "RouteInMemAssign"))
  Mem4Mem = AverRouteLength * RouteInMemAssign;
  auto tmp = ColPool4Mem;
  ColPool4Mem = new int[Mem4Mem];
  copy(tmp, tmp + PoolBeg4Mem, ColPool4Mem);
  delete[]tmp;
}

int CVRP::checkPricingPool() const {
  if (PoolBeg4Pricing >= size_t(CONFIG::FracMemTolerance * double(Mem4Pricing))) {
    return true;
  }
  return false;
}

void CVRP::reallocatePricingPool() {
  cout << "bad behavior! because PricingPool are reallocated!" << endl;
  RouteInPricingAssign *= 2;
  safe_Hyperparameter(checkAssignGreaterMAXINT(RouteInPricingAssign, "RouteInPricingAssign"))
  Mem4Pricing = AverRouteLength * RouteInPricingAssign;
  auto tmp = ColPool4Pricing;
  ColPool4Pricing = new int[Mem4Pricing];
  copy(tmp, tmp + PoolBeg4Pricing, ColPool4Pricing);
  delete[]tmp;
}

void CVRP::reallocateSolverPtr(size_t num) {
  cout << "reallocateSolverPtr!" << endl;
  MaxNonZeroEntry = num * 2;
//  cout << "MaxNonZeroEntry= " << MaxNonZeroEntry << endl;
//  cout << "num= " << num << endl;
  delete[]solver_ind;
  delete[]solver_ind2;
  delete[]solver_val;
  solver_ind = new int[MaxNonZeroEntry];
  solver_ind2 = new int[MaxNonZeroEntry];
  solver_val = new double[MaxNonZeroEntry];
}

int CVRP::checkSolverPtr(size_t num) const {
  if (num >= size_t(CONFIG::FracMemTolerance * double(MaxNonZeroEntry))) {
//    cout << "num= " << num << endl;
//    cout << "MaxNonZeroEntry= " << MaxNonZeroEntry << endl;
//    cout << "CONFIG::FracMemTolerance= " << CONFIG::FracMemTolerance << endl;
//    cout << "size_t(CONFIG::FracMemTolerance * double(MaxNonZeroEntry))= "
//         << size_t(CONFIG::FracMemTolerance * double(MaxNonZeroEntry)) << endl;
    return true;
  }
  return false;
}






