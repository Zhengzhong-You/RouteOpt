
#include "CVRP.hpp"

using namespace std;

int CVRP::checkMaximumNumberCustomers() const {
  if (dim >= MAX_NUM_CUSTOMERS) {
	cerr << "The dim is " << dim << " the solver cannot handle Customers over " << MAX_NUM_CUSTOMERS << endl;
	return -1;
  }
  return 0;
}

int CVRP::checkMaximumNumberR1Cs(int num_r1cs) {
  if (num_r1cs >= MAX_NUM_R1CS) {
	cerr << "num_r1cs= " << num_r1cs << " but MaxNum_SR3= " << MAX_NUM_R1CS << endl;
	return -2;
  }
  return 0;
}

int CVRP::checkCSTLimit() const {
  if (num_row >= CST_LIMIT) {
	cerr << "num_row= " << num_row << " but CST_LIMIT= " << CST_LIMIT << endl;
	return -3;
  }
  return 0;
}

int CVRP::checkIfAssignGreaterMAXINT(const size_t &value, const string &str) {
  if (value >= size_t(MAX_INT)) {
	cerr << str << " value= " << value << " but MAX_INT= " << MAX_INT << endl;
	return -4;
  }
  return 0;
}

void CVRP::reallocateLabel() {
  cout << "bad behavior! because labels are reallocated! and labeling will be rolled back!" << endl;
  delete[]all_label;
  delete[]all_valid_r1cs;
  delete[]all_valid_r1c_multi;
  delete[]all_r1c_multi_mem;

  label_assign *= 2;
  safe_Hyperparameter(checkIfAssignGreaterMAXINT(label_assign, "label_assign"))
  all_label = new Label[label_assign];
  all_valid_r1cs = new int[label_assign * MAX_NUM_R1CS];
  all_valid_r1c_multi = new int[(label_assign) * MAX_NUM_R1C_MULTI];
  all_r1c_multi_mem = new int[(label_assign) * MAX_NUM_R1C_MULTI];
  for (size_t i = 0; i < label_assign; ++i) {
	all_label[i].valid_rank1_cut = all_valid_r1cs + i * MAX_NUM_R1CS;
	all_label[i].valid_rank1_cut_multi = all_valid_r1c_multi + i * MAX_NUM_R1C_MULTI;
	all_label[i].rank1_cut_mem_multi = all_r1c_multi_mem + i * MAX_NUM_R1C_MULTI;
  }
}

int CVRP::checkMemoryPool() const {
  if (pool_beg4_mem >= size_t(Config::FracMemTolerance * double(mem4_mem))) {
	cout << "reallocate!" << endl;
	return true;
  }
  return false;
}

void CVRP::reallocateMemoryPool() {
  cout << "bad behavior! because MemPool are reallocated!" << endl;
  route_in_mem_assign *= 2;
  safe_Hyperparameter(checkIfAssignGreaterMAXINT(route_in_mem_assign, "route_in_mem_assign"))
  mem4_mem = aver_route_length * route_in_mem_assign;
  auto tmp = col_pool4_mem;
  col_pool4_mem = new int[mem4_mem];
  copy(tmp, tmp + pool_beg4_mem, col_pool4_mem);
  delete[]tmp;
}

int CVRP::checkPricingPool() const {
  if (pool_beg4_pricing >= size_t(Config::FracMemTolerance * double(mem4_pricing))) {
	return true;
  }
  return false;
}

void CVRP::reallocatePricingPool() {
  cout << "bad behavior! because PricingPool are reallocated!" << endl;
  route_in_pricing_assign *= 2;
  safe_Hyperparameter(checkIfAssignGreaterMAXINT(route_in_pricing_assign, "route_in_pricing_assign"))
  mem4_pricing = aver_route_length * route_in_pricing_assign;
  auto tmp = col_pool4_pricing;
  col_pool4_pricing = new int[mem4_pricing];
  copy(tmp, tmp + pool_beg4_pricing, col_pool4_pricing);
  delete[]tmp;
}






