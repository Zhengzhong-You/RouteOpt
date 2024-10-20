#include "cvrp.hpp"
#include <new> // For std::bad_alloc

using namespace std;

int CVRP::checkMaximumNumberCustomers() const {
  if (dim >= MAX_NUM_CUSTOMERS) {
	cerr << "The dim is " << dim << " the solver cannot handle Customers over " << MAX_NUM_CUSTOMERS << endl;
	return -1;
  }
  return 0;
}

int CVRP::checkIfAssignGreaterMAXINT(const size_t &value, const string &str) {
  if (value >= size_t(numeric_limits<int>::max())) {
	cerr << str << " value= " << value << " but MAX_INT= " << numeric_limits<int>::max() << endl;
	return -4;
  }
  return 0;
}

void CVRP::reallocateLabel() {
  try {
	delete[] all_label;
	delete[] label_int_space;

	label_assign *= 2;
	safe_Hyperparameter(checkIfAssignGreaterMAXINT(label_assign, "label_assign"))
	all_label = new Label[label_assign];
	assignMemory4Label();
	label_int_space_len = MAX_NUM_R1CS_IN_PRICING + MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX;
	label_int_space = new int[label_int_space_len * label_assign];

	size_t cut_map_beg = 0;
	for (size_t i = 0; i < label_assign; ++i) {
	  auto &r1c = all_label[i].r1c;
	  r1c.cut_map = label_int_space + cut_map_beg;
	  cut_map_beg += MAX_NUM_R1CS_IN_PRICING;
	  r1c.valid_cut_idx = label_int_space + cut_map_beg;
	  cut_map_beg += MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX;
	}
  } catch (const std::bad_alloc &) {
	cerr << "Memory allocation failed during reallocateLabel. Not enough memory." << endl;
	exit(EXIT_FAILURE);
  }
}

int CVRP::checkPricingPool() const {
  if (pool_beg4_pricing >= size_t(Config::FracMemTolerance * double(mem4_pricing))) {
	return true;
  }
  return false;
}

void CVRP::reallocatePricingPool(size_t num) {
  try {
	if (num == 0) {
	  route_in_pricing_assign *= 2;
	} else {
	  if (num < size_t((double)route_in_pricing_assign * aver_route_length)) return;
	  route_in_pricing_assign = size_t(2 * (double)num / aver_route_length);
	}
	cout << "bad behavior! PricingPool is reallocated!" << endl;
	safe_Hyperparameter(checkIfAssignGreaterMAXINT(route_in_pricing_assign, "route_in_pricing_assign"))
	mem4_pricing = size_t(aver_route_length * (double)route_in_pricing_assign);
	auto tmp = col_pool4_pricing;
	col_pool4_pricing = new int[mem4_pricing];
	copy(tmp, tmp + pool_beg4_pricing, col_pool4_pricing);
	delete[] tmp;
  } catch (const std::bad_alloc &) {
	cerr << "Memory allocation failed during reallocatePricingPool. Not enough memory." << endl;
	exit(EXIT_FAILURE);
  }
}
