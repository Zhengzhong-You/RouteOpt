
#include "vrptw.hpp"

using namespace std;
VRPTW::VRPTW(const InstanceData &instanceData) : CVRP(instanceData) {
  cout << "This is the VRPTW solver" << endl;
  cout << "Notice: arc-based rank1 memory may be activated." << endl;
  cout << "In asymmetric cases, arc memory fixes the meet point resource for pricing after adding rank1 cuts." << endl;
  cout
	  << "This approach enforces the same meet point resource for both enumeration and pricing, potentially degrading enumeration performance significantly."
	  << endl;
  transformed_number = 1;
}

void VRPTW::setResourceInBucketGraph() {
  auto &travel_mat = travel_time_matrix.empty() ? cost_mat4_vertex : travel_time_matrix;
#ifdef USE_TWO_RESOURCE
#ifdef CAPACITY_AS_MAIN_RESOURCE
  resource = {roundAndConvertResLong(cap), roundAndConvertResLong(info_vertex[0][5])};
  meet_point_resource_in_bi_dir = double(resource.first_res) / 2;
  resource_across_arcs_in_forward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_forward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_forward_sense[i][j] = {roundAndConvertResLong(info_vertex[i][3]),
													 roundAndConvertResLong(
														 info_vertex[i][6] + travel_mat[i][j])};
	}
  }
  resource_across_arcs_in_backward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_backward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_backward_sense[i][j] = {roundAndConvertResLong(info_vertex[j][3]),
													  roundAndConvertResLong(
														  info_vertex[j][6] + travel_mat[i][j])};
	}
  }
  lb4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	lb4_vertex[i] = {0, roundAndConvertResLong(info_vertex[i][4])};
  }
  ub4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	ub4_vertex[i] = {roundAndConvertResLong(cap), roundAndConvertResLong(info_vertex[i][5])};
  }
#else
  resource = {roundAndConvertResLong(info_vertex[0][5]), roundAndConvertResLong(cap)};
  meet_point_resource_in_bi_dir = double(resource.first_res) / 2;
  resource_across_arcs_in_forward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_forward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_forward_sense[i][j] = {roundAndConvertResLong(info_vertex[i][6] + travel_mat[i][j]),
													 roundAndConvertResLong(info_vertex[i][3])};
	}
  }

  resource_across_arcs_in_backward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_backward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_backward_sense[i][j] =
		  {roundAndConvertResLong(info_vertex[j][6] + travel_mat[i][j]),
		   roundAndConvertResLong(info_vertex[j][3])};
	}
  }
  lb4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	lb4_vertex[i] = {roundAndConvertResLong(info_vertex[i][4]), 0};
  }
  ub4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	ub4_vertex[i] = {roundAndConvertResLong(info_vertex[i][5]), roundAndConvertResLong(cap)};
  }
#endif
#else
  resource.first_res = roundAndConvertResLong(info_vertex[0][5]);
  meet_point_resource_in_bi_dir = double(resource.first_res) / 2;
  resource_across_arcs_in_forward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_forward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_forward_sense[i][j].first_res =
		  roundAndConvertResLong(info_vertex[i][6] + travel_mat[i][j]);
	}
  }

  resource_across_arcs_in_backward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_backward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_backward_sense[i][j].first_res =
		  roundAndConvertResLong(info_vertex[j][6] + travel_mat[i][j]);
	}
  }
  lb4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	lb4_vertex[i].first_res = roundAndConvertResLong(info_vertex[i][4]);
  }
  ub4_vertex.resize(dim);
  ub4_vertex[0].first_res = resource.first_res;
  for (int i = 1; i < dim; ++i) {
	ub4_vertex[i].first_res = roundAndConvertResLong(info_vertex[i][5]);
  }
#endif
}


double VRPTW::transformCost(double x) {
  return (floor(10 * (x + VRPTW_DISTANCE_TOLERANCE))) / 10;
}

void VRPTW::getLowerBoundofMinimumNumberCars() {
  double demand_sum = 0;
  for (int i = 1; i < dim; ++i) {
	demand_sum += info_vertex[i][3];
  }
  num_vehicle = (int)ceil(demand_sum / cap);
  cout << "LBNumVehicle= " << num_vehicle << endl;
}

void VRPTW::checkSolutionFeasibleByCapacity(bool &feasible) {
  feasible = true;
  for (auto &route : ip_opt_sol) {
	double local_cap = 0;
	for (auto &i : route) {
	  local_cap += demand[i];
	}
	if (local_cap > cap + TOLERANCE) {
	  cout << "cap for route ";
	  for (auto &i : route) {
		cout << i << " ";
	  }
	  cout << " is " << local_cap << ", while cap= " << cap << endl;
	  feasible = false;
	  break;
	}
  }
}

void VRPTW::cleanColumnsCapInfeasible(BbNode *node) {
  /**
   * only focus on the lp, the enumeration pool is considered elsewhere
   */
  double res;
  vector<int> col_idx;

  for (int i = 1; i < num_col; ++i) {
	res = 0;
	for (auto current_node : node->getCols()[i].col_seq) {
	  res += demand[current_node];
	}
	if (res > cap + TOLERANCE) {
	  col_idx.emplace_back(i);
	}
  }
  rmLPCols(node, col_idx);
  safe_solver(node->getSolver().reoptimize(SOLVER_PRIMAL_SIMPLEX))
  safe_solver(node->getSolver().getObjVal(&lp_val))
#if VERBOSE_MODE == 1
  cout << "after remove " << col_idx.size() << " res-cap infeasible columns, lpval= " << lp_val << endl;
#endif
}

void VRPTW::popArcGraph(BbNode *node) {
  vector<double> X(num_col);
  safe_solver(node->getSolver().getX(0, num_col, X.data()))

  /**
   * do not update the lp sol of the node
   */

  pair<vector<SequenceInfo>, vector<double>> local_all_sol;

  for (int i = 0; i < num_col; ++i) {
	if (X[i] > TOLERANCE) {
	  local_all_sol.second.emplace_back(X[i]);
	  local_all_sol.first.emplace_back(node->getCols()[i]);
	}
  }

  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  arc_graph[i][j] = 0;
	}
  }

  auto opt_size = (int)local_all_sol.first.size();
  for (int i = 0; i < opt_size; ++i) {
	auto &seq = local_all_sol.first[i].col_seq;
	double val = local_all_sol.second[i];
	int b4 = 0;
	for (auto j : seq) {
	  arc_graph[b4][j] += val;
	  b4 = j;
	}
	arc_graph[b4][0] += val;
  }

  num_edge = 0;

  /**
   * no need to revise arc_graph(since arc_graph_revised is only for branching)
   */
  for (int i = 0; i < dim; ++i) {
	for (int j = i + 1; j < dim; ++j) {
	  arc_graph[i][j] += arc_graph[j][i];
	  if (arc_graph[i][j] > TOLERANCE) {
		++num_edge;
	  }
	}
  }
}

void VRPTW::tryGetTravelTimeMatrix(const std::string &file_name) {
  ifstream file(file_name);
  if (!file.is_open()) {
	throw runtime_error(file_name + " cannot be opened!");
  }

  std::string line;
  bool foundMatrix = false;

  while (getline(file, line)) {
	if (line.find("Travel Time Matrix") != std::string::npos) {
	  foundMatrix = true;
	  break;
	}
  }

  if (!foundMatrix) {
	return;
  }

  travel_time_matrix.resize(dim, vector<double>(dim));

  for (size_t i = 0; i < dim; ++i) {
	getline(file, line);
	std::istringstream iss(line);
	for (size_t j = 0; j < dim; ++j) {
	  iss >> travel_time_matrix[i][j];
	}
  }

  file.close();
}

VRPTW::~VRPTW() = default;
