
#include "VRPTW.hpp"

#ifdef SOLVER_VRPTW
using namespace std;
VRPTW::VRPTW(const InstanceData &instanceData) : CVRP(instanceData) {
  cout << "This is the VRPTW solver" << endl;
  transformed_number = 1;
}

void VRPTW::setResourceInBucketGraph() {
  max_main_resource = info_vertex[0][5];//due date depot time
  meet_point_resource_in_bi_dir = max_main_resource / 2;
  main_resource_across_arcs_in_forward_sense.resize(dim);
  for (auto &vertex : main_resource_across_arcs_in_forward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      main_resource_across_arcs_in_forward_sense[i][j] = info_vertex[i][6] + cost_mat4_vertex[i][j];
    }
  }

  main_resource_across_arcs_in_backward_sense.resize(dim);
  for (auto &vertex : main_resource_across_arcs_in_backward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
	  main_resource_across_arcs_in_backward_sense[i][j] = info_vertex[j][6] + cost_mat4_vertex[i][j];
    }
  }
  lb4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
    lb4_vertex[i] = info_vertex[i][4];
  }
  ub4_vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
    ub4_vertex[i] = info_vertex[i][5];
  }
}

double VRPTW::transformCost(double x) {
  return (floor(10 * (x+ TOLERANCE))) / 10;
}

void VRPTW::getLowerBoundofMinimumNumberCars() {
  vector<pair<int, int>> edges;
  for (int i = 1; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      if (lb4_vertex[i] + main_resource_across_arcs_in_forward_sense[i][j] > ub4_vertex[j]
          && lb4_vertex[j] + main_resource_across_arcs_in_forward_sense[j][i] > ub4_vertex[i]) {
        edges.emplace_back(i, j);
      } else {
        double tmp_res = 0;
        int tag = 0;
        increaseMainResourceConsumption(tmp_res, tmp_res, 0, i);
        if (!increaseMainResourceConsumption(tmp_res, tmp_res, i, j)
            || !increaseMainResourceConsumption(tmp_res, tmp_res, j, 0)) {
          ++tag;
        }
        if (tag) {
          tmp_res = 0;
          increaseMainResourceConsumption(tmp_res, tmp_res, 0, j);
          if (!increaseMainResourceConsumption(tmp_res, tmp_res, j, i) ||
              !increaseMainResourceConsumption(tmp_res, tmp_res, i, 0)) {
            ++tag;
          }
        }
        if (tag == 2 || demand[i] + demand[j] > cap) {
          edges.emplace_back(i, j);
        }
      }
    }
  }
  if (!edges.empty()) {
    vector<vector<int>> incompatibilities_graph(dim, vector<int>(dim));
    for (auto &e : edges) {
      incompatibilities_graph[e.first].emplace_back(e.second);
      incompatibilities_graph[e.second].emplace_back(e.first);
    }
  }
  if (edges.empty()) {
    cout << "edges is empty!" << endl;
  } else {
    cout << "edges is not empty!" << endl;
  }

  int max_clique_k = 0;
  cout << "max_clique_k= " << max_clique_k << endl;

  double sum_demand = accumulate(demand + 1, demand + dim, 0.0);
  int cap_k = ceil(sum_demand / cap);

  cout << "cap_k= " << cap_k << endl;

  int m = 2;
  vector<vector<double>> cost;
  cout << "WARNING: No check time window resource!" << endl;
  goto NEXT;
  cost.resize(dim, vector<double>(dim));
  for (int i = 0; i < dim; ++i) {// can be achieved at 0
    for (int j = 0; j < dim; ++j) {
      cost[i][j] = max(cost_mat4_vertex[i][j], lb4_vertex[j] - ub4_vertex[i] - info_vertex[i][6]);
      if (lb4_vertex[i] + info_vertex[i][6] + cost[i][j] > ub4_vertex[j] || i == j) {
        cost[i][j] = numeric_limits<double>::max();
      }
    }
  }
  m = max(cap_k, max_clique_k);
  for (;; ++m) {
    if (checkFeasibilityByEnergeticReasoning(m, cost)) break;
  }
  NEXT:
  num_vehicle = m;
  cout << "LBNumVehicle= " << num_vehicle << endl;
}

bool VRPTW::checkFeasibilityByEnergeticReasoning(int m, const std::vector<std::vector<double>> &cost) {
  int NumVertex = dim + 2 * m;
  int beg_va = dim + m;
  vector<vector<double>> new_cost(NumVertex, vector<double>(NumVertex));
  for (int i = 1; i < NumVertex; ++i) {
    for (int j = 1; j < NumVertex; ++j) {
      if (i < dim) {
        if (j < dim) {
          new_cost[i][j] = cost[i][j];
        } else if (j >= beg_va && j < NumVertex) {
          new_cost[i][j] = cost[i][0];
        } else if (j >= dim && j < beg_va) {
          new_cost[i][j] = numeric_limits<double>::max();
        }
      } else if (i >= dim && i < beg_va) {
        if (j < dim) {
          new_cost[i][j] = cost[0][j];
        }
      } else if (i >= beg_va && i < NumVertex) {
        new_cost[i][j] = numeric_limits<double>::max();
      }
    }
  }
  vector<vrptwVertexInfo> vrptwInfoVertex(NumVertex);

  for (int i = 1; i < dim; ++i) {
    auto min_it = min_element(new_cost[i].begin() + 1, new_cost[i].end());
    vrptwInfoVertex[i] = {lb4_vertex[i], ub4_vertex[i], *min_it};
  }

  vector<double> m_smallest(dim - 1);
  for (int i = 1; i < dim; ++i) {
    m_smallest[i - 1] = cost[0][i];//yes cost
  }
  nth_element(m_smallest.begin(), m_smallest.begin() + m, m_smallest.end());

  for (int i = dim; i < beg_va; ++i) {
    vrptwInfoVertex[i] = {0, 0, m_smallest[i - dim]};
  }
  for (int i = beg_va; i < NumVertex; ++i) {
    vrptwInfoVertex[i] = {0, max_main_resource, 0};
  }

  for (int i = 1; i < dim; ++i) {
    for (int j = 1; j < NumVertex; ++j) {
      auto min_it = min_element(new_cost[i].begin() + 1, new_cost[i].end());
      new_cost[i][j] -= *min_it;
    }
  }
  for (int i = dim; i < beg_va; ++i) {
    for (int j = 1; j < NumVertex; ++j) {
      new_cost[i][j] = max(0.0, new_cost[i][j] - m_smallest[i - dim]);
    }
  }
  for (int i = 1; i < NumVertex; ++i) {
    vector<double> tmp(NumVertex - 1);
    for (int j = 1; j < NumVertex; ++j) {
      tmp[j - 1] = new_cost[j][i];//from j to i
    }
    auto min_it = min_element(tmp.begin(), tmp.end());
    vrptwInfoVertex[i].s_time += *min_it;
    vrptwInfoVertex[i].e_time = max(0.0, vrptwInfoVertex[i].e_time - *min_it);
    vrptwInfoVertex[i].l_time = max(0.0, vrptwInfoVertex[i].l_time - *min_it);
  }

  set<double, double_same_Tolerance> T_1;
  set<double, double_same_Tolerance> T_2;
  for (int i = 1; i < NumVertex; ++i) {
    T_1.emplace(vrptwInfoVertex[i].e_time);
    T_1.emplace(vrptwInfoVertex[i].l_time);
    double e_s = vrptwInfoVertex[i].e_time + vrptwInfoVertex[i].s_time;
    T_1.emplace(e_s);
    T_2.emplace(e_s);
    T_2.emplace(vrptwInfoVertex[i].l_time + vrptwInfoVertex[i].s_time);
    T_2.emplace(vrptwInfoVertex[i].l_time);
  }
  vector<double> T1(T_1.begin(), T_1.end());
  vector<double> T2(T_2.begin(), T_2.end());

  for (auto t1 : T1) {
    for (auto t2 : T2) {
      if (t1 < t2 - TOLERANCE) {
        double w = 0;
        for (int i = 1; i < NumVertex; ++i) {
          w += calculateW(i, t1, t2, vrptwInfoVertex[i]);
        }
        if (w > m * (t2 - t1) + TOLERANCE) {
          return false;
        }
      }
    }
  }
  return true;
}

double VRPTW::calculateW(int i, double t1, double t2, const vrptwVertexInfo &vrptwInfo) {
  double w_left = min({t2 - t1, vrptwInfo.s_time, max(0.0, vrptwInfo.e_time + vrptwInfo.s_time - t1)});
  double w_right = min({t2 - t1, vrptwInfo.s_time, max(0.0, t2 - vrptwInfo.l_time)});
  return min(w_left, w_right);
}

void VRPTW::checkSolutionFeasibleByCapacity(bool &feasible) {
  feasible= true;
  for (auto &route : ip_opt_sol) {
	double local_cap = 0;
	for (int i = 0; i < route.size() - 1; ++i) {
	  local_cap += demand[route[i]];
	}
	if (local_cap > cap+TOLERANCE) {
	  cout<<"cap for route ";
	  for (auto i : route) {
		cout<<i<<" ";
	  }
	  cout<<" is "<<local_cap<<", while cap= "<<cap<<endl;
	  feasible=false;
	  break;
	}
  }
}

void VRPTW::cleanColumnsNonFeasible(BbNode *node) {
	int len = 0, keep = 1, current_node;
	bool if_break;
	yzzLong pi;
	vector<int> solver_ind(num_col);
	for (int i = keep; i < num_col; ++i) {
	  double local_cap=0;
	  for (size_t j = node->index_columns[i] + 1;; ++j) {
		current_node = col_pool4_pricing[j];
		if (current_node == 0) break;
		local_cap+=demand[current_node];
	  }
	  if(local_cap>cap+TOLERANCE){
		solver_ind[len++] = i;
	  }
	  else{
		node->index_columns[keep++] = node->index_columns[i];
	  }
	}

	if (len) {
	  safe_solver(node->solver.delVars(len, solver_ind.data()))
	  safe_solver(node->solver.reoptimize())
	  safe_solver(node->solver.getNumCol(&num_col))
	}

	safe_solver(node->solver.getObjVal(&lp_val))
	cout << "after clean non-feasible routes lpval= " << lp_val << endl;
}

void VRPTW::popArcGraph(BbNode* node){
  vector<double> X(num_col);
  safe_solver(node->solver.getX(0, num_col, X.data()))

  int parent_column =  node->num_parent_cols;
  vector<pair<size_t, double>> tmp_solindex;

  for (int i = 0; i < parent_column; ++i) {
	if (X[i] > TOLERANCE) {
	  tmp_solindex.emplace_back(node->index_columns[i], X[i]);
	}
  }

  int break_record = (int) tmp_solindex.size();

  for (int i = parent_column; i < num_col; ++i) {
	if (X[i] > TOLERANCE) {
	  tmp_solindex.emplace_back(node->index_columns[i], X[i]);
	}
  }

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      arc_graph[i][j] = 0;
    }
  }

  int beg = 0;
  int end = break_record;
  int *colpool=col_pool4_mem;

  AGAIN:
  for (int i = beg; i < end; ++i) {
    for (size_t j = tmp_solindex[i].first;;) {
      arc_graph[colpool[j]][colpool[j + 1]] += tmp_solindex[i].second;
      if (!colpool[++j]) break;
    }
  }

  if (end != tmp_solindex.size()) {
    beg = end;
    end = (int) tmp_solindex.size();
    colpool = col_pool4_pricing;
    goto AGAIN;
  }


  num_edge = 0;

  for (int i = 0; i < dim; ++i) {
    for (int j = i + 1; j < dim; ++j) {
      arc_graph[i][j] += arc_graph[j][i];
	  if (arc_graph[i][j] > TOLERANCE) {
		++num_edge;
	  }
    }
  }
}

VRPTW::~VRPTW() = default;
#endif
