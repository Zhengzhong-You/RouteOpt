//
// Created by Zhengzhong You on 11/28/22.
//

#include "VRPTW.hpp"

#ifdef SOLVER_VRPTW
using namespace std;
VRPTW::VRPTW(const InstanceData &instanceData) : CVRP(instanceData) {
  //specialize other parameters used in VRPTW
  cout << "This is the VRPTW solver" << endl;
  transformed_number = 1;
}

void VRPTW::setResourceInBucketGraph() {
  //initial DP transit functions
  MaxMainResource = InfoVertex[0][5];//due date depot time
  MeetPointResourceInBiDir = MaxMainResource / 2;
  MainResourceAcrossArcsInForwardSense.resize(Dim);
  for (auto &vertex : MainResourceAcrossArcsInForwardSense) vertex.resize(Dim);
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
      MainResourceAcrossArcsInForwardSense[i][j] = InfoVertex[i][6] + CostMat4Vertex[i][j];
    }
  }

  MainResourceAcrossArcsInBackwardSense.resize(Dim);
  for (auto &vertex : MainResourceAcrossArcsInBackwardSense) vertex.resize(Dim);
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
      MainResourceAcrossArcsInBackwardSense[i][j] = InfoVertex[j][6] + CostMat4Vertex[i][j];
    }
  }
  lb4Vertex.resize(Dim);
  for (int i = 0; i < Dim; ++i) {
    lb4Vertex[i] = InfoVertex[i][4];
  }
  ub4Vertex.resize(Dim);
  for (int i = 0; i < Dim; ++i) {
    ub4Vertex[i] = InfoVertex[i][5];
  }
}

void VRPTW::specialize_MaxLengthEleRoute() {
  int max_len_ele_route;
  //for demand
  vector<double> demand(RealDim);
  for (int i = 1; i < Dim; ++i) {
    demand[i - 1] = InfoVertex[i][3];
  }
  std::stable_sort(demand.begin(), demand.end());
  double acc_cap = 0;
  max_len_ele_route = RealDim + 3;
  for (int i = 0; i < RealDim; ++i) {
    acc_cap += demand[i];
    if (acc_cap > Cap - TOLERANCE) {
      max_len_ele_route = i + 3;
      break;
    }
  }
  MaxLengthEleRoute = max(MaxLengthEleRoute, max_len_ele_route);
  cout << "MaxLengthEleRoute in cap context=  " << max_len_ele_route << endl;
  //for time
  vector<double> time;
  time.reserve(RealDim * RealDim);
  for (int i = 0; i < Dim; ++i) {
    for (int j = 0; j < Dim; ++j) {
      if (i != j)
        time.emplace_back(MainResourceAcrossArcsInForwardSense[i][j]);
    }
  }
  std::stable_sort(time.begin(), time.end());
  double acc_time = 0;
  max_len_ele_route = RealDim + 3;
  for (int i = 0; i < time.size(); ++i) {
    acc_time += time[i];
    if (acc_time >= MaxMainResource) {
      max_len_ele_route = i + 3;
      break;
    }
  }
  MaxLengthEleRoute = min(MaxLengthEleRoute, max_len_ele_route);
  cout << "MaxLengthEleRoute in time context=  " << max_len_ele_route << endl;
}

double VRPTW::transformCost(double x) {
  return (floor(10 * x + ten_FOLAT_TOLERANCE)) / 10;
}

void VRPTW::getLowerBoundofMinimumNumCars() {
  //compute by max clique
  vector<pair<int, int>> edges;
  for (int i = 1; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      if (lb4Vertex[i] + MainResourceAcrossArcsInForwardSense[i][j] > ub4Vertex[j]
          && lb4Vertex[j] + MainResourceAcrossArcsInForwardSense[j][i] > ub4Vertex[i]) {
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
        if (tag == 2 || Demand[i] + Demand[j] > Cap) {
          edges.emplace_back(i, j);
        }
      }
    }
  }
  if (!edges.empty()) {
    //construct incompatibilities graph
    vector<vector<int>> incompatibilities_graph(Dim, vector<int>(Dim));
    for (auto &e : edges) {
      incompatibilities_graph[e.first].emplace_back(e.second);
      incompatibilities_graph[e.second].emplace_back(e.first);
    }
  }
//  cout << "edge= " << endl;
//  for (auto &e : edges) {
//    cout << "i= " << e.first << " j=" << e.second << endl;
//  }
  if (edges.empty()) {
    cout << "edges is empty!" << endl;
  } else {
    cout << "edges is not empty!" << endl;
  }

  ///TODO:maximum clique 还没有用上！
  int max_clique_k = 0;
  cout << "max_clique_k= " << max_clique_k << endl;

  //compute by capacity
  double sum_demand = accumulate(Demand + 1, Demand + Dim, 0.0);
  int cap_k = ceil(sum_demand / Cap);

  cout << "cap_k= " << cap_k << endl;
//  cout << "LBNumVehicle= " << K << endl;

///TODO: 反向的那部分还没有去实现！
//construct new graph
  int m = 2;
  vector<vector<double>> cost;
  cout << "WARNING: No check time window resource!" << endl;
  goto NEXT;
  cost.resize(Dim, vector<double>(Dim));
  for (int i = 0; i < Dim; ++i) {// can be achieved at 0
    for (int j = 0; j < Dim; ++j) {
      cost[i][j] = max(CostMat4Vertex[i][j], lb4Vertex[j] - ub4Vertex[i] - InfoVertex[i][6]);
      if (lb4Vertex[i] + InfoVertex[i][6] + cost[i][j] > ub4Vertex[j] || i == j) {
        cost[i][j] = numeric_limits<double>::max();
      }
    }
  }
  m = max(cap_k, max_clique_k);
  for (;; ++m) {
    if (checkFeasibilityByEnergeticReasoning(m, cost)) break;
  }
  NEXT:
  K = m;
  cout << "LBNumVehicle= " << K << endl;
}

bool VRPTW::checkFeasibilityByEnergeticReasoning(int m, const std::vector<std::vector<double>> &cost) {
  int NumVertex = Dim + 2 * m;
  int beg_va = Dim + m;
  vector<vector<double>> new_cost(NumVertex, vector<double>(NumVertex));
  for (int i = 1; i < NumVertex; ++i) {
    for (int j = 1; j < NumVertex; ++j) {
      if (i < Dim) {
        if (j < Dim) {
          new_cost[i][j] = cost[i][j];
        } else if (j >= beg_va && j < NumVertex) {
          new_cost[i][j] = cost[i][0];
        } else if (j >= Dim && j < beg_va) {
          new_cost[i][j] = numeric_limits<double>::max();
        }
      } else if (i >= Dim && i < beg_va) {
        if (j < Dim) {
          new_cost[i][j] = cost[0][j];
        }
      } else if (i >= beg_va && i < NumVertex) {
        new_cost[i][j] = numeric_limits<double>::max();
      }
    }
  }
  vector<vrptwVertexInfo> vrptwInfoVertex(NumVertex);

  for (int i = 1; i < Dim; ++i) {
    auto min_it = min_element(new_cost[i].begin() + 1, new_cost[i].end());
    vrptwInfoVertex[i] = {lb4Vertex[i], ub4Vertex[i], *min_it};
  }

  vector<double> m_smallest(Dim - 1);
  for (int i = 1; i < Dim; ++i) {
    m_smallest[i - 1] = cost[0][i];//yes cost
  }
  nth_element(m_smallest.begin(), m_smallest.begin() + m, m_smallest.end());

  for (int i = Dim; i < beg_va; ++i) {
    vrptwInfoVertex[i] = {0, 0, m_smallest[i - Dim]};
  }
  for (int i = beg_va; i < NumVertex; ++i) {
    vrptwInfoVertex[i] = {0, MaxMainResource, 0};
  }

  //(18)
  for (int i = 1; i < Dim; ++i) {
    for (int j = 1; j < NumVertex; ++j) {
      auto min_it = min_element(new_cost[i].begin() + 1, new_cost[i].end());
      new_cost[i][j] -= *min_it;
    }
  }
  for (int i = Dim; i < beg_va; ++i) {
    for (int j = 1; j < NumVertex; ++j) {
      new_cost[i][j] = max(0.0, new_cost[i][j] - m_smallest[i - Dim]);
    }
  }
  //(19)-(21)
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

//get data
//only one side, not really useful
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

VRPTW::~VRPTW() = default;
#endif
