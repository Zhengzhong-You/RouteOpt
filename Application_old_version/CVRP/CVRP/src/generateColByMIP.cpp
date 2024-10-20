//
// Created by You, Zhengzhong on 11/5/23.
//

#include "CVRP.hpp"
#include "techniqueControl.hpp"

#ifdef GENERATE_COL_BY_MIP
using namespace std;
void CVRP::buildMIPModel(BbNode *node) {

  node->solver_mip.getEnv(&solver);
  const char *model_name = "BC.lp";
  num_mip_edge = dim * real_dim / 2;
  num_mip_var = num_mip_edge + real_dim;//duplicate ones for y variables

  vector<double> obj(num_mip_var, 0);
  int cnt = 0;
  for (int i = 0; i < dim; ++i) {
	for (int j = i + 1; j < dim; ++j) {
	  map_mip_idx_edge[cnt] = {i, j};
	  map_edge_mip_idx[{i, j}] = cnt;
	  cnt++;
	}
  }
  vector<double> lb_local(num_mip_var, 0);
  vector<double> ub_local(num_mip_var, 1);
  for (int j = 1; j < dim; ++j) ub_local[map_edge_mip_idx[{0, j}]] = 2;
  for (int i = cnt; i < num_mip_var; ++i) ub_local[i] = SOLVER_INFINITY;
  vector<char> x_type(num_mip_var, GRB_INTEGER);
  safe_solver(node->solver_mip.newModel(model_name, num_mip_var, obj.data(), lb_local.data(),
										ub_local.data(), x_type.data(), nullptr))

  vector<int> vbeg(dim + 1);
  vector<char> sense(dim, GRB_EQUAL);
  vector<double> rhs(dim, 0);
  rhs[0] = 2;
  vector<int> vind;
  vector<double> vval;
  vind.reserve(dim);
  vval.reserve(dim);

  // i=0
  vbeg[0] = 0;
  for (int j = 1; j < dim; ++j) {
	vind.emplace_back(map_edge_mip_idx[{0, j}]);
	vval.emplace_back(1);
  }

  //i!=0
  for (int i = 1; i < dim; ++i) {
	vbeg[i] = (int)vind.size();
	for (int j = 0; j < dim; ++j) {
	  if (i == j) continue;
	  if (i < j) vind.emplace_back(map_edge_mip_idx[{i, j}]);
	  else vind.emplace_back(map_edge_mip_idx[{j, i}]);
	  vval.emplace_back(1);
	}
	vind.emplace_back(num_mip_edge + i - 1);
	vval.emplace_back(-2);
  }
  vbeg[dim] = (int)vind.size();
  safe_solver(node->solver_mip.addConstraints(dim, vind.size(), vbeg.data(), vind.data(), vval.data(), sense.data(),
											  rhs.data(), nullptr))
  safe_solver(node->solver_mip.updateModel())
  safe_solver(node->solver_mip.getNumRow(&num_mip_constr))
  //have not checked this model!
}

struct callback_data {
  CVRP *cvrp;
};

//int __stdcall
//CVRP::mycallback(GRBmodel *model,
//				 void *cbdata,
//				 int where,
//				 void *usrdata) {
////  auto data = (struct callback_data *)usrdata;
////  CVRP *cvrp = data->cvrp;
//
//  if (where == GRB_CB_MIPSOL) {
//	vector<double> sol;
//	sol.resize(num_mip_var);
//	GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, sol.data());//only consider one route this time!
//	vector<int> tour;
//
//	//ng-sub-tour elimination
//	findNgSubTour(sol, tour);
//  }
//  return 0;
//}
//
//void CVRP::findNgSubTour(void *cbdata, const std::vector<std::vector<int>> &tours) {
//  if (tours.empty()) return;
//  for (auto &tour : tours) {
//	vector<int> vind;
//	vector<double> vval;
//	for (int i = 0; i < tour.size(); ++i) {
//	  for (int j = i + 1; j < tour.size(); ++j) {
//		tour[i] < tour[j] ? vind.emplace_back(Idx4Var[{tour[i], tour[j]}]) : vind.emplace_back(Idx4Var[{tour[j],
//																										tour[i]}]);
//		vval.emplace_back(1);
//	  }
//	}
//	GRBcblazy(cbdata, (int)vind.size(), vind.data(), vval.data(), GRB_LESS_EQUAL, (int)tour.size() - 1);
//  }
//}

int __stdcall
mycallback(GRBmodel *model,
		   void *cbdata,
		   int where,
		   void *usrdata) {
  return 0;
}

void CVRP::setObjCoeffs(BbNode *node) {
  vector<double> obj(num_mip_edge);
  for (int i = 0; i < dim; ++i) {
	for (int j = i + 1; j < dim; ++j) {
	  obj[map_edge_mip_idx[{i, j}]] = chg_cost_mat4_vertex[i][j];
	}
  }
  safe_solver(node->solver_mip.changeObj(0, num_mip_edge, obj.data()))
}

void CVRP::cbSolve(BbNode *node) {
  double cut_off;
  int lazy_constr;
  safe_solver(node->solver_mip.getEnvCutoff(&cut_off))
  safe_solver(node->solver_mip.setEnvCutoff(-TOLERANCE))

  safe_solver(node->solver_mip.getEnvLazyConstrStat(&lazy_constr))
  safe_solver(node->solver_mip.setEnvLazyConstrStat(1))

  safe_solver(node->solver_mip.setcallbackfunc(mycallback, nullptr))

  setObjCoeffs(node);

  auto beg = chrono::high_resolution_clock::now();
  safe_solver(node->solver_mip.optimize())

  int sol_count;

  safe_solver(node->solver_mip.getSolCount(&sol_count))

  if (sol_count) {
	vector<double> sol;
	sol.resize(num_mip_edge);
	safe_solver(node->solver_mip.getX(0, num_mip_edge, sol.data()))
	for (int i = 0; i < num_mip_edge; ++i) {
	  if (sol[i] > 0.5) {
		cout << map_mip_idx_edge[i].first << " " << map_mip_idx_edge[i].second << ": " << sol[i] << endl;
	  }
	}
	exit(0);
	vector<vector<int>> tours;
  }

  auto end = chrono::high_resolution_clock::now();
  auto eps = chrono::duration<double>(end - beg).count();
  cout << "solve time: " << eps << endl;
  safe_solver(node->solver_mip.setEnvCutoff(cut_off))
  safe_solver(node->solver_mip.setEnvLazyConstrStat(lazy_constr))
}

#endif