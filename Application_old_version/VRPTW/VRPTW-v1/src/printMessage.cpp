
#include "CVRP.hpp"
#include <filesystem>

using namespace std;
using namespace chrono;

#ifdef VERBOSE_MODE

void printSolution(std::ostream &os, const vector<vector<int>> &ip_opt_sol, double ub,
				   double glo_eps, int num_explored_nodes, double global_gap, const string &file_name, double cap) {
  os << "<Solution> " << endl;
  for (auto &r : ip_opt_sol) {
	os << "0";
	for (auto &i : r) {
	  os << "-" << i;
	}
	os << "-0" << endl;
  }
  os << "<Optval  " << ub << ">  <Etime  " << glo_eps
	 << ">  <Nodes_explored  " << num_explored_nodes << ">"
	 << "  <global_gap  " << global_gap * 100 << "%>" << endl;
  os << "<Instance  " << file_name << "  Capacity  " << cap << ">" << endl;
}

void CVRP::printOptIntSol() {
  glo_end = chrono::high_resolution_clock::now();
  glo_eps = duration<double>(glo_end - glo_beg).count();
#ifdef Resolve_Ins_with_Optimal_Val
  if (!Config::Marker4FindLargeGapIns) {
	if (ip_opt_sol.empty()) {
	  cout << "No solution found" << endl;
	  Config::Marker4FindLargeGapIns = -1;
	  return;
	}
	Config::ub = ub;
	return;
  }
#endif
  printSolution(cout, ip_opt_sol, ub, glo_eps, num_explored_nodes, global_gap, file_name, cap);

#ifdef WRITE_SOL_OUT_TO_FILE
  if (!std::filesystem::exists(WRITE_SOL_OUT_TO_FILE)) {
	std::filesystem::create_directories(WRITE_SOL_OUT_TO_FILE);
  }

  string sol_file_name = string(WRITE_SOL_OUT_TO_FILE) + "/" + file_name + ".sol";
  ofstream sol_file(sol_file_name);
  printSolution(sol_file, ip_opt_sol, ub, glo_eps, num_explored_nodes, global_gap, file_name, cap);
  sol_file.close();
#endif
#ifdef writeAccOut2File
  if (!std::experimental::filesystem::exists(writeAccOut2File)) {
	std::experimental::filesystem::create_directories(writeAccOut2File);
  }
  string acc_file_name = string(writeAccOut2File) + "/" + file_name + ".acc";
  ofstream acc_file(acc_file_name);
  if (top_n_percentage_stage1.empty())
	acc_file << "all branching times= 0" << endl;
  else{
	acc_file << "all branching times= " << int(top_n_percentage_stage1.begin()->second.second) << endl;
  acc_file << "stage one:" << endl;
  vector<pair<int, pair<int, double>>>
	  top_n_stage_1_vec(top_n_percentage_stage1.begin(), top_n_percentage_stage1.end());
  sort(top_n_stage_1_vec.begin(), top_n_stage_1_vec.end(), [](const pair<int, pair<int, double>> &a,
															  const pair<int, pair<int, double>> &b) {
	return a.first < b.first;
  });
  for (auto &pr : top_n_stage_1_vec) {
	acc_file << pr.first << "= " << pr.second.first / pr.second.second << endl;
  }
  acc_file << "stage two:" << endl;
  for (int j = 0; j < top_n_percentage_stage2_takeout_n.size(); ++j) {
	acc_file << "take out= " << j + 1 << endl;
	vector<pair<int, pair<int, double>>>
		top_n_stage_2_vec(top_n_percentage_stage2_takeout_n[j].begin(), top_n_percentage_stage2_takeout_n[j].end());
	sort(top_n_stage_2_vec.begin(), top_n_stage_2_vec.end(), [](const pair<int, pair<int, double>> &a,
																const pair<int, pair<int, double>> &b) {
	  return a.first < b.first;
	});
	for (auto &pr : top_n_stage_2_vec) {
	  acc_file << pr.first << "= " << pr.second.first / pr.second.second << endl;
	}
  }
	}
  acc_file.close();
#endif
#ifdef VERBOSE_MODE
  printBrDecisions();
#endif
}

void CVRP::printBrDecisions() {
  cout << BIG_PHASE_SEPARATION
	   << "Branch choice\n";
  int cnt = 0;
  for (auto &i : branch_choice) {
	if (!i.second)continue;
	cout << "<(" << i.first.first << ", " << i.first.second << ")  "
		 << i.second << ">  ";
	if (++cnt % 5 == 0)cout << "\n";
  }
  cout << "\nBranch Times= " << branch_times << endl;
}

void CVRP::printInfoLabeling(int iter, int num_added_col, int num_col, int num_row,
							 double mt, double spt, double et, double lp, double lb_val, double ub) {
  ios init(nullptr);
  init.copyfmt(cout);
  cout << fixed << setprecision(1);
  cout << "it= " << setw(2) << left << iter
	   << "  chgcol= " << setw(2) << left << num_added_col << "  ncol= " << setw(2) << left << num_col
	   << "  ncstr= " << setw(2) << left << num_row
	   << "  mt= " << setw(2)
	   << left << mt << "  spt= " << setw(2) << left << spt << "  et= "
	   << setw(2) << left << et << "  lpval= " << setw(2) << left << lp << "  lb= " << setw(2)
	   << left << lb_val << "  ub= " << setw(2) << left << ub << endl;
  cout.copyfmt(init);
}

#endif

void CVRP::printCutsInformation(BbNode *node) const {
  map<size_t, pair<int, double>> cuts_mem;
  for (auto &r1c : node->r1cs) {
#if LIMITED_MEMORY_TYPE == 1
	cuts_mem[r1c.info_r1c.first.size()].second += (int)r1c.mem.size();
#elif LIMITED_MEMORY_TYPE == 2
	cuts_mem[r1c.info_r1c.first.size()].second += (int)r1c.arc_mem.size();
#endif
	++cuts_mem[r1c.info_r1c.first.size()].first;
  }

  cout << SMALL_PHASE_SEPARATION;
  cout << "rccs= " << node->rccs.size() << endl;
  cout << "r1cs= " << node->r1cs.size() << endl;
  for (auto &it : cuts_mem) {
	cout << "R1c" << it.first << "s= " << it.second.first << " mem: "
		 << it.second.second / it.second.first << "|";
  }
  cout << "num_row= " << num_row << endl;
  cout << SMALL_PHASE_SEPARATION;
}

