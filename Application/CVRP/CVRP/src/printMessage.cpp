//
// Created by Zhengzhong You on 6/3/22.
//

#include "CVRP.hpp"
#include <sys/stat.h>
#include <experimental/filesystem>

using namespace std;
using namespace chrono;

#ifdef PrintOptimalSolution

void printSolution(std::ostream &os, const vector<vector<int>> &IPOptSol, double UB,
				   double GloEps, int NumExploredNodes, double GlobalGap, const string &FileName, double Cap) {
  os << "<Solution> " << endl;
  for (auto &r : IPOptSol) {
	os << "0";
	for (int i = 1; i < r.size(); ++i) {
	  os << "-" << r[i];
	}
	os << endl;
  }
  os << "<Optval  " << UB << ">  <Etime  " << GloEps
	 << ">  <Nodes_explored  " << NumExploredNodes << ">"
	 << "  <GlobalGap  " << GlobalGap * 100 << "%>" << endl;
  os << "<Instance  " << FileName << "  Capacity  " << Cap << ">" << endl;
}

void CVRP::printOptIntSol() {
  GloEnd = chrono::high_resolution_clock::now();
  GloEps = duration<double>(GloEnd - GloBeg).count();
#ifdef Resolve_Ins_with_Optimal_Val
  if (!CONFIG::Marker4FindLargeGapIns) {
	if (IPOptSol.empty()) {
	  cout << "No solution found" << endl;
	  CONFIG::Marker4FindLargeGapIns = -1;
	  return;
	}
	CONFIG::UB = UB;
	return;
  }
#endif
  printSolution(cout, IPOptSol, UB, GloEps, NumExploredNodes, GlobalGap, FileName, Cap);

#ifdef writeSolOut2File
//check if dir exists
  if (!std::experimental::filesystem::exists(writeSolOut2File)) {
	std::experimental::filesystem::create_directories(writeSolOut2File);
  }

  string sol_file_name = string(writeSolOut2File) + "/" + FileName + ".sol";
  ofstream sol_file(sol_file_name);
  printSolution(sol_file, IPOptSol, UB, GloEps, NumExploredNodes, GlobalGap, FileName, Cap);
  sol_file.close();
#endif
#ifdef writeAccOut2File
  if (!std::experimental::filesystem::exists(writeAccOut2File)) {
	std::experimental::filesystem::create_directories(writeAccOut2File);
  }
  string acc_file_name = string(writeAccOut2File) + "/" + FileName + ".acc";
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
#ifdef VERBOSE
  printBrDecisions();
#endif
#ifdef if_draw_BBT_graph
  drawBBTgraph();
#endif
//#ifdef AccuracyTest
//  cout << "average phase 1= " << average_rank_phase1.first / average_rank_phase1.second << endl;
//  cout << "average phase 2= " << average_rank_phase2.first / average_rank_phase2.second << endl;
//#endif
}

#endif

#ifdef VERBOSE

void CVRP::printBrDecisions() {
  cout << BIG_PHASE_SEPARATION
	   << "Branch choice\n";
  int cnt = 0;
  for (auto &i : BranchChoice) {
	if (!i.second)continue;
	cout << "<(" << i.first.first << ", " << i.first.second << ")  "
		 << i.second << ">  ";
	if (++cnt % 5 == 0)cout << "\n";
  }
  cout << "\nBranch Times= " << BranchTimes << endl;
}

void CVRP::printInfoLabeling(int iter, int num_added_col, int num_col, int num_row,
							 double mt, double spt, double et, double lp, double lb_val, double ub) {
  ios init(nullptr);
  init.copyfmt(cout);
  cout << fixed << setprecision(2);
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

void CVRP::printCutsInfo(BBNODE *node) const {
  map<size_t, pair<int, double>> cuts_mem;
  for (auto &r1c : node->R1Cs) {
	cuts_mem[r1c.InfoR1C.size()].second += (int)r1c.Mem.size();
	++cuts_mem[r1c.InfoR1C.size()].first;
  }

  for (auto &r1c : node->R1Cs_multi) {
	cuts_mem[r1c.InfoR1C.first.size()].second += (int)r1c.Mem.size();
	++cuts_mem[r1c.InfoR1C.first.size()].first;
  }

  cout << SMALL_PHASE_SEPARATION;
  cout << "RCCs= " << node->RCCs.size() << endl;
  cout << "R1Cs= " << node->R1Cs.size() + node->R1Cs_multi.size() << endl;
  for (auto &it : cuts_mem) {
	cout << "R1C" << it.first << "s= " << it.second.first << " Mem: "
		 << it.second.second / it.second.first << endl;
  }
  cout << "NumRow= " << NumRow << endl;
  cout << SMALL_PHASE_SEPARATION;
}

