#include "cvrp.hpp"
#include <filesystem>

#include "branching.hpp"
#include "read_node_in.hpp"
using namespace std;
using namespace chrono;

#ifdef VERBOSE_MODE

#ifdef HEURISTIC
#include "heuristic.hpp"
#endif

void printSolution(std::ostream &os, const vector<vector<int> > &ip_opt_sol,
                   int num_explored_nodes,
                   const string &file_name, double cap) {
    os << "<Solution> " << endl;
    for (auto &r: ip_opt_sol) {
        os << "0";
        for (auto &i: r) {
            os << "-" << i;
        }
        os << "-0" << endl;
    }
#ifdef HEURISTIC
    if (Heuristic::if_fail) {
        os << "<Optval  " << BaseBranching::ub << ">  <Etime  " << -Heuristic::heuristic_mip_time
                << ">  <Nodes_explored  " << 0 << ">"
                << "  <global_gap  " << BaseBranching::global_gap * 100 << "%>" << endl;
    } else {
        os << "<Optval  " << BaseBranching::ub << ">  <Etime  " << BaseBranching::glo_eps
                << ">  <Nodes_explored  " << num_explored_nodes << ">"
                << "  <global_gap  " << BaseBranching::global_gap * 100 << "%>" << endl;
    }
#else
    os << "<Optval  " << BaseBranching::ub << ">  <Etime  " << BaseBranching::glo_eps
            << ">  <Nodes_explored  " << num_explored_nodes << ">"
            << "  <global_gap  " << BaseBranching::global_gap * 100 << "%>" << endl;
#endif
    os << "<Instance  " << file_name << "  Capacity  " << cap << ">" << endl;
}

void CVRP::printOptIntSol() {
    BaseBranching::glo_end = chrono::high_resolution_clock::now();
    BaseBranching::glo_eps = duration<double>(BaseBranching::glo_end - BaseBranching::glo_beg).count();
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
    printSolution(cout,
                  ip_opt_sol,
                  BaseBranching::num_explored_nodes,
                  file_name,
                  cap);

#ifdef WRITE_SOL_OUT_TO_FILE
    if (!std::filesystem::exists(WRITE_SOL_OUT_TO_FILE)) {
        std::filesystem::create_directories(WRITE_SOL_OUT_TO_FILE);
    }

    string sol_file_name = string(WRITE_SOL_OUT_TO_FILE) + "/" + file_name + ".sol";
    ofstream sol_file(sol_file_name);
    printSolution(sol_file,
                  ip_opt_sol,
                  BaseBranching::num_explored_nodes,
                  file_name,
                  cap);
    sol_file.close();
#endif
    read_node_in_call(ReadNodeIn::tryUpdateUB())
#ifdef VERBOSE_MODE
    printBrDecisions();
#endif
}

void CVRP::printBrDecisions() {
    cout << BIG_PHASE_SEPARATION
            << "Branch choice\n";
    int cnt = 0;
    for (auto &i: BaseBranching::branching_history.branch_choice) {
        if (!i.second)continue;
        cout << "<(" << i.first.first << ", " << i.first.second << ")  "
                << i.second << ">  ";
        if (++cnt % 5 == 0)cout << "\n";
    }
    cout << "\nBranch Times= " << BaseBranching::num_br << endl;
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
    map<size_t, pair<int, double> > cuts_mem;
    for (auto &r1c: node->r1cs) {
        cuts_mem[r1c.info_r1c.first.size()].second += (int) r1c.arc_mem.size();
        ++cuts_mem[r1c.info_r1c.first.size()].first;
    }

    cout << SMALL_PHASE_SEPARATION;
    cout << "rccs= " << node->rccs.size() << endl;
    cout << "r1cs= " << node->r1cs.size() << endl;
    for (auto &it: cuts_mem) {
        cout << "R1c" << it.first << "s= " << it.second.first << " mem: "
                << it.second.second / it.second.first << "|";
    }
    cout << "num_row= " << num_row << endl;
    cout << SMALL_PHASE_SEPARATION;
}
