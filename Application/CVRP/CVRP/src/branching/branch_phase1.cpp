#include <utility>
#include "branching.hpp"
#include "heuristic.hpp"


#include "machine_learning.hpp"


using namespace std;
using namespace std::chrono;

#define PRINT_LP_SCORE

void BaseBranching::testLP(int num, bool if_writeBranch_pair, bool if_record_sb_scores,
                           bool if_record_LP_improvement) {
    auto &branch_pair = current_branching_info.branch_pair;
    if (branch_pair.size() == 1) {
        cout << "No LP Testing: a single candidate" << endl;
        return;
    }
    time_point<high_resolution_clock> beg, end;
    double eps;
    beg = high_resolution_clock::now();
    double org_val = node->getCurrentNodeVal();
    int BeforeNumRow = cvrp->getNumRow();
    int cnt = 0;
    vector<int> solver_ind;
    vector<double> solver_val;
    bool if_changed = false;

    vector<EdgeScoreInfo> edge_info(branch_pair.size());

    barrier_call(safe_solver(node->getSolver().setEnvCrossOver(SOLVER_CROSSOVER_DOWN)))
    for (auto &edge: branch_pair) {
        double tmp_val;
        cvrp->getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
        if (!if_changed) {
            safe_solver(cvrp->addBranchConstraint(solver_ind,
                solver_val,
                SOLVER_EQUAL, //use equal to avoid 0-1-2-1-2-0, (vio of triangle inequality)
                0,
                nullptr,
                node->getSolver()))
            if_changed = true;
        } else {
            cvrp->changeBranchConstraint(solver_ind,
                                         solver_val,
                                         SOLVER_EQUAL,
                                         0,
                                         BeforeNumRow,
                                         node->getSolver());
        }
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        safe_solver(node->getSolver().getObjVal(&tmp_val))
        auto dif1 = calculateDifference(tmp_val, org_val);
        safe_solver(cvrp->inverseLastBranchConstraint(SOLVER_EQUAL, 1, node->getSolver()))
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        safe_solver(node->getSolver().getObjVal(&tmp_val))
        auto dif2 = calculateDifference(tmp_val, org_val);
        edge_info[cnt].dif1 = dif1;
        edge_info[cnt].dif2 = dif2;
        edge_info[cnt++].edge = edge;
#ifdef PRINT_LP_SCORE
        verbose_call(
            auto product = dif1 * dif2;
            cout << "edge: " << edge.first << "-" << edge.second << " ldf: " << dif1 << " rdf: " << dif2 << " product: "
            << product << endl;)
#endif
    }
    safe_solver(node->getSolver().delConstraints(1, &BeforeNumRow))
    safe_solver(node->getSolver().updateModel())
    BaseBranching::reviseExtremeUnbalancedScore(edge_info, if_record_LP_improvement);
    dynamic_call({
        end = high_resolution_clock::now();
        eps = duration<double>(end - beg).count();
        Dynamics::getAverageT4LPNHeuristic(eps / (int)branch_pair.size() / 2);
        })
    vector<pair<pair<int, int>, double> > branch_val(edge_info.size());
    transform(edge_info.begin(),
              edge_info.end(),
              branch_val.begin(),
              [](const auto &a) {
                  return make_pair(a.edge, a.dif1 * a.dif2);
              });
    if (if_writeBranch_pair) {
        int size = min(num, int(branch_val.size()));
        sort(branch_val.begin(), branch_val.end(), [](const auto &a, const auto &b) {
            return a.second > b.second;
        });
        branch_pair.resize(size);
        transform(branch_val.begin(),
                  branch_val.begin() + size,
                  branch_pair.begin(),
                  [](const auto &a) {
                      return a.first;
                  });
    }

    if (if_record_sb_scores) {
        auto is_sorted = std::is_sorted(branch_val.begin(), branch_val.end(), [](const auto &a, const auto &b) {
            return a.second > b.second;
        });
        if (!is_sorted) {
            sort(branch_val.begin(), branch_val.end(), [](const auto &a, const auto &b) {
                return a.second > b.second;
            });
        }
        current_branching_info.branch_pair_val = branch_val;
    }

    auto &branch_lp = current_branching_info.branch_lp;
    std::unordered_map<std::pair<int, int>, double, PairHasher> branch_LP_tmp;
    branch_LP_tmp.reserve(branch_pair.size());
    for (auto &edge: branch_pair) {
        branch_LP_tmp[edge] = branch_lp[edge];
    }
    branch_lp = branch_LP_tmp;
    barrier_call(safe_solver(node->getSolver().setEnvCrossOver(SOLVER_CROSSOVER_DEFAULT)))
    safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();
    cout << "testLP: ";
    for (auto &edge: branch_pair) {
        cout << edge.first << "-" << edge.second << " " << branch_lp[edge] << ",";
    }
    cout << " time= " << eps << endl;
}
