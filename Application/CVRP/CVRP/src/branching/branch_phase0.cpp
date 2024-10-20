#include <utility>
#include "branching.hpp"

#ifdef USE_M_DYNAMICS
#include "dynamics.hpp"
#endif

using namespace std;
using namespace std::chrono;

#define EDGE_TOLERANCE 1e-4

void printColumns(BbNode *node, const Brc &brc);

void checkFracSet(BbNode *node);


void BaseBranching::initialScreen(
    bool if_record_source,
    bool if_fill_candidates_gap,
    int num,
    double pseudo_frac) {
    cvrp->getEdgeInfo(node, true);
    checkFracSet(node);
    vector<tuple<int, int, double> > fracEdges;
    vector<tuple<int, int, double> > OldBranch;
    int num_all_frac_edge = 0;
    for (int i = 1; i <= node->getNumEdges(); ++i) {
        if (node->getEdgeValue()[i] < 1 - EDGE_TOLERANCE) {
            ++num_all_frac_edge;
            int ai = node->getEdgeTail()[i];
            int aj = node->getEdgeHead()[i];
            double frac_down = node->getEdgeValue()[i];
            double frac_up = 1 - frac_down;
            pair<int, int> edge = {ai, aj};
            std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_up_ptr;
            std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher> *improvement_down_ptr;
            std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator up_iter;
            std::unordered_map<std::pair<int, int>, std::pair<double, int>, PairHasher>::iterator down_iter;
            int loop = 1;
        HERE:
            switch (loop) {
                case 1: {
                    improvement_up_ptr = &BaseBranching::branching_history.exact_improvement_up;
                    improvement_down_ptr = &BaseBranching::branching_history.exact_improvement_down;
                    break;
                }
                case 2: {
                    improvement_up_ptr = &BaseBranching::branching_history.heuristic_improvement_up;
                    improvement_down_ptr = &BaseBranching::branching_history.heuristic_improvement_down;
                    break;
                }
                case 3: {
                    improvement_up_ptr = &BaseBranching::branching_history.lp_testing_improvement_up;
                    improvement_down_ptr = &BaseBranching::branching_history.lp_testing_improvement_down;
                    break;
                }
                default: {
                    goto END;
                }
            }
            up_iter = improvement_up_ptr->find(edge);
            down_iter = improvement_down_ptr->find(edge);
            if (up_iter != improvement_up_ptr->end() && down_iter != improvement_down_ptr->end()) {
                double up = up_iter->second.first / up_iter->second.second;
                double down = down_iter->second.first / down_iter->second.second;
                OldBranch.emplace_back(ai, aj, up * frac_up * down * frac_down);
                continue;
            } else {
                ++loop;
                goto HERE;
            }
        END:
            fracEdges.emplace_back(ai, aj, abs(node->getEdgeValue()[i] - 0.5));
        }
    }
    std::sort(OldBranch.begin(), OldBranch.end(), [](const auto &a, const auto &b) {
        return get<2>(a) > get<2>(b);
    });
    std::sort(fracEdges.begin(), fracEdges.end(), [](const auto &a, const auto &b) {
        return get<2>(a) < get<2>(b);
    });
    dynamic_call(
        {
        Dynamics::evaluateM1();
        ml_call(MASTER_VALVE_ML == 0, dynamic_call(Dynamics::giveDynamicK(num))) //for 3PB
        }
    )
    int all_branch_phase1 = min(num, num_all_frac_edge);
    int sudo_cap = min(int(OldBranch.size()), int(all_branch_phase1 * pseudo_frac));
    int frac_cap = min(int(fracEdges.size()), all_branch_phase1 - sudo_cap);
    if (if_fill_candidates_gap) {
        if (frac_cap < all_branch_phase1 - sudo_cap) {
            sudo_cap = min(int(OldBranch.size()), all_branch_phase1 - frac_cap);
        }
    }
    all_branch_phase1 = sudo_cap + frac_cap;
    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    branch_pair.resize(all_branch_phase1);
    transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, branch_pair.begin(), [](const auto &a) {
        return make_pair(get<0>(a), get<1>(a));
    });
    transform(fracEdges.begin(), fracEdges.begin() + frac_cap, branch_pair.begin() + sudo_cap, [](const auto &a) {
        return make_pair(get<0>(a), get<1>(a));
    });
    if (if_record_source) {
        auto &branch_pair_from_pseudo = BaseBranching::current_branching_info.branch_pair_from_pseudo;
        branch_pair_from_pseudo.resize(sudo_cap);
        transform(OldBranch.begin(), OldBranch.begin() + sudo_cap, branch_pair_from_pseudo.begin(), [](const auto &a) {
            return make_pair(get<0>(a), get<1>(a));
        });
        auto &branch_pair_from_fractional = BaseBranching::current_branching_info.branch_pair_from_fractional;
        branch_pair_from_fractional.resize(frac_cap);
        transform(fracEdges.begin(), fracEdges.begin() + frac_cap, branch_pair_from_fractional.begin(),
                  [](const auto &a) {
                      return make_pair(get<0>(a), get<1>(a));
                  });
        verbose_call({
            cout << "pseudo: " << sudo_cap << endl;
            cout << "frac: " << frac_cap << endl;
            })
    }
    verbose_call(cout << "initialScreen branch_pair size: " << branch_pair.size() << endl)
}



double BaseBranching::calculateDifference(double tmp_val, double prior_val, bool if_consider_ub) {
    if (if_consider_ub) {
        tmp_val = min(tmp_val, BaseBranching::ub);
    }
    return max(tmp_val - prior_val, TOLERANCE);
}


void printColumns(BbNode *node, const Brc &brc) {
    int num_col;
    safe_solver(node->getSolver().getNumCol(&num_col))
    vector<double> x(num_col);
    safe_solver(node->getSolver().getX(0, num_col, x.data()))
    double sum_ = 0;
    for (int i = 0; i < num_col; ++i) {
        if (x[i] > TOLERANCE) {
            bool is = false;
            for (auto &seq: node->getCols()[i].col_seq) {
                if (seq == brc.edge.first || seq == brc.edge.second) {
                    is = true;
                    break;
                }
            }
            if (is) {
                for (auto &seq: node->getCols()[i].col_seq) {
                    cout << seq << " ";
                }
                cout << " , " << x[i] << endl;
                sum_ += x[i];
            }
        }
    }
    cout << "sum: " << sum_ << endl;
}

void checkFracSet(BbNode *node) {
    unordered_map<pair<int, int>, double, PairHasher> edge_set;
    for (int i = 1; i <= node->getNumEdges(); ++i) {
        edge_set[make_pair(node->getEdgeTail()[i], node->getEdgeHead()[i])] = node->getEdgeValue()[i];
    }
    for (auto &brc: node->getBrCs()) {
        if (brc.br_dir) {
            if (edge_set.find(brc.edge) == edge_set.end()) {
                printColumns(node, brc);
                throw runtime_error(
                    "edge: " + to_string(brc.edge.first) + " " + to_string(brc.edge.second) + " not found");
            }
            if (abs(edge_set[brc.edge] - 1) > EDGE_TOLERANCE) {
                printColumns(node, brc);
                throw runtime_error("edge: " + to_string(brc.edge.first) + " " + to_string(brc.edge.second) + " :"
                                    + to_string(edge_set[brc.edge]) + ", not 1");
            }
        } else {
            if (edge_set.find(brc.edge) != edge_set.end()) {
                throw runtime_error("edge: " + to_string(brc.edge.first) + " " + to_string(brc.edge.second) + " :"
                                    + to_string(edge_set[brc.edge]));
            }
        }
    }
}
