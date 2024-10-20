//
// Created by You, Zhengzhong on 7/23/24.
//


#include "robust_control.hpp"

using namespace std;

CVRP *RobustControl::cvrp = nullptr;
BbNode *RobustControl::node = nullptr;
bool RobustControl::rank1_cuts_mode = false;
int RobustControl::initial_max_num_r1c_per_round = 0;
int RobustControl::initial_max_num_r1c3_per_round = 0;
std::vector<R1c> RobustControl::new_generated_cuts;
bool RobustControl::if_fix_resource_point = false;
int RobustControl::pricing_hard_level = 0;
bool RobustControl::if_ever_roll_back = false;
std::vector<yzzLong> RobustControl::solution_arcs;
double RobustControl::val_b4_rank1 = 0;

void selectCutsByMemoryConstraint(BbNode *node, int dim,
                                  vector<R1c> &cuts);


void RobustControl::robustControl(int oldNum, double &prior_nodeVal, int &goto_state) {
    goto_state = 3; //nothing happens;
    if (cvrp->getIfInEnuState()) return;
    if (node->index == 0) {
        if (rank1_cuts_mode == false) {
            if (!cvrp->if_exact_cg_finished) {
                convert2ArcBasedMemory();
                prior_nodeVal = val_b4_rank1;
                cvrp->if_tail_off = false;
                goto_state = 0; //goto pricing;
            }
        } else {
            if (cvrp->if_roll_back) {
                cvrp->rollbackEasyWay(node, oldNum);
                if (RobustControl::restoreCuts()) {
                    goto_state = 0; //goto pricing;
                    return;
                }
                goto_state = 1; //goto quit;
            }
        }
    } else {
        if (cvrp->if_roll_back) {
            cvrp->rollbackEasyWay(node, oldNum);
            goto_state = 1; //goto quit;
        }
    }
}

void RobustControl::convert2ArcBasedMemory() {
    if (cvrp->getIfInEnuState()) return;
    if (rank1_cuts_mode == false) {
        verbose_call(cout << "switch to arc-based memory! the lb might be lower!" << endl;)
        pricing_hard_level = 2;
        rank1_cuts_mode = true;
        vector<vector<int> > solution_set(cvrp->dim);
        for (int i = 1; i < cvrp->dim; ++i) {
            auto &arc = solution_arcs[i];
            auto &set = solution_set[i];
            set.resize(arc.count());
            int cnt = 0;
            for (int j = 1; j < cvrp->dim; ++j) {
                if (arc.test(j)) {
                    set[cnt++] = j;
                }
            }
        }
        for (auto &c: node->r1cs) {
            for (auto &it: c.arc_mem) {
                it.first = solution_set[it.second];
            }
        }
        cvrp->getVCutMapLP(node);
        vector<int> idx(node->r1cs.size());
        iota(idx.begin(), idx.end(), 0);
        cvrp->addR1CAtOnce(node, idx);
#ifdef SOLVER_VRPTW
	if_fix_resource_point = true;
#endif
    }
}

void RobustControl::recordSolutionArcs() {
    if (cvrp->getIfInEnuState()) return;
    if (node->index == 0 && rank1_cuts_mode == false) {
        if (solution_arcs.empty()) solution_arcs.resize(cvrp->dim, 0);
        for (auto &c: node->all_lp_sol.first) {
            if (c.col_seq.size() == 1) continue;
            for (int i = 0; i < c.forward_concatenate_pos; ++i) {
                solution_arcs[c.col_seq[i + 1]].set(c.col_seq[i]); //i goto i+1
            }
            int after_forward = c.forward_concatenate_pos + 1;
            for (int i = (int) c.col_seq.size() - 1; i > after_forward; --i) {
                solution_arcs[c.col_seq[i - 1]].set(c.col_seq[i]); //i goto i-1
            }
            if (c.forward_concatenate_pos >= 0 && after_forward < (int) c.col_seq.size())
                solution_arcs[c.col_seq[after_forward]].set(c.col_seq[c.forward_concatenate_pos]);
        }
    } else {
        solution_arcs.clear();
    }
}

void RobustControl::init(CVRP *pr_cvrp) {
    cvrp = pr_cvrp;
    getInitialMaxNumR1CPerRound();
}

void RobustControl::getInitialMaxNumR1CPerRound() {
    initial_max_num_r1c_per_round = Config::MaxNumR1CPerRound;
    initial_max_num_r1c3_per_round = Config::MaxNumR1C3PerRound;
}

void RobustControl::estimatedNumberOfNonRobustCutsCanBeAdded() {
    if (cvrp->getIfInEnuState()) {
#ifdef DYNAMIC_ADD_CUTS
	ORIGIN:
#endif
        Config::MaxNumR1CPerRound = initial_max_num_r1c_per_round;
        Config::MaxNumR1C3PerRound = initial_max_num_r1c3_per_round;
        return;
    }

    if (cvrp->if_tail_off) {
        Config::MaxNumR1C3PerRound = 0;
        Config::MaxNumR1CPerRound = 0;
    } else {
#ifdef DYNAMIC_ADD_CUTS
	const int ROUND_UNIT = 20;
	double ratio = 1.0 - min(cvrp->max_labeling_time_for_node_so_far / Config::HardTimeThresholdInPricing, 1.0);
	if (ratio > 0.8) goto ORIGIN;
	Config::MaxNumR1CPerRound =
		int((initial_max_num_r1c_per_round * ratio + ROUND_UNIT / 2.) / ROUND_UNIT) * ROUND_UNIT;
	Config::MaxNumR1C3PerRound =
		int((initial_max_num_r1c3_per_round * ratio + ROUND_UNIT / 2.) / ROUND_UNIT) * ROUND_UNIT;
	if (Config::MaxNumR1CPerRound + Config::MaxNumR1C3PerRound < MINIMUM_CUTS_ADDED) {
	  Config::MaxNumR1CPerRound = 0;
	  Config::MaxNumR1C3PerRound = 0;
	}
	verbose_call(cout << "Config::MaxNumR1CPerRound: " << Config::MaxNumR1CPerRound << endl;
					 cout << "Config::MaxNumR1C3PerRound: " << Config::MaxNumR1C3PerRound << " ratio: " << 1.0 - min(
						 cvrp->max_labeling_time_for_node_so_far / Config::HardTimeThresholdInPricing, 1.0)
						  << " max_labeling_time_for_node_so_far: "
						  << cvrp->max_labeling_time_for_node_so_far << endl;)
#else
        Config::MaxNumR1CPerRound = initial_max_num_r1c_per_round;
        Config::MaxNumR1C3PerRound = initial_max_num_r1c3_per_round;
#endif
    }
}

void RobustControl::recordNewGeneratedCuts(int oldNum) {
    if (cvrp->getIfInEnuState()) return;
    new_generated_cuts.clear();
    for (auto &r1c: node->r1cs) {
        if (r1c.idx_r1c >= oldNum) {
            new_generated_cuts.emplace_back();
            auto &new_r1c = new_generated_cuts.back();
            new_r1c = r1c;
            new_r1c.arc_mem.clear();
        }
    }
}

void RobustControl::updateNode(BbNode *pr_node) {
    node = pr_node;
}

bool RobustControl::restoreCuts() {
    if (cvrp->getIfInEnuState()) return true;
    if (cvrp->max_labeling_time_for_node_so_far > Config::CutGenTimeThresholdInPricingInitial) {
        Config::CutGenTimeThresholdInPricingInitial /= Config::SoftTimeDecreaseFactorInPricing;
        new_generated_cuts.clear();
        verbose_call(
            cout << "Config::CutGenTimeThresholdInPricingInitial: " << Config::CutGenTimeThresholdInPricingInitial <<
            endl;
            cout << "no cuts are recommended to add!" << endl;)
    } else new_generated_cuts.resize(int(new_generated_cuts.size() / 2));
    verbose_call(cout << "new_generated_cuts.size():" << new_generated_cuts.size() << endl;)
    if (new_generated_cuts.size() < MINIMUM_CUTS_ADDED) return false;
    Rank1CutsSeparator::findMemory4Cuts(new_generated_cuts, false);
    selectCutsByMemoryConstraint(node, cvrp->dim, new_generated_cuts);
    cvrp->addLimitedMemoryR1Cs(node, new_generated_cuts);
    return true;
}

void selectCutsByMemoryConstraint(BbNode *node, int dim,
                                  vector<R1c> &cuts) {
    vector<int> vertex_related_r1c(dim, MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX - 1);

    for (auto &r1c: node->r1cs) {
        for (auto &it: r1c.info_r1c.first) {
            --vertex_related_r1c[it];
        }
        for (auto &it: r1c.arc_mem) {
            --vertex_related_r1c[it.second];
        }
    }

    for (int i = 0; i < dim; ++i) {
        if (vertex_related_r1c[i] < 0) {
            cout << "vertex_related_r1c[i] < 0" << endl;
            throw runtime_error("vertex_related_r1c[i] < 0");
        }
    }

    vector<R1c> select_cut(cuts.size());
    int num = 0;
    for (auto &c: cuts) {
        bool if_keep = true;
        auto tmp = vertex_related_r1c;
        for (auto j: c.info_r1c.first) {
            if (tmp[j] <= 0) {
                if_keep = false;
                break;
            }
            --tmp[j];
        }
        for (auto &j: c.arc_mem) {
            if (tmp[j.second] <= 0) {
                if_keep = false;
                break;
            }
            --tmp[j.second];
        }
        if (!if_keep) {
            cout << "some vertex has been used!" << endl;
            continue;
        }
        vertex_related_r1c = tmp;
        select_cut[num++] = c;
    }
    select_cut.resize(num);
    cuts = select_cut;
}
