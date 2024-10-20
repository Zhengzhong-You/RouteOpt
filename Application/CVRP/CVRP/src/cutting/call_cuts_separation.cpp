#include "cvrp.hpp"
#ifdef HEURISTIC
#include "heuristic.hpp"
#endif
#include "rank1_cuts_separator.hpp"
#include "branching.hpp"
#include "robust_control.hpp"

using namespace std;

void CVRP::separateHybridCuts(BbNode *&node) {
    int cnt_tail_off = 0;
    double standard;
    int count = 0;
    bool if_opt = false;
    std::chrono::high_resolution_clock::time_point t1, t2;
    double duration;

    bool if_ban_cuts = false;
    setting_2_call(banCutsAtNonRoot(node->index == 0, if_ban_cuts))
    if (if_ban_cuts) return;

    if (if_in_enu_state) goto HYBRID;
    if (!node->index) {
        force_not_rollback = true;
        solver_vrptw_call(augmentNonGuillotineRound(node))
        separateRCCs(node);
        force_not_rollback = false;

        if (!node || node->getIfTerminated()) {
            if_opt = true;
            goto QUIT;
        }
        assignInitialLabelingMemory();
        RobustControl::val_b4_rank1 = node->getCurrentNodeVal();
    }

HYBRID:


    cout << MID_PHASE_SEPARATION;
    cout << "Switch on Hybrid...\n";
    while (true) {
        RobustControl::estimatedNumberOfNonRobustCutsCanBeAdded();

        if (Config::MaxNumR1C3PerRound == 0) {
            cout << "No R1C3s allowed. Skip the separation.\n";
            eliminateArcs(node);
            enumerateMIP(node);
            goto QUIT;
        }

        if (!if_in_enu_state && RobustControl::if_ever_roll_back) {
            verbose_call(cout << "copy model\n";)
            if (rollback_solver.model)rollback_solver.freeModel();
            rollback_solver.model = node->getSolver().copyModel();
            rollback_cols = node->getCols();
        }

        ++count;
        cout << BIG_PHASE_SEPARATION;
        cout << "Hybrid separation round " << count << endl;

        double prior_nodeVal = node->getCurrentNodeVal();

        int oldNum = num_row;

        auto beg = chrono::high_resolution_clock::now();
        if_in_enu_state ? generateRCCsInEnum(node) : generateRCCs(node);

        separateNAddRank1Cuts(node);

        auto end = chrono::high_resolution_clock::now();
        verbose_call(cout << "sep time: " << chrono::duration<double>(end - beg).count() << endl;)
        int cuts_sum = num_row - oldNum;
        if (!cuts_sum) {
            cout << "sep fails due to zero cuts" << endl;
            if (!if_in_enu_state && Rank1CutsSeparator::tellIfCutsAreChangedByMem(node->r1cs)) {
                cout << "fill mem! resolve the cg!" << endl;
                cnt_tail_off = Config::CutsTolerance;
                goto SOLVE_CG;
            }
            goto QUIT;
        }

        t1 = chrono::high_resolution_clock::now();
        bool if_give_up_this_cutting_round;
        deleteNonActiveCutsSafely(node, oldNum, if_give_up_this_cutting_round); //must be kept!
        t2 = chrono::high_resolution_clock::now();
        verbose_call(printCutsInformation(node))
        duration = chrono::duration<double>(t2 - t1).count();
        verbose_call(cout << "add time: " << duration << endl;)
        if (if_give_up_this_cutting_round) goto QUIT;

        RobustControl::recordNewGeneratedCuts(oldNum);
    SOLVE_CG:
        if (if_in_enu_state) {
            solveLPByInspection(node, false, false, true);
            if (!node->index) BaseBranching::updateLowerBound();
            if (node->getIfTerminated())goto QUIT;
            if (node->size_enumeration_col_pool + num_col <= MaxNumRoute4Mip) goto QUIT;
        } else {
        PRICING:
            getVCutMapLP(node);
            solveLPInLabeling(node);

            if (node->getIfTerminated()) {
                delete node;
                node = nullptr;
                goto QUIT;
            }

            int goto_state;
            RobustControl::robustControl(oldNum, prior_nodeVal, goto_state);
            if (goto_state == 0) goto PRICING;
            else if (goto_state == 1) goto QUIT;

            if (node->index == 0) BaseBranching::updateLowerBound();

            eliminateArcs(node);
            enumerateMIP(node);
            if (!node) goto QUIT;

            cleanIndexColForNode(node, true);
            findNonActiveCuts(node); //cannot be freely deleted in enu since change matrix is very expensive!
        }

        standard = calculateGapImprovement(node->getCurrentNodeVal(), prior_nodeVal);

        verbose_call(cout << "local gap= " << (double(BaseBranching::ub - node->getCurrentNodeVal()) / BaseBranching::ub > TOLERANCE ?
                double(BaseBranching::ub - node->getCurrentNodeVal()) / BaseBranching::ub * 100 : 0)
            << endl;
            cout << "gap improved= " << standard << endl;
            cout << SMALL_PHASE_SEPARATION;)
        if (BaseBranching::ub <= BaseBranching::lb_transformed) {
            if_opt = true;
            goto QUIT;
        }
        if (node->index != 0) {
            verbose_call(cout << "br_value_improved= " << node->br_value_improved << endl;)
            if (node->getCurrentNodeVal() - prior_nodeVal < CUTTING_BRANCHING_RATIO * node->br_value_improved) goto QUIT;
        }
#ifdef CONTROL_ROOT_GAP
	if (1 - node->getCurrentNodeVal() / BaseBranching::ub < CONTROL_ROOT_GAP) {
	  cout << "gap control\n";
	  goto QUIT;
	}
#endif
        if (standard < TOLERANCE) goto QUIT;
        else if (standard < Config::CutsTailOff) {
            ++cnt_tail_off;
            cout << "cnt_tail_off= " << cnt_tail_off << endl;
            if (cnt_tail_off >= Config::CutsTolerance)goto QUIT;
        }
    }

QUIT:
    if (!if_in_enu_state) {
        if (rollback_solver.model)
            rollback_solver.
                    freeModel();
        if (node && !if_opt) {
            if (node->index == 0) {
                tryCallHeuristicB4Branching(node);
                if (!node) {
                    goto
                            OUT;
                }
            }

            getVCutMapLP(node);
            cleanIndexColForNode(node,
                                 false);
            verbose_call(cout << "we clean the col pool!" << endl;
                cout
                << "Hybrid tail off detected. Halting cut separation. Columns remaining: " << num_col
                << endl;)
            cout << BIG_PHASE_SEPARATION;
            heuristic_call(Heuristic::heuristicMIP(this, node);)
        } else {
        OUT:
            cout << "Terminate the node!\n" << "Stop cut-separation!" <<
                    endl;
            cout << SMALL_PHASE_SEPARATION;
        }
    } else {
        verbose_call(
            cout << "Hybrid tail off detected. Halting cut separation. Columns remaining: " << num_col << endl;)
        cout << BIG_PHASE_SEPARATION;
    }
    if (node) {
        recordOptimalColumn(node); //very essential
        if (!node->index)
            BaseBranching::updateLowerBound();
    }
}






void CVRP::separateRCCs(BbNode *&node) {
    if (if_in_enu_state) throw runtime_error("this separator for rccs cannot be called in enu state");
    cout << BIG_PHASE_SEPARATION;
    cout << "Begin Rcc separation...\n";
    CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

    CMGR_CreateCMgr(&MyCutsCMP, 100);
    CMGR_CreateCMgr(&MyOldCutsCMP, 100);
    vector<int> solver_ind;
    vector<double> solver_val;
    int oldNum = num_row;

    double old_val = node->getCurrentNodeVal();
    while (true) {
        updateEdgeColMap(node, false);
        getEdgeInfo(node, false);

        for (int i = 1; i <= node->getNumEdges(); ++i) {
            if (!node->getEdgeTail()[i]) {
                node->getEdgeTail()[i] = dim;
            } else break;
        }

        CAPSEP_SeparateCapCuts(real_dim, demand, cap, node->getNumEdges(), node->getEdgeTail().data(),
                               node->getEdgeHead().data(), node->getEdgeValue().data(), MyOldCutsCMP,
                               MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                               &if_int_n_feasible, &max_vio, MyCutsCMP);

        for (int i = 1; i <= node->getNumEdges(); ++i) {
            if (node->getEdgeTail()[i] == dim) {
                node->getEdgeTail()[i] = 0;
            } else break;
        }

        if (!MyCutsCMP->Size) break;
        int cnt = 0;
        for (int i = 0; i < MyCutsCMP->Size; ++i) {
            Rcc rcc;
            auto &tmp_customerInfo = rcc.info_rcc_customer;
            for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
                tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
            }

            if (tmp_customerInfo.size() <= dim / 2) {
                rcc.form_rcc = true;
                rcc.rhs = MyCutsCMP->CPL[i]->RHS;
            } else {
                rcc.form_rcc = false;
                rcc.rhs = real_dim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
                auto &tmp_NoCustomerInfo = rcc.info_rcc_outside_customer;
                vector<bool> tmp(dim, false);
                for (int j: tmp_customerInfo) {
                    tmp[j] = true;
                }
                for (int j = 1; j < dim; ++j) {
                    if (!tmp[j]) {
                        tmp_NoCustomerInfo.emplace_back(j);
                    }
                }
            }

            if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) {
                continue;
            }

            ++cnt;
            rcc.idx_rcc = num_row;
            node->rccs.emplace_back(rcc);
            getCoefficientRCC(node, rcc, solver_ind, solver_val);
            safe_solver(
                node->getSolver().addConstraint(solver_ind.size(),
                    solver_ind.data(),
                    solver_val.data(),
                    SOLVER_LESS_EQUAL,
                    rcc.rhs,
                    nullptr))
            safe_solver(node->getSolver().updateModel())
            safe_solver(node->getSolver().getNumRow(&num_row))
        }

        if (!cnt) break;

        for (int i = 0; i < MyCutsCMP->Size; ++i) {
            CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
        }

#ifdef VERBOSE_MODE
        cout << "add rcc= " << cnt << endl;
#endif
        MyCutsCMP->Size = 0;

        getVCutMapLP(node);
        solveLPInLabeling(node);

        if (node->getIfTerminated()) {
            delete node;
            node = nullptr;
            goto QUIT;
        }
        eliminateArcs(node);
        enumerateMIP(node);

        if (!node) goto QUIT;
        cleanIndexColForNode(node, true);
        findNonActiveCuts(node); //cannot be freely

        BaseBranching::lb = node->getCurrentNodeVal();
        BaseBranching::lb_transformed = ceilTransformedNumberRelated(BaseBranching::lb - TOLERANCE);
        if (ceilTransformedNumberRelated(node->getCurrentNodeVal() - TOLERANCE) + TOLERANCE >= BaseBranching::ub) {
            node->getIfTerminated() = true;
            cout << TERMINATED_MESSAGE_SEP_RCC;
            break;
        }

        if (abs(node->getCurrentNodeVal() - old_val) < 0.1) break;
        old_val = node->getCurrentNodeVal();

#ifdef VERBOSE_MODE
        cout << "gap= " << (double(BaseBranching::ub - BaseBranching::lb) / BaseBranching::ub > TOLERANCE
                                ? double(BaseBranching::ub - BaseBranching::lb) / BaseBranching::ub * 100
                                : 0) << "%\n";
        cout << SMALL_PHASE_SEPARATION;
#endif
    }

QUIT:
    CMGR_FreeMemCMgr(&MyOldCutsCMP);
    CMGR_FreeMemCMgr(&MyCutsCMP);
}

void CVRP::generateRCCs(BbNode *node) {
    int cnt = 0;

    CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
    CMGR_CreateCMgr(&MyCutsCMP, 100);
    CMGR_CreateCMgr(&MyOldCutsCMP, 100);
    vector<int> solver_ind;
    vector<double> solver_val;
#if SOLVER_VRPTW == 1
  if (if_force_keep_rcc) popArcGraph(node);
#endif
    updateEdgeColMap(node, false);
    getEdgeInfo(node, false);
    for (int i = 1; i <= node->getNumEdges(); ++i) {
        if (!node->getEdgeTail()[i]) {
            node->getEdgeTail()[i] = dim;
        } else break;
    }

    CAPSEP_SeparateCapCuts(real_dim, demand, cap, node->getNumEdges(), node->getEdgeTail().data(),
                           node->getEdgeHead().data(), node->getEdgeValue().data(), MyOldCutsCMP,
                           MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                           &if_int_n_feasible, &max_vio, MyCutsCMP);

    for (int i = 1; i <= node->getNumEdges(); ++i) {
        if (node->getEdgeTail()[i] == dim) {
            node->getEdgeTail()[i] = 0;
        } else break;
    }

    if (!MyCutsCMP->Size) goto QUIT;
    cnt = 0;
    for (int i = 0; i < MyCutsCMP->Size; ++i) {
        Rcc rcc;
        auto &tmp_customerInfo = rcc.info_rcc_customer;
        for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
            tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
        }

        if (tmp_customerInfo.size() <= dim / 2) {
            rcc.form_rcc = true;
            rcc.rhs = MyCutsCMP->CPL[i]->RHS;
        } else {
            rcc.form_rcc = false;
            rcc.rhs = real_dim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;
            auto &tmp_NoCustomerInfo = rcc.info_rcc_outside_customer;
            vector<bool> tmp(dim, false);
            for (int j: tmp_customerInfo) tmp[j] = true;
            for (int j = 1; j < dim; ++j) {
                if (!tmp[j]) {
                    tmp_NoCustomerInfo.emplace_back(j);
                }
            }
        }

        if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) {
            continue;
        }

        ++cnt;

        rcc.idx_rcc = num_row;
#if SOLVER_VRPTW == 1
	if (if_force_keep_rcc) rcc.if_keep = true;
#endif
        node->rccs.emplace_back(rcc);

        getCoefficientRCC(node, rcc, solver_ind, solver_val);
        safe_solver(node->getSolver().addConstraint(solver_ind.size(),
            solver_ind.data(),
            solver_val.data(),
            SOLVER_LESS_EQUAL,
            rcc.rhs,
            nullptr))
        safe_solver(node->getSolver().updateModel())
        safe_solver(node->getSolver().getNumRow(&num_row))
    }
QUIT:
    for (int i = 0; i < MyCutsCMP->Size; ++i) {
        CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
    }
    MyCutsCMP->Size = 0;

    CMGR_FreeMemCMgr(&MyOldCutsCMP);
    CMGR_FreeMemCMgr(&MyCutsCMP);
}


void CVRP::separateNAddRank1Cuts(BbNode *node) {
    vector<RouteInfo> sol(node->only_frac_sol.first.size());
    for (int i = 0; i < sol.size(); ++i) {
        auto &s = sol[i];
        const auto &node_s = node->only_frac_sol.first[i];
        s.col_seq = node_s.col_seq;
        s.forward_concatenate_pos = node_s.forward_concatenate_pos;
        s.frac_x = node->only_frac_sol.second[i];
    }
    Rank1CutsSeparator::updateInfo(if_in_enu_state
                                       ? NO_MEMORY
                                       : (RobustControl::rank1_cuts_mode == 0 ? NODE_MEMORY : ARC_MEMORY),
                                   RobustControl::pricing_hard_level,
                                   sol, node->r1cs);
    Rank1CutsSeparator::separateRank1Cuts();
    if (!if_in_enu_state) {
        Rank1CutsSeparator::selectR1CsByVioNMemory();
    }
    vector<R1c> new_cuts;
    Rank1CutsSeparator::getSeparatedCuts(new_cuts);

    if (if_in_enu_state) {
        addR1CAtOnceInEnum(node, new_cuts);
    } else {
        RobustControl::rank1_cuts_mode == 0
            ? addLimitedMemoryR1CsNodeBased(node, new_cuts)
            : addLimitedMemoryR1Cs(node, new_cuts);
    }
}
