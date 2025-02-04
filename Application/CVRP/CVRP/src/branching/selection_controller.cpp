#include "branching.hpp"

#include "train.hpp"
#include "predict.hpp"


using namespace std;
using namespace std::chrono;

void BaseBranching::controlBrSelection(pair<int, int> &info) {
    if (node->getIfTerminated()) return;
    dynamic_call(bool if_cal_b4 = false;
        if (cvrp->getIfInEnuState() && Dynamics::n_nodes_enu == 0) {
        if_cal_b4 = true;
        Dynamics::c_enu = Dynamics::solve_a_node_time;
        ++Dynamics::n_nodes_enu;
        }
    )

    cvrp->getIfAllowChangeCol() = false;
    cvrp->updateEdgeColMap(node,
                           true);
    current_branching_info.branch_pair.clear();

    ml_call(MASTER_VALVE_ML == ML_GET_DATA_1, GetTrainingData::generateModelPhase1())
    ml_call(MASTER_VALVE_ML == ML_GET_DATA_2, GetTrainingData::generateModelPhase2())
    ml_call(MASTER_VALVE_ML == ML_USE_MODEL, Predict::useMLInGeneralFramework())
    ml_call(MASTER_VALVE_ML == ML_USE_MODEL_1, Predict::useML1_N_LPTesting())
    ml_call(MASTER_VALVE_ML == 0, useDefaultSB())

    ml_call(MASTER_VALVE_ML, MachineLearning::cleanLastData())
    info = current_branching_info.branch_pair[0];
    cout << SMALL_PHASE_SEPARATION;
    cout << "brc= " << "( " << info.first << " , " << info.second << " )\n";
    ++branching_history.branch_choice[info];
    ++num_br;
    cvrp->getIfAllowChangeCol() = true;

    dynamic_call({
        double eps = Dynamics::solve_a_node_time + (Dynamics::solve_cg_time / 2);
        Dynamics::solve_cg_time = 0;
        if (cvrp->getIfInEnuState()) {
        if (if_cal_b4) {
        Dynamics::c_enu = 0;
        Dynamics::n_nodes_enu = 0;
        }
        Dynamics::updateStateAverage(eps, Dynamics::c_enu, Dynamics::n_nodes_enu);
        ++Dynamics::n_nodes_enu;
        } else {
        Dynamics::updateStateAverage(eps, Dynamics::c_b4, Dynamics::n_nodes_b4);
        ++Dynamics::n_nodes_b4;
        }
        }
    )
}

void BaseBranching::useDefaultSB() {
    int num1, num2;

    if (cvrp->getIfInEnuState()) {
        num1 = Config::BranPhase1InEnu;
        num2 = Config::BranPhase2InEnu;
    } else {
        num1 = Config::BranPhase1;
        num2 = Config::BranPhase2;
    }

    initialScreen(true, true, num1, Config::Frac4sudoCostBranPhase0);
    dynamic_call(Dynamics::giveDynamicK(num2, false))
    testLP(num2, true); //analyze what the m is!
    testCG(false, false, true, true);
}
