//
// Created by You, Zhengzhong on 4/26/24.
//

#include "best_bound_first_branching.hpp"
#include "read_node_in.hpp"
#ifdef WRITE_ENUMERATION_TREES
#include "write_enumeration_tree.hpp"
#endif

using namespace std;
using namespace std::chrono;

QueueTree BestBoundFirstBranching::bbt{};
QueueTree BestBoundFirstBranching::sub_bbt{};

void BestBoundFirstBranching::solve(QueueTree &tree) {
#ifdef KILL_ENUMERATION_TREES
  if (cvrp->getIfInEnuState()) {
	node = tree.top();
	tree.pop();
	freeNode();
	cvrp->getIfInEnuState() = false;
  	++num_explored_nodes;
	return;
  }
#endif
    while (!tree.empty()) {
        glo_end = high_resolution_clock::now();
        glo_eps = duration<double>(glo_end - glo_beg).count();

        node = tree.top();
        updateLowerBound();
        read_node_in_call(ReadNodeIn::tryUpdateUB())
        tree.pop();
        ++num_explored_nodes;
        nd_rmn = static_cast<int>(tree.size());

        tellIfTerminateTree();
        if (if_terminate_tree) {
            if_terminate_tree = false;
            goto QUIT;
        }

        resetEnv();

        write_node_out_call(
            if(nd_rmn>=1 && !cvrp->getIfInEnuState()) {
            WriteNodeOut::writeNodeOut(node);
            if (!node) {
            --num_explored_nodes; //since it has been written out
            continue;
            }
            }
        )

        solveNode();


        if (!node) continue;
        pair<int, int> info;
        controlBrSelection(info);


        if (cvrp->getIfInEnuState()) {
            cvrp->addBranchCutToUnsolvedInEnu(node, info);
        } else cvrp->addBranchCutToUnsolved(node, info);
    }

    if (!cvrp->getIfInEnuState())lb_transformed = ub;
QUIT:
    while (!tree.empty()) {
        const BbNode *extra_node = tree.top();
        tree.pop();
        ++num_explored_nodes;
        delete extra_node;
    }
    if (cvrp->getIfInEnuState()) {
        cvrp->getIfInEnuState() = false;
    } else {
        global_gap = (ub - lb_transformed) / ub;
        cvrp->printOptIntSol();
        write_enumeration_trees_call(WriteEnumerationTree::writeEnuCols())
        read_node_in_call(ReadNodeIn::rmNodeFile())
    }
}

void BestBoundFirstBranching::updateLowerBound() {
    double tmp_lb;
    if (!bbt.empty()) {
        tmp_lb = bbt.top()->value;
    } else {
        if (!sub_bbt.empty()) {
            tmp_lb = sub_bbt.top()->value;
        } else {
            tmp_lb = node->getCurrentNodeVal();
        }
    }
    tmp_lb = std::min(tmp_lb, node->getCurrentNodeVal());
    if (lb != tmp_lb) {
        lb = tmp_lb;
        lb_transformed = cvrp->ceilTransformedNumberRelated(lb - TOLERANCE);
    }
    if (!sub_bbt.empty()) {
        tmp_lb = sub_bbt.top()->value;
    } else {
        tmp_lb = node->getCurrentNodeVal();
    }
    tmp_lb = std::min(tmp_lb, node->getCurrentNodeVal());
    if (sub_lb != tmp_lb) {
        sub_lb = tmp_lb;
        sub_lb_transformed = cvrp->ceilTransformedNumberRelated(sub_lb - TOLERANCE);
    }
}
