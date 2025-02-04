/**
 * searchInDepthFashion: search depth-first
 */

#include "depth_first_branching.hpp"
#include "read_node_in.hpp"
#ifdef WRITE_ENUMERATION_TREES
#include "write_enumeration_tree.hpp"
#endif

using namespace std;
using namespace chrono;

StackTree DepthFirstBranching::bbt{};
StackTree DepthFirstBranching::sub_bbt{};
std::list<BbNode *> DepthFirstBranching::bbt_set{};
std::list<BbNode *> DepthFirstBranching::sub_bbt_set{};

void DepthFirstBranching::solveOne(StackTree &tree) {
    if (tree.empty()) {
        node = nullptr;
        return;
    }
    takeNodeOut(tree);
    updateLowerBound();

    resetEnv();
    write_node_out_call(
        if(nd_rmn>=1 && !cvrp->getIfInEnuState()) {
        WriteNodeOut::writeNodeOut(node);
        if (!node) {
        --num_explored_nodes; //since it has been written out
        return;
        }
        }
    )
    solveNode();
}

void DepthFirstBranching::solve(StackTree &tree) {
BEGIN:
    solveOne(tree);
    double node_val1 = node ? node->getCurrentNodeVal() : numeric_limits<double>::max();
    auto node1 = node;
    solveOne(tree);
    double node_val2 = node ? node->getCurrentNodeVal() : numeric_limits<double>::max();
    auto node2 = node;

    read_node_in_call(ReadNodeIn::tryUpdateUB())

    if (node1) {
        if (node2) {
            if (node_val1 < node_val2) {
                addNodeIn(tree, node1);
                node = node2;
            } else {
                addNodeIn(tree, node2);
                node = node1;
            }
        } else node = node1;
    } else node = node2;

    if (!node) {
        if (tree.empty()) {
            goto QUIT;
        } else {
            takeNodeOut(tree);
            updateLowerBound();
        }
    }

    if (node) {
        glo_end = high_resolution_clock::now();
        glo_eps = duration<double>(glo_end - glo_beg).count();

        tellIfTerminateTree();
        if (if_terminate_tree) {
            if_terminate_tree = false;
            goto QUIT;
        }

        resetEnv();
        cvrp->recordOptimalColumn(node, true);

        pair<int, int> info;
        controlBrSelection(info);
        if (cvrp->getIfInEnuState()) {
            cvrp->addBranchCutToUnsolvedInEnu(node, info);
        } else cvrp->addBranchCutToUnsolved(node, info);
        goto BEGIN;
    }
    if (!cvrp->getIfInEnuState()) lb_transformed = ub;
QUIT:
    while (!tree.empty()) {
        takeNodeOut(tree);
        freeNode();
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

void DepthFirstBranching::updateLowerBound() {
    double tmp_lb;
    if (!bbt_set.empty()) {
        tmp_lb = (*min_element(bbt_set.begin(), bbt_set.end(), [](BbNode *a, BbNode *b) {
            return a->getCurrentNodeVal() < b->getCurrentNodeVal();
        }))->getCurrentNodeVal();
    } else {
        tmp_lb = node->getCurrentNodeVal();
    }
    tmp_lb = std::min(tmp_lb, node->getCurrentNodeVal());
    if (lb != tmp_lb) {
        lb = tmp_lb;
        lb_transformed = cvrp->ceilTransformedNumberRelated(lb - TOLERANCE);
    }
    if (!sub_bbt_set.empty()) {
        tmp_lb = (*min_element(sub_bbt_set.begin(), sub_bbt_set.end(), [](BbNode *a, BbNode *b) {
            return a->getCurrentNodeVal() < b->getCurrentNodeVal();
        }))->getCurrentNodeVal();
    } else {
        tmp_lb = node->getCurrentNodeVal();
    }
    tmp_lb = std::min(tmp_lb, node->getCurrentNodeVal());
    if (sub_lb != tmp_lb) {
        sub_lb = tmp_lb;
        sub_lb_transformed = cvrp->ceilTransformedNumberRelated(sub_lb - TOLERANCE);
    }
}

void DepthFirstBranching::takeNodeOut(StackTree &tree) {
    node = tree.top();
    tree.pop();
    if (&tree == &bbt) {
        bbt_set.remove(node);
    } else if (&tree == &sub_bbt) {
        sub_bbt_set.remove(node);
    }
    ++num_explored_nodes;
    nd_rmn = static_cast<int>(tree.size());
}

void DepthFirstBranching::addNodeIn(StackTree &tree, BbNode *pr_node) {
    tree.push(pr_node);
    if (&tree == &bbt)
        bbt_set.emplace_front(pr_node);
    else if (&tree == &sub_bbt)
        sub_bbt_set.emplace_front(pr_node);
}
