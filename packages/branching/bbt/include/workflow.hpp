/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file workflow.hpp
 * @brief Implementation of the branch-and-bound tree (BBT) workflow for solving optimization problems.
 *
 * This file contains the implementation of the BBTController class, which manages the workflow of
 * solving a branch-and-bound tree.
 */

#ifndef ROUTE_OPT_WORKFLOW_HPP
#define ROUTE_OPT_WORKFLOW_HPP

#include <set>
#include "candidate_selector_controller.hpp"

namespace RouteOpt::Branching::BBT {
    /**
     * @brief Solves the branch-and-bound tree (BBT).
     *
     * The solve function performs a branch-and-bound tree search to solve the optimization problem.
     * It processes nodes from the tree, performs pricing and cutting on each node, and uses BKF
     * controllers and branching testers to decide on the best branching candidate.
     *
     * @tparam Node       Type representing a node in the branch-and-bound tree.
     * @tparam BrCType    Type representing the branching candidate.
     * @tparam Hasher     Type for hashing branching candidates.
     * @param root_node   Pointer to the root node of the branch-and-bound tree.
     */
    template<typename Node, typename BrCType, typename Hasher>
    void BBTController<Node, BrCType, Hasher>::solve(Node *root_node) {
        // Add the root node to the tree.
        addNode(root_node);
        // Optionally, read additional information for the root node if the callback is defined.
        if (tryReadNodeIn != nullptr)
            tryReadNodeIn(root_node, branching_history, bkf_data_shared);

        // Continue processing nodes until the tree is empty.
        while (!tree.empty()) {
            // Select and remove the first node from the tree.
            auto node = *tree.begin();
            tree.erase(node);

            // Print current tree details for debugging/monitoring.
            printTreeDetail();

            // Extract the current node's value and branching candidate details.
            auto val = valueExtractor(node);
            auto [edge, dir, tree_level] = lastBrcExtractor(node);

            // Measure the time taken for the pricing operation on the node.
            auto eps = TimeSetter::measure([&]() {
                pricing(node);
            });

            // If the node is not terminated and is not the root (tree_level != 0), record branching information.
            if (!isTerminate(node) && tree_level != 0) {
                branching_history.recordExactPerScore(edge, val, valueExtractor(node), dir, tree_level);
                brcValueOutput(node) = branching_history.tellBranchingIncreaseVal(tree_level);
            }

            // Measure additional time taken for the cutting operation.
            eps += TimeSetter::measure([&]() {
                cutting(node);
            });

            // If the node remains active (not terminated), proceed with branching.
            if (!isTerminate(node)) {
                // For non-root nodes, update the r_star value via LP testing if BKF controllers are available.
                if (tree_level != 0) {
                    if (!bkf_controllers.empty()) {
                        bkf_data_shared.calculateRStar(valueExtractor(node) - val, tree_level, dir, idxExtractor(node),
                                                       bkf_controllers[0]); // Estimate r_star using LP testing.
                    }
                }

                // Determine the best branching candidate.
                BrCType brc;
                std::vector<int> bst_ks(bkf_controllers.size());
                for (int i = 0; i < bst_ks.size(); ++i) {
                    // Set the BKF controller state based on the current enumeration state.
                    bkf_controllers[i].setIfB4(!enuStateExtractor(node));
                    // Retrieve the best 'K' value from the BKF controller.
                    bst_ks[i] = bkf_controllers[i].getBestK(bkf_data_shared, ub_ref.get(), valueExtractor(node));

                    // Record the best K value for each phase (supports up to 3 phases).
                    if (i == 0)
                        branching_tester.setNumPhase0(bst_ks[i]);
                    else if (i == 1)
                        branching_tester.setNumPhase1(bst_ks[i]);
                    else if (i == 2)
                        branching_tester.setNumPhase2(bst_ks[i]);
                    else
                        THROW_RUNTIME_ERROR("BKFController only supports 3 phases, but got " + std::to_string(i));
                }

                // Select the branching candidate using either a self-defined function or the default tester.
                if (getBestCandidateBySelfDefined) {
                    brc = getBestCandidateBySelfDefined(node, branching_history, branching_data_shared,
                                                        branching_tester);
                } else {
                    brc = branching_tester.getBestCandidate(node, branching_history, branching_data_shared,
                                                            getBranchingCandidates(node));
                }

                // Update BKF timing information and state for each BKF controller.
                if (!bkf_controllers.empty()) {
                    branching_tester.updateBKFtime(eps, bkf_controllers);
                    for (auto &bkf: bkf_controllers) {
                        bkf.updateTimeMeasure(bkf_data_shared);
                    }
                    // Update node counters based on enumeration state.
                    enuStateExtractor(node) ? bkf_data_shared.increaseNodeAf() : bkf_data_shared.increaseNodeB4();
                }

                // Increment the branch choice count for the selected candidate.
                ++branching_history.branch_choice[brc];

                // Impose the branching decision on the current node and generate its children.
                std::vector<Node *> children;
                imposeBranching(node, brc, children);

                // Process each child node.
                for (auto &child: children) {
                    // If an output callback is defined, attempt to write the node out.
                    if (child != node && tryWriteNodeOut != nullptr) {
                        tryWriteNodeOut(child, branching_history, bkf_data_shared);
                        // If the child node is terminated after writing, delete it.
                        if (isTerminate(child)) {
                            delete child;
                            continue;
                        }
                    }
                    // Add the child node to the tree for further exploration.
                    addNode(child);
                    // Update the optimal K values in BKF controllers for the new child node.
                    for (int i = 0; i < bst_ks.size(); ++i) {
                        bkf_controllers[i].updateOptK(bst_ks[i], idxExtractor(child));
                    }
                }
            } else {
                // If the node is terminated, update BKF shared data and delete the node.
                if (!bkf_controllers.empty()) {
                    bkf_data_shared.updateF(ub_ref.get() - valueExtractor(node));
                }
                delete node;
            }
            // Increment the counter for the number of nodes explored.
            ++num_nodes_explored;
            // Update the global bounds based on the current tree state.
            updateBounds();
        }
    }
} // namespace RouteOpt::Branching::BBT

#endif // ROUTE_OPT_WORKFLOW_HPP
