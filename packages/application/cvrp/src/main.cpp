/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file main.cpp
 * @brief Main entry point for the CVRP application.
 *
 * This file contains the main function that initializes and runs the CVRP, or VRPTW solver.
 * It sets up the branch-and-bound tree controller and handles the overall solving process.
 */

#include "cvrp.hpp"                // CVRP solver header.
#include "vrptw.hpp"               // VRPTW solver header.
#include "cvrp_macro.hpp"          // Macros specific to CVRP.
#include "bbt_controller.hpp"      // Branch-and-bound tree controller.
using namespace std;
using namespace RouteOpt;
using namespace Application::CVRP;

int main(int argc, char *argv[]) {
    // Create an instance of the solver based on the application type.
    // If the application type is CVRP, create a CVRPSolver; otherwise, create a VRPTW solver.

    std::cout << "application type: " <<
            (app_type == APPLICATION_TYPE::CVRP ? "CVRP" : "VRPTW") << std::endl;


    CVRPSolver *cvrp = app_type == APPLICATION_TYPE::CVRP
                           ? new CVRPSolver(argc, argv)
                           : new VRPTW(argc, argv);
    // Create a new branch-and-bound node (BbNode) as the root node.
    auto node = new BbNode();

    // Define a candidate selection function pointer for machine learning based candidate selection.
    // Initialize to nullptr by default.
    CandidateSelectionFuncType ml_candidate_selection = nullptr;
    // If machine learning is enabled (ml_type is not ML_NO_USE), set the candidate selection function.
    if constexpr (ml_type != ML_TYPE::ML_NO_USE) {
        ml_candidate_selection = [cvrp](auto arg1, auto &arg2, auto &arg3, auto &arg4) -> std::pair<int, int> {
            return cvrp->callMLCandidateSelection(arg1, arg2, arg3, arg4);
        };
    }

    // Define an output function pointer for writing nodes.
    OutNodeFuncType node_out_func = nullptr;
    // If the configuration enables writing node outputs, set the function accordingly.
    if constexpr (IF_WRITE_NODE_OUT) {
        node_out_func = [cvrp](auto arg1, auto &arg2, auto &arg3) -> void {
            cvrp->callWriteNodeOut(arg1, arg2, arg3);
        };
    }

    // Define a function for reading node information.
    auto node_in_func = [cvrp](auto arg1, auto &arg2, auto &arg3) -> void {
        cvrp->callReadNodeIn(arg1, arg2, arg3);
    };

    // Create an instance of the Branch-and-Bound Tree (BBT) Controller.
    // The controller is templated with:
    //   - BbNode: the type for branch-and-bound nodes.
    //   - std::pair<int, int>: the type for branching candidates.
    //   - PairHasher: hash function for branching candidates.
    Branching::BBT::BBTController<BbNode, std::pair<int, int>, PairHasher> bbt_controller{
        cvrp->getDim(), // Problem/Instance dimension.
        cvrp->refUB(), // Reference to the current upper bound.
        static_cast<int>(NUM_TESTING::PHASE0), // Number of output candidates for phase 0.
        static_cast<int>(NUM_TESTING::PHASE1), // Number of output candidates for phase 1.
        static_cast<int>(NUM_TESTING::PHASE2), // Number of output candidates for phase 2.
        static_cast<int>(NUM_TESTING::PHASE3), // Number of output candidates for phase 3.
        IF_BKF
            ? std::vector<std::pair<int, int> >{
                {static_cast<int>(BKF_TYPE::M_LP), static_cast<int>(BKF_TYPE::N_LP)},
                {static_cast<int>(BKF_TYPE::M_HEUR), static_cast<int>(BKF_TYPE::N_HEUR)}
            }
            : std::vector<std::pair<int, int> >{}, // BKF parameters if BKF is enabled.
        BbNode::defineBetterNode, // Function to compare nodes.
        BbNode::getNodeValue, // Function to extract node value.
        BbNode::getNodeIdx, // Function to extract node index.
        BbNode::getNodeIfInEnumState, // Function to check node enumeration state.
        BbNode::getLastBrc, // Function to get the last branching candidate.
        BbNode::refNodeBrIncreaseVal, // Function to reference node branching increase value.
        BbNode::getNodeIfTerminate, // Function to check if node is terminated.
        BbNode::obtainSolEdgeMap, // Function to obtain solution edge map.
        // LP testing function for processing LP testing.
        [cvrp](auto arg1, auto &arg2, auto &arg3, auto &arg4) -> decltype(auto) {
            return cvrp->processLPTesting(arg1, arg2, arg3, arg4);
        },
        // Heuristic testing function.
        [cvrp](auto arg1, auto &arg2, auto &arg3, auto &arg4) -> decltype(auto) {
            return cvrp->processCGTesting<false>(arg1, arg2, arg3, arg4);
        },
        // Exact testing function.
        [cvrp](auto arg1, auto &arg2, auto &arg3, auto &arg4) -> decltype(auto) {
            return cvrp->processCGTesting<true>(arg1, arg2, arg3, arg4);
        },
        // Function to perform pricing at the beginning of processing a node.
        [cvrp](auto arg1) -> decltype(auto) {
            return cvrp->callPricingAtBeg(arg1);
        },
        // Function to perform cutting on a node.
        [cvrp](auto arg1) -> decltype(auto) {
            return cvrp->callCutting(arg1);
        },
        // Function to impose a branching decision on a node, generating child nodes.
        [cvrp](auto arg1, auto arg2, auto &arg3) -> decltype(auto) {
            return cvrp->imposeBranching(arg1, arg2, arg3);
        },
        // Self-Defined Branching Selection Function (e.g., 2LBB).
        ml_candidate_selection,
        // Node output function for writing nodes to external storage.
        node_out_func,
        // Node input function for reading node information.
        node_in_func
    };

    // Solve the branch-and-bound tree starting from the root node.
    bbt_controller.solve(node, TIME_LIMIT);

    // After solving, print the optimal solution along with statistics:
    // number of nodes explored and the lower bound of the solution.
    cvrp->printOptSol(std::cout, bbt_controller.getNumNodesExplored(), bbt_controller.getLowerBound());

    // Clean up: delete the solver instance.
    delete cvrp;

    return 0;
}
