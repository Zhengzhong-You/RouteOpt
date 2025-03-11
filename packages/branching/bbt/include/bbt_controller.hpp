/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file bbt_controller.hpp
 * @brief Header file for the BBTController class template.
 *
 * This header defines the BBTController class template, which manages the branch-and-bound
 * tree for the optimization solver using a branching strategy. It integrates various
 * components such as BKF controllers, candidate selectors, and testing functions.
 */

#ifndef ROUTE_OPT_BBT_CONTROLLER_HPP
#define ROUTE_OPT_BBT_CONTROLLER_HPP

#include <set>
#include "bkf_controller.hpp"
#include "candidate_selector_controller.hpp"

namespace RouteOpt::Branching::BBT {
    /**
     * @brief BBTController class template for branch-and-bound tree management.
     *
     * This class template manages the branch-and-bound tree in a branch-and-bound algorithm.
     * It coordinates various components (like BKF controllers and candidate selectors) to
     * facilitate the branching selection process and maintain the tree structure.
     *
     * @tparam Node The node type representing a branch-and-bound tree node.
     * @tparam BrCType The branching candidate type.
     * @tparam Hasher A hash functor for the branching candidate type.
     */
    template<typename Node, typename BrCType, typename Hasher>
    class BBTController {
    public:
        /**
         * @brief Constructs a BBTController with the provided parameters.
         *
         * Initializes the controller with instance dimensions, phases, BKF sizes, and various callback
         * functions used for node comparisons, value extraction, branching candidate extraction,
         * and testing functions.
         *
         * @param dim Instance Dimension.
         * @param ub Reference to the upper bound value.
         * @param num_phase0 Number of phase 0 outputs, also the number of LP testings.
         * @param num_phase1 Number of phase 1 outputs, also the number of heuristic testings.
         * @param num_phase2 Number of phase 2 outputs, also the number of exact testings.
         * @param num_phase3 Number of phase 3 outputs, typically be 1.
         * @param bkf_size Vector of pairs indicating sizes for each BKF controller.
         * @param nodeComparator Function to compare two nodes.
         * @param valueExtractor Function to extract local lower bound from a node.
         * @param idxExtractor Function to extract an index from a node.
         * @param enuStateExtractor Function to extract enumeration state flag from a node.
         * @param lastBrcExtractor Function to extract last branching constraint info (BrC, Branching direction , Local tree size) from a node.
         * @param brcValueOutput Function to output a reference of exact improvement of LP value from the node to its parent node.
         * @param isTerminate Function to check if the node has been marked as terminated.
         * @param getBranchingCandidates Function to retrieve branching candidates from a node.
         * @param processLPTestingFunction Function to process LP testing.
         * @param processHeurTestingFunction Function to process heuristic testing.
         * @param processExactTestingFunction Function to process exact testing.
         * @param pricing Function to perform pricing operations.
         * @param cutting Function to perform cutting operations.
         * @param imposeBranching Function to impose branching decisions on a node.
         * @param getBestCandidateBySelfDefined (Optional) Custom function to choose the best branching candidate.
         * @param tryWriteNodeOut (Optional) Function to output node details externally.
         * @param tryReadNodeIn (Optional) Function to read node details from an external source.
         */
        BBTController(
            int dim,
            double &ub,
            int num_phase0,
            int num_phase1,
            int num_phase2,
            int num_phase3,
            const std::vector<std::pair<int, int> > &bkf_size,
            const std::function<bool(Node *, Node *)> &nodeComparator,
            const std::function<double(Node *)> &valueExtractor,
            const std::function<int(Node *)> &idxExtractor,
            const std::function<bool(Node *)> &enuStateExtractor,
            const std::function<std::tuple<BrCType, bool, int>(Node *)> &lastBrcExtractor,
            std::function<double&(Node *)> brcValueOutput,
            const std::function<bool(Node *)> &isTerminate,
            const std::function<std::unordered_map<BrCType, double, Hasher>(Node *)> &getBranchingCandidates,
            const std::function<void(Node *, const BrCType &, double &, double &)> &processLPTestingFunction,
            const std::function<void(Node *, const BrCType &, double &, double &)> &processHeurTestingFunction,
            const std::function<void(Node *, const BrCType &, double &, double &)> &processExactTestingFunction,
            const std::function<void(Node *)> &pricing,
            const std::function<void(Node *)> &cutting,
            const std::function<void(Node *, const BrCType &, std::vector<Node *> &)> &imposeBranching,
            const std::function<BrCType(Node *, BranchingHistory<BrCType, Hasher> &,
                                        BranchingDataShared<BrCType, Hasher> &,
                                        CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &)> &
                    getBestCandidateBySelfDefined = nullptr,
            const std::function<void(Node *,
                                     const BranchingHistory<BrCType, Hasher> &,
                                     const BKF::BKFDataShared &)> &tryWriteNodeOut = nullptr,
            const std::function<void(Node *,
                                     BranchingHistory<BrCType, Hasher> &,
                                     BKF::BKFDataShared &)> &tryReadNodeIn = nullptr)
            : branching_data_shared(dim),
              branching_tester(num_phase0, num_phase1, num_phase2, num_phase3,
                               processLPTestingFunction,
                               processHeurTestingFunction,
                               processExactTestingFunction),
              pricing(pricing),
              cutting(cutting),
              imposeBranching(imposeBranching),
              getBranchingCandidates(getBranchingCandidates),
              getBestCandidateBySelfDefined(getBestCandidateBySelfDefined),
              tryWriteNodeOut(tryWriteNodeOut),
              tryReadNodeIn(tryReadNodeIn),
              isTerminate(isTerminate),
              valueExtractor(valueExtractor),
              idxExtractor(idxExtractor),
              enuStateExtractor(enuStateExtractor),
              lastBrcExtractor(lastBrcExtractor),
              brcValueOutput(brcValueOutput),
              tree(nodeComparator),
              ub_ref(ub) {
            // Initialize BKF controllers based on provided sizes.
            bkf_controllers.resize(bkf_size.size());
            for (int i = 0; i < bkf_size.size(); ++i) {
                bkf_controllers[i].setMN(bkf_size[i].first, bkf_size[i].second);
            }
        }

        /**
         * @brief Solve the branch-and-bound tree starting from the root node.
         *
         * This method initiates the solving process by exploring nodes, updating bounds, and
         * applying branching decisions until termination criteria are met.
         *
         * @param root_node Pointer to the root node of the branch-and-bound tree.
         */
        void solve(Node *root_node);

        /**
         * @brief Get the total number of nodes explored during the solving process.
         *
         * @return Number of nodes explored.
         */
        [[nodiscard]] int getNumNodesExplored() const {
            return num_nodes_explored;
        }

        /**
         * @brief Get the current number of nodes in the branch-and-bound tree.
         *
         * @return Total number of nodes in the tree.
         */
        [[nodiscard]] int getNumNodes() const {
            return static_cast<int>(tree.size());
        }

        /**
         * @brief Retrieve the current lower bound on the objective function.
         *
         * @return Current lower bound value.
         */
        [[nodiscard]] double getLowerBound() const {
            return lb;
        }

        // Default constructor.
        BBTController() = default;

        // Default destructor.
        ~BBTController() = default;

    private:
        int num_nodes_explored{}; ///< Counter for the number of nodes explored.

        // Shared data for branching operations.
        BranchingDataShared<BrCType, Hasher> branching_data_shared;

        // History of branching decisions.
        BranchingHistory<BrCType, Hasher> branching_history{};

        // Tester for branching candidate selection.
        CandidateSelector::BranchingTesting<Node, BrCType, Hasher> branching_tester;

        // Shared data for the BKF component.
        BKF::BKFDataShared bkf_data_shared{};

        // Collection of BKF controllers for managing BKF-specific tasks.
        std::vector<BKF::BKFController> bkf_controllers{};

        // Function for performing pricing operations.
        std::function<void(Node *)> pricing{};

        // Function for performing cutting operations.
        std::function<void(Node *)> cutting{};

        // Function to impose branching decisions on a node.
        std::function<void(Node *, const BrCType &, std::vector<Node *> &)> imposeBranching{};

        // Function to get branching candidates from a node.
        std::function<std::unordered_map<BrCType, double, Hasher>(Node *)> getBranchingCandidates{};

        // Optional custom function to select the best branching candidate.
        std::function<BrCType(Node *, BranchingHistory<BrCType, Hasher> &, BranchingDataShared<BrCType, Hasher> &,
                              CandidateSelector::BranchingTesting<Node, BrCType, Hasher> &)>
        getBestCandidateBySelfDefined{};

        // Optional function to write node details to external storage.
        std::function<void(Node *,
                           const BranchingHistory<BrCType, Hasher> &,
                           const BKF::BKFDataShared &)> tryWriteNodeOut{};

        // Optional function to read node details from external storage.
        std::function<void(Node *,
                           BranchingHistory<BrCType, Hasher> &,
                           BKF::BKFDataShared &)> tryReadNodeIn{};

        // Function to check if termination criteria have been met.
        std::function<bool(Node *)> isTerminate;

        // Function to extract a value from a node (e.g., objective value).
        std::function<double(Node *)> valueExtractor;

        // Function to extract an index from a node.
        std::function<int(Node *)> idxExtractor;

        // Function to extract the enumeration state from a node.
        std::function<bool(Node *)> enuStateExtractor;

        // Function to extract branching candidate information from a node.
        std::function<std::tuple<BrCType, bool, int>(Node *)> lastBrcExtractor;

        // Function to obtain a reference to a branching candidate value from a node.
        std::function<double&(Node *)> brcValueOutput;

        // Set representing the branch-and-bound tree, ordered by the custom comparator.
        std::set<Node *, std::function<bool(Node *, Node *)> > tree;

        double lb{}; ///< Lower bound on the objective function.
        std::reference_wrapper<double> ub_ref; ///< Reference to the upper bound value.

        /**
         * @brief Add a node to the branch-and-bound tree.
         *
         * @param node Pointer to the node to be added.
         */
        void addNode(Node *node) {
            tree.insert(node);
        }

        /**
         * @brief Update the lower bound based on the nodes present in the tree.
         *
         * If the tree is empty, the lower bound is set to the upper bound; otherwise, it is updated
         * to the minimum value extracted from the nodes in the tree.
         */
        void updateBounds() {
            if (tree.empty()) {
                lb = ub_ref.get();
            } else {
                double minVal = std::numeric_limits<double>::infinity();
                for (auto &nodePtr: tree) {
                    double val = valueExtractor(nodePtr);
                    if (val < minVal) {
                        minVal = val;
                    }
                }
                lb = minVal;
            }
        }

        /**
         * @brief Print detailed information about the current branch-and-bound tree.
         *
         * This function outputs the current tree size, lower bound, and upper bound.
         */
        void printTreeDetail() {
            std::cout << BIG_PHASE_SEPARATION;
            std::cout << "tree size= " << tree.size() << " lb= " << lb << " ub= " << ub_ref.get()
                    << std::endl;
            std::cout << BIG_PHASE_SEPARATION;
        }
    };
} // namespace RouteOpt::Branching::BBT

#include "workflow.hpp"
#endif // ROUTE_OPT_BBT_CONTROLLER_HPP
