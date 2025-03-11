/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rank1_separation_controller.hpp
 * @brief Rank-1 Separation Controller for handling the generation, updating, and processing of Rank-1 cuts.
 *
 * This header defines the Rank1SeparationController class which coordinates the separation process for Rank-1 cuts.
 * It integrates shared Rank-1 cuts data, cut generation, memory finding, and cut selection modules,
 * and interacts with the solver and the vertex cost matrix.
 */

#ifndef ROUTE_OPT_RANK1_SEPARATION_CONTROLLER_HPP
#define ROUTE_OPT_RANK1_SEPARATION_CONTROLLER_HPP

#include "rank1_data_shared.hpp"
#include "rank1_cuts_generator.hpp"
#include "rank1_memory_finder.hpp"
#include "rank1_cuts_selector.hpp"
#include "solver.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    /**
     * @brief Controller for the Rank-1 separation process.
     *
     * The Rank1SeparationController class manages the separation of Rank-1 cuts.
     * It utilizes several modules including shared Rank-1 cuts data, cut generation,
     * memory finding, and cut selection. The controller updates internal data,
     * generates new cuts, converts cuts to an arc memory representation, and interfaces
     * with the underlying solver.
     */
    class Rank1SeparationController {
    public:
        /**
         * @brief Constructs a Rank1SeparationController.
         *
         * @param rank1CutsDataShared Reference to the shared Rank-1 cuts data (used for retrieving multiplier plans).
         * @param max_row_rank1 Maximum number of Rank-1 rows (dimension) allowed.
         * @param max_num_r1c3_per_round Maximum number of R1C3 cuts to be generated per round.
         * @param max_num_r1c_per_round Maximum number of R1C cuts (dimension!=3) to be generated per round.
         * @param solver Reference to the IP solver.
         * @param cost_mat4_vertex Cost matrix associated with vertices.
         */
        Rank1SeparationController(const Rank1CutsDataShared &rank1CutsDataShared,
                                  int max_row_rank1,
                                  int max_num_r1c3_per_round,
                                  int max_num_r1c_per_round,
                                  Solver &solver,
                                  const std::vector<std::vector<double> > &cost_mat4_vertex);

        /**
         * @brief Updates the controller's internal information.
         *
         * This method updates the internal state with new information including limited memory type,
         * pricing difficulty level, solution collection flag, current solution routes, and existing cuts.
         *
         * @param limited_memory_type Memory type limitation used during processing.
         * @param pricing_hard_level Current pricing difficulty level.
         * @param if_collect_sol Flag indicating whether solutions should be collected.
         * @param sol Vector of RouteInfo objects representing current solution routes.
         * @param old_cuts Vector of existing Rank-1 cuts.
         */
        void updateInfo(MemoryType limited_memory_type,
                        PRICING_HARD_LEVEL pricing_hard_level,
                        bool if_collect_sol,
                        const std::vector<RouteInfo> &sol,
                        const std::vector<R1c> &old_cuts);

        /**
         * @brief Performs the separation of Rank-1 cuts.
         *
         * Generates new Rank-1 cuts based on current data. Optionally applies memory searching
         * and cut selection strategies.
         *
         * @param cuts [in,out] Vector of Rank-1 cuts to be updated or augmented.
         * @param if_mem Flag indicating whether to apply memory searching.
         * @param if_select_cuts Flag indicating whether to perform cut selection.
         */
        void separateRank1Cuts(std::vector<R1c> &cuts,
                               bool if_mem = true, bool if_select_cuts = true);

        /**
         * @brief Converts existing Rank-1 cuts with node memory representation into an arc memory representation.
         *
         * This method processes the provided cuts and converts them into an arc-based format,
         * facilitating pricing problem.
         *
         * @param existing_cuts Vector of existing Rank-1 cuts to be converted.
         */
        void convert2ArcMemory(std::vector<R1c> &existing_cuts);

        /**
         * @brief Retrieves the flag indicating whether the "no symmetry memory" has been used.
         *
         * @return auto Boolean flag value.
         */
        [[nodiscard]] auto getIfOnceUseNoSymmetryMem() const {
            return if_once_use_no_symmetry_mem;
        }

        // Delete default constructor to enforce the use of the parameterized constructor.
        Rank1SeparationController() = delete;

        // Default destructor.
        ~Rank1SeparationController() = default;

    private:
        // Vector storing solution routes, each represented as a vector of RouteInfo.
        std::vector<std::vector<RouteInfo> > sol_vec{};

        // Reference to shared Rank-1 cuts data (used for retrieving multiplier plans).
        std::reference_wrapper<const Rank1CutsDataShared> rank1CutsDataShared_ref;
        // Shared data structure for managing common information in the separation process.
        DataShared sharedData;
        // Module responsible for generating Rank-1 cuts.
        CutGenerator cutGen;
        // Module for finding memory.
        MemGenerator memGen;
        // Module for selecting the most effective cuts from the generated cuts.
        CutSelector cutSelector;

        // Flag indicating whether the "no symmetry memory" has been used once.
        bool if_once_use_no_symmetry_mem{false};

        /**
         * @brief Cleans internal data in preparation for a new round of processing.
         *
         * This private method resets the state of shared data, cut generation,
         * memory finding, and cut selection modules.
         */
        void cleanData() {
            sharedData.cleanData();
            cutGen.cleanData();
            memGen.cleanData();
            cutSelector.cleanData();
        }

        /**
         * @brief Sets the right-hand side (rhs) values for a vector of Rank-1 cuts.
         *
         * This method updates the rhs of each cut based on current internal computations
         * and data.
         *
         * @param cuts Vector of Rank-1 cuts whose rhs values are to be updated.
         */
        void setRhs(std::vector<R1c> &cuts);
    };
} // namespace RouteOpt::Rank1Cuts::Separation

#endif // ROUTE_OPT_RANK1_SEPARATION_CONTROLLER_HPP
