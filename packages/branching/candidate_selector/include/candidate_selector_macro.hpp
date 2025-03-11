/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file candidate_selector_macro.hpp
 * @brief Macros and constants for candidate selection in RouteOpt.
 *
 * This header defines constants, macros, and an enumeration class used in the candidate selection process.
 * It includes settings for invalid indices, score adjustment ratios, and console output formatting.
 */

#ifndef ROUTE_OPT_CANDIDATE_SELECTOR_MACRO_HPP
#define ROUTE_OPT_CANDIDATE_SELECTOR_MACRO_HPP

namespace RouteOpt::Branching::CandidateSelector {
    /**
     * @brief Indicator for an invalid index.
     *
     * This constant is used to denote an index that is invalid or uninitialized.
     */
    constexpr int INVALID_INDEX = -1;

    /**
     * @brief Ratio used for adjusting scores in candidate selection.
     *
     * This constant defines the multiplier applied when adjusting candidate scores.
     * It helps scale the scores appropriately during comparisons.
     */
    constexpr double ADJUSTMENT_SCORE_RATIO = 3;

    /**
     * @brief Print column width for console output.
     *
     * This constant specifies the width of each column when candidate data is
     * printed in a tabular format.
     */
    constexpr int PRINT_COL_WIDTH = 10;

    /**
     * @brief Number of columns to print in console output.
     *
     * This constant defines the number of columns displayed per row in the printed output.
     */
    constexpr int PRINT_NUM_COLS = 5;

    /**
     * @brief Enumeration of testing phases for candidate selection.
     *
     * This enumeration defines the different phases used during the candidate selection process:
     * - LP: Testing using a linear programming (LP) approach.
     * - Heuristic: Testing based on LP + heuristic Column Generation (CG).
     * - Exact: Testing based on LP + exact CG.
     */
    enum class TestingPhase {
        LP,
        Heuristic,
        Exact
    };

    // Define CANDIDATE_SELECTOR_VERBOSE to enable verbose debugging for candidate selection.
#define CANDIDATE_SELECTOR_VERBOSE

#ifdef CANDIDATE_SELECTOR_VERBOSE
    /**
     * @brief Macro to execute verbose candidate selector code.
     *
     * When CANDIDATE_SELECTOR_VERBOSE is defined, the code passed to this macro is executed.
     * This is useful for inserting debugging statements or additional logging.
     */
#define CANDIDATE_SELECTOR_VERBOSE_EXEC(...) __VA_ARGS__;
#else
    // If verbose mode is not enabled, the macro does nothing.
    CANDIDATE_SELECTOR_VERBOSE_EXEC();
#endif
} // namespace RouteOpt::Branching::CandidateSelector

#endif // ROUTE_OPT_CANDIDATE_SELECTOR_MACRO_HPP
