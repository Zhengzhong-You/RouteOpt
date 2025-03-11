/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
*/

/*
 * @file rank1_separation_macro.hpp
 * @brief Macro definitions and type aliases for Rank-1 separation in RouteOpt.
 *
 * This header defines various constants, type aliases, and helper macros that are used
 * in the Rank-1 separation process. These definitions include parameters for tolerances,
 * memory settings, and other generator parameters, as well as custom types based on the Eigen library.
 */

#ifndef ROUTE_OPT_RANK1_SEPARATION_MACRO_HPP
#define ROUTE_OPT_RANK1_SEPARATION_MACRO_HPP

#include <bitset>
#include <Eigen/Sparse>
#include <vector>
#include "route_opt_macro.hpp"

namespace RouteOpt::Rank1Cuts::Separation {
    // Verbosity flag for Rank-1 separation debugging.
    constexpr bool RANK1_SEPARATION_VERBOSE = false;

    // Tolerance for solution x values in Rank-1 separation.
    constexpr double SOL_X_RANK1_TOLERANCE = 1e-3;
    // Factor to control the maximum memory size for cut generation.
    constexpr double MAX_CUT_MEM_FACTOR = 0.15;
    // Initial size for the pool of Rank-1 multi-labels.
    constexpr int INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE = 2048;
    // Threshold value to decide whether to use enumeration or MIP (Mixed-Integer Programming) for memory finding.
    constexpr int FIND_MEM_USE_ENUMERATION_OR_MIP = 1000;
    // Maximum number of segments allowed for a single column in the separation process.
    constexpr int MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN = 16;
    // Maximum number of undetermined arcs.
    constexpr int MAX_UNDETERMINED_ARC_NUMBER = 128;
    // Maximum number of labels allowed.
    constexpr int MAX_LABELS = 1024;
    // Maximum number of neighbor customers for separating Rank-1 cuts.
    constexpr int MAX_HEURISTIC_SEP_ROW_RANK1 = 16;
    // Minimum number of cuts to be added per separation round.
    constexpr int MIN_CUTS_ADDED_PER_ROUND = 50;
    // Time limit (in seconds) for MIP-based memory finding.
    constexpr double TIME_LIMIT_FOR_MIP_FIND_MEM = 0.2;
    // Violation factor for cuts.
    constexpr double CUT_VIO_FACTOR = 0.1;
    // Density estimation for Rank-1 cuts.
    constexpr double RANK1_ROW_DENSITY = 0.1;
    // Limit for keeping solutions in memory.
    constexpr int SOL_KEEPER_LIMIT = 20;

    /**
     * @brief Enumeration for memory type used in Rank-1 separation.
     *
     * This enum defines the memory mode for the separation process:
     * - NO_MEMORY: No memory is used.
     * - NODE_MEMORY: Node-based memory is utilized.
     * - ARC_MEMORY: Arc-based memory is applied.
     */
    enum class MemoryType {
        NO_MEMORY = 0,
        NODE_MEMORY = 1,
        ARC_MEMORY = 2
    };

    // Forward declarations for types used in the separation process.
    struct Rank1MultiLabel;
    struct RouteInfo;

    /**
     * @brief Custom hash functor for a vector of integers.
     *
     * This structure provides a hash function for std::vector<int>,
     * which can be used in unordered containers.
     */
    struct VectorHashInRank1 {
        size_t operator()(const std::vector<int> &v) const {
            std::size_t hash = 0;
            for (const int num: v) {
                hash ^= std::hash<int>{}(num) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };

    // Type alias for a bitset representing undetermined arcs.
    using arcBit = std::bitset<MAX_UNDETERMINED_ARC_NUMBER>;
    // Alias for the long bitset type used in cut representations.
    using cutLong = routeOptLong;
    // Alias for a row-major sparse matrix of integers using the Eigen library.
    using sparseRowMatrixXI = Eigen::SparseMatrix<int, Eigen::RowMajor>;
    // Alias for a row-major sparse vector of integers using the Eigen library.
    using sparseRowVectorXI = Eigen::SparseVector<int, Eigen::RowMajor>;

    /**
     * @brief Structure containing route information.
     *
     * The RouteInfo structure encapsulates key information for a route in the context
     * of Rank-1 separation. It stores a fractional value, the sequence of column indices,
     * and the position for forward concatenation.
     */
    struct RouteInfo {
        // Fractional value associated with the route (e.g., LP solution value).
        double frac_x{};
        // Sequence of customer indices representing the route, excluding 0.
        std::vector<int> col_seq{};
        // Position for forward concatenation in the route, e.g., col_seq= 1,2,3 and forward_concatenate_pos=1, means
        // this route is concated by 0-1-2 (forward) and 0-3 (backward).
        int forward_concatenate_pos{};
    };

    /**
     * @brief Macro for conditionally executing a statement based on verbosity.
     *
     * This macro executes the provided statement only if RANK1_SEPARATION_VERBOSE is set to true.
     *
     * @param stmt The statement to execute conditionally.
     */
#define RANK1_VERBOSE_EXEC(stmt) do { \
    if (RANK1_SEPARATION_VERBOSE) { \
        stmt; \
    } \
} while (0);
} // namespace RouteOpt::Rank1Cuts::Separation

#endif // ROUTE_OPT_RANK1_SEPARATION_MACRO_HPP
