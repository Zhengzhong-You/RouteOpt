/*
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file route_opt_macro.hpp
 * @brief Macros, constants, and utility functions for RouteOpt.
 *
 * This header provides definitions for constants (such as maximum customer count and tolerances),
 * utility functions (e.g., for directory creation and timing), custom data structures (e.g., SequenceInfo),
 * and logging macros for the RouteOpt project.
 */

#ifndef ROUTE_OPT_ROUTE_OPT_MACRO_HPP
#define ROUTE_OPT_ROUTE_OPT_MACRO_HPP

// Standard library headers for bit manipulation, math functions, I/O, and filesystem operations.
#include <bitset>
#include <cmath>
#include <vector>
#include <string>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace RouteOpt {
    // Maximum number of customers. This is used to size bitsets and define limits in the application.
    constexpr int MAX_NUM_CUSTOMERS = 1002;

    // Tolerance values for numerical comparisons.
    constexpr double TOLERANCE = 1e-6;
    constexpr double SOL_X_TOLERANCE = 1e-8;
    constexpr double RC_TOLERANCE = -1e-5;
    constexpr double DUAL_TOLERANCE = 1e-6;

    // Define a bitset type for representing customer routes with a fixed maximum size.
    using routeOptLong = std::bitset<MAX_NUM_CUSTOMERS>;

    // Separation lines for formatted console output.
    constexpr auto SMALL_PHASE_SEPARATION =
            "-----------------------------------------------------------------------------------------\n";
    constexpr auto MID_PHASE_SEPARATION =
            "*****************************************************************************************\n";
    constexpr auto BIG_PHASE_SEPARATION =
            "#########################################################################################\n";

    /**
     * @brief Custom hash functor for a pair of integers.
     *
     * This structure provides a hash function for std::pair<int, int> objects,
     * which is useful when using pairs as keys in unordered containers.
     */
    struct PairHasher {
        size_t operator()(const std::pair<int, int> &V) const {
            return V.first * MAX_NUM_CUSTOMERS + V.second;
        }
    };

    /**
     * @brief Structure representing sequence information.
     *
     * Contains a vector of integers representing a forward sequence (excluding 0)
     * and a position indicator for the last processed element in a forward partial label.
     */
    struct SequenceInfo {
        std::vector<int> col_seq; // Forward sequence of customer indices (excludes 0).
        /*
         * Records the last position of the forward partial label.
         * For example, if the route is [0, 2, 3, 4], and forward partial label is [0, 2], then forward_concatenate_pos = 1.
         */
        int forward_concatenate_pos{};
        // Equality operator comparing both the sequence and the concatenation position.
        bool operator==(const SequenceInfo &other) const {
            return col_seq == other.col_seq &
                   forward_concatenate_pos == other.forward_concatenate_pos;
        }
    };

    /**
     * @brief Compares two double values within a specified tolerance.
     *
     * @param x First value.
     * @param y Second value.
     * @param tolerance Acceptable difference between x and y (default is TOLERANCE).
     * @return true if the difference is less than the tolerance, false otherwise.
     */
    constexpr bool equalFloat(double x, double y, double tolerance = TOLERANCE) {
        return std::fabs(x - y) < tolerance;
    }

    /**
     * @brief Checks if a number has a fractional component.
     *
     * @param x The number to check.
     * @return true if x is not an integer (has a fractional part), false otherwise.
     */
    constexpr bool checkFrac(double x, double tolerance = TOLERANCE) {
        return !equalFloat(std::fmod(x, 1.), 0, tolerance);
    }

    /**
     * @brief Creates a directory if it does not exist.
     *
     * This function checks if a directory exists at the given path; if not, it creates the directory.
     *
     * @param path1 The path of the directory.
     */
    inline void mkDir(const std::string &path1) {
        namespace fs = std::filesystem;
        if (fs::exists(path1) && fs::is_directory(path1)) {
            std::cout << path1 + " already exists" << std::endl;
        } else {
            fs::create_directory(path1);
            std::cout << path1 + " created" << std::endl;
        }
    }

    /**
     * @brief Utility class to measure the execution time of functions.
     *
     * Provides a templated static function to time any callable object's execution.
     */
    struct TimeSetter {
        template<typename Func, typename... Args>
        static double measure(Func f, Args &&... args) {
            auto start = std::chrono::high_resolution_clock::now();
            f(std::forward<Args>(args)...);
            auto end = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double>(end - start).count();
        }
    };

    /**
     * @brief Macro to safely execute an Eigen function.
     *
     * Wraps a function call in a try-catch block to catch and report any exceptions,
     * then exits the program in case of an error.
     */
#define SAFE_EIGEN(func) \
try { \
    func; \
} catch (const std::exception &e) { \
    std::ostringstream oss; \
    oss << "\x1b[91mError in file " << __FILE__ \
        << " at line " << __LINE__ << ": " << e.what() << "\x1b[0m"; \
    std::cerr << oss.str() << std::endl; \
    std::exit(EXIT_FAILURE); \
}

    /**
     * @brief Prints a formatted message with the elapsed time.
     *
     * @param message A custom message to print.
     * @param time The elapsed time in seconds.
     */
    inline void printTimeMessage(const std::string &message, double time) {
        std::ios init(nullptr);
        init.copyfmt(std::cout);
        std::cout << message << std::fixed << std::setprecision(2) << " time= " << time << " s." << std::endl;
        std::cout.copyfmt(init);
    }

    /**
     * @brief Prints a headline message with formatting.
     *
     * @param message The headline text.
     */
    inline void printHeadLines(const std::string &message) {
        std::cout << " ****** " << message << " ****** " << std::endl;
    }

    /**
     * @brief Global timer for tracking elapsed time.
     *
     * The GlobTimer starts upon instantiation and provides methods to report or retrieve
     * the elapsed time.
     */
    struct GlobTimer {
        std::chrono::time_point<std::chrono::high_resolution_clock> start;

        // Constructor initializes the timer.
        GlobTimer() {
            start = std::chrono::high_resolution_clock::now();
        }

        /**
         * @brief Reports the elapsed time to standard output.
         */
        void report() const {
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - start).count();

            std::cout << SMALL_PHASE_SEPARATION;
            printTimeMessage("elapsed", elapsed);
            std::cout << SMALL_PHASE_SEPARATION;
        }

        /**
         * @brief Returns the elapsed time in seconds.
         *
         * @return Elapsed time since the timer started.
         */
        [[nodiscard]] double getTime() const {
            return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
        }
    };

    // A global instance of GlobTimer for tracking runtime across the application.
    inline GlobTimer glob_timer{};

    /**
     * @brief Rounds a double value to a specified number of decimal places.
     *
     * @param value The value to round.
     * @param decimals Number of decimal places (default is 1).
     * @return The rounded value.
     */
    inline double roundTo(double value, int decimals = 1) {
        double factor = std::pow(10.0, decimals);
        return std::round(value * factor) / factor;
    }

    /**
     * @brief Enumeration for different error and logging types.
     */
    enum class ErrorType { ERROR_EXIT, WARNING, DEBUG, REMIND };

    /**
     * @brief Logs a message with detailed context information.
     *
     * Depending on the error type, this function prints a formatted message including
     * the file name, line number, and function name, and may exit the program.
     *
     * @tparam type The logging type (ERROR_EXIT, WARNING, DEBUG, REMIND).
     * @param msg The message to log.
     * @param file The source file name.
     * @param line The line number.
     * @param func The function name.
     */
    template<ErrorType type>
    void logMessage(const std::string &msg, const char *file, int line, const char *func) {
        std::ostringstream oss;
        oss << msg << " in function " << func << " at " << file << ":" << line;

        if constexpr (type == ErrorType::ERROR_EXIT) {
            std::cerr << "\x1b[91m[FATAL ERROR] " << oss.str() << "\x1b[0m" << std::endl;
            std::exit(EXIT_FAILURE);
        } else if constexpr (type == ErrorType::WARNING) {
            std::cerr << "\x1b[1;35m[WARNING] " << oss.str() << "\x1b[0m" << std::endl;
        } else if constexpr (type == ErrorType::DEBUG) {
            std::cerr << "\x1b[94m[DEBUG] " << oss.str() << "\x1b[0m" << std::endl;
        } else if constexpr (type == ErrorType::REMIND) {
            std::cout << "\x1b[1;36m[REMIND] " << oss.str() << "\x1b[0m" << std::endl;
        }
    }

    // Convenience macros for logging messages with contextual information.
#define THROW_RUNTIME_ERROR(msg) logMessage<ErrorType::ERROR_EXIT>((msg), __FILE__, __LINE__, __func__);
#define PRINT_WARNING(msg) logMessage<ErrorType::WARNING>((msg), __FILE__, __LINE__, __func__);
#define PRINT_DEBUG(msg) logMessage<ErrorType::DEBUG>((msg), __FILE__, __LINE__, __func__);
#define PRINT_REMIND(msg) logMessage<ErrorType::REMIND>((msg), __FILE__, __LINE__, __func__);
} // namespace RouteOpt

#endif // ROUTE_OPT_ROUTE_OPT_MACRO_HPP
