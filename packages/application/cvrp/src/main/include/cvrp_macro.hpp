/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_CVRP_MACRO_HPP
#define ROUTE_OPT_CVRP_MACRO_HPP
#include <Eigen/Sparse>
#include "route_opt_macro.hpp"

namespace RouteOpt::Application::CVRP {
    constexpr double TIME_LIMIT{std::numeric_limits<float>::max()};

    enum class APPLICATION_TYPE {
        CVRP,
        VRPTW,
    };

    enum class VRPTW_TYPE {
        SINGLE_RESOURCE,
        CAP_MAIN_RESOURCE,
        TW_MAIN_RESOURCE,
    };

    //learning section
    enum class ML_TYPE {
        ML_NO_USE,
        ML_GET_DATA_1,
        ML_GET_DATA_2,
        ML_USE_MODEL
    };

    constexpr APPLICATION_TYPE app_type{APPLICATION_TYPE::CVRP};
    constexpr VRPTW_TYPE vrptw_type{VRPTW_TYPE::SINGLE_RESOURCE};
    constexpr ML_TYPE ml_type{ML_TYPE::ML_USE_MODEL};

    constexpr bool IF_WRITE_NODE_OUT{false};

    enum class BKF_TYPE {
        M_LP = ml_type == ML_TYPE::ML_NO_USE ? 60 : 37,
        M_HEUR = 3,
        N_LP = 100,
        N_HEUR = 4,
    };

    constexpr bool IF_BKF = true && (ml_type != ML_TYPE::ML_GET_DATA_1 && ml_type != ML_TYPE::ML_GET_DATA_2);

    constexpr int get_transformed_number(APPLICATION_TYPE app) {
        switch (app) {
            case APPLICATION_TYPE::CVRP:
                return 0;
            case APPLICATION_TYPE::VRPTW:
                return 1;
        }
        return 0;
    }

    constexpr double pow10_int(int exp) {
        double value = 1.0;
        for (int i = 0; i < exp; ++i) {
            value *= 10.0;
        }
        return value;
    }

    constexpr int transformed_number = get_transformed_number(app_type);
    constexpr double round_up_tolerance = -1.0 / pow10_int(transformed_number) + TOLERANCE; // std::pow not constexpr in C++20

    inline double ceilTransformedNumberRelated(double x) {
        for (int i = 0; i < transformed_number; ++i) {
            x *= 10;
        }
        x = std::ceil(x);
        for (int i = 0; i < transformed_number; ++i) {
            x /= 10;
        }
        return x;
    }

    constexpr bool get_if_tw_type_match(VRPTW_TYPE tw) {
        return app_type == APPLICATION_TYPE::VRPTW && (tw == vrptw_type);
    }

    constexpr int CapResourceIdx = (app_type == APPLICATION_TYPE::CVRP || get_if_tw_type_match(
                                        VRPTW_TYPE::CAP_MAIN_RESOURCE))
                                       ? 0
                                       : 1;

    constexpr int NUM_RESOURCE = app_type == APPLICATION_TYPE::CVRP ? 1 : 2;
    constexpr int TWResourceIdx = 1 - CapResourceIdx;
    constexpr int MarkerSameNodeSetResourceIdx = CapResourceIdx;
    constexpr int NUM_RESOURCE_CMP_IN_DOMINANCE = (app_type == APPLICATION_TYPE::CVRP || get_if_tw_type_match(
                                                       VRPTW_TYPE::SINGLE_RESOURCE))
                                                      ? 1
                                                      : 2;

    constexpr int NODE_MEMORY_ROUTE_LENGTH{20};

    enum class NUM_TESTING {
        PHASE0 = 30,
        PHASE1 = 3,
        PHASE2 = 1,
        PHASE3 = 1
    };


    constexpr bool IF_HGS_HEURISTIC = false;

    constexpr bool if_symmetry_prohibit(APPLICATION_TYPE app) {
        switch (app) {
            case APPLICATION_TYPE::CVRP:
                return false;
            case APPLICATION_TYPE::VRPTW:
                return true;
        }
        return false;
    }

    constexpr bool IF_SYMMETRY_PROHIBIT = if_symmetry_prohibit(app_type);
    constexpr int MaxNumRoutesInExactInspection = 50;
    constexpr int MaxNumRoutesInExactPricingLow = 100;
    constexpr int MaxNumRoutesInExactPricingHigh = 1000;
    constexpr int MaxNumRoutesInHeavierHeur = 30;
    constexpr int MaxNumRoutesInLighterHeur = 30;
    constexpr double CUTTING_BRANCHING_RATIO = 0.2;
    constexpr double CUTTING_BRANCHING_RATIO_LOW = 0.1;
    constexpr double FracMemTolerance = 0.8;
    constexpr int NUM_THREADS_LP = 1;
    constexpr double MIP_GAP_TOLERANCE = 1e-9;
    constexpr double FeasibilityTol = 1e-7;
    constexpr int LP_COL_FINAL_LIMIT = 10000;
    constexpr double COL_KEEP_FRAC = 0.67;
    constexpr double EDGE_IF_ONE_TOLERANCE = 1e-4;


    constexpr double HEURISTIC_LIGHT_TESTING_MAX_COLUMN_RATIO = 1.5;
    constexpr double HEURISTIC_HEAVY_TESTING_MAX_COLUMN_RATIO = 2;


    //cuts
    constexpr int MAX_ROW_RANK1{5};
    constexpr int MAX_NUM_R1C3_PER_ROUND = 150;
    constexpr int MAX_NUM_R1C_PER_ROUND = 100;
    constexpr int CutsTolerance = 3;
    constexpr double NOMINAL_IMPROVE_THRESHOLD = 1e-1;
    constexpr double RELATIVE_IMPROVE_THRESHOLD = 0.04;
    constexpr double LB_RELATIVE_IMPROVE_THRESHOLD = 0.00015;
    constexpr double GAP_TOO_CLOSE_TO_ONE = 0.1;

    constexpr std::string_view model_1 = app_type == APPLICATION_TYPE::CVRP
                                             ? "model/cvrp_model_1.2.bin"
                                             : "model/vrptw_model_1.1.bin";
    constexpr std::string_view model_2 = app_type == APPLICATION_TYPE::CVRP
                                             ? "model/cvrp_model_2.2.bin"
                                             : "model/vrptw_model_2.1.bin";

    constexpr const char *MODEL_NAME = "CVRP";


    enum class PRICING_LEVEL {
        EXACT,
        HEAVY,
        LIGHT
    };

    using sparseRowMatrixXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using sparseColMatrixXd = Eigen::SparseMatrix<double, Eigen::ColMajor>;
    using sparseRowMatrixXI = Eigen::SparseMatrix<int, Eigen::RowMajor>;
    using denseRowMatrixXi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using denseColMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using RowVectorXT = Eigen::Matrix<size_t, 1, Eigen::Dynamic>;
    using RowVectorXd = Eigen::Matrix<double, 1, Eigen::Dynamic>;
}

#endif // ROUTE_OPT_CVRP_MACRO_HPP
