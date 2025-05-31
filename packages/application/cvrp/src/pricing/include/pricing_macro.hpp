/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_PRICING_MACRO_HPP
#define ROUTE_OPT_PRICING_MACRO_HPP

namespace RouteOpt::Application::CVRP {
#ifdef VALGRIND_MEM_CHECK
    constexpr double CutGenTimeThresholdInPricingInitial = std::numeric_limits<float>::max();
    constexpr double TailOffHardTime = std::numeric_limits<float>::max();
    constexpr double HardTimeThresholdInPricing = std::numeric_limits<float>::max();
    constexpr double HardTimeThresholdInArcEliminationLastHalf = std::numeric_limits<float>::max();
    constexpr double HardTimeThresholdInArcEliminationMidConcatenate = std::numeric_limits<float>::max();
    constexpr double ArcEliminationTimeIncreased = std::numeric_limits<float>::max();
    constexpr double HardTimeThresholdInAllEnumeration = std::numeric_limits<float>::max();
    constexpr int LABEL_ASSIGN = 80000;
#else
    constexpr double CutGenTimeThresholdInPricingInitial = 6;
    constexpr double TailOffHardTime = 60;
    constexpr double HardTimeThresholdInPricing = 10;
    constexpr double HardTimeThresholdInArcEliminationLastHalf = 50;
    constexpr double HardTimeThresholdInArcEliminationMidConcatenate = 5;
    constexpr double ArcEliminationTimeIncreased = 20;
    constexpr double HardTimeThresholdInAllEnumeration = 50;
    constexpr int LABEL_ASSIGN = 8000000;
#endif


    constexpr bool CHECK_PRICING_LABELS = false;
    constexpr bool INSPECT_COLUMN_FEASIBILITY = false;
    constexpr double MeetPointFactor{0.2};
    constexpr double MeetPointFactorDecayFactor{2.};
    constexpr double NumberOfOverLabelsInMeetPoint = 0.2;
    constexpr double BucketResizeFactorRatioDominanceChecksNonDominant = 800;
    constexpr double BucketResizeFactorNumBucketArcsPerVertex = 10000;
    constexpr double EXPECTED_AVER_ROUTE_LENGTH = 10;
    constexpr int InitialNumBuckets = 25;
    constexpr int InitialNGSize = 8;
    constexpr int MAX_NG_SIZE = 16;
    constexpr int CYCLE_SIZE = 8;
    constexpr double NGAugTailOff = 0.5;
    constexpr double NGAugTimeSoftThresholdFactor = 1.;
    constexpr double NGAugTimeHardThresholdFactor = 2.;
    constexpr int RESOURCE_FACTOR = 100;
    constexpr int MAX_ROUTE_PRICING = 6400000;
    constexpr int MAX_NUM_REGENERATE_BUCKET = 4;
    constexpr int LABEL_LIMIT_PER_BIN = 128;
    constexpr int PRINT_LABELING_STEP_SIZE{30};
    constexpr double InitGapTolerance4ArcEliminationNEnumeration = 1e-4;
    constexpr double GapImproved4ArcEliminationNEnumeration = 1.1;
    constexpr double GapBarIncreased4ArcEliminationNEnumeration = 0.05;
    constexpr double InitialMaxGap2TryEnumeration = 0.008;
    constexpr double EnumerationFailFactor = 0.8;
    constexpr double min_enumeration_exploration_added = 0.001;
    constexpr int MaxNumLabelInEnumeration{500000};
    constexpr int MaxNumLabelInEnumerationForce{1000000};
    constexpr int NumCheckLabelInEnumeration = 50000;
    constexpr int LabelsCheckBinInEnumeration{20};
    constexpr double MaxNumLabelsCheckBinFactorInEnumeration{4.};
    constexpr double ToleranceFactorInEnumerationLabelChecking{1.3};
    constexpr int MaxNumRouteInEnumeration{5000000};
    constexpr int MaxNumRouteInEnumerationForce{10000000};
    constexpr int MaxNumRoute4Mip{10000}; //often adjustable parameter
    constexpr double PricingWarningThreshold = 0.9;
    constexpr double SOL_NG_X_TOLERANCE = 1e-4;
    constexpr int MaxNumColsInNGAug = 200;

    // stabilization
    constexpr bool IF_USE_STAB{true};
    constexpr int MIN_STAB_DELTA = 1;
    constexpr double MIN_STAB_GAMMA = 1e-4;
    constexpr int STAB_ARTI_COLUMN_IDX = -2;
    constexpr double STAB_DELTA_DECAY_FACTOR = 2.;
    constexpr double STAB_GAMMA_DECAY_FACTOR = 10.;
    constexpr double INITIAL_STAB_GAMMA{0.1};
    constexpr double INITIAL_STAB_KAPA{1};
    constexpr double VERY_INITIAL_STAB_KAPA{10.};
    constexpr int POSITIVE_ARTI_COLUMN_COL_SEQ = 1;
    constexpr int NEGATIVE_ARTI_COLUMN_COL_SEQ = -1;

    enum class PRICING_STATE {
        NORMAL,
        OUT_OF_MEMORY
    };

    enum class CONCATENATE_STATE: int {
        FAIL = -1,
        TOTALLY_FAIL = -2
    };

    enum class ENUMERATION_STATE {
        NORMAL = 0,
        OUT_OF_MEMORY = 1,
        TIME_LIMIT = 2,
        LABEL_FOR_LIMIT = 3,
        ROUTE_LIMIT = 4,
        LABEL_BACK_LIMIT = 5
    };

    enum class dominanceCoreInEnumeration_STATE {
        NO_DOMINANCE,
        KI_DOMINATE_KJ,
        KJ_DOMINATE_KI,
    };
}
#endif // ROUTE_OPT_PRICING_MACRO_HPP
