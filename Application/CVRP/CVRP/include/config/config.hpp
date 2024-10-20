#ifndef CONFIG_HPP
#define CONFIG_HPP

#include<vector>
#include<string>
#include "macro.hpp"

class Config {
public:
    static int InitialNumBuckets;
    static double ub;

    static int MaxNumRoutesInExact;
    static int MaxNumRoutesInLighterHeur;
    static int MaxNumRoutesInHeavierHeur;

    static double GapImproved4ArcEliminationNEnumeration;
    static double GapBarIncreased4ArcEliminationNEnumeration;
    static double ArcEliminationTimeIncreased;
    static double MaxGap2TryEnumeration;
    static double EnumerationFailFactor;
    static double min_enumeration_exploration_added;
    static int MaxNumLabelInEnumeration;
    static int MinNumRoute4MIP;
    static int MaxNumRouteInEnumeration_half;
    static int MaxNumRouteInEnumeration;
    static int MaxNumRoute4Mip;
    static double InitGapTolerance4ArcEliminationNEnumeration;
    static double MIPInEnumerationTimeLimit;
    static double MIPInEnumerationTimeLimitChgFactor;

    static int MaxRowRank1;
    static int MaxNumR1C3PerRound;
    static int MaxNumR1CPerRound;
    static int CutsTolerance;
    static double CutsTailOff;

    static double RatioGapImprovedVSTimeIncreased_UB;
    static double RatioGapImprovedVSTimeIncreased_LB;
    static double CutGenTimeThresholdInPricingInitial;
    static double SoftTimeDecreaseFactorInPricing;
    static double HardTimeThresholdInPricing;
    static double Ratio_Adjusted_GapImprovedVSTimeIncreased;
    static double HardTimeThresholdInAllEnumeration;
    static double HardTimeThresholdInArcEliminationMidConcatenate;
    static double HardTimeThresholdInArcEliminationLastHalf;

    static int CycleSize;
    static int InitialNGSize;
    static int MaxNGSize;
    static int MaxNumColsInNGAug;
    static double NGAugTailOff;
    static double NGAugTimeSoftThresholdFactor;
    static double NGAugTimeHardThresholdFactor;

    static int ML_BranchPhase0;
    static int BranPhase1;
    static int BranPhase2;
    static int BranPhase1InEnu;
    static int BranPhase2InEnu;
    static double Frac4sudoCostBranPhase0;

    static double FracMemTolerance;

    static double BucketResizeFactorRatioDominanceChecksNonDominant;
    static double BucketResizeFactorNumBucketArcsPerVertex;

    static double LeftThresholdRCFixing4EnumerationPool;

    static double MeetPointFactor;
    static double NumberOfOverLabelsInMeetPoint;

    static std::string tree_path;
    static std::string col_pool_path;
};

#endif //CONFIG_HPP
