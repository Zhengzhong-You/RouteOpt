
#ifndef CVRP_CONFIG_HPP
#define CVRP_CONFIG_HPP

#include<vector>
#include<string>
#include "MACRO.hpp"

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
  static int MaxNumLabelInEnumeration;
  static int MinNumRoute4MIP;
  static int MaxNumRouteInEnumeration_half;
  static int MaxNumRouteInEnumeration;
  static int max_num_route4_mip;
  static double InitGapTolerance4ArcEliminationNEnumeration;
  static double MIPInEnumerationTimeLimit;
  static double MIPInEnumerationTimeLimitChgFactor;

  static int MaxRowRank1;
  static int MaxHeurSepMem4RowRank1;
  static int MaxHeurInitialCsetSize4RowRank1C;
  static int MaxNumR1C3PerRound;
  static int MaxNumR1CPerRound;
  static double CutVioFactor;
  static int CutsTolerance;
  static double MaxCutMemFactor;
  static double CutsTailOff;

  static double RatioGapImprovedVSTimeIncreased_UB;
  static double RatioGapImprovedVSTimeIncreased_LB;
  static double CutGenTimeThresholdInPricingInitial;
  static double SoftTimeDecreaseFactorInPricing;
  static double HardTimeThresholdInPricing;
  static double Ratio_Adjusted_GapImprovedVSTimeIncreased;
  static int CounterTryHardPricingButFailed;
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

  static int MaxNumLabelReCalculateRC;


  static void checkParameters();
#ifdef READ_ENUMERATION_TREES
  static std::string tree_path;
  static std::string col_pool_path;
#endif

#ifdef EXTRA_ARC_MEMORY
  static double extra_arc_memory_percentage;
  static int neighborhood_indicator_size;
#endif

};

#endif //CVRP_CONFIG_HPP
