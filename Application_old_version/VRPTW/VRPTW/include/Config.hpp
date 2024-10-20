
#ifndef CVRP_CONFIG_HPP
#define CVRP_CONFIG_HPP

#include<vector>
#include<string>
#include "MACRO.hpp"

class Config {
 public:
  static double ChangeDualByTimeRatio;

  static int InitialNumBuckets;
  static double ub;

  static int MaxNumRoutesInExact;
  static int MaxNumRoutesInLighterHeur;
  static int MaxNumRoutesInHeavierHeur;

  static double GapImproved4ArcEliminationNEnumeration;
  static double MaxGap2TryEnumeration;
  static double EnumerationFailFactor;
  static int MaxNumLabelInEnumeration;
  static int MinNumRoute4MIP;
  static int MaxNumRouteInEnumeration_half;
  static int MaxNumRouteInEnumeration;
  static int max_num_route4_mip;
  static double InitGapTolerance4ArcEliminationNEnumeration;
  static int InitialTolerance4tryEnumerationWhenArcEliminationFails;
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
  static std::vector<double> CutsTailOff;

  static double CutGenTimeThresholdInPricingInitial;
  static double CutGenTimeThresholdInPricingReduceFactor;
  static double HardTimeThresholdInPricing;
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

  static bool If_ExactIsMixedStyle;

  static double MeetPointFactor;
  static double NumberOfOverLabelsInMeetPoint;

  static void checkParameters();
#ifdef READ_ENUMERATION_TREES
  static std::string tree_path;
  static std::string col_pool_path;
#endif

};

#endif //CVRP_CONFIG_HPP
