//
// Created by Zhengzhong You on 7/8/22.
//

#ifndef CVRP_CONFIG_HPP
#define CVRP_CONFIG_HPP

#include<vector>
#include<string>
#include "MACRO.hpp"
#include "ML2.hpp"

class CONFIG {
 public:

  static double convert_dual_pricingTime_vs_LPTime;

#ifdef readEnumerationTrees
  static std::string tree_path;
  static std::string colPool_path;
#endif

  static int Marker4FindLargeGapIns;

  //Buckets
  static int InitialNumBuckets;

  // NG Mem

  static double UB;
  //Labeling
  static int MaxNumRoutesInExact;
  static int MaxNumRoutesInLighterHeur;
  static int MaxNumRoutesInHeavierHeur;
  static int MaxNumCols;

  //MIP
  static double GapImproved4ArcEliminationNEnumeration;
  static double MaxGap2TryEnumeration;
  static double EnumerationFailFactor;
  static int MaxNumLabelInEnumeration;
  static double MIPInEnumerationTimeLimit;
  static double MIPInEnumerationTimeLimitChgFactor;
  static double MIPKeepPercentageRowsByDensity;
  static int MinNumRoute4MIP;
  static int MaxNumRouteInEnumeration_half;
  static int MaxNumRouteInEnumeration;
  static int MaxNumRoute4MIP;
  static int MaxNumRoute4RootMIP;
  static double InitGapTolerance4ArcEliminationNEnumeration;
  static int InitialTolerance4tryEnumerationWhenArcEliminationFails;

  static void readParameters();

  static void checkParameters();

  //Cuts
  static int MaxRowRank1;
  static int MaxHeurSepMem4RowRank1;
  static int MaxHeurInitialCsetSize4RowRank1C;
  static int MaxNumR1C3PerRound;
  static int MaxNumR1CPerRound;
  static int MaxNumR1C3PerRoundSecondSelect;
  static int MaxNumR1CPerRoundSecondSelect;
  static double CutVioFactor;
  static double VioThresholdInEnu;
  static std::vector<double> CutsTailOff;
  static int ExtraCuts4Root;
  static int CutsTolerance;
  static int CutsToleranceInEnu;

  //ng memory
  static int CycleSize;
  static int InitialNGSize;
  static int MaxNGSize;
  static int MaxNumColsInNGAug;
  static double NGAugTailOff;
  static double NGAugTimeSoftThresholdFactor;
  static double NGAugTimeHardThresholdFactor;

  //Branch
  static double Frac4sudoCostBranPhase0;//0.5 & sudo-cost
  static int ML_BranchPhase0;
  static int BranPhase1;//lp
  static double estimateFactorPhase1;//
  static int BranPhase2;//heur test
  static double estimateFactorPhase2;//
  static int BranPhase1InEnu;
  static int BranPhase2InEnu;
  static int UseModelBranRootStageModelNormal;
  static int UseModelBranRootStageModelPseudo;
  static int UseModelBranNormalStage4ModelNormal;
  static int UseModelBranNormalStage4ModelPseudo;
  static int MaxUseModelBranPhase2;
  static int MaxUseModelBranPhase1;
  static std::vector<std::pair<int, int>> TreeLevel_Num_Vec;
  static std::vector<std::pair<int, int>> TreeLevel_Num_Vec_enu;

  //other
  static double FracMemTolerance;
  static double GapRequirementReduction4Arc;

  //Cuts robust control
  static double CutGenTimeThresholdInPricingInitial;
  static double CutGenTimeThresholdInPricing_reduce_factor;
  static double HardTimeThresholdInPricing;
  static double HardTimeThresholdInAllEnumeration;
  static double HardTimeThresholdInArcElimination_Mid_concatenate;
  static double HardTimeThresholdInArcElimination_last_half;

  //Mem control
  static double MemFactor;
  static double CutsKeptFactor;
  static double MaxCutMemFactor;

  //bucket resize
  static double BucketResizeFactor_Ratio_DominanceChecks_NonDominant;
  static double BucketResizeFactor_Num_BucketArcsPerVertex;

  //initial gap
  static double InitialGapInPricing;

  //meetpoint control in non-symmetry case
  static double MeetPointFactor;
  static double NumberOfOverLabelsInMeetPoint;

  //early heuristic enumeration to search for better solutions
  static double EarlyHeurEnumerationGuessGap;
  static double EarlyHeurEnumerationGuessGapChgFactor;
  static double controlGuessGap;
  static double MinimumGapChgUnit;
  static double HardTimeThresholdInHeurEnumeration;
  static double HardTimeThresholdInHeurArcElimination_last_half;
  static int MaxNumRouteInHeurEnumeration_half;

  //enumeration
  static double LeftThresholdRCFixing4EnumerationPool;
};

#endif //CVRP_CONFIG_HPP
