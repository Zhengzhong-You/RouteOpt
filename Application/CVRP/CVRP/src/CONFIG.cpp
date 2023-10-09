//
// Created by Zhengzhong You on 7/8/22.
//

#include "CONFIG.hpp"
#include "MACRO.hpp"
#include <iostream>
#include <limits>

using namespace std;

#ifdef readEnumerationTrees
std::string CONFIG::tree_path{};
std::string CONFIG::colPool_path{};
#endif

//decide if open technique to convert dual
#ifdef changeModel
double CONFIG::convert_dual_pricingTime_vs_LPTime{2};
#elif defined (changeModel_mustUse)
double CONFIG::convert_dual_pricingTime_vs_LPTime{0};
#else
double CONFIG::convert_dual_pricingTime_vs_LPTime{numeric_limits<double>::max()};
#endif

#ifdef Resolve_Ins_with_Optimal_Val
int CONFIG::Marker4FindLargeGapIns{0};
#endif
//Buckets
int CONFIG::InitialNumBuckets{25};
double CONFIG::UB{1000000};

int CONFIG::MaxNumRoutesInExact{1000};
int CONFIG::MaxNumRoutesInLighterHeur{20};
int CONFIG::MaxNumRoutesInHeavierHeur{20};

//MIP
#ifdef DynamicAdjustmentNotAllowed
double CONFIG::GapImproved4ArcEliminationNEnumeration{1e9};
double CONFIG::EnumerationFailFactor{1};
#else
double CONFIG::GapImproved4ArcEliminationNEnumeration{1.1};
double CONFIG::EnumerationFailFactor{0.8};
#endif
double CONFIG::MaxGap2TryEnumeration{0.1};

#ifdef BranchInEnuNotAllowed
int CONFIG::MaxNumLabelInEnumeration{50000};//old 500000
int CONFIG::MaxNumRouteInEnumeration_half{10000};//old 50000
int CONFIG::MaxNumRouteInEnumeration{10000};//old 5000000
double CONFIG::HardTimeThresholdInAllEnumeration{25};
#elif defined(BranchInDefNotAllowed)
int CONFIG::MaxNumLabelInEnumeration{500000};//old 500000
int CONFIG::MaxNumRouteInEnumeration_half{50000};//old 50000
int CONFIG::MaxNumRouteInEnumeration{5000000};//old 5000000
double CONFIG::HardTimeThresholdInAllEnumeration{250};
#else
int CONFIG::MaxNumLabelInEnumeration{500000};//old 500000
int CONFIG::MaxNumRouteInEnumeration_half{50000};//old 50000
int CONFIG::MaxNumRouteInEnumeration{5000000};//old 5000000
double CONFIG::HardTimeThresholdInAllEnumeration{25};
#endif

#ifdef NoAdjustmentMIPTimeLimit
double CONFIG::MIPInEnumerationTimeLimit{10000000};//old 20
#else
double CONFIG::MIPInEnumerationTimeLimit{20};//old 20
#endif

#ifdef DynamicAdjustmentNotAllowed
double CONFIG::MIPInEnumerationTimeLimitChgFactor{1};
#else
double CONFIG::MIPInEnumerationTimeLimitChgFactor{0.8};
#endif
double CONFIG::MIPKeepPercentageRowsByDensity{1};
int CONFIG::MinNumRoute4MIP{3000};
#ifdef noSolveMIPDirectly
int CONFIG::MaxNumRoute4MIP{1000};
#else
int CONFIG::MaxNumRoute4MIP{10000};
#endif
int CONFIG::MaxNumRoute4RootMIP{10000};// this is not for general mip
double CONFIG::InitGapTolerance4ArcEliminationNEnumeration{1e-4};
int CONFIG::InitialTolerance4tryEnumerationWhenArcEliminationFails{2};
int CONFIG::MaxNumCols
    {max(MaxNumColsBasic, CONFIG::MaxNumRouteInEnumeration) + LPCol_FINAL_LIMIT + LIMIT_NEW_ADDED_COLS};

//Cuts

#ifdef SOLVER_VRPTW
int CONFIG::MaxRowRank1{5};
#else
int CONFIG::MaxRowRank1{5};
#endif
double CONFIG::CutVioFactor{0.1};
int CONFIG::MaxHeurSepMem4RowRank1{16};
int CONFIG::MaxHeurInitialCsetSize4RowRank1C{CONFIG::MaxRowRank1 + 1};// need to be lower than MaxHeurSepMem4RowRank1
double CONFIG::VioThresholdInEnu{1e-6};
#ifdef NewWay2AddCuts
int CONFIG::MaxNumR1C3PerRound{150};
int CONFIG::MaxNumR1CPerRound{100};
#else
int CONFIG::MaxNumR1C3PerRound{1000};
int CONFIG::MaxNumR1CPerRound{1000};
#endif
int CONFIG::MaxNumR1C3PerRoundSecondSelect{150};
int CONFIG::MaxNumR1CPerRoundSecondSelect{100};
vector<double> CONFIG::CutsTailOff{0.015, 0.015, 0.015};
int CONFIG::CutsTolerance{3};
int CONFIG::ExtraCuts4Root{0};
int CONFIG::CutsToleranceInEnu{3};

//ng memory
int CONFIG::CycleSize{8};
#ifdef SOLVER_VRPTW
int CONFIG::InitialNGSize{8};
#else
int CONFIG::InitialNGSize{8};
#endif
int CONFIG::MaxNGSize{16};
int CONFIG::MaxNumColsInNGAug{200};
double CONFIG::NGAugTailOff{0.5};
double  CONFIG::NGAugTimeSoftThresholdFactor{1};
double CONFIG::NGAugTimeHardThresholdFactor{2};

//Branch
double CONFIG::Frac4sudoCostBranPhase0{0.5};//0.5 & sudo-cost

auto branchValue = [](const std::vector<int> &vec_val) {
  return vec_val[branch_mode - 1];
};

int CONFIG::ML_BranchPhase0 = 100;
int CONFIG::BranPhase1 = branchValue({100, 15, 100});
int CONFIG::BranPhase2 = branchValue({3, 2, 2});
int CONFIG::BranPhase1InEnu = branchValue({15, 15, 100});
int CONFIG::BranPhase2InEnu = branchValue({2, 2, 2});

double CONFIG::FracMemTolerance{0.8};
#ifdef DynamicAdjustmentNotAllowed
double CONFIG::GapRequirementReduction4Arc{1e9};
#else
double CONFIG::GapRequirementReduction4Arc{0.9};
#endif

//Cuts robust control
#ifdef if_constraint_on_small_ins
double CONFIG::CutGenTimeThresholdInPricingInitial{0.5};
#else
//double CONFIG::CutGenTimeThresholdInPricingInitial{1.5};
#endif

#ifdef rollbackOff
double CONFIG::CutGenTimeThresholdInPricingInitial{1};//2000
double CONFIG::CutGenTimeThresholdInPricing_reduce_factor{0.9};
double CONFIG::HardTimeThresholdInPricing{5};//1500000
//double CONFIG::HardTimeThresholdInAllEnumeration{20};
double CONFIG::HardTimeThresholdInArcElimination_Mid_concatenate{2};
double CONFIG::HardTimeThresholdInArcElimination_last_half{20};
#else
double CONFIG::CutGenTimeThresholdInPricingInitial{3};//2000
double CONFIG::CutGenTimeThresholdInPricing_reduce_factor{0.9};
double CONFIG::HardTimeThresholdInPricing{10};//1500000
double CONFIG::HardTimeThresholdInArcElimination_Mid_concatenate{5};
double CONFIG::HardTimeThresholdInArcElimination_last_half{25};
#endif

//Mem control
double CONFIG::MemFactor{2};
double CONFIG::CutsKeptFactor{0.9};
double CONFIG::MaxCutMemFactor{0.2};

//bucket resize
double CONFIG::BucketResizeFactor_Ratio_DominanceChecks_NonDominant{500};
double CONFIG::BucketResizeFactor_Num_BucketArcsPerVertex{10000};

//initial gap
double CONFIG::InitialGapInPricing{2};

//meetpoint control
#ifdef SYMMETRY_PROHIBIT
double CONFIG::MeetPointFactor{0.05};
double CONFIG::NumberOfOverLabelsInMeetPoint{0.2};
#endif

//early heuristic enumeration to search for better solutions
#ifdef Resolve_Ins_with_Optimal_Val
double CONFIG::EarlyHeurEnumerationGuessGap{0.007};
#else
double CONFIG::EarlyHeurEnumerationGuessGap
    {0.003};//0.003 is enough for 180, but not for 200, 200 needs 0.006 to have a fair guess
#endif
double CONFIG::EarlyHeurEnumerationGuessGapChgFactor{0.001};
double CONFIG::controlGuessGap{0.002};
double CONFIG::MinimumGapChgUnit{1};
int CONFIG::MaxNumRouteInHeurEnumeration_half{10000};
double CONFIG::HardTimeThresholdInHeurEnumeration{10};
double CONFIG::HardTimeThresholdInHeurArcElimination_last_half{10};

//enumeration
double CONFIG::LeftThresholdRCFixing4EnumerationPool{0.6};

void CONFIG::readParameters() {

}

void CONFIG::checkParameters() {
  cout << "No time limit for cutting separation! only by its tail off condition!" << endl;
  if (MaxNumRouteInEnumeration + 2 * LPCol_FINAL_LIMIT > MAX_ROUTE_PRICING) {
    cerr << "MaxNumRouteInEnumeration too big!" << endl;
    exit(0);
  }
  if (MaxNumRouteInEnumeration > MaxNumCols) {
    cerr << "MaxNumCols too small!" << endl;
    exit(0);
  }

}
//
//#define TGap 0.02///Standard for tailing-off
//#define TNum 4///Max times for reaching the limits
