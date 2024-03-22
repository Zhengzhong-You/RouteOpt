
#include "Config.hpp"
#include "MACRO.hpp"
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

int Config::InitialNumBuckets{25};
double Config::ub{1000000};

int Config::MaxNumRoutesInExact{1000};
int Config::MaxNumRoutesInLighterHeur{50};
int Config::MaxNumRoutesInHeavierHeur{30};

double Config::GapImproved4ArcEliminationNEnumeration{1.1};
double Config::GapBarIncreased4ArcEliminationNEnumeration{0.05};//if arc fail, gap required is increased.
double Config::ArcEliminationTimeIncreased{20};
double Config::EnumerationFailFactor{0.8};
double Config::MaxGap2TryEnumeration{0.008};

int Config::MaxNumLabelInEnumeration{500000};
int Config::MaxNumRouteInEnumeration_half{50000};
int Config::MaxNumRouteInEnumeration{5000000};
double Config::MIPInEnumerationTimeLimitChgFactor{0.8};
int Config::max_num_route4_mip{10000};
double Config::InitGapTolerance4ArcEliminationNEnumeration{1e-4};

int Config::MaxRowRank1{5};
double Config::CutVioFactor{0.1};//the vio standard for cut
int Config::MaxHeurSepMem4RowRank1{16};
int Config::MaxHeurInitialCsetSize4RowRank1C{Config::MaxRowRank1 + 1};
int Config::MaxNumR1C3PerRound{150};
int Config::MaxNumR1CPerRound{100};
int Config::CutsTolerance{3};
double Config::MaxCutMemFactor{0.15};//if the memory is larger than dim*factor, then abandon the cut

#if LIMITED_MEMORY_TYPE == 1
double Config::CutsTailOff{0.015};
#elif LIMITED_MEMORY_TYPE == 2
double Config::CutsTailOff{0.010};
#endif

int Config::CycleSize{8};
int Config::MaxNumColsInNGAug{200};
double Config::NGAugTailOff{0.5};
double  Config::NGAugTimeSoftThresholdFactor{1};
double Config::NGAugTimeHardThresholdFactor{2};

double Config::Frac4sudoCostBranPhase0{0.5};//0.5 & sudo-cost
auto branchValue = [](const std::vector<int> &vec_val) {
  return vec_val[BRANCH_CANDIDATES_TYPE - 1];
};
int Config::ML_BranchPhase0 = 100;
int Config::BranPhase1 = branchValue({100, 15, 100, 100});
int Config::BranPhase2 = branchValue({3, 2, 2, 2});
int Config::BranPhase1InEnu = branchValue({15, 15, 100, 100});
int Config::BranPhase2InEnu = branchValue({2, 2, 2, 2});

double Config::FracMemTolerance{0.8};

double Config::BucketResizeFactorRatioDominanceChecksNonDominant{800};
double Config::BucketResizeFactorNumBucketArcsPerVertex{10000};

double Config::MeetPointFactor{0.05};//meet point can only be changed before adding rank1 cuts!
double Config::NumberOfOverLabelsInMeetPoint{0.2};
double Config::LeftThresholdRCFixing4EnumerationPool{0.6};

int Config::CounterTryHardPricingButFailed{3};
double Config::Ratio_Adjusted_GapImprovedVSTimeIncreased{1.1};
double Config::SoftTimeDecreaseFactorInPricing{1.1};
#ifdef VALGRIND_CHECK
double Config::RatioGapImprovedVSTimeIncreased_UB{0};
double Config::RatioGapImprovedVSTimeIncreased_LB{0};
double Config::HardTimeThresholdInAllEnumeration{10000000};
double Config::CutGenTimeThresholdInPricingInitial{10000000};
double Config::HardTimeThresholdInPricing{10000000};
double Config::HardTimeThresholdInArcEliminationMidConcatenate{10000000};
double Config::HardTimeThresholdInArcEliminationLastHalf{10000000};
#else
double Config::RatioGapImprovedVSTimeIncreased_UB{0.05};
double Config::RatioGapImprovedVSTimeIncreased_LB{0.005};
#ifdef PRICING_HARD
double Config::HardTimeThresholdInAllEnumeration{100};
double Config::CutGenTimeThresholdInPricingInitial{3};
double Config::HardTimeThresholdInPricing{20};
double Config::HardTimeThresholdInArcEliminationMidConcatenate{5};
double Config::HardTimeThresholdInArcEliminationLastHalf{25};
#else
double Config::HardTimeThresholdInAllEnumeration{20};
double Config::CutGenTimeThresholdInPricingInitial{20};//2
double Config::HardTimeThresholdInPricing{100};//5
double Config::HardTimeThresholdInArcEliminationMidConcatenate{5};
double Config::HardTimeThresholdInArcEliminationLastHalf{25};
#endif

#endif
int Config::InitialNGSize{8};
int Config::MaxNGSize{12};

int Config::MaxNumLabelReCalculateRC{10000};

void Config::checkParameters() {
  if (MaxHeurInitialCsetSize4RowRank1C >= MaxHeurSepMem4RowRank1) {
	throw runtime_error("MaxHeurInitialCsetSize4RowRank1C >= MaxHeurSepMem4RowRank1");
  }
}

#ifdef EXTRA_ARC_MEMORY
double Config::extra_arc_memory_percentage{0.5};
int Config::neighborhood_indicator_size{6};
#endif

#ifdef READ_ENUMERATION_TREES
std::string Config::tree_path{};
std::string Config::col_pool_path{};
#endif

#if ADJUST_MIP_TIME_LIMIT == 1
double Config::MIPInEnumerationTimeLimit{10000000};
int Config::MinNumRoute4MIP{10000};
#elif ADJUST_MIP_TIME_LIMIT == 0
double Config::MIPInEnumerationTimeLimit{20};
int Config::MinNumRoute4MIP{3000};
#else
throw std::runtime_error("ADJUST_MIP_TIME_LIMIT is not well defined!");
#endif