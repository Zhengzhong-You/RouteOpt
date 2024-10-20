#include "config.hpp"
#include "macro.hpp"
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
double Config::GapBarIncreased4ArcEliminationNEnumeration{0.05};
double Config::ArcEliminationTimeIncreased{20};
double Config::EnumerationFailFactor{0.8};
double Config::min_enumeration_exploration_added{0.001};

#ifdef SOLVER_RVRPSTW
double Config::MaxGap2TryEnumeration{0.05};
#else
double Config::MaxGap2TryEnumeration{0.008};
#endif
int Config::MaxNumLabelInEnumeration{500000};
int Config::MaxNumRouteInEnumeration_half{10000};
int Config::MaxNumRouteInEnumeration{5000000};
double Config::MIPInEnumerationTimeLimitChgFactor{0.8};
int Config::MaxNumRoute4Mip{10000};
double Config::InitGapTolerance4ArcEliminationNEnumeration{1e-4};

int Config::MaxRowRank1{5};
#ifdef LESS_CUTS
int Config::MaxNumR1C3PerRound{15};
int Config::MaxNumR1CPerRound{10};
#else
int Config::MaxNumR1C3PerRound{150};
int Config::MaxNumR1CPerRound{100};
#endif
int Config::CutsTolerance{3};
double Config::CutsTailOff{0.015};
int Config::CycleSize{8};
int Config::MaxNumColsInNGAug{200};
double Config::NGAugTailOff{0.5};
double Config::NGAugTimeSoftThresholdFactor{1};
double Config::NGAugTimeHardThresholdFactor{2};

double Config::Frac4sudoCostBranPhase0{0.5}; //0.5 & sudo-cost
auto branchValue = [](const std::vector<int> &vec_val) {
    return vec_val[BRANCH_CANDIDATES_TYPE - 1];
};
int Config::ML_BranchPhase0 = 100;
int Config::BranPhase1 = branchValue({100, 10, 100, 30, 60, 10, 10, 30, 30, 5, 10, 15, 10, 40});
int Config::BranPhase2 = branchValue({3, 2, 2, 3, 6, 10, 1, 2, 1, 1, 1, 1, 3, 1});
int Config::BranPhase1InEnu = branchValue({10, 10, 100, 30, 60, 10, 10, 30, 30, 5, 10, 15, 10, 40});
int Config::BranPhase2InEnu = branchValue({2, 2, 2, 3, 6, 10, 1, 2, 1, 1, 1, 1, 3, 1});

double Config::FracMemTolerance{0.8};

double Config::BucketResizeFactorRatioDominanceChecksNonDominant{800};
double Config::BucketResizeFactorNumBucketArcsPerVertex{10000};

double Config::MeetPointFactor{0.05}; //meet point can only be changed before adding rank1 cuts!
double Config::NumberOfOverLabelsInMeetPoint{0.2};
double Config::LeftThresholdRCFixing4EnumerationPool{0.6};

double Config::RatioGapImprovedVSTimeIncreased_LB{0.005};
double Config::Ratio_Adjusted_GapImprovedVSTimeIncreased{1.1};
double Config::SoftTimeDecreaseFactorInPricing{1.1};
#ifdef VALGRIND_CHECK
double Config::RatioGapImprovedVSTimeIncreased_UB{-1};
double Config::HardTimeThresholdInAllEnumeration{10000000};
double Config::CutGenTimeThresholdInPricingInitial{10000000};
double Config::HardTimeThresholdInPricing{10000000};
double Config::HardTimeThresholdInArcEliminationMidConcatenate{10000000};
double Config::HardTimeThresholdInArcEliminationLastHalf{10000000};
#else
double Config::RatioGapImprovedVSTimeIncreased_UB{0.05};
double Config::HardTimeThresholdInAllEnumeration{60};
double Config::CutGenTimeThresholdInPricingInitial{6};
#ifdef SOLVER_RVRPSTW
double Config::HardTimeThresholdInPricing{25};
#else
double Config::HardTimeThresholdInPricing{10};
#endif
double Config::HardTimeThresholdInArcEliminationMidConcatenate{5};
double Config::HardTimeThresholdInArcEliminationLastHalf{50};

#endif
int Config::InitialNGSize{8};
int Config::MaxNGSize{12};

std::string Config::tree_path{};
std::string Config::col_pool_path{};

#if ADJUST_MIP_TIME_LIMIT == 1
double Config::MIPInEnumerationTimeLimit{10000000};
int Config::MinNumRoute4MIP{10000};
#elif ADJUST_MIP_TIME_LIMIT == 0
double Config::MIPInEnumerationTimeLimit{20};
int Config::MinNumRoute4MIP{3000};
#else
throw std::runtime_error("ADJUST_MIP_TIME_LIMIT is not well defined!");
#endif
