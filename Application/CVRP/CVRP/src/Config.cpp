
#include "Config.hpp"
#include "MACRO.hpp"
#include <iostream>
#include <limits>

using namespace std;

int Config::InitialNumBuckets{25};
double Config::ub{1000000};

int Config::MaxNumRoutesInExact{1000};
int Config::MaxNumRoutesInLighterHeur{20};
int Config::MaxNumRoutesInHeavierHeur{20};

double Config::GapImproved4ArcEliminationNEnumeration{1.1};
double Config::EnumerationFailFactor{0.8};
double Config::MaxGap2TryEnumeration{0.1};

int Config::MaxNumLabelInEnumeration{500000};
int Config::MaxNumRouteInEnumeration_half{50000};
int Config::MaxNumRouteInEnumeration{5000000};
double Config::HardTimeThresholdInAllEnumeration{25};
double Config::MIPInEnumerationTimeLimitChgFactor{0.8};
int Config::max_num_route4_mip{10000};
double Config::InitGapTolerance4ArcEliminationNEnumeration{1e-4};
int Config::InitialTolerance4tryEnumerationWhenArcEliminationFails{2};

int Config::MaxRowRank1{5};
double Config::CutVioFactor{0.1};
int Config::MaxHeurSepMem4RowRank1{16};
int Config::MaxHeurInitialCsetSize4RowRank1C{Config::MaxRowRank1 + 1};
int Config::MaxNumR1C3PerRound{150};
int Config::MaxNumR1CPerRound{100};
int Config::CutsTolerance{3};
double Config::MaxCutMemFactor{0.15};//if the memory is larger than dim*factor, then abandon the cut
vector<double> Config::CutsTailOff{0.015, 0.015, 0.015};

int Config::CycleSize{8};
int Config::InitialNGSize{8};
int Config::MaxNGSize{16};
int Config::MaxNumColsInNGAug{200};
double Config::NGAugTailOff{0.5};
double  Config::NGAugTimeSoftThresholdFactor{1};
double Config::NGAugTimeHardThresholdFactor{2};

double Config::Frac4sudoCostBranPhase0{0.5};//0.5 & sudo-cost
auto branchValue = [](const std::vector<int> &vec_val) {
  return vec_val[BRANCH_CANDIDATES_TYPE - 1];
};
int Config::ML_BranchPhase0 = 100;
int Config::BranPhase1 = branchValue({100, 15, 100});
int Config::BranPhase2 = branchValue({3, 2, 2});
int Config::BranPhase1InEnu = branchValue({15, 15, 100});
int Config::BranPhase2InEnu = branchValue({2, 2, 2});

double Config::FracMemTolerance{0.8};

double Config::CutGenTimeThresholdInPricingInitial{3};
double Config::CutGenTimeThresholdInPricingReduceFactor{0.9};
double Config::HardTimeThresholdInPricing{10};
double Config::HardTimeThresholdInArcEliminationMidConcatenate{5};
double Config::HardTimeThresholdInArcEliminationLastHalf{25};

double Config::BucketResizeFactorRatioDominanceChecksNonDominant{500};
double Config::BucketResizeFactorNumBucketArcsPerVertex{10000};

bool Config::If_ExactIsMixedStyle{true};

double Config::MeetPointFactor{0.05};
double Config::NumberOfOverLabelsInMeetPoint{0.2};

double Config::LeftThresholdRCFixing4EnumerationPool{0.6};

void Config::checkParameters() {
  if (MaxHeurInitialCsetSize4RowRank1C >= MaxHeurSepMem4RowRank1) {
	cout << "MaxHeurInitialCsetSize4RowRank1C >= MaxHeurSepMem4RowRank1" << endl;
	exit(0);
  }
}

#ifdef READ_ENUMERATION_TREES
std::string Config::tree_path{};
std::string Config::col_pool_path{};
#endif

#if CHANGE_DUAL == 0
double Config::ChangeDualByTimeRatio{2};
#elif CHANGE_DUAL == 1
double Config::ChangeDualByTimeRatio{1};
#elif CHANGE_DUAL == 2
double Config::ChangeDualByTimeRatio{numeric_limits<double>::max()};
#else
throw std::runtime_error("ChangeDual is not well defined!");
#endif

#if ADJUST_MIP_TIME_LIMIT == 0
double Config::MIPInEnumerationTimeLimit{10000000};
int Config::MinNumRoute4MIP{10000};
#elif ADJUST_MIP_TIME_LIMIT == 1
double Config::MIPInEnumerationTimeLimit{20};
int Config::MinNumRoute4MIP{3000};
#else
throw std::runtime_error("ADJUST_MIP_TIME_LIMIT is not well defined!");
#endif