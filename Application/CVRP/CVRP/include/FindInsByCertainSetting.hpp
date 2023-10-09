//
// Created by Zhengzhong You on 3/28/23.
//

#ifndef FINDINSBYCERTAINSETTING_HPP_
#define FINDINSBYCERTAINSETTING_HPP_

#define NoAdjustmentMIPTimeLimit

//#define Setting_NoBranching_In_Enumeration
//#define Setting_Enumeration_At_Root

#ifdef Setting_NoBranching_In_Enumeration
#define BranchInEnuNotAllowed
#define CutRollBackNotAllowed
#define CutAtNonRootNotAllowed
#endif

#ifdef Setting_Enumeration_At_Root
//#define noSolveMIPDirectly //do not use this!
#define openCutsAndEnumerationAtEachEndNode
#define BranchInDefNotAllowed
#define CutRollBackNotAllowed
#define CutAtNonRootNotAllowed
#endif

#if defined(Setting_NoBranching_In_Enumeration) && defined(Setting_Enumeration_At_Root)
#error "Setting_NoBranching_In_Enumeration and Setting_Enumeration_At_Root cannot be defined at the same time"
#endif

#endif //FINDINSBYCERTAINSETTING_HPP_
