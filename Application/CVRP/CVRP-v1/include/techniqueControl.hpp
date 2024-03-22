
#ifndef FINDINSBYCERTAINSETTING_HPP_
#define FINDINSBYCERTAINSETTING_HPP_

//#define WRITE_ENUMERATION_TREES
//#define READ_ENUMERATION_TREES

//#define FIND_INDICATOR
//#define FIND_DIF_LP

#ifdef FIND_INDICATOR
#define VALGRIND_CHECK
#endif

//#define VALGRIND_CHECK

//#define MIX_STYLE
//#define CHANGE_DUAL

#define LIMITED_MEMORY_TYPE 1
/**
 * 1. vertex_based
 * 2. arc_based
 * 3. hybrid_based
 */

#define BRANCH_CANDIDATES_TYPE 2
/**
 * 1. default 100, 3; 15, 1
 * 2. ml 100, 15, 2
 * 3. 3pb 100, 2
 * 4. no model but dynamically change, all 100 as initial, please turn off the ML!
 */

#define MASTER_VALVE_ML

#ifdef MASTER_VALVE_ML
#define ML_STATE 4

/**¬
 * 1. get data 1
 * 2. get data 2
 * 3. use model
 * 4. use it in 3 but with \hat \theta changing dynamically
 */
#endif

//#define SOLVER_VRPTW 3
/**
 * do not use 0, since 0 means undefined
 * 1. by default: solve type-2 instances, single resource
 * 2. consider two resource with cap being main resource
 * 3. consider two resource with tw being main resource
 */

#define ADJUST_MIP_TIME_LIMIT 1
/**
 * 0. by default (not fully implemented yet)
 * 1. forbid to use
 */

#define SETTING 1
/**
 * 1. default: all function works
 * 2. forbid enumeration
 * 	  forbid cutting at non-root
 * 	  forbid roll back
 */

//#define HGS_APPLIED

//#define DELUXING_APPLIED 10
/**
 * DELUXING_ROUND
 */

//#define VALGRIND_CHECK

#define VERBOSE_MODE 2

/**
 * 1. all information
 * 2. highly reduced information
 */


#endif //FINDINSBYCERTAINSETTING_HPP_
