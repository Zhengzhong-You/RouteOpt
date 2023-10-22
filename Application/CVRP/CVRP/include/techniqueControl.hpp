
#ifndef FINDINSBYCERTAINSETTING_HPP_
#define FINDINSBYCERTAINSETTING_HPP_

//#define WRITE_ENUMERATION_TREES
//#define READ_ENUMERATION_TREES

#define BRANCH_CANDIDATES_TYPE 2
/**
 * 1. default 100, 3; 15, 1
 * 2. ml 100, 15, 2
 * 3. 3pb 100, 2
 */

#define MASTER_VALVE_ML

#ifdef MASTER_VALVE_ML
#define ML_STATE 4

/**
 * 1. get data 1
 * 2. get data 2
 * 3. use model
 * 4. use it in 3 but with \hat \theta changing dynamically
 */
#endif

//#define SOLVER_VRPTW

#define CHANGE_DUAL 2
/**
 * 0. by default
 * 1. force to use
 * 2. forbid to use
 */

#define ADJUST_MIP_TIME_LIMIT 1
/**
 * 0. by default
 * 1. force to use
 */

#define SETTING 1
/**
 * 1. default: all function works
 * 2.forbid enumeration
 * 	 forbid cutting at non-root
 * 	 forbid roll back
 */

//#define HGS_APPLIED

//#define DELUXING_APPLIED 10
/**
 * DELUXING_ROUND
 */


#endif //FINDINSBYCERTAINSETTING_HPP_
