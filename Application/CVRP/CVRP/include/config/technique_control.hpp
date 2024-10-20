#ifndef FINDINSBYCERTAINSETTING_HPP_
#define FINDINSBYCERTAINSETTING_HPP_

//#define WRITE_ENUMERATION_TREES
//#define READ_ENUMERATION_TREES

//#define FASTER_DUAL
//#define DUAL_SMOOTHING
//#define COMBINE_FASTER_SMOOTHING
#define CUTTING_BRANCHING_RATIO 0.5
#define INITIAL_CUTTING_BRANCHING_RATIO 0.2

//#define DYNAMIC_ADD_CUTS
//#define LESS_CUTS

//#define HEURISTIC

//#define CONTROL_ROOT_GAP 0.005

#define SOLVER_TYPE 0
//#define VRPSOLVER_CPLEX_PARAMETER

//#define KILL_ENUMERATION_TREES

#define USE_M_DYNAMICS
#define FIX_M 60 //3PB parameter
#define R_DISCOUNT 0.72 //for tw benchmark 0.5 would be better!
#define M_BEGIN 30
#define M_DECAY_SPEED 1// beg(exp(kx/(x-1.01))+1)

#define BRANCH_CANDIDATES_TYPE 13
/**
 * 1. default 100, 3; 10, 1
 * 2. ml 100, 10, 2
 * 3. 3pb 100, 2
 * 4. 3pb 30, 3
 * 5. 3pb 60, 6
 * 6. 3pb 10, 10
 * 7. 3pb 10, 1
 * 8. 3pb 30, 1
 * 9. 3pb 30, 2
 * 10. 5， 1
 * 11. 10, 1
 * 12. 15, 1
 * 13. 10, 3
 * 14. 40, 1
 */

#define MASTER_VALVE_ML 3
/**
 * 0. forbid to use
 * 1. get data 1
 * 2. get data 2
 * 3. use model
 * 4. use model 1 + lp testing
 */

#if MASTER_VALVE_ML != 0
#undef FIX_M
#define FIX_M 37 //for tw benchmark 45 would be better!
#endif

#if MASTER_VALVE_ML == 2
#undef BRANCH_CANDIDATES_TYPE
#define BRANCH_CANDIDATES_TYPE 14
#endif

#define SOLVER_VRPTW 0
/**
 * 0: no define
 * 1. by default: solve type-2 instances, single resource
 * 2. consider two resource with cap being main resource
 * 3. consider two resource with tw being main resource
 */


//#define SOLVER_RVRPSTW
#if defined(SOLVER_RVRPSTW) && SOLVER_VRPTW != 0
#error "Only one solver can be defined"
#endif

#ifdef SOLVER_RVRPSTW
#define ONLY_ROOT_NODE
#endif

#define ADJUST_MIP_TIME_LIMIT 1
/**
 * 0. by default (not fully implemented yet)
 * 1. forbid to use
 */

#define SETTING 1
/**
 * 1. default: all function works
 * 2. all cuts in root, no cuts/arc elimination/no cuts + columns(rc) delete/ at non-root,
 * only implement for the lower bound prioritized
 */

//#define HGS_APPLIED
//#define FIND_ALL_SOLUTIONS

#define DELUXING_APPLIED 20
/**
 * HEURISTIC_DELUXING_ROUND
 */

//#define VALGRIND_CHECK

#define VERBOSE_MODE 1

/**
 * 1. all information
 * 2. highly reduced information
 */


#define SOLUTION_TYPE 1
/**
 * 1. lower bound
 * 2. depth
 */

//#define DYNAMICS_TYPE 0
/**
 * 0, linear to gap
 * 1. exponential to gap (also, linear to tree size)
 */
#if DYNAMICS_TYPE == 1
#define DYNAMICS_COEFFICIENT 0.2
#elif DYNAMICS_TYPE == 0
#define DYNAMICS_COEFFICIENT 1.23
#endif

#endif //FINDINSBYCERTAINSETTING_HPP_
