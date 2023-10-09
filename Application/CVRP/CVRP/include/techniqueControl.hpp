//
// Created by You, Zhengzhong on 6/17/23.
//

#ifndef VRPTW_INCLUDE_TECHNIQUECONTROL_HPP_
#define VRPTW_INCLUDE_TECHNIQUECONTROL_HPP_

//#define NominalBranchingInEnu
/**
 * this is the nominal branching in the enumeration. branching by half optGap!
 * 1. record the col cannot be used over twice!
 * 2. branching on using these cols or not!
 * 3. calculate the fixed (\sum x=1) lp node first, and if need branching then put it into the queue!
 */

//#define SolveMIPNotFixAll
/**
 * do not change variables to 0 or 1, about columns whose rc is larger than 1/2 gap
 */


#endif //VRPTW_INCLUDE_TECHNIQUECONTROL_HPP_


