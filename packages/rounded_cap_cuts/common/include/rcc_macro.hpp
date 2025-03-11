/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file rcc_macro.hpp
 * @brief Definitions for Rounded Capacity Cuts (RCCs) and related data structures in RouteOpt.
 *
 * This header defines enumerations of RCC form, structures of RCC, and utility operators used for managing
 * RCCs in RouteOpt. RCCs are used to strengthen the LP relaxation of formulations for capacitated vehicle routing problems.
 */


#ifndef ROUTE_OPT_RCC_MACRO_HPP
#define ROUTE_OPT_RCC_MACRO_HPP

#include <vector>

namespace RouteOpt::RCCs {
    // Enumeration to represent different forms of RCCs
    // Each enumerator corresponds to a specific formulation for the RCC.
    enum class RCCForm {
        RCC_FORM_1 = 1, // Standard RCC formulation based on internal customer grouping.
        RCC_FORM_2 = 2, // Alternate RCC formulation; typically uses internal and outside customer information.
        RCC_FORM_3 = 3 // A strengthened version of capacity cuts, supposed to be used only in enumeration state.
    };

    // Structure representing a RCC.
    struct Rcc {
        /**
  * The following formulas describe the general form of the Route Cut Constraint (RCC):
  *
  * 1. Basic Constraint for a Set S:
  *      x(S:S) ≤ |S| - k(S)
  *
  *    - x(S:S): Represents the total "flow" that start and end within the customer set S.
  *    - |S|: The total number of customers in set S.
  *    - k(S): The lower bound on the number of vehicles required to serve the customers in S.
  *
  * 2. Constraint for the Complement Set (𝑆̅):
  *      x(𝑆̅:𝑆̅) + ½·x({0}:𝑆̅) - ½·x({0}:S) ≤ |𝑆̅| - k(S)
  *
  *    - x(𝑆̅:𝑆̅): Total flow among vertices within the complement set (𝑆̅).
  *    - x({0}:𝑆̅): Flow from the depot to vertices in 𝑆̅.
  *    - x({0}:S): Flow from the depot to vertices in S.
  *    - ½ factors: These adjust the contribution of the depot flows, balancing the inflow and outflow between S and 𝑆̅.
  *
  * 3. Strengthened Constraint in the Enumeration Phase:
  *      ∑ₚ∈P I{(|p ∩ S| ≥ 1)} xₚ ≥ ⎡(1/Q) ∑ᵢ∈S qᵢ⎤
  *
  *    - ∑ₚ∈P: Summation over a set P, where each p represents a possible path or subset.
  *    - I{(|p ∩ S| ≥ 1)}: An indicator function that equals 1 if the path p includes at least one customer from S; otherwise, it is 0.
  *    - xₚ: The decision variable associated with path p.
  *    - qᵢ: A quantity associated with customer i (for instance, demand or weight).
  *    - Q: Capacity of the vehicle.
  *    - ⎡ ⎤: The ceiling function, which rounds up to the nearest integer.
  */
        int form_rcc{}; // 'form_rcc' determines which RCC formulation is used:

        std::vector<int> info_rcc_customer{}; // Customer set S.

        std::vector<int> info_rcc_outside_customer{}; // Customer set outside S.

        double rhs{}; // The right-hand side (rhs) value of the RCC.

        int idx_rcc{}; // Row index of the RCC in the constraint matrix.

        bool if_keep{}; // If this cut is kept throughout the process. Default is false.
    };

    // Overloaded equality operator to compare two RCC instances.
    inline bool operator==(const Rcc &lhs, const Rcc &rhs) {
        // Check if the right-hand side values are equal.
        if (lhs.rhs != rhs.rhs) return false;
        // Check if both RCCs have the same form.
        if (lhs.form_rcc != rhs.form_rcc) return false;
        // For RCC_FORM_1 and RCC_FORM_3, compare the info_rcc_customer vector.
        if (lhs.form_rcc == static_cast<int>(RCCForm::RCC_FORM_1) ||
            lhs.form_rcc == static_cast<int>(RCCForm::RCC_FORM_3)) {
            if (lhs.info_rcc_customer != rhs.info_rcc_customer) {
                return false;
            }
        } else {
            // For RCC_FORM_2, compare the info_rcc_outside_customer vector.
            if (lhs.info_rcc_outside_customer != rhs.info_rcc_outside_customer) {
                return false;
            }
        }
        // If all checks pass, the two RCC instances are considered equal.
        return true;
    }
}

#endif // ROUTE_OPT_RCC_MACRO_HPP
