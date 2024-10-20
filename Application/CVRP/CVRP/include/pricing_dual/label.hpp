//
// Created by You, Zhengzhong on 6/1/24.
//

#ifndef INCLUDE_LABEL_HPP_
#define INCLUDE_LABEL_HPP_

#include "rank1_cuts_separator.hpp"
using R1CINDEX = std::bitset<MAX_NUM_R1CS_IN_PRICING>;


struct R1CPricingStat {
    int num{};
    int *valid_cut_idx{};
    int *cut_map{};

    void copyFrom(const R1CPricingStat &other) {
        num = other.num;
        std::copy(other.valid_cut_idx, other.valid_cut_idx + num, valid_cut_idx);
        std::copy(other.cut_map, other.cut_map + MAX_NUM_R1CS_IN_PRICING, cut_map);
    }
};

struct R1CUseStates {
    R1CINDEX union_map{}; // >0 and =0 is 1
    std::vector<int> v_union_mem; // >0 and =0 at least for one node to j
    std::vector<std::pair<int, int> > sparse{}; // >0 idx set
    R1CINDEX sparse_map{}; // >0 is 1 in this bitmap
};

struct ResTuple {
    res_int first_res{};
#ifdef USE_TWO_RESOURCE
  res_int second_res{};
#endif
    ResTuple operator+(const ResTuple &other) const {
        ResTuple result;
        result.first_res = this->first_res + other.first_res;
#ifdef USE_TWO_RESOURCE
	result.second_res = this->second_res + other.second_res;
#endif
        return result;
    }

    ResTuple operator-(const ResTuple &other) const {
        ResTuple result;
        result.first_res = this->first_res - other.first_res;
#ifdef USE_TWO_RESOURCE
	result.second_res = this->second_res - other.second_res;
#endif
        return result;
    }
};

template<typename DerivedLabel>
struct LabelBase {
    bool is_extended{};
    int end_vertex{};
    ResTuple res{};
    double rc{};
    double cost{};
    yzzLong pi{};
    R1CPricingStat r1c{};
    DerivedLabel *p_label{};

    virtual ~LabelBase() = default;
};

struct CVRPLabel : public LabelBase<CVRPLabel> {
};

struct RVRPSTWLabel : public LabelBase<RVRPSTWLabel> {
    int w{};
    double return_2_depot_cost{};
    double extra_term_in_backward_rc{};
    double *b_service_time{};
};

#endif //INCLUDE_LABEL_HPP_
