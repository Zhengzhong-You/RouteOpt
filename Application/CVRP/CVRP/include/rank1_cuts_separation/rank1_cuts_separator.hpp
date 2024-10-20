//
// Created by Ricky You on 10/11/24.
//

#ifndef RANK1_CUTS_INFO_HPP
#define RANK1_CUTS_INFO_HPP

#include "macro.hpp"
#include "solver.hpp"
#include "solver_grb.hpp"
#include "solver_cplex.hpp"
#include <Eigen/Sparse>
#include <unordered_set>
#include <vector>


#define PRINT_DETAILED_CUT_INFO
#define INITIAL_RANK_1_MULTI_LABEL_POOL_SIZE 2048
#define FIND_MEM_USE_ENUMERATION_OR_MIP 1000 //when to use MIP to find the least memory
#define NO_MEMORY 0
#define NODE_MEMORY 1
#define ARC_MEMORY 2
#define INITIAL_IDX_R1C (-1)
#define MAX_NUMBER_SEGMENT_FOR_ONE_COLUMN 16
#define MAX_UNDETERMINED_ARC_NUMBER 128
#define MAX_LABELS 1024
#define MAX_NUM_R1CS_IN_PRICING 2048
#define MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX 128
using ARCBIT = std::bitset<MAX_UNDETERMINED_ARC_NUMBER>;

#ifdef FIND_MEM_USE_ENUMERATION_OR_MIP
#define TIME_LIMIT_FOR_MIP_FIND_MEM 0.2
#endif

using cutLong = yzzLong;
using sparseRowMatrixXI = Eigen::SparseMatrix<int, Eigen::RowMajor>;
using sparseRowVectorXI = Eigen::SparseVector<int, Eigen::RowMajor>;

struct VectorHashInRank1;

struct R1c {
    std::pair<std::vector<int>, int> info_r1c{}; // cut and plan
    int idx_r1c{INITIAL_IDX_R1C};
    int rhs{}; // get<2>map[cut.size()][plan_idx]
    std::vector<std::pair<std::vector<int>, int> > arc_mem{};
};

struct Rank1MultiLabel {
    std::vector<int> c;
    std::vector<int> w_no_c;
    int plan_idx{};
    double vio{};
    char search_dir{};

    Rank1MultiLabel(std::vector<int> c, std::vector<int> w_no_c, int plan_idx, double vio,
                    char search_dir) : c(std::move(c)),
                                       w_no_c(std::move(w_no_c)),
                                       plan_idx(plan_idx),
                                       vio(vio),
                                       search_dir(search_dir) {
    }

    Rank1MultiLabel() = default;
};

struct RouteInfo {
    double frac_x{};
    std::vector<int> col_seq{};
    int forward_concatenate_pos{};
};

class Rank1CutsSeparator {
    static int max_row_rank1;
    static double cut_vio_factor;
    static int max_heuristic_sep_mem4_row_rank1;
    static int max_heuristic_initial_seed_set_size_row_rank1c; //usually max_row_rank1+1
    static int max_num_r1c3_per_round;
    static int max_num_r1c_per_round;
    static double max_cut_mem_factor;
    static const double tolerance;
    static std::vector<RouteInfo> sol;
    static std::vector<R1c> old_cuts; //only need for limited memory mode
    static std::vector<R1c> cuts; //current generated cuts
    static int limited_memory_type; //0: no, 1: node, 2: arc
    static int dim; //number of customers
    static int pricing_hard_level;
    static Solver solver; // only need t pass env
    static std::unordered_map<int, std::vector<std::tuple<std::vector<int>, int, int> > >
    map_rank1_multiplier; //states & denominator & rhs
    static std::unordered_map<int, std::unordered_set<int> > map_rank1_multiplier_dominance;
    static std::vector<std::vector<std::vector<std::vector<int> > > > record_map_rank1_combinations;
    static std::unordered_map<std::vector<int>, std::unordered_set<int>, VectorHashInRank1> cut_record;
    static std::vector<std::unordered_map<int, int> > v_r_map;
    static std::vector<std::pair<std::vector<int>, std::vector<int> > > c_N_noC;
    static std::vector<yzzLong> rank1_sep_heur_mem4_vertex;
    static std::vector<std::vector<double> > cost_mat4_vertex;

    static std::unordered_map<cutLong, std::vector<std::pair<std::vector<int>, double> > > map_cut_plan_vio;
    static std::unordered_map<int, std::vector<std::tuple<cutLong, int, double> > > generated_rank1_multi_pool;
    static std::vector<Rank1MultiLabel> rank1_multi_label_pool;
    static int num_label;
    static std::unordered_map<std::vector<int>, std::vector<std::vector<int> >, VectorHashInRank1>
    rank1_multi_mem_plan_map;

public:
    static void setInitialInfo(
        int max_row_rank1,
        int max_num_r1c3_per_round,
        int max_num_r1c_per_round,
        int dim,
        Solver &solver,
        const std::vector<std::vector<double> > &cost_mat4_vertex);

    static void updateInfo(int limited_memory_type,
                           int pricing_hard_level,
                           const std::vector<RouteInfo> &sol,
                           const std::vector<R1c> &old_cuts);

    static void generateOptimalMultiplier();

    static void generateSepHeurMem4Vertex();

    static void separateRank1Cuts();

    static void getRank1CutMemory();

    static void fillMemory();


    static void generatePermutations(std::unordered_map<int, int> &count_map,
                                     std::vector<int> &result,
                                     std::vector<std::vector<int> > &results,
                                     int remaining);

    static void generateR1C1();

    static void generateR1C3();

    static void chooseCuts(const std::vector<R1c> &tmp_cuts,
                           std::vector<R1c> &chosen_cuts,
                           int numCuts);

    static void getHighDimCuts();

    static void constructVRMapAndSeedCrazy();

    static void initialSupportVector();

    static void startSeedCrazy();

    static void exactFindBestPermutationForOnePlan(std::vector<int> &cut, int plan_idx, double &vio);

    static void addSearchCrazy(int plan_idx,
                               const std::vector<int> &c,
                               const std::vector<int> &w_no_c,
                               double &new_vio,
                               int &add_j);

    static void removeSearchCrazy(int plan_idx,
                                  const std::vector<int> &c,
                                  double &new_vio,
                                  int &remove_j);

    static void swapSearchCrazy(int plan_idx,
                                const std::vector<int> &c,
                                const std::vector<int> &w_no_c,
                                double &new_vio,
                                std::pair<int, int> &swap_i_j);

    static void operationsCrazy(Rank1MultiLabel &label, int &i);

    static void constructCutsCrazy();

    static void constructMemoryVertexBased();


    static void combinationUtil(const std::vector<int> &arr,
                                std::vector<int> &tmp,
                                std::vector<std::vector<int> > &data,
                                int start,
                                int end,
                                int index,
                                int r);

    static void findMemAggressively(const std::vector<std::vector<std::vector<int> > > &array,
                                    const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                                    std::unordered_set<int> &mem);

    static void combinations(const std::vector<std::vector<std::vector<int> > > &array,
                             const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                             int i,
                             const std::vector<int> &accum,
                             const std::unordered_set<int> &mem,
                             int &record_min,
                             std::unordered_set<int> &new_mem);

    static void getMemoryByMIP(const std::vector<std::vector<std::vector<int> > > &array,
                               const std::vector<std::vector<std::unordered_set<int> > > &vec_segment,
                               std::unordered_set<int> &mem, bool &if_suc);

    static void findMemoryForR1CsMode3ByEnumerationAndMIP(const cutLong &v_comb,
                                                          std::unordered_set<int> &mem,
                                                          bool &if_suc);

    static void findPlanForRank1Multi(const std::vector<int> &vis, int denominator, cutLong &mem,
                                      std::vector<std::unordered_set<int> > &segment,
                                      std::vector<std::vector<int> > &plan);

    static void findMemoryForRank1Multi(
        const std::pair<std::vector<int>, int> &cut_pair,
        std::unordered_set<int> &mem,
        bool &if_suc);

    static void constructMemoryArcBased();

    static void findLeastMemoryArcBased(const sparseRowMatrixXI &sol_matrix, R1c &cut, bool &if_suc);

    static void selectR1CsByVioNMemory();

    static void getSeparatedCuts(std::vector<R1c> &cuts);

    static void getMapPlanInfo(std::vector<int> &states, int &denominator, int &rhs, int cut_dim, int plan_idx);

    static void getMapRhs(int &rhs, int cut_dim, int plan_idx);

    static void getMapDenominator(int &denominator, int cut_dim, int plan_idx);

    static bool tellIfCutsAreChangedByMem(const std::vector<R1c> &current_cuts);

    static void findMemory4Cuts(std::vector<R1c> &cuts, bool is_node_memory);
};


struct VectorHashInRank1 {
    size_t operator()(const std::vector<int> &v) const {
        std::size_t hash = 0;
        for (const int num: v) {
            hash ^= std::hash<int>{}(num) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

#ifdef PRINT_DETAILED_CUT_INFO
#define print_cuts(...) __VA_ARGS__
#else
#define print_cuts(...)
#endif

#endif //RANK1_CUTS_INFO_HPP
