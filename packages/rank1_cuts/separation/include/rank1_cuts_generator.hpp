/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#ifndef ROUTE_OPT_RANK1_CUTS_GENERATOR_HPP
#define ROUTE_OPT_RANK1_CUTS_GENERATOR_HPP
#include <vector>
#include "rank1_data_shared.hpp"
#include "rank1_macro.hpp"


namespace RouteOpt::Rank1Cuts::Separation {
    class CutGenerator {
    public:
        CutGenerator(const Rank1CutsDataShared &rank1CutsDataShared, DataShared &sharedData) : rank1CutsDataShared_ref(
                std::ref(rank1CutsDataShared)),
            sharedData_ref(std::ref(sharedData)) {
        }

        void generateRank1Cuts() {
            generateR1C1();
            generateR1C3();
            getHighDimCuts();
        }

        void preprocess() {
            generateRecordMapRank1Combinations();
            generateMapRank1MultiplierDominance();
            generateSepHeurMem4Vertex();
        }

        void cleanData() {
            num_label = 0;
            c_N_noC.clear();
            map_cut_plan_vio.clear();
            generated_rank1_multi_pool.clear();
            rank1_multi_label_pool.clear();
        }

        //setters
        void setMaxNumR1c3PerRound(int value) { max_num_r1c3_per_round = value; }
        void setMaxNumR1cPerRound(int value) { max_num_r1c_per_round = value; }
        void setMaxRowRank1(int value) { max_row_rank1 = value; }

        void setMaxHeuristicInitialSeedSetSizeRowRank1c(int value) {
            max_heuristic_initial_seed_set_size_row_rank1c = value;
            if (max_heuristic_initial_seed_set_size_row_rank1c > MAX_RANK_ROW)
                THROW_RUNTIME_ERROR("max_heuristic_initial_seed_set_size_row_rank1c: " + std::to_string(
                    max_heuristic_initial_seed_set_size_row_rank1c) +
                " should be no larger than MAX_RANK_ROW_IN_MEM: " +
                std::to_string(MAX_RANK_ROW));
        }


        //getters
        int getMaxNumR1c3PerRound() const { return max_num_r1c3_per_round; }
        int getMaxNumR1cPerRound() const { return max_num_r1c_per_round; }
        int getMaxRowRank1() const { return max_row_rank1; }

        CutGenerator() = delete;

        ~CutGenerator() = default;

    private:
        //need to keep
        int max_num_r1c3_per_round{};
        int max_num_r1c_per_round{};
        int max_row_rank1{};
        int max_heuristic_initial_seed_set_size_row_rank1c{};
        std::vector<cutLong> rank1_sep_heur_mem4_vertex{};
        std::vector<std::vector<std::vector<std::vector<int> > > > record_map_rank1_combinations{};
        std::unordered_map<int, std::unordered_set<int> > map_rank1_multiplier_dominance{};
        //need to clean
        int num_label{};
        std::vector<std::pair<std::vector<int>, std::vector<int> > > c_N_noC{};
        std::unordered_map<cutLong, std::vector<std::pair<std::vector<int>, double> > > map_cut_plan_vio{};
        std::unordered_map<int, std::vector<std::tuple<cutLong, int, double> > > generated_rank1_multi_pool{};
        std::vector<Rank1MultiLabel> rank1_multi_label_pool{};
        //ref
        const std::reference_wrapper<const Rank1CutsDataShared> rank1CutsDataShared_ref;
        std::reference_wrapper<DataShared> sharedData_ref;

        void generateRecordMapRank1Combinations();

        void generateMapRank1MultiplierDominance();

        void generateSepHeurMem4Vertex();

        void generateR1C1();

        void generateR1C3();

        void getHighDimCuts();

        void constructCuts();

        void constructVRMapNSeed();

        void startSeed();

        void exactFindBestPermutationForOnePlan(std::vector<int> &cut, int plan_idx,
                                                double &vio);

        void addSearch(int plan_idx,
                       const std::vector<int> &c,
                       const std::vector<int> &w_no_c,
                       double &new_vio,
                       int &add_j);

        void removeSearch(int plan_idx,
                          const std::vector<int> &c,
                          double &new_vio,
                          int &remove_j);

        void swapSearch(int plan_idx,
                        const std::vector<int> &c,
                        const std::vector<int> &w_no_c,

                        double &new_vio,
                        std::pair<int, int> &swap_i_j);

        void operationsControl(
            Rank1MultiLabel &label, int &i);
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
} // namespace RouteOpt ::Rank1Cuts

#endif // ROUTE_OPT_RANK1_CUTS_GENERATOR_HPP
