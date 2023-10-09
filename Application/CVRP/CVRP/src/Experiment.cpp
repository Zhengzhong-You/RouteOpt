//
// Created by Zhengzhong You on 2/14/23.
//


#include "Experiment.hpp"

//using namespace std;

#ifdef OpenTest
#ifdef TestIfPhase0IsWorthy
std::unordered_set<std::pair<int, int>, PairHasher> Experiment::Phase0Candidates;
std::unordered_map<int, std::tuple<double, double, int>>
    Experiment::Phase0TopNAcc; // top n & acc (sum_acc & sum_take_up & all)
void Experiment::printPhase0TopNAcc() {
  cout << "Phase0TopNAcc" << endl;
  vector<pair<int, tuple<double, double, int>>> tmp(Phase0TopNAcc.begin(), Phase0TopNAcc.end());
  sort(tmp.begin(),
       tmp.end(),
       [](const pair<int, tuple<double, double, int>> &a, const pair<int, tuple<double, double, int>> &b) {
         return a.first < b.first;
       });
  for (auto &p : tmp) {
    cout << "top_" << p.first << "= " << get<0>(p.second) / get<2>(p.second) << ", take_up= "
         << get<1>(p.second) / get<2>(p.second)
         << endl;
  }
}
#endif

#ifdef TestIfMLIsWorthy
std::pair<int, int> Experiment::ML_1TopNAcc;
std::pair<int, int> Experiment::ML_2TopNAcc;
void Experiment::printMLTopNAcc() {
  cout << "ML_1TopNAcc= " << double(ML_1TopNAcc.first) / ML_1TopNAcc.second << endl;
  cout << "ML_2TopNAcc= " << double(ML_2TopNAcc.first) / ML_2TopNAcc.second << endl;
}

void Experiment::mappingfunctor(const std::vector<std::pair<std::pair<int, int>, double> > &Branch_Val,
                                std::unordered_map<std::pair<int, int>, int, PairHasher> &rank) {}
#endif

#endif

