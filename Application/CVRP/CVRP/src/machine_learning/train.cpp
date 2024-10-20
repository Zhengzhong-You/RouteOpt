//
// Created by You, Zhengzhong on 5/5/24.
//

#include "train.hpp"
#include "predict.hpp"
#include "branching.hpp"

using namespace std;

int GetTrainingData::qid{};
std::string GetTrainingData::lp_output_path{};
std::string GetTrainingData::exact_output_path{};
std::unordered_map<std::pair<int, int>, std::pair<double, double>, PairHasher> GetTrainingData::edge_cg_change{};

void GetTrainingData::initOutputPath() {
    self_mkdir("train_lp");
    self_mkdir("train_exact");
    lp_output_path = "train_lp/" + cvrp->file_name + ".txt";
    exact_output_path = "train_exact/" + cvrp->file_name + ".txt";
}

void GetTrainingData::generateModelPhase1() {
    int num = Config::ML_BranchPhase0;
    BaseBranching::initialScreen(true, true, num, Config::Frac4sudoCostBranPhase0);
    getFeatureDataPhase1();
    simulateWriteLPPseudoCost();
    BaseBranching::testCG(true, true, false, true);
    writeTrainingLPFile();
}

void GetTrainingData::generateModelPhase2() {
    int num = Config::ML_BranchPhase0;
    BaseBranching::initialScreen(true, true, num, Config::Frac4sudoCostBranPhase0);
    num = cvrp->getIfInEnuState() ? Config::BranPhase1InEnu : Config::BranPhase1;
    num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
    Predict::useModelPhase1(num);
    getFeatureDataPhase2();
    num = 3; //fix to 3
    num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
    chooseBestNCandidate(num);
    BaseBranching::testCG(true, true, false, true);
    writeTrainingExactFile();
}

void GetTrainingData::simulateWriteLPPseudoCost() {
    auto tmp_pair = BaseBranching::current_branching_info.branch_pair;
    BaseBranching::current_branching_info.branch_pair.clear();
    auto &branch_pair_from_pseudo = BaseBranching::current_branching_info.branch_pair_from_pseudo;
    auto &branch_pair_from_fractional = BaseBranching::current_branching_info.branch_pair_from_fractional;
    auto ratio_pseudo =
    ((double) branch_pair_from_pseudo.size()
     / (double) (branch_pair_from_pseudo.size() + branch_pair_from_fractional.size()));
    int num = cvrp->getIfInEnuState() ? Config::BranPhase1InEnu : Config::BranPhase1;
    num = min(num, int(tmp_pair.size()));
    int pseudo_size = (int) (num * ratio_pseudo);
    int frac_size = num - pseudo_size;
    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    transform(branch_pair_from_pseudo.begin(), branch_pair_from_pseudo.begin() + pseudo_size,
              back_inserter(branch_pair),
              [](pair<int, int> &p) { return make_pair(p.first, p.second); });
    transform(branch_pair_from_fractional.begin(),
              branch_pair_from_fractional.begin() + frac_size,
              back_inserter(branch_pair),
              [](pair<int, int> &p) { return make_pair(p.first, p.second); });
    BaseBranching::testLP(0, false);
    branch_pair = tmp_pair;
}

void GetTrainingData::writeTrainingLPFile() {
    for (auto &edge: BaseBranching::current_branching_info.branch_pair_val) {
        edge_tmp_info[edge.first].sb_scores = edge.second;
    }
    auto &path = lp_output_path;
    ofstream trainingData;
    trainingData.open(path, ios::app);
    vector<pair<int, int> > record;
    record.reserve(edge_tmp_info.size());
    for (auto &tmp_info: edge_tmp_info) {
        auto find = find_if(tmp_info.second.basic_features.begin(), tmp_info.second.basic_features.end(),
                            [](const pair<string, double> &p) {
                                return p.first == PseudoMark + "ever_lp_find";
                            });
        if (abs(find->second - 1) < TOLERANCE) {
            trainingData << tmp_info.second.sb_scores;
            trainingData << " qid:" << qid;
            int cnt = 0;
            for (auto &feature: tmp_info.second.basic_features) {
                trainingData << " " << cnt << ":" << (float) feature.second;
                debug_input_data_call(trainingData << " (" << feature.first << ") ";)
                ++cnt;
            }
            trainingData << endl;
        } else {
            record.emplace_back(tmp_info.first);
        }
    }

    ++qid;

    for (auto &pr: record) {
        auto &tmp_info = edge_tmp_info[pr];
        trainingData << tmp_info.sb_scores;
        trainingData << " qid:" << qid;
        int cnt = 0;
        for (auto &feature: tmp_info.basic_features) {
            trainingData << " " << cnt << ":" << (float) feature.second;
            debug_input_data_call(trainingData << " (" << feature.first << ") ";)
            ++cnt;
        }
        trainingData << endl;
    }

    ++qid;

    trainingData.close();
}

void GetTrainingData::findCGScore4OneEdge(const std::pair<int, int> &edge, double dif, bool if_left) {
    if_left ? edge_cg_change[edge].first = dif : edge_cg_change[edge].second = dif;
}

void GetTrainingData::writeOneEdgeInfo(const std::pair<int, int> &edge, ofstream &trainingData, bool if_left) {
    auto &tmp_info = edge_tmp_info[edge];
    if (tmp_info.resolving_lp_features.empty())return;
    if (tmp_info.resolving_lp_features.back().first != "product") {
        throw runtime_error("Error: the last feature is not product but "
                            + tmp_info.resolving_lp_features.back().first);
    }
    trainingData << (if_left
                         ? (edge_lp_change[edge].first / edge_cg_change[edge].first)
                         : (edge_lp_change[edge].second
                            / edge_cg_change[edge].second)); //this is gap! (this time, we try if ratio would work!)
    trainingData << " qid:" << qid;
    int cnt = 0;
    for (auto &feature: tmp_info.basic_features) {
        trainingData << " " << cnt << ":" << (float) feature.second;
        debug_input_data_call(cout << " " << cnt << " " << feature.first << " " << feature.second << " ";)
        ++cnt;
    }
    auto &extra = if_left ? tmp_info.extra_features_edge0 : tmp_info.extra_features_edge1;
    for (auto &feature: extra) {
        trainingData << " " << cnt << ":" << (float) feature.second;
        debug_input_data_call(cout << " " << cnt << " " << feature.first << " " << feature.second << " ";)
        ++cnt;
    }
    for (auto &feature: tmp_info.resolving_lp_features) {
        trainingData << " " << cnt << ":" << (float) feature.second;
        debug_input_data_call(cout << " " << cnt << " " << feature.first << " " << feature.second << " ";)
        ++cnt;
    }
    trainingData << endl;
    debug_input_data_call(cout << endl;)
}


void GetTrainingData::writeTrainingExactFile() {
    for (auto &edge: BaseBranching::current_branching_info.branch_pair_val) {
        edge_tmp_info[edge.first].sb_scores = edge.second;
    }
    auto &path = exact_output_path;
    ofstream trainingData;
    trainingData.open(path, ios::app);
    for (auto &edge: BaseBranching::current_branching_info.branch_pair_val) {
        writeOneEdgeInfo(edge.first, trainingData, true);
        writeOneEdgeInfo(edge.first, trainingData, false);
    }

    ++qid;

    trainingData.close();
    edge_cg_change.clear();
    debug_input_data_call(printFeatures();)
}
