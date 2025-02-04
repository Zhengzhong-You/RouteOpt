//
// Created by You, Zhengzhong on 5/5/24.
//


#include "predict.hpp"
#include "branching.hpp"

using namespace std;
using namespace chrono;

int Predict::num_deep_testing = 0;
BoosterHandle Predict::booster_1{};
BoosterHandle Predict::booster_2{};
std::unordered_map<std::pair<int, int>,
    std::pair<LPNPreInfo, LPNPreInfo>, PairHasher> Predict::edge_lp_pre{};
std::vector<std::pair<std::pair<LPNPreInfo, LPNPreInfo>, std::pair<int, int> > > Predict::his_recording;
std::vector<std::pair<std::pair<LPNPreInfo, LPNPreInfo>, std::pair<int, int> > > Predict::his_recording_after_testing;
std::vector<std::pair<double, int> > Predict::aver_decrease_estimator;

void checkMLInputData(int num_row, float *data, int num_features);

void Predict::loadModel(int phase) {
    const std::string model_path = "model";
    string m1 = "cvrp_model_1.bin";
    solver_vrptw_call(m1 = "vrptw_model_1.bin")
    string m2 = "cvrp_model_2.bin";
    solver_vrptw_call(m2 = "vrptw_model_2.bin")

    if (phase == 1) {
        auto path1 = model_path + "/" + m1;
        safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_1))
        safe_xgboost(XGBoosterLoadModel(booster_1, path1.c_str()))
        cout << m1 + " loaded successfully" << endl;
    } else if (phase == 2) {
        auto path1 = model_path + "/" + m2;
        safe_xgboost(XGBoosterCreate(nullptr, 0, &booster_2))
        safe_xgboost(XGBoosterLoadModel(booster_2, path1.c_str()))
        cout << m2 + " loaded successfully" << endl;
    } else {
        throw std::runtime_error("phase is not valid");
    }
}

void Predict::predict_model_1(std::vector<std::pair<std::pair<int, int>, double> > &branch_val) {
    if (branch_val.empty()) return;

    DMatrixHandle test;
    int numFeatures;
    BoosterHandle booster;
    bst_ulong output_length;
    const float *output_result;
    output_result = nullptr;

    auto &edge = edge_tmp_info[branch_val[0].first];
    numFeatures = (int) (edge.basic_features.size());
    booster = booster_1;
    auto data = new float[branch_val.size() * numFeatures];
    for (int i = 0; i < branch_val.size(); i++) {
        auto &tmp_edge = edge_tmp_info[branch_val[i].first];
        int j = 0;
        for (auto &fs: tmp_edge.basic_features) {
            debug_input_data_call(debugInputData(fs))
            data[i * numFeatures + j] = (float) fs.second;
            ++j;
        }
    }

    debug_input_data_call(checkMLInputData((int)branch_val.size(), data, numFeatures))

    safe_xgboost(XGDMatrixCreateFromMat(data,
        (int)branch_val.size(),
        numFeatures,
        numeric_limits<float>::quiet_NaN(),
        &test))

    safe_xgboost(XGBoosterPredict(booster, test, 1, 0, 0, &output_length, &output_result))

    for (unsigned int i = 0; i < output_length; i++) {
        branch_val[i].second = output_result[i];
    }
    sort(branch_val.begin(), branch_val.end(), [](const auto &a, const auto &b) {
        return a.second > b.second;
    });
    safe_xgboost(XGDMatrixFree(test))
    delete[] data;
}

void Predict::predict_model_2() {
    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    if (branch_pair.empty()) return;
    edge_lp_pre.clear();
    his_recording.clear();
    his_recording_after_testing.clear();
    for (auto &edge: branch_pair) {
        edge_lp_pre[edge].first.lp = edge_lp_change[edge].first;
        edge_lp_pre[edge].second.lp = edge_lp_change[edge].second;
    }
    DMatrixHandle test;
    int numFeatures;
    BoosterHandle booster;
    bst_ulong output_length;
    const float *output_result;
    output_result = nullptr;

    auto &edge = edge_tmp_info[branch_pair[0]];
    numFeatures =
            (int) (edge.basic_features.size() + edge.extra_features_edge0.size() + edge.resolving_lp_features.size());
    booster = booster_2;
    int num_row = 2 * (int) branch_pair.size();
    auto data = new float[num_row * numFeatures];

    for (int i = 0; i < num_row; i++) {
        auto &tmp_edge = edge_tmp_info[branch_pair[i / 2]];
        int j = 0;
        for (auto &fs: tmp_edge.basic_features) {
            debug_input_data_call(debugInputData(fs))
            data[i * numFeatures + j] = (float) fs.second;
            ++j;
        }
        auto &extra = i % 2 == 0 ? tmp_edge.extra_features_edge0 : tmp_edge.extra_features_edge1;
        for (auto &fs: extra) {
            debug_input_data_call(debugInputData(fs))
            data[i * numFeatures + j] = (float) fs.second;
            ++j;
        }
        for (auto &fs: tmp_edge.resolving_lp_features) {
            debug_input_data_call(debugInputData(fs))
            data[i * numFeatures + j] = (float) fs.second;
            ++j;
        }
    }

    debug_input_data_call(checkMLInputData(num_row, data, numFeatures))

    safe_xgboost(XGDMatrixCreateFromMat(data,
        num_row,
        numFeatures,
        numeric_limits<float>::quiet_NaN(),
        &test))

    safe_xgboost(XGBoosterPredict(booster, test, 1, 0, 0, &output_length, &output_result))

    float min_val = *std::min_element(output_result, output_result + output_length);
    float max_val = *std::max_element(output_result, output_result + output_length);
    for (unsigned int i = 0; i < output_length; i++) {
        auto &e = edge_lp_pre[branch_pair[i / 2]];
        auto val = max(int((output_result[i] - min_val) / (max_val - min_val) * (MAX_R_SECOND_STAGE + 1) - TOLERANCE),
                       0);
        i % 2 == 0 ? e.first.pre = val : e.second.pre = val;
    }

    ml_call(MASTER_VALVE_ML == ML_USE_MODEL, dynamic_call(
                bool if_stricter =
                cvrp->getIfInEnuState() ? (int(branch_pair.size()) > Config::BranPhase2InEnu) : (int(branch_pair.size())
                    > Config::BranPhase2);
                if (if_stricter) {
                for (unsigned int i = 0; i < output_length; i++) {
                auto &e = edge_lp_pre[branch_pair[i / 2]];
                i % 2 == 0 ? e.first.pre = MAX_R_SECOND_STAGE : e.second.pre = MAX_R_SECOND_STAGE;
                }
                }
            ))

    aver_decrease_estimator.assign(MAX_RANK_DEVIATION, make_pair(0.0, 0));
    safe_xgboost(XGDMatrixFree(test))
    delete[] data;
}

void Predict::useModelPhase1(int num) {
    if (BaseBranching::current_branching_info.branch_pair.size() == 1) {
        cout << "No Model Phase1: a single candidate" << endl;
        return;
    }
    getFeatureDataPhase1();
    auto local_ratio_pseudo =
    ((double) BaseBranching::current_branching_info.branch_pair_from_pseudo.size()
     / (double) (BaseBranching::current_branching_info.branch_pair_from_pseudo.size()
                 + BaseBranching::current_branching_info.branch_pair_from_fractional.size()));
    int pseudo_size, frac_size;
    if (num == 1) {
        pseudo_size = BaseBranching::current_branching_info.branch_pair_from_pseudo.empty() ? 0 : 1;
        frac_size = 1 - pseudo_size;
    } else {
        pseudo_size =
                min((int) (num * local_ratio_pseudo),
                    (int) BaseBranching::current_branching_info.branch_pair_from_pseudo.size());
        frac_size = min(num - pseudo_size,
                        (int) BaseBranching::current_branching_info.branch_pair_from_fractional.size());
    }
    std::vector<std::pair<std::pair<int, int>, double> >
            Branch_Val_model_pseudo(BaseBranching::current_branching_info.branch_pair_from_pseudo.size());
    transform(BaseBranching::current_branching_info.branch_pair_from_pseudo.begin(),
              BaseBranching::current_branching_info.branch_pair_from_pseudo.end(),
              Branch_Val_model_pseudo.begin(),
              [](const auto &edge) {
                  return make_pair(edge, 0.0);
              });
    predict_model_1(Branch_Val_model_pseudo);
    transform(Branch_Val_model_pseudo.begin(),
              Branch_Val_model_pseudo.begin() + pseudo_size,
              BaseBranching::current_branching_info.branch_pair.begin(),
              [](const auto &a) {
                  return a.first;
              });
    std::vector<std::pair<std::pair<int, int>, double> > Branch_Val_model_fractional(
        BaseBranching::current_branching_info.branch_pair_from_fractional.size());
    transform(BaseBranching::current_branching_info.branch_pair_from_fractional.begin(),
              BaseBranching::current_branching_info.branch_pair_from_fractional.end(),
              Branch_Val_model_fractional.begin(),
              [](const auto &edge) {
                  return make_pair(edge, 0.0);
              });
    predict_model_1(Branch_Val_model_fractional);
    transform(Branch_Val_model_fractional.begin(),
              Branch_Val_model_fractional.begin() + frac_size,
              BaseBranching::current_branching_info.branch_pair.begin() + pseudo_size,
              [](const auto &a) {
                  return a.first;
              });
    BaseBranching::current_branching_info.branch_pair.resize(num);
    verbose_call(cout << "candidates in phase1: " << BaseBranching::current_branching_info.branch_pair.size() << endl;
    )
}

double Predict::getDecreaseEstimator(int pre, bool if_over_estimate) {
    double red = -1;
    if (if_over_estimate) {
        if (pre >= TRUST_PRE_BAR) ++pre;
        for (int i = pre; i < MAX_RANK_DEVIATION; ++i) {
            if (aver_decrease_estimator[i].second == 0) continue;
            red = aver_decrease_estimator[i].first / aver_decrease_estimator[i].second;
            break;
        }
        if (red == -1)red = 0;
    } else {
        if (pre >= TRUST_PRE_BAR) --pre;
        for (int i = pre; i >= 0; --i) {
            if (aver_decrease_estimator[i].second == 0) continue;
            red = aver_decrease_estimator[i].first / aver_decrease_estimator[i].second;
            break;
        }
        if (red == -1)red = 1; //if not find, then must be 0!
    }
    return red;
}

bool Predict::compareTwoCandidates(const std::pair<int, int> &edge1, const std::pair<int, int> &edge2) {
    auto &e1 = edge_lp_pre[edge1];
    auto &e2 = edge_lp_pre[edge2];
    if (e1.first.pre <= e2.first.pre) {
        if ((e1.first.pre == e2.first.pre) && (e1.first.pre >= TRUST_PRE_BAR)) { return false; }
        if (e1.second.pre <= e2.second.pre) {
            if (e1.second.pre == e2.second.pre && e1.second.pre >= TRUST_PRE_BAR) return false;
            return e1.first.lp * e1.second.lp >= e2.first.lp * e2.second.lp;
        } else {
            double red = getDecreaseEstimator(e1.second.pre, true);
            double red2 = getDecreaseEstimator(e2.second.pre, false);
            return e1.first.lp * max(e1.second.lp * red, 0.) >= e2.first.lp * max(e2.second.lp * red2, 0.) + TOLERANCE;
        }
    } else {
        if (e1.second.pre <= e2.second.pre) {
            if (e1.second.pre == e2.second.pre && e1.second.pre >= TRUST_PRE_BAR) return false;
            double red = getDecreaseEstimator(e1.first.pre, true);
            double red2 = getDecreaseEstimator(e2.first.pre, false);
            return max(e1.first.lp * red, 0.) * e1.second.lp >= max(e2.first.lp * red2, 0.) * e2.second.lp + TOLERANCE;
        } else {
            double red = getDecreaseEstimator(e1.first.pre, true);
            double red3 = getDecreaseEstimator(e1.second.pre, true);
            double red2 = getDecreaseEstimator(e2.first.pre, false);
            double red4 = getDecreaseEstimator(e2.second.pre, false);
            return max(e1.first.lp * red, 0.) * max(e1.second.lp * red3, 0.)
                   >= max(e2.first.lp * red2, 0.) * max(e2.second.lp * red4, 0.) + TOLERANCE;
        }
    }
}


bool Predict::isDominated(const std::pair<int, int> &edge2,
                          const std::vector<std::pair<std::pair<LPNPreInfo, LPNPreInfo>,
                              std::pair<int, int> > > &recordings) {
    for (auto &info: recordings) {
        auto &e = info.second;
        if (e == edge2) continue;
        auto &e1_left = info.first.first;
        auto &e1_right = info.first.second;

        if (e1_left.pre >= TRUST_PRE_BAR || e1_right.pre >= TRUST_PRE_BAR) continue;

        double min_e_pre, max_e_pre, min_edge2_pre, max_edge2_pre;
        if (e1_left.pre <= e1_right.pre) {
            min_e_pre = e1_left.pre;
            max_e_pre = e1_right.pre;
        } else {
            min_e_pre = e1_right.pre;
            max_e_pre = e1_left.pre;
        }
        if (edge_lp_pre[edge2].first.pre <= edge_lp_pre[edge2].second.pre) {
            min_edge2_pre = edge_lp_pre[edge2].first.pre;
            max_edge2_pre = edge_lp_pre[edge2].second.pre;
        } else {
            min_edge2_pre = edge_lp_pre[edge2].second.pre;
            max_edge2_pre = edge_lp_pre[edge2].first.pre;
        }
        if (min_e_pre <= min_edge2_pre && max_e_pre <= max_edge2_pre) {
            double prod = e1_left.lp * e1_right.lp;
            if (prod > edge_lp_pre[edge2].first.lp * edge_lp_pre[edge2].second.lp - TOLERANCE) {
                return true;
            }
        }
    }
    return false;
}

bool Predict::checkIfComparisonNecessary(const std::pair<int, int> &edge2) {
    if (isDominated(edge2, his_recording)) {
        verbose_call(cout << "(" << edge2.first << "," << edge2.second << ") is dominated before testing" << endl;)
        return false;
    }

    if (isDominated(edge2, his_recording_after_testing)) {
        verbose_call(cout << "(" << edge2.first << "," << edge2.second << ") is dominated after testing" << endl;)
        return false;
    }

    return true;
}

void Predict::compareOnePair(std::pair<int, int> &current_bst, const std::pair<int, int> &edge2) {
    auto &e1 = edge_lp_pre[current_bst];
    auto &e2 = edge_lp_pre[edge2];
    if (!checkIfComparisonNecessary(edge2)) return;
    double prod1, prod2;
    pair<pair<int, int>, bool> e_info;
HERE:
    prod1 = e1.first.lp * e1.second.lp;
    prod2 = e2.first.lp * e2.second.lp;
    if (prod1 < prod2) {
        if (compareTwoCandidates(edge2, current_bst)) {
            current_bst = edge2;
            goto OUT;
        }
        e_info.first = edge2;
        if (e2.first.pre >= TRUST_PRE_BAR) {
            e_info.second = e2.first.pre >= e1.first.pre;
        } else if (e2.second.pre >= TRUST_PRE_BAR) {
            e_info.second = !(e2.second.pre >= e1.second.pre);
        } else e_info.second = e2.first.pre > e1.first.pre;
        testCGOneSide(e_info);
        goto HERE;
    } else {
        if (compareTwoCandidates(current_bst, edge2)) {
            goto OUT;
        }
        e_info.first = current_bst;
        if (e1.first.pre >= TRUST_PRE_BAR) {
            e_info.second = e1.first.pre >= e2.first.pre;
        } else if (e1.second.pre >= TRUST_PRE_BAR) {
            e_info.second = !(e1.second.pre >= e2.second.pre);
        } else e_info.second = e1.first.pre > e2.first.pre;
        testCGOneSide(e_info);
        goto HERE;
    }
OUT:;
}

void Predict::testCGOneSide(std::pair<std::pair<int, int>, bool> &info) {
    ++num_deep_testing;
    auto &num_row = cvrp->getNumRow();
    auto &num_col = cvrp->getNumCol();
    int BeforeNumRow = num_row;
    bool if_left = info.second;
    auto &edge = info.first;
    if (if_left ? (edge_lp_pre[edge].first.pre == BEST_PRE) : (edge_lp_pre[edge].second.pre == BEST_PRE)) return;
    auto org_val = node->getCurrentNodeVal();

    Brc bf;
    bf.idx_brc = num_row;
    node->getBrCs().emplace_back(bf);
    int col_start = num_col;

    if (cvrp->getIfInEnuState()) {
        sparseRowMatrixXd mat(1, node->size_enumeration_col_pool);
        node->matrix_in_enumeration.push_back(std::move(mat));
    }
    vector<int> solver_ind;
    vector<double> solver_val;
    vector<int> const_for_branching(num_col);
    iota(const_for_branching.begin(), const_for_branching.end(), num_col);

    ++num_row;
    verbose_call(
        cout << "Evaluate on ( " << edge.first << " , " << edge.second << " ) " << (if_left ? "left" : "right") << endl
        ;)
    cvrp->getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
    safe_solver(cvrp->addBranchConstraint(solver_ind,
        solver_val,
        SOLVER_EQUAL,
        if_left ? 0 : 1,
        nullptr,
        node->getSolver()))
    safe_solver(node->getSolver().updateModel())
    node->getBrCs().back().edge = edge;
    node->getBrCs().back().br_dir = !if_left;

    if (cvrp->getIfInEnuState()) {
        cvrp->addBranchConstraint2ColPoolInEnumByColMap(node, edge);
        cvrp->solveLPByInspection(node, true, true, false);
    } else {
        BaseBranching::if_test_cg = true;
        cvrp->solveLPInLabeling(node, true, false, false);
        BaseBranching::if_test_cg = false;
    }
    verbose_call(cout << SMALL_PHASE_SEPARATION;)
    node->getIfTerminated() = false;
    double temp_val;
    safe_solver(node->getSolver().getObjVal(&temp_val))
    auto dif1 = BaseBranching::calculateDifference(temp_val, org_val);
    recordDiscrepancyLongInfo(edge, dif1, if_left);
    auto dif1_consider_ub = BaseBranching::calculateDifference(temp_val, org_val, true);
    verbose_call(
        cout << "dif=" << dif1 << ", dif_consider_ub=" << dif1_consider_ub << endl;)
    int len = num_col - col_start;
    if (const_for_branching.size() < len) {
        int arr_beg = (int) const_for_branching.size();
        int arr_val = const_for_branching.back() + 1;
        const_for_branching.resize(len);
        iota(const_for_branching.begin() + arr_beg, const_for_branching.end(), arr_val);
    }
    cvrp->rmLPCols(node, vector<int>(const_for_branching.begin(), const_for_branching.begin() + len));
    auto &pseudo = if_left
                       ? BaseBranching::branching_history.heuristic_improvement_down
                       : BaseBranching::branching_history.heuristic_improvement_up;

    (pseudo)[edge].first += dif1;
    ++(pseudo)[edge].second;

    safe_solver(node->getSolver().delConstraints(1, &BeforeNumRow))
    safe_solver(node->getSolver().updateModel())
    --num_row;

    if (cvrp->getIfInEnuState()) {
        node->matrix_in_enumeration.pop_back();
    }
    node->getBrCs().pop_back();
    node->getIfTerminated() = false;
    node->getCurrentNodeVal() = org_val;
    /**
     * change edge_lp_pre
     */
    auto &e_info = if_left ? edge_lp_pre[edge].first : edge_lp_pre[edge].second;
    aver_decrease_estimator[e_info.pre].first +=
            dif1 / (if_left ? edge_lp_change[edge].first : edge_lp_change[edge].second);
    ++aver_decrease_estimator[e_info.pre].second;
    e_info.lp = dif1_consider_ub;
    e_info.pre = BEST_PRE;
    his_recording_after_testing.emplace_back(edge_lp_pre[edge], edge);
}

void Predict::deeperTesting(int num) {
    num_deep_testing = 0;
    double val = max(BaseBranching::ub - node->getCurrentNodeVal(), 0.);
    for (auto &edge: edge_lp_pre) {
        auto &e = edge.second;
        e.first.lp = e.first.lp > val ? val : e.first.lp;
        e.second.lp = e.second.lp > val ? val : e.second.lp;
    }
    auto &branching_info = BaseBranching::current_branching_info.branch_pair;
    std::sort(branching_info.begin(), branching_info.end(), [](const auto &a, const auto &b) {
        auto &e1 = edge_lp_pre[a];
        auto &e2 = edge_lp_pre[b];
        if (abs(e1.first.lp * e1.second.lp - e2.first.lp * e2.second.lp) < TOLERANCE) {
            return pow(2, e1.first.pre) * pow(2, e1.second.pre) < pow(2, e2.first.pre) * pow(2, e2.second.pre);
        }
        return e1.first.lp * e1.second.lp > e2.first.lp * e2.second.lp;
    });

    verbose_call(for (auto &edge : branching_info) {
        cout << "edge: " << edge.first << "-" << edge.second << " ";
        cout << "left: " << edge_lp_pre[edge].first.lp << " right: " << edge_lp_pre[edge].second.lp << " prod: "
        << edge_lp_pre[edge].first.lp * edge_lp_pre[edge].second.lp << " " <<
        "left: " << edge_lp_pre[edge].first.pre << " right: " << edge_lp_pre[edge].second.pre << endl;
        his_recording.emplace_back(edge_lp_pre[edge], edge);
        })

    std::pair<int, int> current_bst = branching_info[0];
    if (num > 1) {
        for (int i = 1; i < branching_info.size(); ++i) {
            compareOnePair(current_bst, branching_info[i]);
        }
    }
    cout << "the best is (" << current_bst.first << "," << current_bst.second << ")" << endl;
    BaseBranching::current_branching_info.branch_pair = {current_bst};
}


void Predict::useModelPhase2(int num) {
    getFeatureDataPhase2();
    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    num = min(num, (int) branch_pair.size());
    chooseBestNCandidate(num);
    predict_model_2();
    auto beg = high_resolution_clock::now();
    deeperTesting(num);
    auto end = high_resolution_clock::now();
    if (num > 1) {
        verbose_call(cout << "time for deeper testing: " << duration<double>(end - beg).count() << endl;)
        dynamic_call(if (num_deep_testing)
            Dynamics::getAverageT4LPNHeuristic(duration<double>(end - beg).count() / num_deep_testing,
                false);)
    }
    dynamic_call(node->getDynamicKNode() = Dynamics::opt_k)
}

void Predict::useMLInGeneralFramework() {
    int num = Config::ML_BranchPhase0;
    BaseBranching::initialScreen(true, true, num, Config::Frac4sudoCostBranPhase0);
    num = cvrp->getIfInEnuState() ? Config::BranPhase1InEnu : Config::BranPhase1;
    num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
    dynamic_call(Dynamics::giveDynamicK(num))
    useModelPhase1(num);
    num = cvrp->getIfInEnuState() ? Config::BranPhase2InEnu : Config::BranPhase2;
    dynamic_call(Dynamics::giveDynamicK(num, false))
    useModelPhase2(num);
}

void Predict::useML1_N_LPTesting() {
    int num = Config::ML_BranchPhase0;
    BaseBranching::initialScreen(true, true, num, Config::Frac4sudoCostBranPhase0);
    num = cvrp->getIfInEnuState() ? Config::BranPhase1InEnu : Config::BranPhase1;
    num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
    dynamic_call(Dynamics::giveDynamicK(num))
    useModelPhase1(num);
    num = cvrp->getIfInEnuState() ? Config::BranPhase2InEnu : Config::BranPhase2;
    BaseBranching::testLP(num, true);
    BaseBranching::testCG(false, false, true, true);
}

void Predict::freeModel() {
    if (booster_1) {
        safe_xgboost(XGBoosterFree(booster_1))
        booster_1 = nullptr;
    }
    if (booster_2) {
        safe_xgboost(XGBoosterFree(booster_2))
        booster_2 = nullptr;
    }
}


void checkMLInputData(int num_row, float *data, int num_features) {
    cout << "check in MLInputData" << endl;
    size_t tmp_len = num_row * num_features;
    for (int i = 0; i < tmp_len; i++) {
        if (isnan(data[i]))
            throw std::runtime_error("nan in data");
        else if (isinf(data[i]))
            throw std::runtime_error("inf in data");
    }
}
