#include "machine_learning.hpp"
#include "cvrp.hpp"
#include "branching.hpp"

using namespace std;
using namespace std::chrono;

CVRP *MachineLearning::cvrp{};
BbNode *MachineLearning::node{};
double MachineLearning::max_edge_cost{};
double MachineLearning::max_mid_point_edge_cord_2_depot{};
double MachineLearning::cluster_coeff{};
double MachineLearning::depot_2_center{};
std::pair<int, double> MachineLearning::average_route_length;
std::vector<std::vector<std::pair<double, double> > > MachineLearning::mid_point_edge_cord{};
std::vector<std::vector<double> > MachineLearning::mid_point_edge_cord_2_depot{};
std::vector<std::vector<std::vector<int> > > MachineLearning::node_density_in_std_dis_vec_form{};
std::vector<std::vector<double> > MachineLearning::edge_2_other_convert_dis{};
std::unordered_map<std::pair<int, int>, TmpEdgeRelatedData, PairHasher> MachineLearning::edge_tmp_info{};
std::unordered_map<std::pair<int, int>, LongEdgeRelatedData, PairHasher> MachineLearning::edge_long_info{};
std::vector<double> MachineLearning::is_in_solution{};
std::unordered_map<std::pair<int, int>, double, PairHasher> MachineLearning::edge_val{};

void checkOptimalStatus(BbNode *node) {
    int status;
    safe_solver(node->getSolver().getStatus(&status))
    if (status != SOLVER_OPTIMAL && status != SOLVER_SUBOPTIMAL) {
        if (status == SOLVER_NUMERIC) {
            barrier_call(int crossover;
                safe_solver(node->getSolver().getCrossOver(&crossover))
                if (crossover == SOLVER_CROSSOVER_DOWN) {
                safe_solver(node->getSolver().setEnvCrossOver(SOLVER_CROSSOVER_DEFAULT))
                }
                safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
                safe_solver(node->getSolver().getStatus(&status))
                if (crossover == SOLVER_CROSSOVER_DOWN) {
                safe_solver(node->getSolver().setEnvCrossOver(SOLVER_CROSSOVER_DOWN))
                })
            if (status != SOLVER_OPTIMAL && status != SOLVER_SUBOPTIMAL) {
                cout << "WARNING: status is not optimal or suboptimal again: " << status << endl;
            }
        } else {
            cout << "before status: " << status << endl;
            safe_solver(node->getSolver().setEnvOutputFlag(1, false))
            safe_solver(node->getSolver().reoptimize())
            safe_solver(node->getSolver().setEnvOutputFlag(0, false))
            safe_solver(node->getSolver().getStatus(&status))
            throw runtime_error("status is not optimal or suboptimal:" + to_string(status));
        }
    }
}

void MachineLearning::init(CVRP *pr_cvrp) {
    cvrp = pr_cvrp;
}

void MachineLearning::updateNode(BbNode *pr_node) {
    node = pr_node;
}

void MachineLearning::calculatePrerequisites() {
    auto dim = cvrp->dim;
    auto &info_vertex = cvrp->info_vertex;
    auto &cost_mat4_vertex = cvrp->cost_mat4_vertex;
    max_edge_cost = 0;
    int real_dim = dim - 1;
    for (int i = 0; i < real_dim; ++i) {
        double pre_cost = *max_element(cost_mat4_vertex[i].begin() + i + 1, cost_mat4_vertex[i].end());
        max_edge_cost = max(max_edge_cost, pre_cost);
    }
    mid_point_edge_cord.resize(dim, vector<pair<double, double> >(dim));
    for (int i = 0; i < dim; ++i) {
        for (int j = i; j < dim; ++j) {
            mid_point_edge_cord[i][j].first = (info_vertex[i][1] + info_vertex[j][1]) / 2;
            mid_point_edge_cord[i][j].second = (info_vertex[i][2] + info_vertex[j][2]) / 2;
            mid_point_edge_cord[j][i] = mid_point_edge_cord[i][j];
        }
    }
    mid_point_edge_cord_2_depot.resize(dim, vector<double>(dim));
    for (int i = 0; i < dim; ++i) {
        for (int j = i; j < dim; ++j) {
            mid_point_edge_cord_2_depot[i][j] =
                    sqrt(float(
                        (mid_point_edge_cord[i][j].first - info_vertex[0][1])
                        * (mid_point_edge_cord[i][j].first - info_vertex[0][1]) +
                        (mid_point_edge_cord[i][j].second - info_vertex[0][2])
                        * (mid_point_edge_cord[i][j].second - info_vertex[0][2])));
            mid_point_edge_cord_2_depot[j][i] = mid_point_edge_cord_2_depot[i][j];
        }
    }
    max_mid_point_edge_cord_2_depot = 0;
    for (int i = 0; i < real_dim; ++i) {
        double max_dis = *max_element(mid_point_edge_cord_2_depot[i].begin() + i + 1,
                                      mid_point_edge_cord_2_depot[i].end());
        max_mid_point_edge_cord_2_depot = max(max_mid_point_edge_cord_2_depot, max_dis);
    }
    double geo_dis = exp(accumulate(cost_mat4_vertex[0].begin() + 1,
                                    cost_mat4_vertex[0].end(),
                                    0.0,
                                    [](double a, double b) { return a + log(b); }) / (dim - 1));
    vector<yzzLong> v_neighbor(dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            if (cost_mat4_vertex[i][j] <= geo_dis) {
                v_neighbor[i].set(j);
            }
        }
    }
    vector<vector<yzzLong> > density_std_dis(dim, vector<yzzLong>(dim));
    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            density_std_dis[i][j] = v_neighbor[i] & v_neighbor[j];
            density_std_dis[j][i] = density_std_dis[i][j];
        }
    }
    node_density_in_std_dis_vec_form.resize(dim, vector<vector<int> >(dim));
    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
                if (density_std_dis[i][j].test(k)) {
                    node_density_in_std_dis_vec_form[i][j].emplace_back(k);
                }
            }
            node_density_in_std_dis_vec_form[j][i] = node_density_in_std_dis_vec_form[i][j];
        }
    }
    edge_2_other_convert_dis.resize(dim, vector<double>(dim));
    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            double aver_dis = 0;
            for (int k = 0; k < node_density_in_std_dis_vec_form[i][j].size(); ++k) {
                for (int l = k + 1; l < node_density_in_std_dis_vec_form[i][j].size(); ++l) {
                    double dif_x =
                            mid_point_edge_cord[node_density_in_std_dis_vec_form[i][j][k]][
                                node_density_in_std_dis_vec_form[i][j][l]].first
                            - mid_point_edge_cord[i][j].first;
                    double dif_y =
                            mid_point_edge_cord[node_density_in_std_dis_vec_form[i][j][k]][
                                node_density_in_std_dis_vec_form[i][j][l]].second
                            - mid_point_edge_cord[i][j].second;
                    aver_dis += sqrt(float(dif_x * dif_x + dif_y * dif_y));
                }
            }
            if (node_density_in_std_dis_vec_form[i][j].empty())
                aver_dis = 0;
            else
                aver_dis /=
                        (double) pow(node_density_in_std_dis_vec_form[i][j].size(), 2) / 2;
            edge_2_other_convert_dis[i][j] = aver_dis / geo_dis;
            edge_2_other_convert_dis[j][i] = edge_2_other_convert_dis[i][j];
        }
    }
    calculateClusteringCoefficient();
    calculateDisDepot2Center();
}

void MachineLearning::recordEdgeLongInfo(int i, int j) {
    auto &tmp = edge_long_info[{i, j}].aver_edge_lp;
    tmp.first += cvrp->arc_graph_revised[i][j];
    ++tmp.second;
}

void MachineLearning::recordDiscrepancyLongInfo(const std::pair<int, int> &edge, double cg_change, bool if_left) {
    auto &e =
            if_left
                ? edge_long_info[edge].aver_exact_lp_discrepancy_down
                : edge_long_info[edge].aver_exact_lp_discrepancy_up;
    double lp = if_left ? edge_lp_change[edge].first : edge_lp_change[edge].second;
    double dis;
    if (lp < TOLERANCE) {
        cout << "lp change is near zero!" << endl;
        dis = 0;
    } else {
        dis = 1 - cg_change / lp;
    }
    e.first += dis;
    ++e.second;
}

void MachineLearning::cleanLastData() {
    edge_tmp_info.clear();
    edge_lp_change.clear();
}

void MachineLearning::getFeatureDataPhase1() {
    if (BaseBranching::current_branching_info.branch_pair.size() == 1) {
        cout << "No Training Data Phase1: a single candidate" << endl;
        return;
    }
    double org_val = node->getCurrentNodeVal();
    int BeforeNumRow = cvrp->getNumRow();
    collectEdgeRelatedFeatures(org_val);
    collectStaticFeatures();
    vector<int> solver_ind(cvrp->getNumCol());
    vector<double> solver_val(cvrp->getNumCol());
    for (auto &edge: BaseBranching::current_branching_info.branch_pair) {
        cvrp->getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
        collectVariableRelatedFeatures(edge, solver_ind.data(), BeforeNumRow, (int) solver_ind.size(), org_val);
    }
}


void MachineLearning::getFeatureDataPhase2() {
    if (BaseBranching::current_branching_info.branch_pair.size() == 1) {
        cout << "getFeatureDataPhase2: only one edge, no need to get this data!" << endl;
        return;
    }
    auto beg = high_resolution_clock::now();
    collectOneSideEdgeFeatures();
    int BeforeNumRow = cvrp->getNumRow();
    auto org_val = node->getCurrentNodeVal();

    Brc bf;
    bf.idx_brc = cvrp->getNumRow();
    node->getBrCs().emplace_back(bf);

    vector<int> solver_ind;
    vector<double> solver_val;
    bool if_changed = false;

    barrier_call(safe_solver(node->getSolver().setEnvCrossOver(SOLVER_CROSSOVER_DOWN)))

    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    vector<EdgeScoreInfo> edge_info(branch_pair.size());
    vector<DualRC> dual_rc(branch_pair.size());

    ++cvrp->getNumRow();
    int cnt = 0;
    for (auto &edge: branch_pair) {
        double tmp_val;
        cvrp->getNewConstraintCoefficientByEdge(node, edge, solver_ind, solver_val);
        if (!if_changed) {
            safe_solver(cvrp->addBranchConstraint(solver_ind,
                solver_val,
                SOLVER_EQUAL,
                0,
                nullptr,
                node->getSolver()))
            if_changed = true;
        } else {
            cvrp->changeBranchConstraint(solver_ind,
                                         solver_val,
                                         SOLVER_EQUAL,
                                         0,
                                         BeforeNumRow,
                                         node->getSolver());
        }
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        safe_solver(node->getSolver().getObjVal(&tmp_val))
        auto dif1 = BaseBranching::calculateDifference(tmp_val, org_val);
        collectResolvingDualRC(dual_rc[cnt], edge, BeforeNumRow, true);
        safe_solver(cvrp->inverseLastBranchConstraint(SOLVER_EQUAL, 1, node->getSolver()))
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        safe_solver(node->getSolver().getObjVal(&tmp_val))
        auto dif2 = BaseBranching::calculateDifference(tmp_val, org_val);
        collectResolvingDualRC(dual_rc[cnt], edge, BeforeNumRow, false);
        edge_info[cnt].dif1 = dif1;
        edge_info[cnt].dif2 = dif2;
        edge_info[cnt++].edge = edge;
        cout << "edge: (" << edge.first << ", " << edge.second << ")" << " dif1: " << dif1 << " dif2: " << dif2 <<
                " prod: "
                << dif1 * dif2 << endl;
    }
    --cvrp->getNumRow();
    barrier_call(safe_solver(node->getSolver().setEnvCrossOver(SOLVER_CROSSOVER_DEFAULT)))
    safe_solver(node->getSolver().delConstraints(1, &BeforeNumRow))
    safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
    node->getBrCs().pop_back();

    BaseBranching::reviseExtremeUnbalancedScore(edge_info, true);

    collectResolvingFeatures(edge_info, dual_rc);

    debug_input_data_call(
        for (auto &edge : branch_pair) {
        auto &e = edge_tmp_info[edge];
        cout << "edge: " << edge.first << " " << edge.second << endl;
        for (auto &feature : e.resolving_lp_features) {
        cout << feature.first << " " << feature.second << endl;
        }
        cout << BIG_PHASE_SEPARATION;
        }
    )
    auto end = high_resolution_clock::now();
    auto eps = chrono::duration<double>(end - beg).count();
    dynamic_call({

        Dynamics::getAverageT4LPNHeuristic(
            eps / (int)branch_pair.size() / 2);
        }
    )
    verbose_call(cout << "testLP time: " << eps << endl;)
}

void MachineLearning::collectStaticFeatures() {
    calculateAverageRouteLength();
    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    for (auto &edge: branch_pair) {
        auto &e = edge_tmp_info[edge];
        e.basic_features.emplace_back("cluster_coeff", cluster_coeff);
        e.basic_features.emplace_back("depot_2_center", depot_2_center);
        e.basic_features.emplace_back("average_route_length", average_route_length.second / average_route_length.first);
    }
}

void MachineLearning::collectEdgeRelatedFeatures(double org_val) {
    auto num_col = cvrp->getNumCol();
    auto dim = cvrp->dim;
    auto &branch_pair = BaseBranching::current_branching_info.branch_pair;
    auto &cost_mat4_vertex = cvrp->cost_mat4_vertex;
    auto max_main_resource = cvrp->resource.first_res;
    auto &resource_across_arcs_in_forward_sense = cvrp->resource_across_arcs_in_forward_sense;
    if (!edge_tmp_info.empty()) throw std::runtime_error("edge_tmp_info is not empty");

    edge_val.clear();
    for (int i = 1; i <= node->getNumEdges(); ++i) {
        edge_val[make_pair(node->getEdgeTail()[i], node->getEdgeHead()[i])] = node->getEdgeValue()[i];
    }

    for (auto &edge: branch_pair) {
        auto &e = edge_tmp_info[edge];
        e.basic_features.emplace_back("tree_level", node->getTreeLevel());
        e.basic_features.emplace_back("edge_cost", cost_mat4_vertex[edge.first][edge.second] / max_edge_cost);

        double back_res = -1;
        symmetry_prohibit_call(
            back_res = cvrp->resource_across_arcs_in_backward_sense[edge.first][edge.second].first_res)
        double res;
        if (back_res != -1) {
            res = double(back_res + resource_across_arcs_in_forward_sense[edge.first][edge.second].first_res) / 2;
        } else {
            res = double(resource_across_arcs_in_forward_sense[edge.first][edge.second].first_res);
        }
        e.basic_features.emplace_back("edge_res",
                                      res / max_main_resource);

        e.basic_features.emplace_back("edge_dis_2_depot",
                                      mid_point_edge_cord_2_depot[edge.first][edge.second]
                                      / max_mid_point_edge_cord_2_depot);

        e.basic_features.emplace_back("node_density_in_std_dis_vec_form",
                                      (double) pow(node_density_in_std_dis_vec_form[edge.first][edge.second].size(),
                                                   2) / dim
                                      / dim);
        e.basic_features.emplace_back("edge_2_other_convert_dis", edge_2_other_convert_dis[edge.first][edge.second]);
    }
    is_in_solution.resize(num_col);
    safe_solver(node->getSolver().getX(0, num_col, is_in_solution.data()))
}

void MachineLearning::calculateAverageRouteLength() {
    int num_col = cvrp->getNumCol();
    vector<double> x(num_col);
    safe_solver(node->getSolver().getX(0, num_col, x.data()))
    double sum_length = 0;
    int cnt = 0;
    for (int i = 0; i < num_col; ++i) {
        if (abs(x[i]) > TOLERANCE) {
            sum_length += (double) node->getCols()[i].col_seq.size();
            ++cnt;
        }
    }
    sum_length /= cnt;
    average_route_length.second += sum_length;
    ++average_route_length.first;
}

void MachineLearning::collectVariableRelatedFeatures(
    const pair<int, int> &edge,
    const int *solver_ind,
    int BeforeNumRow,
    int numnz,
    double org_val) {
    auto &e = edge_tmp_info[edge];
    vector<double> frac_up;
    frac_up.reserve(numnz);
    for (int i = 0; i < numnz; ++i) {
        int col_idx = solver_ind[i];
        if (abs(is_in_solution[col_idx]) > TOLERANCE) {
            frac_up.emplace_back(1 - is_in_solution[col_idx]);
        }
    }
    if (frac_up.empty()) throw runtime_error("frac_up is empty");
    e.basic_features.emplace_back("mean_frac_up",
                                  accumulate(frac_up.begin(), frac_up.end(), 0.0) / (int) frac_up.size());
    e.basic_features.emplace_back("min_frac_up", *min_element(frac_up.begin(), frac_up.end()));
    e.basic_features.emplace_back("max_frac_up", *max_element(frac_up.begin(), frac_up.end()));
    e.basic_features.emplace_back("frac_edge", edge_val[edge]);
    auto &exact_improvement_up = BaseBranching::branching_history.exact_improvement_up;
    auto &exact_improvement_down = BaseBranching::branching_history.exact_improvement_down;
    auto frac_edge_down = edge_val[edge];
    auto frac_edge_up = 1 - frac_edge_down;
    auto if_up = exact_improvement_up.find(edge) == exact_improvement_up.end() ? false : true;
    auto if_down = exact_improvement_down.find(edge) == exact_improvement_down.end() ? false : true;
    double improvement_up = if_up ? (exact_improvement_up[edge].first / exact_improvement_up[edge].second) : 0;
    double improvement_down = if_down ? (exact_improvement_down[edge].first / exact_improvement_down[edge].second) : 0;
    double pseudo_cost_up = max(improvement_up * frac_edge_up, 0.0);
    double pseudo_cost_down = max(improvement_down * frac_edge_down, 0.0);
    double pseudo_cost_mean = sqrt(pseudo_cost_up * pseudo_cost_down);
    e.basic_features.emplace_back("pseudo_cost_geomean_ratio", pseudo_cost_mean / org_val);
    e.basic_features.emplace_back("ever_geomean", if_up && if_down);

    auto &lp_testing_improvement_up = BaseBranching::branching_history.lp_testing_improvement_up;
    auto &lp_testing_improvement_down = BaseBranching::branching_history.lp_testing_improvement_down;
    bool if_find = lp_testing_improvement_down.find(edge) == lp_testing_improvement_down.end() ? false : true;
    double improvement_lp_up, improvement_lp_down;
    if (if_find) {
        improvement_lp_up = lp_testing_improvement_up[edge].first / lp_testing_improvement_up[edge].second;
        improvement_lp_down = lp_testing_improvement_down[edge].first / lp_testing_improvement_down[edge].second;
    } else {
        improvement_lp_up = improvement_lp_down = 0;
    }
    double pseudo_cost_lp_up = max(improvement_lp_up * frac_edge_up, 0.0);
    double pseudo_cost_lp_down = max(improvement_lp_down * frac_edge_down, 0.0);
    double pseudo_cost_lp_mean = sqrt(pseudo_cost_lp_up * pseudo_cost_lp_down);
    e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_up_ratio", pseudo_cost_lp_up / org_val);
    e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_down_ratio", pseudo_cost_lp_down / org_val);
    e.basic_features.emplace_back(PseudoMark + "pseudo_cost_lp_geomean_ratio", pseudo_cost_lp_mean / org_val);
    e.basic_features.emplace_back(PseudoMark + "improvement_lp_up_ratio", improvement_lp_up / org_val);
    e.basic_features.emplace_back(PseudoMark + "improvement_lp_down_ratio", improvement_lp_down / org_val);
    e.basic_features.emplace_back(PseudoMark + "ever_lp_find", if_find);
    e.basic_features.emplace_back("branch_times", BaseBranching::branching_history.branch_choice[edge]);
    auto &lp = edge_long_info[edge].aver_edge_lp;
    e.basic_features.emplace_back("aver_edge_lp", lp.first / lp.second);
}


void MachineLearning::findDiscrepancyResolvingFeatures(const std::pair<int, int> &edge, bool if_left) {
    auto &long_ =
            if_left
                ? edge_long_info[edge].aver_exact_lp_discrepancy_down
                : edge_long_info[edge].aver_exact_lp_discrepancy_up;
    auto &e = if_left ? edge_tmp_info[edge].extra_features_edge0 : edge_tmp_info[edge].extra_features_edge1;
    auto string1 = if_left ? "aver_exact_lp_discrepancy_down" : "aver_exact_lp_discrepancy_up";
    if (long_.second == 0) {
        e.emplace_back("dis_times", 0);
        e.emplace_back(string1, 0);
    } else {
        e.emplace_back("dis_times", long_.second);
        e.emplace_back(string1, long_.first / long_.second);
    }
}

void MachineLearning::collectResolvingFeatures(const std::vector<EdgeScoreInfo> &edge_info,
                                               const std::vector<DualRC> &dual_rc) {
    for (int i = 0; i < edge_info.size(); ++i) {
        auto &e_info = edge_info[i];
        auto &d_info = dual_rc[i];
        auto &e = edge_tmp_info[edge_info[i].edge];
        auto &edge = e_info.edge;
        if (e_info.edge != d_info.edge) throw runtime_error("edge is not equal!");
        e.extra_features_edge0.emplace_back("rc_edge", d_info.rc1);
        e.extra_features_edge0.emplace_back("dual", d_info.dual1);
        e.extra_features_edge0.emplace_back("dif", e_info.dif1);
        updateOneSideLPChange(edge, e_info.dif1, true);
        findDiscrepancyResolvingFeatures(edge, true);
        e.extra_features_edge1.emplace_back("rc_edge", d_info.rc2);
        e.extra_features_edge1.emplace_back("dual", d_info.dual2);
        e.extra_features_edge1.emplace_back("dif", e_info.dif2);
        updateOneSideLPChange(edge, e_info.dif2, false);
        findDiscrepancyResolvingFeatures(edge, false);
        double pro = edge_lp_change.at(edge).first * edge_lp_change.at(edge).second;
        e.resolving_lp_features.emplace_back("product", pro);
    }
}

void MachineLearning::collectResolvingDualRC(DualRC &dual_rc, const std::pair<int, int> &edge,
                                             int BeforeNumRow, bool if_left) {
    vector<double> duals(3, 0);
    checkOptimalStatus(node);
    if (edge.first) safe_solver(node->getSolver().getDual(edge.first - 1, 1, duals.data()))
    safe_solver(node->getSolver().getDual(edge.second - 1, 1, duals.data() + 1))
    safe_solver(node->getSolver().getDual(BeforeNumRow, 1, duals.data() + 2))
    double rc_edge = 1 - (0.5 * (duals[0] + duals[1]) + duals[2]) / cvrp->cost_mat4_vertex[edge.first][edge.second];
    if (if_left) {
        dual_rc.edge = edge;
        dual_rc.dual1 = duals[2];
        dual_rc.rc1 = rc_edge;
    } else {
        dual_rc.dual2 = duals[2];
        dual_rc.rc2 = rc_edge;
    }
}


void MachineLearning::debugInputData(const std::pair<std::string, double> &fs) {
    if (isnan(fs.second)) {
        cout << fs.first << " is nan" << endl;
    } else if (isinf(fs.second)) {
        cout << fs.first << " is inf" << endl;
    }
}

void MachineLearning::printFeatures() {
    for (auto &tmp_info: edge_tmp_info) {
        int cnt = 0;
        for (auto &feature: tmp_info.second.basic_features) {
            cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
            ++cnt;
        }
        for (auto &feature: tmp_info.second.resolving_lp_features) {
            cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
            ++cnt;
        }
        cout << BIG_PHASE_SEPARATION;
    }
    auto &edge = *edge_tmp_info.begin();
    int cnt = 0;
    vector<int> f_set;
    for (auto &feature: edge.second.basic_features) {
        if (feature.first.find(PseudoMark) != string::npos) {
            cout << " " << cnt << ":" << feature.first << " " << feature.second << endl;
            f_set.emplace_back(cnt);
        }
        ++cnt;
    }
    cout << "f_set: " << endl;
    cout << "[";
    for (auto &tmp: f_set) {
        cout << "," << tmp;
    }
    cout << "]";
}

void MachineLearning::calculateDisDepot2Center() {
    auto dim = cvrp->dim;
    auto &info_vertex = cvrp->info_vertex;
    double center_x = 0;
    double center_y = 0;
    for (int i = 1; i < dim; ++i) {
        center_x += info_vertex[i][1];
        center_y += info_vertex[i][2];
    }
    center_x /= dim - 1;
    center_y /= dim - 1;
    depot_2_center = sqrt(float((center_x - info_vertex[0][1]) * (center_x - info_vertex[0][1]) +
                                (center_y - info_vertex[0][2]) * (center_y - info_vertex[0][2])));
    depot_2_center /= max_edge_cost;
}

void MachineLearning::chooseBestNCandidate(int num) {
    num = min(num, (int) BaseBranching::current_branching_info.branch_pair.size());
    double val = max(BaseBranching::ub - node->getCurrentNodeVal(), 0.);
    vector<tuple<pair<int, int>, double, double> > edge_lp_prod(
        BaseBranching::current_branching_info.branch_pair.size());
    int cnt = 0;
    for (auto &edge: BaseBranching::current_branching_info.branch_pair) {
        auto &e = edge_lp_prod[cnt];
        get<0>(e) = edge;
        get<1>(e) = min(edge_lp_change[edge].first, val) * min(edge_lp_change[edge].second, val);
        get<2>(e) = edge_lp_change[edge].first * edge_lp_change[edge].second;
        ++cnt;
    }
    std::sort(edge_lp_prod.begin(), edge_lp_prod.end(),
              [](const auto &a, const auto &b) {
                  if (std::get<1>(a) == std::get<1>(b)) {
                      return std::get<2>(a) > std::get<2>(b); // Sort by get<2> if get<1> is equal
                  }
                  return std::get<1>(a) > std::get<1>(b); // Otherwise, sort by get<1>
              });
    double max_value = get<1>(edge_lp_prod[0]);
    auto num_max = (int) count_if(edge_lp_prod.begin(), edge_lp_prod.end(), [max_value](const auto &a) {
        return abs(get<1>(a) - max_value) < TOLERANCE;
    });
    num = max(num, num_max);
    edge_lp_prod.resize(num);
    BaseBranching::current_branching_info.branch_pair.resize(num);
    transform(edge_lp_prod.begin(), edge_lp_prod.end(), BaseBranching::current_branching_info.branch_pair.begin(),
              [](const auto &a) {
                  return get<0>(a);
              });
}
