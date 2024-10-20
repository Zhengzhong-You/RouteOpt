#include "cvrp.hpp"
#include "cutting/get_rank1_matrix.hpp"

#ifdef HEURISTIC
#include "heuristic.hpp"
#endif

#ifdef FASTER_DUAL
#include "faster_dual.hpp"
#endif

#ifdef DUAL_SMOOTHING
#include "dual_smoothing.hpp"
#endif

#ifdef COMBINE_FASTER_SMOOTHING
#include "combine_faster_smoothing.hpp"
#endif


using namespace std;
using namespace chrono;


void inspectColumns(CVRP *cvrp, BbNode *node, const ResTuple &resource, const RowVectorXd &rc, int ccnt, int num_col);

void checkRC(BbNode *node,
             const std::vector<std::tuple<Label *, Label *, double> > &negative_rc_label_tuple,
             const std::vector<SequenceInfo> &new_cols,
             const RowVectorXd &rc,
             const RowVectorXd &pi4_labeling,
             const sparseColMatrixXd &mat,
             int num_row);

void CVRP::calculateColumnCoefficientsB4Enumeration(BbNode *node, sparseColMatrixXd &mat,
                                                    Eigen::RowVectorXd &cost,
                                                    unordered_map<pair<int, int>, vector<int>, PairHasher> &edge_map) {
    int ccnt_cnt;
    int past_node;
    double cost_sum;

    int ccnt = (int) new_cols.size();
    edge_map.clear();
    edge_map.reserve(dim * dim);
    mat.resize(num_row, ccnt);
    mat.setZero();
    cost.resize(ccnt);

    unordered_map<pair<int, int>, double, PairHasher> coeff_map;
    coeff_map.reserve(dim * dim);

    ccnt_cnt = 0;
    for (auto &c: new_cols) {
        auto &col = c.col_seq;
        past_node = 0;
        cost_sum = 0;
        for (int curr_node: col) {
            cost_sum += cost_mat4_vertex[past_node][curr_node];
            ++coeff_map[{curr_node - 1, ccnt_cnt}];
            auto pr = past_node < curr_node ? make_pair(past_node, curr_node) : make_pair(curr_node, past_node);
            edge_map[pr].emplace_back(ccnt_cnt);
            past_node = curr_node;
        }
        edge_map[make_pair(0, past_node)].emplace_back(ccnt_cnt);
        cost_sum += cost_mat4_vertex[past_node][0];
        cost(ccnt_cnt) = cost_sum;
        ++ccnt_cnt;
    }

    for (auto &rcc: node->rccs) {
        if (rcc.form_rcc) {
            auto &info = rcc.info_rcc_customer;
            int idx = rcc.idx_rcc;
            for (auto it = info.begin(); it != info.end(); ++it) {
                int ai = *it;
                auto it_inner = it;
                ++it_inner;
                for (; it_inner != info.end(); ++it_inner) {
                    int aj = *it_inner;
                    auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
                    for (auto col: edge_map[pr]) ++coeff_map[{idx, col}];
                }
            }
        } else {
            auto &customer_info = rcc.info_rcc_customer;
            auto &outside_customer_info = rcc.info_rcc_outside_customer;
            int idx = rcc.idx_rcc;
            for (auto it = outside_customer_info.begin(); it != outside_customer_info.end(); ++it) {
                int ai = *it;
                auto it_inner = it;
                ++it_inner;
                for (; it_inner != outside_customer_info.end(); ++it_inner) {
                    int aj = *it_inner;
                    auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
                    for (auto it_map: edge_map[pr]) ++coeff_map[{idx, it_map}];
                }
            }
            for (int aj: outside_customer_info) {
                for (auto it_map: edge_map[make_pair(0, aj)]) coeff_map[{idx, it_map}] += 0.5;
            }
            for (int aj: customer_info) {
                for (auto it_map: edge_map[make_pair(0, aj)]) coeff_map[{idx, it_map}] -= 0.5;
            }
        }
    }

    for (auto &br: node->getBrCs()) {
        int idx = br.idx_br_c;
        if (idx == -1) continue;
        int ai = br.edge.first;
        int aj = br.edge.second;
        auto pr = ai < aj ? make_pair(ai, aj) : make_pair(aj, ai);
        for (auto it_map: edge_map[pr]) ++coeff_map[{idx, it_map}];
    }

    sparseColMatrixXd mat_r1c;
    getLimitedR1CCoeffs(new_cols, mat_r1c);

    size_t n = coeff_map.size() + mat_r1c.nonZeros() + ccnt;
    vector<Eigen::Triplet<double> > triplet(n);
    n = 0;

    for (auto &it: coeff_map) {
        triplet[n++] = {it.first.first, it.first.second, it.second};
    }

    for (int i = 0; i < mat_r1c.outerSize(); ++i) {
        for (sparseColMatrixXd::InnerIterator it(mat_r1c, i); it; ++it) {
            triplet[n++] = {lp_r1c_map[it.row()], (int) it.col(), it.value()};
        }
    }
    //  remove_vehicle_constraint_call()
    for (int i = 0; i < ccnt; ++i) {
        triplet[n++] = {real_dim, i, 1};
    }
    triplet.resize(n);

    safe_eigen(mat.setFromTriplets(triplet.begin(), triplet.end());)
}

void CVRP::addColumns(BbNode *node, int &ccnt, bool if_check_rc) {
    if (ccnt == 0) return;

    heuristic_call(Heuristic::addHeuristicMIPCols(this, node);)

    sparseColMatrixXd mat;
    Eigen::RowVectorXd cost;
    unordered_map<pair<int, int>, vector<int>, PairHasher> edge_map;

    calculateColumnCoefficientsB4Enumeration(node, mat, cost, edge_map);

    RowVectorXd rc;
    if (if_check_rc) {
        Eigen::Map<RowVectorXd> local_pi(pi4_labeling.data(), num_row);
        safe_eigen(rc = cost - local_pi * mat;)
#ifdef CHECK_RC_EVERY_COLUMN
	checkRC(node, negative_rc_label_tuple, new_cols, rc, local_pi, mat, num_row);
#endif
    } else {
        rc.resize(ccnt);
        for (int i = 0; i < ccnt; ++i) rc(i) = -1;
    }
    int ccnt_cnt = 0;
    size_t nzcnt = 0;
    vector<size_t> solver_beg(ccnt + 1);
    vector<int> solver_ind(num_row * ccnt);
    vector<double> solver_val(num_row * ccnt);
    vector<double> solver_obj(ccnt);
    for (int i = 0; i < ccnt; ++i) {
        if (rc(i) < RC_TOLERANCE) {
            solver_beg[ccnt_cnt] = nzcnt;
            solver_obj[ccnt_cnt] = cost(i);
            ++ccnt_cnt;
            for (sparseColMatrixXd::InnerIterator it(mat, i); it; ++it) {
                solver_ind[nzcnt] = (int) it.row();
                solver_val[nzcnt++] = it.value();
            }
            node->getCols().emplace_back(new_cols[i]);
        } else if (rc(i) > -MAX_POSSIBLE_NUM_R1CS_FOR_VERTEX * dim * RC_TOLERANCE) {
            auto &col = new_cols[i];
            for (auto j: col.col_seq) cout << j << " ";
            cout << " | " << col.forward_concatenate_pos << endl;
            cout << "col " << i << " is not allowed! The rc= " << rc(i) << endl;
            throw runtime_error("wrong pricing column");
        }
#ifdef CHECK_PRICING_LABELS
	if (abs(seq_rc[new_cols[i].col_seq] - rc(i)) > TOLERANCE) {
	  cout << "seq_rc= " << seq_rc[new_cols[ i].col_seq] << " rc= " << rc(i) << endl;
	  throw runtime_error("seq_rc is not equal to rc");
	}
#endif
    }

#ifdef CHECK_PRICING_LABELS
  seq_rc.clear();
#endif

    solver_beg[ccnt_cnt] = nzcnt;

    ccnt = ccnt_cnt;
    if (!ccnt) return;

    safe_solver(node->getSolver().XaddVars(ccnt_cnt,
        nzcnt,
        solver_beg.data(),
        solver_ind.data(),
        solver_val.data(),
        solver_obj.data(),
        nullptr,
        nullptr,
        nullptr,
        nullptr))
    safe_solver(node->getSolver().updateModel())
    safe_solver(node->getSolver().getNumCol(&num_col))

    faster_dual_call(FasterDual::addCols2SubPricingModel(ccnt_cnt,
        nzcnt,
        solver_beg,
        solver_ind,
        solver_val,
        solver_obj);)
    dual_smoothing_call(DualSmoothing::addCol2LPMatrix(ccnt_cnt, solver_beg, solver_ind, solver_val, solver_obj))

    combine_faster_call({
        FasterDual::addCols2SubPricingModel(ccnt_cnt,
            nzcnt,
            solver_beg,
            solver_ind,
            solver_val,
            solver_obj);
        DualSmoothing::addCol2LPMatrix(ccnt_cnt, solver_beg, solver_ind, solver_val, solver_obj);
        })
#ifdef INSPECT_COLUMNS
  inspectColumns(this, node, resource, rc, ccnt, num_col);
#endif
}

void cleanNegativeTuple(vector<tuple<Label *, Label *, double> > &negative_rc_label_tuple) {
    if (negative_rc_label_tuple.empty()) return;
    std::sort(negative_rc_label_tuple.begin(), negative_rc_label_tuple.end(),
              [](const std::tuple<Label *, Label *, double> &a, const std::tuple<Label *, Label *, double> &b) {
                  return get<2>(a) < get<2>(b);
              });
    int cnt = 0;
    for (auto it = negative_rc_label_tuple.begin() + 1; it != negative_rc_label_tuple.end();) {
        if (abs(get<2>(*it) - get<2>(*(it - 1))) < TOLERANCE) {
            it = negative_rc_label_tuple.erase(it);
        } else {
            ++it;
            ++cnt;
            if (cnt >= Config::MaxNumRoutesInExact) break;
        }
    }
    if (negative_rc_label_tuple.size() > Config::MaxNumRoutesInExact)
        negative_rc_label_tuple.resize(Config::MaxNumRoutesInExact);
}

void CVRP::writeColumnsInPricingPool() {
    Label *p;
    if (checkPricingPool()) reallocatePricingPool();

    cleanNegativeTuple(negative_rc_label_tuple);
    new_cols.clear();
    new_cols.resize(negative_rc_label_tuple.size());
    int col_idx = 0;

    for (auto &i: negative_rc_label_tuple) {
        auto &ki = get<0>(i);
        auto &kj = get<1>(i);
        auto &col = new_cols[col_idx].col_seq;
        auto &res = new_cols[col_idx].main_res;
        col.reserve(dim);
        res.reserve(dim);
        p = ki;
        while (p && p->end_vertex) {
            col.emplace_back(p->end_vertex);
            res.emplace_back(p->res.first_res);
            p = p->p_label;
        }
        reverse(col.begin(), col.end());
        reverse(res.begin(), res.end());
        new_cols[col_idx].forward_concatenate_pos = (int) col.size() - 1;

        if (kj) {
            p = kj;
            while (p && p->end_vertex) {
                col.emplace_back(p->end_vertex);
                res.emplace_back(p->res.first_res);
                p = p->p_label;
            }
        }
#ifdef CHECK_PRICING_LABELS
	seq_rc[col] = get<2>(i);
#endif
        ++col_idx;
    }
}

void CVRP::initializeLabels(BbNode *node,
                            int mode,
                            bool if_resetLabelPoint,
                            tuple<bool, int, bool> control_cleanAllPtr) {
    if (if_resetLabelPoint) {
#ifdef SYMMETRY_PROHIBIT
	idx_glo = 2 * dim - 1;
#else
        idx_glo = dim;
#endif
        rc_std = RC_TOLERANCE;
        dual_smoothing_call(DualSmoothing::changeRCStd())

        read_dual_sol_call(ReadDualSol::changeRCStd(node, rc_std, pi4_labeling);)
        num_dominance_checks = 0;
        negative_rc_label_tuple.clear();
    }

    if (get<0>(control_cleanAllPtr)) {
        cleanAllPointers(node, get<1>(control_cleanAllPtr), get<2>(control_cleanAllPtr));
    }

    can_leave_depot_forward.reset();
    for (auto j: node->all_forward_buckets[0][0].bucket_arcs) can_leave_depot_forward.set(j);

    symmetry_prohibit_call(
        can_leave_depot_backward.reset();
        for (auto j : node->all_backward_buckets[0][0].bucket_arcs) can_leave_depot_backward.set(j);)

    if (mode == 1) {
        for (int i = 1; i < dim; ++i) {
            if (!can_leave_depot_forward.test(i)) continue;
            auto &new_label = all_label[i];
            new_label.rc = chg_cost_mat4_vertex[0][i];
            new_label.is_extended = false;
            updateR1CStates(new_label.rc, new_label.r1c, all_label->r1c, 0, i);
            int bin = int(new_label.res.first_res / step_size);
            auto &bucket = label_array_in_forward_sense[i][bin];
            bucket.push_front(all_label + i);
            auto &bucket2 = if_exist_extra_labels_in_forward_sense[i][bin];
            bucket2.first[bucket2.second++] = all_label + i;
            initializeExtraLabels(all_label + i, true);
        }
    }
    symmetry_prohibit_call(else if (mode == 2) {
        int max_num = 2 * dim - 1;
        for (int i = dim; i < max_num; ++i) {
        int point = i - dim + 1;
        if (!can_leave_depot_backward.test(point)) continue;
        auto &new_label = all_label[i];
        new_label.rc = chg_cost_mat4_vertex[0][point];
        new_label.is_extended = false;
        updateR1CStates(new_label.rc, new_label.r1c, all_label->r1c, 0, point);
        int bin = int(new_label.res.first_res / step_size);
        auto &bucket = label_array_in_backward_sense[point][bin];
        bucket.push_front(all_label + i);
        auto &bucket2 = if_exist_extra_labels_in_backward_sense[point][bin];
        bucket2.first[bucket2.second++] = all_label + i;
        initializeExtraLabels(all_label + i, false);
        }
        })
}

void CVRP::addPathByRC(double path_rc, Label *ki, Label *kj, int num) {
    if (path_rc < rc_std) {
        negative_rc_label_tuple.emplace_back(ki, kj, path_rc);
        if (negative_rc_label_tuple.size() >= num)
            rc_std = get<2>(negative_rc_label_tuple[negative_rc_label_tuple.size() - num]);
    }
}

void inspectColumns(CVRP *cvrp, BbNode *node, const ResTuple &resource, const RowVectorXd &rc, int ccnt, int num_col) {
    for (int i = 0; i < ccnt; ++i) {
        auto &col = cvrp->new_cols[i];
        auto &col_seq = col.col_seq;
        auto res = ResTuple{};
        int b4 = 0;
        for (int j: col_seq) {
            if (!cvrp->increaseMainResourceConsumption(res, res, b4, j)) {
                throw runtime_error("infeasible solution in extension");
            }
            b4 = j;
        }
        if (res.first_res > resource.first_res
        ) {
            throw runtime_error("infeasible solution");
        }
        if (!cvrp->can_leave_depot_forward.test(col_seq.front())) {
            cout << "col " << i << " is not allowed! The rc= " << rc(i) << endl;
            for (auto j: col_seq) cout << j << " ";
            cout << " | " << col.forward_concatenate_pos << endl;
            throw runtime_error("infeasible solution");
        }
        if (!cvrp->can_leave_depot_backward.test(col_seq.back())) {
            cout << "col " << i << " is not allowed! The rc= " << rc(i) << endl;
            for (auto j: col_seq) cout << j << " ";
            cout << " | " << col.forward_concatenate_pos << endl;
            throw runtime_error("infeasible solution");
        }
    }
}

void checkRC(BbNode *node,
             const std::vector<std::tuple<Label *, Label *, double> > &negative_rc_label_tuple,
             const std::vector<SequenceInfo> &new_cols,
             const RowVectorXd &rc,
             const RowVectorXd &pi4_labeling,
             const sparseColMatrixXd &mat,
             int num_row) {
    for (int i = 0; i < negative_rc_label_tuple.size(); ++i) {
        auto &label_tuple = negative_rc_label_tuple[i];
        auto rc_n = get<2>(label_tuple);
        if (abs(rc_n - rc(i)) > 0.01) {
            cout << "rc_n= " << rc_n << " rc= " << rc(i) << endl;
            auto rc_dif = rc_n - rc(i);
            for (int j = 0; j < num_row; ++j) {
                if (abs(pi4_labeling[j] - rc_dif) < 0.01) {
                    cout << "pi4_labeling[j]= " << pi4_labeling[j] << " rc_dif= " << rc_dif << endl;
                    cout << "j= " << j << endl;
                    for (auto &r1c: node->r1cs) {
                        if (r1c.idx_r1c == j) {
                            cout << " arc_mem= ";
                            for (auto &it: r1c.arc_mem) {
                                for (auto &it2: it.first) cout << it2 << " ";
                                cout << " to " << it.second << endl;
                            }
                            cout << " info= ";
                            for (auto &it: r1c.info_r1c.first)cout << it << " ";
                            cout << endl;
                            cout << " col seq= ";
                            for (auto &it: new_cols[i].col_seq)cout << it << " ";
                            cout << endl;
                            cout << " col res= ";
                            for (auto &it: new_cols[i].main_res)cout << it << " ";
                            cout << endl;
                            cout << "col forward_concatenate_pos= " << new_cols[i].forward_concatenate_pos << endl;
                            cout << " rc col coeff= " << mat.coeff(j, i) << endl;
                        }
                    }
                }
            }
            throw runtime_error("wrong pricing column");
        }
    }
}
