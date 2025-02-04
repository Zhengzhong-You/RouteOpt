#include <iostream>
#include "cvrp.hpp"
#include "template_functors.hpp"
#include "write_node_out.hpp"
#if SOLUTION_TYPE == 1
#include "best_bound_first_branching.hpp"
#elif SOLUTION_TYPE == 2
#include "depth_first_branching.hpp"
#endif

using namespace std;
using namespace chrono;

void CVRP::solveMIP(BbNode *const node, bool if_inEnu) {
    auto beg = high_resolution_clock::now();

    vector<char> xtype(num_col, SOLVER_BINARY);
    safe_solver(node->getSolver().setEnvCutoff(BaseBranching::ub + round_up_tolerance))
    safe_solver(node->getSolver().setVTypeArray(0, num_col, xtype.data()))

    if (if_inEnu) {
        if_mip_enumeration_suc = true;
        if (num_col > Config::MinNumRoute4MIP) {
            safe_solver(node->getSolver().setEnvTimeLimit(Config::MIPInEnumerationTimeLimit))
        }
    }

    safe_solver(node->getSolver().mipOptimize())

    int status;
    safe_solver(node->getSolver().getStatus(&status))
    if (status == SOLVER_MIP_INFEASIBLE) {
        cout << "Model is infeasible after pre_solving!" << endl;
        cout << "status: " << status << endl;
        lp_val = SOLVER_INFINITY;
    } else {
        safe_solver(node->getSolver().getObjVal(&lp_val))
    }
    auto end = high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    cout << "MIP used " << eps << " seconds!" << endl;

    safe_solver(node->getSolver().setEnvCutoff(1e100))

    if (if_inEnu) {
        if (status == SOLVER_TIME_LIMIT) {
            if_mip_enumeration_suc = false;
            cout << "MIP time limit reached!" << endl;
            MaxNumRoute4Mip = max(int(Config::MIPInEnumerationTimeLimitChgFactor * num_col), Config::MinNumRoute4MIP);
            cout << "MaxNumRoute4Mip = " << MaxNumRoute4Mip << endl;
        }
        safe_solver(node->getSolver().setEnvTimeLimit(MAX_TIME_LIMIT_FOR_MIP))
        if (BaseBranching::ub > lp_val + TOLERANCE) {
            vector<double> X(num_col);
            safe_solver(node->getSolver().getX(0, num_col, X.data()))
            BaseBranching::ub = lp_val;
            cout << "solve MIP get " << BaseBranching::ub << endl;
            updateIPOptSol(node, X);
        }
        if (status == SOLVER_TIME_LIMIT) {
            fill_n(xtype.begin(), num_col, SOLVER_CONTINUOUS);
            safe_solver(node->getSolver().setVTypeArray(0, num_col, xtype.data()))
            safe_solver(node->getSolver().reoptimize())
        }
    } else {
        if (BaseBranching::ub > lp_val + TOLERANCE) {
            vector<double> X(num_col);
            safe_solver(node->getSolver().getX(0, num_col, X.data()))
            BaseBranching::ub = lp_val;
            cout << "solve MIP get " << BaseBranching::ub << endl;
            updateIPOptSol(node, X);
        }
        fill_n(xtype.begin(), num_col, SOLVER_CONTINUOUS);
        safe_solver(node->getSolver().setVTypeArray(0, num_col, xtype.data()))
        safe_solver(node->getSolver().reoptimize())
    }
}

void CVRP::enumerateMIP(BbNode *&node) {
    bool if_ban_enu = false;
    setting_2_call(banEnumeration(if_ban_enu))
    if (if_ban_enu)return;

    determineIfEnumeration(node);
    if (!final_decision_4_enumeration) return;
    final_decision_4_enumeration = false;

    cout << BIG_PHASE_SEPARATION;
    cout << "try EnumerateRoutes..." << endl;

    if (abs(meet_point_resource_in_bi_dir_enu) < TOLERANCE)
        meet_point_resource_in_bi_dir_enu = meet_point_resource_in_bi_dir;
    bool if_succeed;
    double time_labeling;

    priceLabeling(node, optimal_dual_vector);
    getRank1DualsInCG(node, optimal_dual_vector);

    auto beg = high_resolution_clock::now();

    if_succeed = enumerateRoutes(node);

    auto end = high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    time_labeling = eps;

    if (if_succeed) {
        cout << "enumeration time= " << time_labeling << " and succeed!" << endl;
    } else {
        if (if_short_memory) cout << "short memory!" << endl;
        if (if_roll_back) cout << "roll back!" << endl;
        cout << "enumeration time= " << time_labeling << " but failed!" << endl;
    }
    adjustEnumerationStdBucketArcs(node, if_succeed);
    double gap = (BaseBranching::ub - node->getCurrentNodeVal()) / BaseBranching::ub;
    if (if_succeed) {
        if (gap > gap_tolerance4_arc_elimination_n_enumeration) {
            gap_tolerance4_arc_elimination_n_enumeration = gap;
        }

#ifdef SYMMETRY_PROHIBIT
	double dif = abs(num_forward_labels_in_enu - num_backward_labels_in_enu);
	double over = dif / min(num_forward_labels_in_enu, num_backward_labels_in_enu);
	if (RobustControl::if_fix_resource_point) {
	  meet_point_resource_in_bi_dir_enu = meet_point_resource_in_bi_dir;
	} else {
	  if (over > Config::NumberOfOverLabelsInMeetPoint) {
		if (num_forward_labels_in_enu > num_backward_labels_in_enu) {
		  meet_point_resource_in_bi_dir_enu *= (1 - Config::MeetPointFactor);
		} else {
		  meet_point_resource_in_bi_dir_enu *= (1 + Config::MeetPointFactor);
		}
		cout << "meet_point_resource_in_bi_dir_enu= " << meet_point_resource_in_bi_dir_enu << endl;
	  }
	}
#endif

        if_enumeration_suc = true;
        write_node_out_call(
            cout<<"num_explored_nodes= "<<BaseBranching::num_explored_nodes<<endl;
            if(BaseBranching::num_explored_nodes>1) {
            WriteNodeOut::writeNodeOut(node, LARGE_MEMORY_USE);
            --BaseBranching::num_explored_nodes;
            goto PASS;
            })
        beg = high_resolution_clock::now();
#if SOLUTION_TYPE == 1
        BaseBranching::prepareEnuTree(BestBoundFirstBranching::sub_bbt);
        BestBoundFirstBranching::solve(BestBoundFirstBranching::sub_bbt);
#elif SOLUTION_TYPE == 2
	BaseBranching::prepareEnuTree(DepthFirstBranching::sub_bbt);
	DepthFirstBranching::solve(DepthFirstBranching::sub_bbt);
#endif
        end = high_resolution_clock::now();
    PASS:;
    } else {
        if_enumeration_suc = false;
    }
    adjustEnumerationStdGap(gap, if_succeed);
    if (if_succeed) {
        dynamic_call(Dynamics::calculateF(duration<double>(end - beg).count()))
    }
    cout << BIG_PHASE_SEPARATION;
}

bool CVRP::enumerateRoutes(BbNode *const node) {
    if_roll_back = false;
    if_short_memory = false;
    int tmp_labels = Config::MaxNumLabelInEnumeration;
    int tmp_routes = Config::MaxNumRouteInEnumeration;
    double time_limit = Config::HardTimeThresholdInAllEnumeration;
    if (if_force_enumeration_suc) {
        Config::MaxNumLabelInEnumeration = numeric_limits<int>::max();
        Config::MaxNumRouteInEnumeration = numeric_limits<int>::max();
        Config::HardTimeThresholdInAllEnumeration = numeric_limits<double>::max();
    }
    int num_routes_now = 0;
    auto &cost_m = node->cost_for_columns_in_enumeration_column_pool;
    auto &ptr = node->index_columns_in_enumeration_column_pool;

    int index;
    opt_gap = calculateOptimalGap(node);
    unordered_map<yzzLong, tuple<Label *, Label *, double> > Tags;
    Tags.reserve(Config::MaxNumRouteInEnumeration);

    auto copy_Forward_bucket = new vector<Label *> *[dim];
#ifdef SYMMETRY_PROHIBIT
  auto copy_Backward_bucket = new vector<Label *> *[dim];//these memory will be freed inside the forward function
#endif

    for (int i = 0; i < dim; ++i) {
        copy_Forward_bucket[i] = new vector<Label *>[num_buckets_per_vertex];
#ifdef SYMMETRY_PROHIBIT
	copy_Backward_bucket[i] = new vector<Label *>[num_buckets_per_vertex];
#endif
        for (int b = 0; b < num_buckets_per_vertex; ++b) {
            copy_Forward_bucket[i][b].assign(label_array_in_forward_sense[i][b].begin(),
                                             label_array_in_forward_sense[i][b].end());
#ifdef SYMMETRY_PROHIBIT
	  copy_Backward_bucket[i][b].assign(label_array_in_backward_sense[i][b].begin(),
										label_array_in_backward_sense[i][b].end());
#endif
        }
    }

    auto beg = high_resolution_clock::now();

    int status =
#ifdef SYMMETRY_PROHIBIT
	  enumerateHalfwardRoutes<true, false>(node, Tags,
										   copy_Backward_bucket, num_routes_now);
#else
            enumerateHalfwardRoutes<true, true>(node, Tags, copy_Forward_bucket, num_routes_now);
#endif

    auto end = high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
#if VERBOSE_MODE == 1
    cout << "Half Forward time= " << eps << "s" << endl;
#endif

    Config::MaxNumLabelInEnumeration = tmp_labels;
    Config::MaxNumRouteInEnumeration = tmp_routes;
    Config::HardTimeThresholdInAllEnumeration = time_limit;

    if (status || if_short_memory || if_roll_back) {
        if (status == 1)
            cout << "the number of labels in Forward reached its limit!" << endl;
        return false;
    }

    write_node_out_call(if(BaseBranching::num_explored_nodes>1)return true)

#ifdef SYMMETRY_PROHIBIT
  beg = high_resolution_clock::now();

  status = enumerateHalfwardRoutes<false, false>(node,
												 Tags,
												 copy_Forward_bucket,
												 num_routes_now);

  end = high_resolution_clock::now();
  eps = duration<double>(end - beg).count();
#if VERBOSE_MODE == 1
  cout << "Half Backward time= " << eps << "s" << endl;
#endif

  if (status || if_short_memory || if_roll_back) {
	if (status == 1)
	  cout << "the number of labels in Backward reached its limit!" << endl;
	return false;
  }
#endif

    beg = high_resolution_clock::now();

    status = concatenateRoutesPriorForwardInEnumeration(node, Tags, num_routes_now);

    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();
#if VERBOSE_MODE == 1
    cout << "Concatenate time= " << eps << "s" << endl;
#endif

    if (status) {
        cout << "the number of routes reached its limit!" << endl;
        return false;
    }

    if_in_enu_state = true;
    cost_m.resize(num_routes_now);
    ptr.resize(num_routes_now);

    auto num = size_t(num_routes_now * aver_route_length);
    reallocatePricingPool(num);
    auto PricingWarning = (size_t) (0.9 * (double) mem4_pricing);
    pool_beg4_pricing = 0;
    index = 0;
    Label *ki, *li, *p;
    vector<int> seq(dim + 1);
    for (auto &tag: Tags) {
        ki = get<0>(tag.second);
        li = get<1>(tag.second);

        cost_m[index] = get<2>(tag.second);
        ptr[index++] = pool_beg4_pricing;

        int cnt = 0;
        p = ki;
        while (p) {
            seq[cnt++] = p->end_vertex;
            p = p->p_label;
        }
        for (int k = 0; k < cnt; ++k) {
            col_pool4_pricing[pool_beg4_pricing++] = seq[cnt - 1 - k];
        }
        if (li) {
            p = li;
            while (p) {
                col_pool4_pricing[pool_beg4_pricing++] = p->end_vertex;
                p = p->p_label;
            }
        } else {
            col_pool4_pricing[pool_beg4_pricing++] = 0;
        }
#if SOLVER_VRPTW == 1
	double local_cap = 0;
	for (size_t i = ptr[index - 1] + 1; i < pool_beg4_pricing; ++i) {
	  local_cap += demand[col_pool4_pricing[i]];
	}
	if (local_cap > cap + TOLERANCE) {
	  --index;
	  pool_beg4_pricing = ptr[index];
	  continue;
	}
#endif
        resizePoolWarning(PricingWarning);
    }

    num_routes_now = index;
    cost_m.conservativeResize(num_routes_now);
    ptr.conservativeResize(num_routes_now);
    node->size_enumeration_col_pool = num_routes_now;
    node->valid_size = num_routes_now;
    cout << "number of routes are " << num_routes_now << endl;

    deleteColumnByNGMemory(node, 1, true);
#if SOLVER_VRPTW == 1
  for (auto &rcc : node->rccs)rcc.if_keep = false;
  cleanColumnsCapInfeasible(node);
#endif
    deleteBranchCutsAndR1C1s(node);
    recoverR1CsInEnum(node);
    recoverRCCsInEnum(node);
    generateVertex2IndexColsAndEdge2IndexCols(node);
    return true;
}

int CVRP::concatenateRoutesPriorForwardInEnumeration(BbNode *node,
                                                     unordered_map<yzzLong, tuple<Label *, Label *, double> > &Tags,
                                                     int &num_routes_now) {
    int status = 0;
    double path_rc, path_cost;
    yzzLong tmp_PI;
#ifdef SYMMETRY_PROHIBIT
  populateRC2TillThisBinNRC2Bin<false>();
#else
    populateRC2TillThisBinNRC2Bin<true>(); //use function here!
#endif
    for (auto &label_list: concatenate_labels_in_forward_cg) {
        int i = label_list.first.first;
        int j = label_list.first.second;
        auto &label_vec = label_list.second;
        for (auto &pr: label_vec) {
            auto &ki = pr.first;
            auto &tmp_Resource = pr.second;
            double tmp_rc = ki->rc + chg_cost_mat4_vertex[i][j];
#ifdef SYMMETRY_PROHIBIT
	  int arr_bj = int((tmp_Resource.first_res) / step_size);
	  if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap)
		continue;

	  if (rc2_bin_in_backward_sense[j][arr_bj] + tmp_rc < opt_gap &&
		  (std::find(node->all_backward_buckets[j][arr_bj].bucket_arcs.begin(),
					 node->all_backward_buckets[j][arr_bj].bucket_arcs.end(), i)
			  != node->all_backward_buckets[j][arr_bj].bucket_arcs.end())) {
		auto &label_arr = label_array_in_backward_sense[j][arr_bj];
#else
            int arr_bj = int((resource.first_res - tmp_Resource.first_res) / step_size);
            if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
                continue;

            if (rc2_bin_in_forward_sense[j][arr_bj] + tmp_rc < opt_gap
            ) {
                auto &label_arr = label_array_in_forward_sense[j][arr_bj];
#endif
            }
#ifdef SYMMETRY_PROHIBIT
	  for (; arr_bj < num_buckets_per_vertex; ++arr_bj) {
		if (rc2_till_this_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap)
		  break;

		if (rc2_bin_in_backward_sense[j][arr_bj] + tmp_rc > opt_gap ||
			(std::find(node->all_backward_buckets[j][arr_bj].bucket_arcs.begin(),
					   node->all_backward_buckets[j][arr_bj].bucket_arcs.end(), i)
				== node->all_backward_buckets[j][arr_bj].bucket_arcs.end()))
		  continue;

		auto &label_arr = label_array_in_backward_sense[j][arr_bj];
#else
            for (; arr_bj >= 0; --arr_bj) {
                if (rc2_till_this_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap)
                    break;

                if (rc2_bin_in_forward_sense[j][arr_bj] + tmp_rc > opt_gap
                )
                    continue;

                auto &label_arr = label_array_in_forward_sense[j][arr_bj];
#endif
                for (auto &kj: label_arr) {
                    path_rc = kj->rc + tmp_rc;
                    if (path_rc > opt_gap) break;

#ifdef SYMMETRY_PROHIBIT
		  if (tellResTupleRelations<'p'>(tmp_Resource, kj->res))continue;
#else
                    if (tellResTupleRelations<'c'>(tmp_Resource, kj->res))continue;
#endif

                    if ((ki->pi & kj->pi).any()) continue;

                    if (!concatenateR1CStates(path_rc, opt_gap, ki->r1c, kj->r1c, i, j))goto here2;

                    path_cost = ki->cost + cost_mat4_vertex[i][j] + kj->cost;
                    tmp_PI = ki->pi | kj->pi;
                    if (Tags.find(tmp_PI) == Tags.end()) {
                        Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
                        ++num_routes_now;
                        if (num_routes_now > Config::MaxNumRouteInEnumeration) {
                            status = 1;
                            goto QUIT;
                        }
                    } else if (get<2>(Tags[tmp_PI]) > path_cost) {
                        Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
                    }
                here2:;
                }
            }
        }
    }
QUIT:
#if VERBOSE_MODE == 1
    cout << "num_routes_now= " << num_routes_now << endl;
#endif
    return status;
}

int CVRP::generateColumnsByInspection(BbNode *node,
                                      bool if_only_need_value) {
    if (node->size_enumeration_col_pool == 0) return 0;
    int size_pool = node->size_enumeration_col_pool;
    auto &deleted_columns_in_enumeration_pool = node->deleted_columns_in_enumeration_pool;
    auto &mat = node->matrix_in_enumeration;

    vector<double> pi(num_row);
    safe_solver(node->getSolver().getDual(0, num_row, pi.data()))
    safe_solver(node->getSolver().getObjVal(&lp_val))

    RowVectorXd rc = node->cost_for_columns_in_enumeration_column_pool;

    int num = 0;
    for (auto &it: mat) {
        RowVectorXd dual(it.rows());
        for (int i = 0; i < it.rows(); ++i)dual(i) = pi[num++];
        rc -= dual * it;
    }

    vector<int> col_to_be_added;
    col_to_be_added.reserve(size_pool);

    int cnt = 0;
    if (if_only_need_value) {
        for (int i = 0; i < size_pool; ++i) {
            if (deleted_columns_in_enumeration_pool[i]) continue;
            if (rc(i) < RC_TOLERANCE) {
                col_to_be_added.emplace_back(i);
                ++cnt;
                if (cnt == Config::MaxNumRoutesInExact) break;
            }
        }
    } else {
        for (int i = 0; i < size_pool; ++i) {
            if (deleted_columns_in_enumeration_pool[i]) continue;
            if (rc(i) < RC_TOLERANCE) {
                col_to_be_added.emplace_back(i);
                deleted_columns_in_enumeration_pool[i] = true;
                ++cnt;
                if (cnt == Config::MaxNumRoutesInExact) break;
            }
        }
    }

    if (col_to_be_added.empty()) {
        if (!if_only_need_value) {
            node->getCurrentNodeVal() = lp_val;
            opt_gap = calculateOptimalGap(node);
            for (int i = 0; i < size_pool; ++i) {
                if (rc(i) > opt_gap) {
                    deleted_columns_in_enumeration_pool[i] = true;
                }
            }
            optimal_dual_vector = pi;
            cleanIndexColForNode(node, true);
            regenerateEnumMat(node, nullptr);
        }
    } else {
        addColumnsByInspection(node, col_to_be_added);
    }
    return cnt;
}

void CVRP::generateVertex2IndexColsAndEdge2IndexCols(BbNode *node) {
    int size_pool = node->size_enumeration_col_pool, curr_node;
    if (size_pool == 0) return;
    auto &ptr = node->index_columns_in_enumeration_column_pool;
    sparseRowMatrixXd mat(num_row, size_pool);

    vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(size_t(double(num_row) * double(size_pool) * 0.1));

    auto &deleted_columns_in_enumeration_pool = node->deleted_columns_in_enumeration_pool;
    deleted_columns_in_enumeration_pool.resize(size_pool, false);

    for (int i = 0; i < size_pool; ++i) {
        for (auto j = ptr[i] + 1;; ++j) {
            curr_node = col_pool4_pricing[j];
            if (!curr_node) break;
            triplets.emplace_back(curr_node - 1, i, 1);
        }
    }

    node->basic_matrix.resize(real_dim, size_pool);
    node->basic_matrix.setFromTriplets(triplets.begin(), triplets.end());
    row_basic_matrix.resize(real_dim, size_pool);
    row_basic_matrix.setFromTriplets(triplets.begin(), triplets.end());

    for (int i = 0; i < size_pool; ++i) triplets.emplace_back(real_dim, i, 1);

    auto beg = high_resolution_clock::now();
    buildRCCInEnuMatrix(node, triplets);
    auto end = high_resolution_clock::now();
    auto eps = duration<double>(end - beg).count();
    verbose_call(cout << "buildRCCInEnuMatrix used " << eps << " seconds!" << endl;)

    beg = high_resolution_clock::now();
    buildAllR1CInEnuMatrix(node, triplets);
    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();
    verbose_call(cout << "buildAllR1CInEnuMatrix used " << eps << " seconds!" << endl;)

    beg = high_resolution_clock::now();
    mat.setFromTriplets(triplets.begin(), triplets.end());
    node->matrix_in_enumeration.push_back(std::move(mat));
    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();
    verbose_call(cout << "setFromTriplets used " << eps << " seconds!" << endl;)

    beg = high_resolution_clock::now();
    rmColByBranchInEnuMatrix(node, node->deleted_columns_in_enumeration_pool, false, node->getBrCs());
    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();
    verbose_call(cout << "rmColByBranchInEnuMatrix used " << eps << " seconds!" << endl;)

    beg = high_resolution_clock::now();
    regenerateEnumMat(node, nullptr, false);
    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();
    verbose_call(cout << "regenerateEnumMat used " << eps << " seconds!" << endl;)
    if (size_pool != node->size_enumeration_col_pool) cout << "now size= " << node->size_enumeration_col_pool << endl;
}

void CVRP::addColumnsByInspection(BbNode *const node, const vector<int> &Col_added) {
    if (Col_added.empty()) return;
    auto &ptr = node->index_columns_in_enumeration_column_pool;
    auto &ptr_cost = node->cost_for_columns_in_enumeration_column_pool;
    auto &mat = node->matrix_in_enumeration;

    int idx = num_col;
    node->getCols().resize(num_col + Col_added.size());
    for (auto col: Col_added) {
        auto &seq = node->getCols()[idx++].col_seq; //no other information is needed!
        for (auto j = ptr[col] + 1;; ++j) {
            int curr_node = col_pool4_pricing[j];
            if (!curr_node)break;
            seq.emplace_back(curr_node);
        }
    }

    vector<int> if_col_added(node->size_enumeration_col_pool, -1);
    int cnt = 0;
    for (auto col: Col_added) {
        if_col_added[col] = cnt++;
    }

    double nonZeros = 0;
    Eigen::SparseMatrix<double, Eigen::ColMajor> tmp_mat(num_row, (int) Col_added.size());
    vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(size_t((double) Col_added.size() * num_row * 0.1));
    size_t num = 0;
    for (auto &it: mat) {
        nonZeros += (double) it.nonZeros();
        for (int i = 0; i < it.rows(); ++i) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_it(it, i); inner_it; ++inner_it) {
                if (if_col_added[inner_it.col()] != -1) {
                    triplets.emplace_back(num, if_col_added[inner_it.col()], inner_it.value());
                }
            }
            ++num;
        }
    }

    safe_eigen(tmp_mat.setFromTriplets(triplets.begin(), triplets.end());)
    vector<double> solver_obj(Col_added.size()), solver_val;
    vector<size_t> solver_beg(Col_added.size() + 1);
    vector<int> solver_ind;
    auto numReserve = size_t(double(Col_added.size()) * aver_route_length);
    solver_ind.reserve(numReserve);
    solver_val.reserve(numReserve);
    int ccnt = 0;
    for (auto col: Col_added) {
        solver_obj[ccnt] = ptr_cost[col];
        solver_beg[ccnt] = solver_ind.size();
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(tmp_mat, ccnt); it; ++it) {
            solver_ind.emplace_back((int) it.row());
            solver_val.emplace_back(it.value());
        }
        ++ccnt;
    }
    solver_beg[ccnt] = solver_ind.size();

    safe_solver(node->getSolver().XaddVars(ccnt,
        solver_ind.size(),
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
}

void CVRP::regenerateEnumMat(BbNode *node, BbNode *node2, bool if_force) {
    auto &deleted_columns_in_enumeration_pool =
            node2 ? node2->deleted_columns_in_enumeration_pool : node->deleted_columns_in_enumeration_pool;
    int size_pool = (int) node->size_enumeration_col_pool;
    if (!size_pool) {
        node->deleted_columns_in_enumeration_pool.clear();
        return;
    }
    BbNode *out_node;
    int del_size = 0;
    if (!node2) {
        out_node = node;
        for (int i = 0; i < size_pool; ++i) {
            if (deleted_columns_in_enumeration_pool[i]) ++del_size;
        }
        if (del_size == 0) return;
        node->valid_size = size_pool - del_size;
        if (!if_force) {
            auto left = double(node->size_enumeration_col_pool - del_size) / node->size_enumeration_col_pool;
            if (left > Config::LeftThresholdRCFixing4EnumerationPool) {
                cout << "columns pending deletion= " << del_size << endl;
                return;
            } else {
                cout << "delete columns= " << del_size << endl;
            }
        }
    } else {
        out_node = node2;
        for (int i = 0; i < size_pool; ++i) {
            if (deleted_columns_in_enumeration_pool[i]) ++del_size;
        }
        node2->valid_size = size_pool - del_size;
    }

    int colIndex = 0;
    vector<int> new_col_map(size_pool, -1);
    for (int i = 0; i < size_pool; ++i) {
        if (!deleted_columns_in_enumeration_pool[i]) {
            new_col_map[i] = colIndex++;
        }
    }

    if (cstr_index.empty()) {
        findNonActiveCuts(node, false);
        if (cstr_index.empty()) {
            cstr_index.resize(num_row);
            iota(cstr_index.begin(), cstr_index.end(), 0);
        }
    }

    int new_size_pool = size_pool - del_size;
    sparseRowMatrixXd tmpMatrix(num_row, new_size_pool);

    int rowIndex = 0;
    vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(size_t(double(num_row) * double(new_size_pool) * 0.1));

    for (auto &it: node->matrix_in_enumeration) {
        for (int i = 0; i < it.rows(); ++i) {
            if (cstr_index[rowIndex] != -1) {
                int j = cstr_index[rowIndex];
                if (j >= num_row)
                    throw runtime_error(
                        "row index out of range!" + to_string(j) + " " + to_string(num_row));
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_it(it, i); inner_it; ++
                     inner_it) {
                    if (new_col_map[inner_it.col()] != -1) {
                        triplets.emplace_back(j, new_col_map[inner_it.col()], inner_it.value());
                    }
                }
            }
            ++rowIndex;
        }
    }

    safe_eigen(tmpMatrix.setFromTriplets(triplets.begin(), triplets.end());)

    out_node->matrix_in_enumeration.clear();
    out_node->matrix_in_enumeration.push_back(std::move(tmpMatrix));

    colIndex = 0;
    for (int i = 0; i < out_node->size_enumeration_col_pool; ++i) {
        if (!out_node->deleted_columns_in_enumeration_pool[i]) {
            out_node->cost_for_columns_in_enumeration_column_pool[colIndex] =
                    out_node->cost_for_columns_in_enumeration_column_pool[i];
            out_node->index_columns_in_enumeration_column_pool[colIndex++] =
                    out_node->index_columns_in_enumeration_column_pool[i];
        }
    }
    out_node->size_enumeration_col_pool = new_size_pool;
    out_node->cost_for_columns_in_enumeration_column_pool.conservativeResize(out_node->size_enumeration_col_pool);
    out_node->index_columns_in_enumeration_column_pool.conservativeResize(out_node->size_enumeration_col_pool);
    out_node->deleted_columns_in_enumeration_pool.assign(out_node->size_enumeration_col_pool, false);

    createBasicMatrix(out_node);

    cstr_index.clear();
}

void CVRP::solveLPByInspection(BbNode *const node, bool if_only_need_value,
                               bool if_heuristic, bool if_record_sol) {
    int ccnt = 0, iter_exact = 0, old_ncol = num_col;

    time_point<high_resolution_clock, duration<long long, ratio<1L, 1000000000LL> > > mt_beg,
            mt_end, spt_beg, spt_end;
    duration<long long int, ratio<1LL, 1000000000LL> > mt_elap{0}, spt_elap{0};
    spt_beg = high_resolution_clock::now();

    optimizeLPForOneIterationInEnum(node);

    int env_method;
    bool if_changed = false;
    safe_solver(node->getSolver().getEnvMethod(&env_method))
    if (env_method != SOLVER_PRIMAL_SIMPLEX) {
        safe_solver(node->getSolver().setEnvMethod(SOLVER_PRIMAL_SIMPLEX))
        if_changed = true;
    }

#if VERBOSE_MODE == 1
    if (if_heuristic)
        cout << "Heuristic phase begin...\n";
    else
        cout << "Exact phase begin...\n";
#endif

    while (true) {
        ++iter_exact;

        spt_beg = high_resolution_clock::now();

        ccnt = generateColumnsByInspection(node, if_only_need_value);

        if (!ccnt) {
            --iter_exact;
            break;
        }

        spt_end = high_resolution_clock::now();
        spt_elap += duration_cast<milliseconds>(spt_end - spt_beg);
        mt_beg = high_resolution_clock::now();

        optimizeLPForOneIterationInEnum(node);

        mt_end = high_resolution_clock::now();
        mt_elap += duration_cast<milliseconds>(mt_end - mt_beg);

#if VERBOSE_MODE == 1
        if (!(iter_exact % PRINT_LABELING_STEP_SIZE)) {
            BaseBranching::glo_end = high_resolution_clock::now();
            BaseBranching::glo_eps = duration<double>(BaseBranching::glo_end - BaseBranching::glo_beg).count();
            printInfoLabeling(iter_exact, num_col - old_ncol, num_col, num_row, double(mt_elap.count()) * 1e-9,
                              double(spt_elap.count()) * 1e-9, BaseBranching::glo_eps,
                              lp_val, BaseBranching::lb, BaseBranching::ub);
            if (!if_only_need_value) {
                cout << "col_pool= " << node->size_enumeration_col_pool << "  remain "
                        << (max_num_enu_col_pool
                                ? (double(node->size_enumeration_col_pool) / max_num_enu_col_pool * 100)
                                : 0)
                        << "%\n";
            }
            spt_elap = duration_values<milliseconds>::zero();
            mt_elap = duration_values<milliseconds>::zero();
            old_ncol = num_col;
        }
#endif
    }

    if (if_record_sol) {
        recordOptimalColumn(node);
    }

    if (if_changed) safe_solver(node->getSolver().setEnvMethod(env_method))
    if (iter_exact % PRINT_LABELING_STEP_SIZE || !iter_exact) {
        BaseBranching::glo_end = high_resolution_clock::now();
        BaseBranching::glo_eps = duration<double>(BaseBranching::glo_end - BaseBranching::glo_beg).count();
        printInfoLabeling(iter_exact, num_col - old_ncol, num_col, num_row, double(mt_elap.count()) * 1e-9,
                          double(spt_elap.count()) * 1e-9, BaseBranching::glo_eps,
                          lp_val, BaseBranching::lb, BaseBranching::ub);
        if (!if_only_need_value) {
            cout << "col_pool= " << node->size_enumeration_col_pool << "  remain "
                    << (max_num_enu_col_pool
                            ? (double(node->size_enumeration_col_pool) / max_num_enu_col_pool * 100)
                            : 0)
                    << "%\n";
        }
    }
}

void CVRP::deleteBranchCutsAndR1C1s(BbNode *const node) {
    if (node->getBrCs().empty() && node->r1cs.empty()) return;
    if (!node->getBrCs().empty()) {
        vector<int> solver_ind(num_col);
        vector<size_t> solver_beg(num_col + 1);
        vector<double> solver_val(num_col);
        set<int> delete_col;
        int *ai_col = new int[num_col];
        int *aj_col = new int[num_col];
        int tmp;
        size_t numnzP;
        int ai, aj;
        for (auto &brc: node->getBrCs()) {
            ai = brc.edge.first;
            aj = brc.edge.second;
            tmp = ai * dim + aj;
            if (brc.br_dir) {
                for (int i = 0; i < num_col; ++i) {
                    ai_col[i] = 0;
                    aj_col[i] = 0;
                }
                if (ai) {
                    numnzP = num_col;
                    safe_solver(node->getSolver().XgetConstraints(&numnzP,
                        solver_beg.data(),
                        solver_ind.data(),
                        solver_val.data(),
                        ai - 1,
                        1))
                    for (size_t i = 0; i < numnzP; ++i)ai_col[solver_ind[i]] = 1;
                }
                numnzP = num_col;
                safe_solver(node->getSolver().XgetConstraints(&numnzP,
                    solver_beg.data(),
                    solver_ind.data(),
                    solver_val.data(),
                    aj - 1,
                    1))
                for (size_t i = 0; i < numnzP; ++i)aj_col[solver_ind[i]] = 1;
                numnzP = num_col;
                safe_solver(node->getSolver().XgetConstraints(&numnzP,
                    solver_beg.data(),
                    solver_ind.data(),
                    solver_val.data(),
                    brc.idx_brc,
                    1))
                for (int j = 0; j < solver_ind[0]; ++j)if (aj_col[j] || ai_col[j]) delete_col.insert(j);
                for (size_t i = 1; i < numnzP; ++i)
                    for (int j = solver_ind[i - 1] + 1; j < solver_ind[i]; ++j)
                        if (aj_col[j] || ai_col[j])delete_col.insert(j);
                for (int j = solver_ind[numnzP - 1] + 1; j < num_col; ++j)
                    if (aj_col[j] || ai_col[j])
                        delete_col.
                                insert(j);
            } else {
            }
        }
        delete_col.erase(0);
        vector<int> col_idx(delete_col.begin(), delete_col.end());
        rmLPCols(node, col_idx);
        safe_solver(node->getSolver().reoptimize(SOLVER_PRIMAL_SIMPLEX))
        delete[]ai_col;
        delete[]aj_col;
    }

    vector<int> solver_ind(num_col);
    int len = 0, keep;
    auto local_cstr_index = new int[num_row];
    vector<int> deleted_cstrs;
    deleted_cstrs.reserve(num_row);
    iota(local_cstr_index, local_cstr_index + num_row, 0);
    for (auto &brc: node->getBrCs()) {
        if (brc.idx_brc == -1) continue;
        keep = brc.idx_brc;
        solver_ind[len++] = keep;
        local_cstr_index[keep] = -1;
        deleted_cstrs.emplace_back(keep);
    }
    for (auto &r1c: node->r1cs) {
        if (r1c.info_r1c.first.size() == 1) {
            keep = r1c.idx_r1c;
            solver_ind[len++] = keep;
            local_cstr_index[keep] = -1;
            deleted_cstrs.emplace_back(keep);
        }
    }
    if (deleted_cstrs.empty()) {
        delete[]local_cstr_index;
        return;
    }
    std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
    int delta = 0;
    auto stop_sign = deleted_cstrs.end() - 1;
    for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
        ++delta;
        for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
    }
    ++delta;
    for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;

    for (auto &rcc: node->rccs) rcc.idx_rcc = local_cstr_index[rcc.idx_rcc];

    for (auto i = node->r1cs.begin(); i < node->r1cs.end();) {
        if (local_cstr_index[i->idx_r1c] == -1) {
            i = node->r1cs.erase(i);
        } else {
            i->idx_r1c = local_cstr_index[i->idx_r1c];
            ++i;
        }
    }

    for (auto &brc: node->getBrCs()) brc.idx_brc = -1;

    safe_solver(node->getSolver().delConstraints(len, solver_ind.data()))
    safe_solver(node->getSolver().reoptimize())
    safe_solver(node->getSolver().getNumCol(&num_col))
    safe_solver(node->getSolver().getNumRow(&num_row))
    safe_solver(node->getSolver().getObjVal(&lp_val))
    delete[]local_cstr_index;
}

void CVRP::recoverR1CsInEnum(BbNode *const node) {
    if (node->r1cs.empty()) return;
    int *ai_col = new int[num_col];
    auto *sum_col = new double[num_col];
    vector<double> val_(num_col);
    int index;
    size_t numnzP;
    vector<size_t> solver_beg(2);
    vector<int> solver_ind, solver_ind2;
    vector<double> solver_val;
    auto numReserve = size_t(num_col * aver_route_length);
    solver_ind.reserve(numReserve);
    solver_ind2.reserve(numReserve);
    solver_val.reserve(numReserve);

    for (auto &r1c: node->r1cs) {
        index = r1c.idx_r1c;
        memset(sum_col, 0, sizeof(double) * num_col);
        vector<int> multi;
        int denominator, rhs;
        Rank1CutsSeparator::getMapPlanInfo(multi, denominator, rhs, r1c.info_r1c.first.size(), r1c.info_r1c.second);
        int count = 0;
        for (auto &i: r1c.info_r1c.first) {
            numnzP = num_col;
            safe_solver(node->getSolver().XgetConstraints(
                &numnzP, solver_beg.data(), ai_col, val_.data(), i - 1, 1))
            for (size_t j = 0; j < numnzP; ++j) sum_col[ai_col[j]] += multi[count];
            ++count;
        }
        for (int i = 0; i < num_col; ++i) {
            int val = int(sum_col[i] / denominator + TOLERANCE);
            if (val) {
                solver_ind.emplace_back(index);
                solver_ind2.emplace_back(i);
                solver_val.emplace_back(val);
            }
        }
    }


    safe_solver(
        node->getSolver().XchangeCoeffs(solver_ind.size(), solver_ind.data(), solver_ind2.data(), solver_val.data()))
    safe_solver(node->getSolver().reoptimize())
    safe_solver(node->getSolver().getObjVal(&lp_val))
#if VERBOSE_MODE == 1
    cout << "after recover rank1c lpval= " << lp_val << endl;
#endif
    delete[]ai_col;
    delete[]sum_col;
}

void CVRP::recoverRCCsInEnum(BbNode *const node) {
    if (node->rccs.empty()) return;

    vector<int> new_rhs(node->rccs.size()); //row idx & rhs
    vector<Eigen::RowVectorXi> v_row_idx_coeff(dim, Eigen::RowVectorXi::Zero((int) node->rccs.size()));
    for (int c = 0; c < (int) node->rccs.size(); ++c) {
        auto &rcc = node->rccs[c];
        int idx = rcc.idx_rcc;
        double k;
        if (rcc.form_rcc) {
            k = (int) rcc.info_rcc_customer.size() - rcc.rhs;
        } else {
            k = (int) rcc.info_rcc_outside_customer.size() - rcc.rhs;
        }
        for (auto i: rcc.info_rcc_customer) {
            v_row_idx_coeff[i][c] = 1;
        }
        new_rhs[c] = (int) k;
        auto sense = SOLVER_GREATER_EQUAL;
        safe_solver(node->getSolver().setRhs(idx, 1, &sense, &k))
    }

    auto size = num_col * node->rccs.size();
    vector<int> cind(size), vind(size);
    vector<double> val(size);
    size = 0;
    for (int c = 0; c < node->rccs.size(); ++c) {
        cind[size] = node->rccs[c].idx_rcc;
        vind[size] = 0;
        val[size] = new_rhs[c];
        ++size;
    }
    for (int i = 1; i < num_col; ++i) {
        Eigen::RowVectorXi sum = Eigen::RowVectorXi::Zero((int) node->rccs.size());
        for (auto j: node->getCols()[i].col_seq) {
            sum += v_row_idx_coeff[j];
        }
        for (int c = 0; c < node->rccs.size(); ++c) {
            vind[size] = i;
            cind[size] = node->rccs[c].idx_rcc;
            if (sum[c]) {
                val[size] = 1;
            } else {
                val[size] = 0;
            }
            ++size;
        }
    }
    safe_solver(node->getSolver().XchangeCoeffs(cind.size(), cind.data(), vind.data(), val.data()))
    safe_solver(node->getSolver().reoptimize())
    safe_solver(node->getSolver().getObjVal(&lp_val))
    verbose_call(cout << "after recover rcc lpval= " << lp_val << endl;)
}

void CVRP::optimizeLPForOneIterationInEnum(BbNode *const node) {
    safe_solver(node->getSolver().reoptimize(SOLVER_PRIMAL_SIMPLEX))
    safe_solver(node->getSolver().getNumCol(&num_col))
    safe_solver(node->getSolver().getObjVal(&lp_val))
    vector<double> X(num_col);
    safe_solver(node->getSolver().getX(0, num_col, X.data()))

    if (tellIfInt(node, X)) {
        if (lp_val + TOLERANCE < BaseBranching::ub) {
            BaseBranching::ub = lp_val;

            updateIPOptSol(node, X);
        }
    }
}

void CVRP::terminateByMIP(BbNode *node) {
#if VERBOSE_MODE == 1
    cout << "b4_col= " << num_col << " ";
#endif
    regenerateEnumMat(node, nullptr, true);
    auto &deleted_columns_in_enumeration_pool = node->deleted_columns_in_enumeration_pool;
    vector<int> added_cols;
    added_cols.reserve(node->size_enumeration_col_pool);
    for (int i = 0; i < node->size_enumeration_col_pool; ++i) {
        if (deleted_columns_in_enumeration_pool[i])continue;
        added_cols.emplace_back(i);
    }
    addColumnsByInspection(node, added_cols);
#if VERBOSE_MODE == 1
    cout << "after_cols= " << num_col << endl;
#endif

    auto beg = high_resolution_clock::now();
    double node_val = node->getCurrentNodeVal();
    solveMIP(node, true);
    if (if_mip_enumeration_suc) {
        auto end = high_resolution_clock::now();
        dynamic_call(Dynamics::calculateF(duration<double>(end - beg).count(), node_val))
        node->getIfTerminated() = true;
    } else {
        cout << "MIP enumeration failed, seek for branch!\n";
        node->getIfTerminated() = false;
        node->size_enumeration_col_pool = 0;
        node->matrix_in_enumeration.clear();
        node->cost_for_columns_in_enumeration_column_pool.resize(0);
        node->index_columns_in_enumeration_column_pool.resize(0);
    }
}

void
CVRP::addBranchCutToUnsolvedInEnu(BbNode *const node, const std::pair<int, int> &info) {
    if (node->getIfTerminated()) return;
    node->if_just_enter_enu = false;
    Brc bf;
    bf.edge = info;
    ++node->getTreeLevel();

    bf.br_dir = true;
    bf.idx_brc = -1;

    cstr_index.resize(num_row);
    iota(cstr_index.begin(), cstr_index.end(), 0);

    auto node2 = new BbNode(node, idx_node + 2, bf);
    reviseEnumColInfoByBrC(node, node2, bf);
    safe_solver(node->getSolver().getNumCol(&num_col))
    bf.br_dir = false;
    node->getBrCs().emplace_back(bf);

    cstr_index.resize(num_row);
    iota(cstr_index.begin(), cstr_index.end(), 0);
    reviseEnumColInfoByBrC(node, node, bf);
    node->index = ++idx_node;
    ++idx_node;

#if SOLUTION_TYPE == 1
    BestBoundFirstBranching::sub_bbt.push(node2);
    BestBoundFirstBranching::sub_bbt.push(node);
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::addNodeIn(DepthFirstBranching::sub_bbt, node2);
  DepthFirstBranching::addNodeIn(DepthFirstBranching::sub_bbt, node);
#endif
}

void CVRP::reviseEnumColInfoByBrC(BbNode *node, BbNode *out_node, const Brc &bf) {
    cleanColumnsInLPByNewBranchCut(out_node, bf);
    if (!node->size_enumeration_col_pool) return;
    rmColByBranchInEnuMatrix(node, out_node->deleted_columns_in_enumeration_pool, true, {bf});
    regenerateEnumMat(node, out_node);
}

void CVRP::cleanColumnsInLPByNewBranchCut(BbNode *const node, const Brc &bf) {
    vector<int> delete_col;
    int *sum_col = new int[num_col]();
    int col;
    int numnzP;
    bool if_keep;
    int ai = bf.edge.first, aj = bf.edge.second;
    if (ai >= aj) throw runtime_error("Wrong in cleanColumnsInLPByNewBranchCut!");
    vector<int> solver_beg(2);
    vector<int> solver_ind(num_col);
    vector<double> solver_val(num_col);
    if (!ai) {
        numnzP = num_col;
        safe_solver(node->getSolver().getConstraints(&numnzP,
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            aj - 1,
            1))
        if (bf.br_dir) {
            for (size_t i = 0; i < numnzP; ++i) {
                col = solver_ind[i];
                if (node->getCols()[col].col_seq.front() == aj || node->getCols()[col].col_seq.back() == aj)continue;
                delete_col.emplace_back(col);
            }
        } else {
            for (size_t i = 0; i < numnzP; ++i) {
                col = solver_ind[i];
                if (node->getCols()[col].col_seq.front() == aj || node->getCols()[col].col_seq.back() == aj) {
                    delete_col.emplace_back(col);
                }
            }
        }
    } else {
        numnzP = num_col;
        safe_solver(node->getSolver().getConstraints(&numnzP,
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            ai - 1,
            1))
        for (size_t i = 0; i < numnzP; ++i)++sum_col[solver_ind[i]];
        numnzP = num_col;
        safe_solver(node->getSolver().getConstraints(&numnzP,
            solver_beg.data(),
            solver_ind.data(),
            solver_val.data(),
            aj - 1,
            1))
        for (size_t i = 0; i < numnzP; ++i)++sum_col[solver_ind[i]];
        if (bf.br_dir) {
            for (int i = 1; i < num_col; ++i) {
                if (sum_col[i]) {
                    if (sum_col[i] == 1) {
                        delete_col.emplace_back(i);
                    } else {
                        auto end = (int) node->getCols()[i].col_seq.size() - 1;
                        for (int j = 0; j < end; ++j) {
                            int current_node = node->getCols()[i].col_seq[j];
                            if (current_node == ai) {
                                if (node->getCols()[i].col_seq[j + 1] != aj) {
                                    delete_col.emplace_back(i);
                                }
                                break;
                            } else if (current_node == aj) {
                                if (node->getCols()[i].col_seq[j + 1] != ai) {
                                    delete_col.emplace_back(i);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        } else {
            for (int i = 1; i < num_col; ++i) {
                if (sum_col[i] == 2) {
                    int end = (int) node->getCols()[i].col_seq.size() - 1;
                    for (int j = 0; j < end; ++j) {
                        int current_node = node->getCols()[i].col_seq[j];
                        if (current_node == ai) {
                            if (node->getCols()[i].col_seq[j + 1] == aj) {
                                delete_col.emplace_back(i);
                            }
                            break;
                        } else if (current_node == aj) {
                            if (node->getCols()[i].col_seq[j + 1] == ai) {
                                delete_col.emplace_back(i);
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    if (std::find(delete_col.begin(), delete_col.end(), 0) != delete_col.end())
        delete_col.erase(std::find(delete_col.begin(), delete_col.end(), 0));

    rmLPCols(node, delete_col);
    safe_solver(node->getSolver().reoptimize())
    delete[]sum_col;
}

void CVRP::generateRCCsInEnum(BbNode *const node) {
    int numnz, cnt = 0;
    auto &ptr_col = node->index_columns_in_enumeration_column_pool;
    int curr_node;
    set<int> tmp_cutinfo;
    vector<int> solver_ind(num_col); //one time thing
    vector<double> solver_val(num_col);

    CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

    CMGR_CreateCMgr(&MyCutsCMP, 100);
    CMGR_CreateCMgr(&MyOldCutsCMP, 100);

    getEdgeInfo(node, false);
    for (int i = 1; i <= node->getNumEdges(); ++i)
        if (!node->getEdgeTail()[i]) node->getEdgeTail()[i] = dim;
        else break;

    CAPSEP_SeparateCapCuts(real_dim, demand, cap, node->getNumEdges(), node->getEdgeTail().data(),
                           node->getEdgeHead().data(), node->getEdgeValue().data(), MyOldCutsCMP,
                           MAX_NUM_OF_CUTS, TOLERANCE, TOLERANCE,
                           &if_int_n_feasible, &max_vio, MyCutsCMP);
    for (int i = 1; i <= node->getNumEdges(); ++i)
        if (node->getEdgeTail()[i] == dim) node->getEdgeTail()[i] = 0;
        else break;

    int old_num_rcc = (int) node->rccs.size();
    vector<Eigen::RowVectorXi> v_col_map;
    Eigen::RowVectorXi tmp;
    if (!MyCutsCMP->Size) goto QUIT;

    v_col_map.resize(dim, Eigen::RowVectorXi::Zero(num_col));
    for (int col = 0; col < num_col; ++col) {
        for (auto j: node->getCols()[col].col_seq) {
            v_col_map[j][col] = 1;
        }
    }
    tmp.resize(num_col);

    for (int i = 0; i < MyCutsCMP->Size; ++i) {
        Rcc rcc;
        tmp.setZero();
        auto &tmp_customerInfo = rcc.info_rcc_customer;
        for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
            tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
            tmp += v_col_map[MyCutsCMP->CPL[i]->IntList[j]];
        }

        rcc.form_rcc = true;
        rcc.rhs = MyCutsCMP->CPL[i]->RHS;

        double k = (int) tmp_customerInfo.size() - rcc.rhs;

        if (std::find(node->rccs.begin(), node->rccs.end(), rcc) != node->rccs.end()) continue;

        ++cnt;

        numnz = 0;
        solver_ind[numnz] = 0;
        solver_val[numnz++] = k;
        for (int j = 1; j < num_col; ++j) {
            if (tmp[j]) {
                solver_ind[numnz] = j;
                solver_val[numnz++] = 1;
            }
        }
        rcc.idx_rcc = num_row;
        node->rccs.emplace_back(rcc);
        safe_solver(node->getSolver().addConstraint(numnz,
            solver_ind.data(),
            solver_val.data(),
            SOLVER_GREATER_EQUAL,
            k,
            nullptr))
        safe_solver(node->getSolver().updateModel())
        safe_solver(node->getSolver().getNumRow(&num_row))
    }

QUIT:
    for (int i = 0; i < MyCutsCMP->Size; ++i) CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
    MyCutsCMP->Size = 0;

    CMGR_FreeMemCMgr(&MyOldCutsCMP);
    CMGR_FreeMemCMgr(&MyCutsCMP);
}

void CVRP::deleteNonActiveCutsSafely(BbNode *const node, int old_num, bool &if_give_up_this_cutting_round) {
    /**
     * since this mode only delete cuts by slack and only new cuts will be tested,
     * therefore, the rcc can be safely deleted
     */
    if_give_up_this_cutting_round = false;
    safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
    if (node->index && !if_in_enu_state) {
        safe_solver(node->getSolver().getObjVal(&lp_val))
        if ((lp_val - node->getCurrentNodeVal() < 2 * CUTTING_BRANCHING_RATIO * node->br_value_improved)
            &&
            (max_labeling_time_last_round_for_node_so_far > Config::CutGenTimeThresholdInPricingInitial)) {
            if_give_up_this_cutting_round = true;
        }
    }
    int cnt = 0;
    vector<double> slack(num_row);
    safe_solver(node->getSolver().getSlack(0, num_row, slack.data()))
    vector<int> local_cstr_index(num_row);
    iota(local_cstr_index.begin(), local_cstr_index.end(), 0);
    int delta = 0;
    vector<int> deleted_cstrs;
    decltype(deleted_cstrs.begin()) stop_sign;
    if (if_give_up_this_cutting_round) {
        deleted_cstrs.resize(num_row - old_num);
        iota(deleted_cstrs.begin(), deleted_cstrs.end(), old_num);
    } else {
        for (int i = old_num; i < num_row; ++i) {
            if (abs(slack[i]) > TOLERANCE) {
                deleted_cstrs.emplace_back(i);
            }
        }
    }
    if (deleted_cstrs.empty()) {
        goto QUIT;
    }

    std::sort(deleted_cstrs.begin(), deleted_cstrs.end());
    for (auto &i: deleted_cstrs) local_cstr_index[i] = -1;
    stop_sign = deleted_cstrs.end() - 1;
    for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
        ++delta;
        for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
    }
    ++delta;
    for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;

    safe_solver(node->getSolver().delConstraints((int)deleted_cstrs.size(), deleted_cstrs.data()))
    safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
    safe_solver(node->getSolver().getNumRow(&num_row)) //cannot be deleted

    for (auto i = node->rccs.begin(); i < node->rccs.end();) {
        if (local_cstr_index[i->idx_rcc] == -1) {
            i = node->rccs.erase(i);
        } else {
            i->idx_rcc = local_cstr_index[i->idx_rcc];
            ++i;
        }
    }

    for (auto i = node->r1cs.begin(); i < node->r1cs.end();) {
        if (local_cstr_index[i->idx_r1c] == -1) {
            i = node->r1cs.erase(i);
        } else {
            i->idx_r1c = local_cstr_index[i->idx_r1c];
            ++i;
        }
    }

QUIT:
    if (if_in_enu_state) {
        if (node->size_enumeration_col_pool) {
            changeEnumMatByCuts(node);
        }
    }
    verbose_call(cout << deleted_cstrs.size() << " cuts are deleted by slack\n")
}

void CVRP::changeEnumMatByCuts(BbNode *node) {
    if (!node->size_enumeration_col_pool) return;

    auto &mat0 = node->matrix_in_enumeration.front();
    int oldNum = 0;
    for (auto &it: node->matrix_in_enumeration) oldNum += (int) it.rows();

    sparseRowMatrixXd mat(num_row - oldNum, node->size_enumeration_col_pool);
    vector<Eigen::Triplet<double> >
            triplets;
    triplets.reserve(size_t(double(num_row) * double(node->size_enumeration_col_pool) * 0.1));

    buildRCCInEnuMatrix(node, triplets, oldNum);
    buildAllR1CInEnuMatrix(node, triplets, oldNum);

    mat.setFromTriplets(triplets.begin(), triplets.end());
    node->matrix_in_enumeration.push_back(std::move(mat));
}

void CVRP::findNonActiveCuts(BbNode *node, bool if_optimal_dual) {
    vector<int> nonactive_cuts;
    vector<double> current_dual_vec;
    nonactive_cuts.reserve(num_row);
    auto &dual_vec = if_optimal_dual ? optimal_dual_vector : current_dual_vec;
    if (if_optimal_dual && optimal_dual_vector.size() != num_row) {
        throw runtime_error("optimal dual vector size is not equal to num_row");
    }
    if (!if_optimal_dual) {
        safe_solver(node->getSolver().reoptimize(barrier_call(SOLVER_BARRIER)))
        current_dual_vec.resize(num_row);
        safe_solver(node->getSolver().getDual(0, num_row, current_dual_vec.data()))
    }

    for (auto &rcc: node->rccs) {
#if SOLVER_VRPTW == 1
	if (rcc.if_keep) continue;
#endif
        int idx = rcc.idx_rcc;
        if (abs(dual_vec[idx]) < TOLERANCE) {
            nonactive_cuts.emplace_back(idx);
        }
    }
    for (auto &r1c: node->r1cs) {
        int idx = r1c.idx_r1c;
        if (abs(dual_vec[idx]) < TOLERANCE) {
            nonactive_cuts.emplace_back(idx);
        }
    }

    deleteNonactiveCuts(node, nonactive_cuts);
}

void CVRP::deleteNonactiveCuts(BbNode *node, std::vector<int> &nonactive_cuts) {
    if (nonactive_cuts.empty()) return;
    int old_num_row = num_row;
    sort(nonactive_cuts.begin(), nonactive_cuts.end());
    vector<int> local_cstr_index(num_row);
    iota(local_cstr_index.begin(), local_cstr_index.end(), 0);
    int cnt = 0;
    vector<int> solver_ind(num_row);
    for (auto &cut: nonactive_cuts) {
        solver_ind[cnt++] = cut;
        local_cstr_index[cut] = -1;
    }
    int delta = 0;
    auto stop_sign = nonactive_cuts.end() - 1;
    for (auto i = nonactive_cuts.begin(); i < stop_sign; ++i) {
        ++delta;
        for (int j = *i + 1; j < *(i + 1); ++j) local_cstr_index[j] = j - delta;
    }
    ++delta;
    for (int j = *stop_sign + 1; j < num_row; ++j) local_cstr_index[j] = j - delta;
    safe_solver(node->getSolver().delConstraints(cnt, solver_ind.data()))
    safe_solver(node->getSolver().reoptimize())
    safe_solver(node->getSolver().getNumRow(&num_row))

    for (auto i = node->rccs.begin(); i < node->rccs.end();) {
        if (local_cstr_index[i->idx_rcc] == -1) {
            i = node->rccs.erase(i);
        } else {
            i->idx_rcc = local_cstr_index[i->idx_rcc];
            ++i;
        }
    }

    for (auto i = node->r1cs.begin(); i < node->r1cs.end();) {
        if (local_cstr_index[i->idx_r1c] == -1) {
            i = node->r1cs.erase(i);
        } else {
            i->idx_r1c = local_cstr_index[i->idx_r1c];
            ++i;
        }
    }
    if (!optimal_dual_vector.empty()) {
        int ccnt = 0;
        for (int i = 0; i < optimal_dual_vector.size(); ++i) {
            if (local_cstr_index[i] == -1) continue;
            optimal_dual_vector[ccnt++] = optimal_dual_vector[i];
        }
        optimal_dual_vector.resize(ccnt);
    }
    if (if_in_enu_state) {
        cstr_index = std::move(local_cstr_index);
    } else {
        for (auto &i: node->getBrCs()) {
            if (i.idx_brc != -1) {
                i.idx_brc = local_cstr_index[i.idx_brc];
            }
        }
    }
}

double CVRP::getGapStdTryEnumeration() const {
    if (success_enumeration_gap.second == 0) {
        return Config::MaxGap2TryEnumeration;
    }
    double mean = success_enumeration_gap.first / success_enumeration_gap.second;
    mean = max(mean / Config::EnumerationFailFactor, mean + Config::min_enumeration_exploration_added);
    mean = (mean + success_enumeration_gap.first) / (success_enumeration_gap.second + 1);
    mean = max(max_enumeration_success_gap, mean);
    mean = min(min_enumeration_fail_gap - TOLERANCE, mean);
    return mean;
}

void CVRP::adjustEnumerationStdBucketArcs(BbNode *node, bool if_suc) {
    if (if_suc) {
        max_bucket_arc_suc_enumeration.first =
                max(max_bucket_arc_suc_enumeration.first, double(node->num_forward_bucket_arcs));
        symmetry_prohibit_call(
            max_bucket_arc_suc_enumeration.second =
            max(max_bucket_arc_suc_enumeration.second, double(node->num_backward_bucket_arcs));)
    } else {
        min_bucket_arc_fail_enumeration.first =
                min(min_bucket_arc_fail_enumeration.first, double(node->num_forward_bucket_arcs));
        symmetry_prohibit_call(
            min_bucket_arc_fail_enumeration.second =
            min(min_bucket_arc_fail_enumeration.second, double(node->num_backward_bucket_arcs));)
        min_bucket_arc_fail_enumeration.first =
                max(max_bucket_arc_suc_enumeration.first + 1, min_bucket_arc_fail_enumeration.first);
        symmetry_prohibit_call(
            min_bucket_arc_fail_enumeration.second =
            max(max_bucket_arc_suc_enumeration.second + 1,
                min_bucket_arc_fail_enumeration.second);)
    }
}

void CVRP::adjustEnumerationStdGap(double gap, bool if_suc) {
    if (if_suc) {
        max_enumeration_success_gap = max(max_enumeration_success_gap, gap);
        double mean = 0;
        if (success_enumeration_gap.second) {
            mean = success_enumeration_gap.first / success_enumeration_gap.second;
        }
        if (gap > mean) {
            success_enumeration_gap.first = gap;
            success_enumeration_gap.second = 1;
        } else {
            success_enumeration_gap.first = getGapStdTryEnumeration() * (success_enumeration_gap.second + 1);
            ++success_enumeration_gap.second;
        }
        Config::MaxGap2TryEnumeration = getGapStdTryEnumeration();
    } else {
        min_enumeration_fail_gap = min(min_enumeration_fail_gap, gap);
        min_enumeration_fail_gap = max(max_enumeration_success_gap + TOLERANCE, min_enumeration_fail_gap);
        if (success_enumeration_gap.second) {
            double mean = success_enumeration_gap.first / success_enumeration_gap.second;
            mean *= Config::EnumerationFailFactor; //reduce the mean
            success_enumeration_gap.first += mean;
            ++success_enumeration_gap.second; //increase the denominator, keep average the same
        } else {
            Config::MaxGap2TryEnumeration *= Config::EnumerationFailFactor;
        }
    }
    verbose_call(
        if (success_enumeration_gap.second)
        cout << "average gap of enumeration: " << success_enumeration_gap.first / success_enumeration_gap.second
        << endl;
        cout << "MaxGap2TryEnumeration: " << Config::MaxGap2TryEnumeration << " next try node value: "
        << BaseBranching::ub * (1 - getGapStdTryEnumeration()) << endl;)
}
