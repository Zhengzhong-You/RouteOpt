//
// Created by Zhengzhong You on 7/16/22.
//
#include <iostream>

#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace chrono;
using namespace Eigen;

#ifdef find_missed_solution//no zero
vector<vector<int>> sol_ =
    {
        {161, 71, 141, 182},
        {123, 77, 121, 73, 7, 128, 193},
        {167, 79, 100, 145},
        {68, 90, 105, 53, 181},
        {56, 151, 164, 80, 44, 143, 140},
        {34, 129, 78, 84, 126, 67},
        {106, 24, 183, 195, 8},
        {43, 13, 70, 2, 39},
        {49, 98, 5, 115, 52},
        {12, 180, 185, 134, 163, 87, 188},
        {21, 118, 60, 66, 104, 35, 64, 168},
        {107, 137, 32, 194, 36},
        {22, 184, 171, 114, 142, 166},
        {127, 65, 92, 97, 91, 103, 148, 31, 191},
        {42, 197, 122, 99, 28},
        {47, 138, 116, 95, 86, 146},
        {18, 169, 63, 26, 173, 196},
        {11, 179, 3, 62, 9},
        {6, 50, 10, 29, 102, 132},
        {23, 165, 74, 110, 33, 192, 136},
        {45, 108, 156, 153, 120},
        {61, 170, 119, 55},
        {93, 1, 178, 94, 69, 147, 41, 198},
        {112, 4, 113, 111, 48, 51, 17},
        {38, 46, 150, 162, 37, 85},
        {133, 177, 81, 72, 190, 176},
        {59, 57, 131, 109, 96},
        {155, 172, 130, 83, 25},
        {175, 159, 89, 144, 117, 54},
        {20, 139, 199, 14, 160, 186},
        {154, 189, 135, 19, 30, 75, 15},
        {101, 174, 152, 27, 187},
        {125, 124, 76, 157, 40},
        {16, 149, 88, 58, 82, 158}
    };
void CVRP::findWhySolutionDisappear(BBNODE *node,
                                    CVRP *cvrp,
                                    vector<int> &data,
                                    bool &if_opt,
                                    bool if_just_check_opt_node) {
  cout << "findWhySolutionDisappear!" << endl;
//  if (cvrp->if_findBetterUsingHeurEnumeration_suc) return;
  if_opt = true;
  double total_cost = 0;
  for (auto &i : sol_) {
    double cost = 0;
    int past = 0;
    for (auto &j : i) {
      cost += cvrp->CostMat4Vertex[past][j];
      past = j;
    }
    cost += cvrp->CostMat4Vertex[past][0];
    total_cost += cost;
  }
  cout << "total cost: " << total_cost << endl;
  set<pair<int, int>> sol_set;
  for (auto &i : sol_) {
    sol_set.insert({0, i[0]});
    sol_set.insert({i[0], 0});
    for (int j = 0; j < i.size() - 1; ++j) {
      sol_set.insert({i[j], i[j + 1]});
      sol_set.insert({i[j + 1], i[j]});
    }
    sol_set.insert({i[i.size() - 1], 0});
    sol_set.insert({0, i[i.size() - 1]});
  }
  for (auto &brc : node->BrCs) {
    if (brc.BrDir) {
      if (sol_set.find(brc.Edge) == sol_set.end()) {
        if_opt = false;
        return;
      }
    } else {
      if (sol_set.find(brc.Edge) != sol_set.end()) {
        if_opt = false;
        return;
      }
    }
  }
  cout << "this is the optimal node! suppose to have solution " << total_cost << endl;
  cout << "check optimal solution !" << endl;

  if (if_just_check_opt_node) {
    data.emplace_back(-1);
    return;
  }

  bool if_wrong = false;
  double tol_cost = 0;
//  {
//    unordered_map<yzzLong, size_t> idx;
//    int NumCol = cvrp->NumCol;
//    auto colpool = cvrp->ColPool4Pricing;
//    for (int i = 0; i < NumCol; ++i) {
//      yzzLong mem = 0;
//      for (auto j = node->IdxCols[i] + 1;; ++j) {
//        int cur_node = colpool[j];
//        if (!cur_node) break;
//        mem.set(cur_node);
//      }
//      idx[mem] = node->IdxCols[i];
//    }
//    int size_enu = node->SizeEnuColPool;
//    for (int i = 0; i < size_enu; ++i) {
//      yzzLong mem = 0;
//      for (auto j = node->IdxColsInEnuColPool[i] + 1;; ++j) {
//        int cur_node = colpool[j];
//        if (!cur_node) break;
//        mem.set(cur_node);
//      }
//      idx[mem] = node->IdxColsInEnuColPool[i];
//    }
//    for (auto &i : sol_) {
//      yzzLong mem = 0;
//      for (int j : i) {
//        mem.set(j);
//      }
//      if (idx.find(mem) == idx.end()) {
//        cout << "solution disappear!" << endl;
//        for (int j : i) {
//          cout << j << " ";
//        }
//        if_wrong = true;
//        cout << endl;
//      } else {
//        vector<int> seq;
//        double cost = 0;
//        int past_node = 0;
//        for (auto j = idx[mem] + 1;; ++j) {
//          int cur_node = colpool[j];
//          cost += cvrp->CostMat4Vertex[past_node][cur_node];
//          if (!cur_node) break;
//          seq.emplace_back(cur_node);
//          past_node = cur_node;
//        }
//        data.emplace_back(cost);
//        cout << "data is cost meaning!" << endl;
//        tol_cost += cost;
//        vector<int> re_seq(seq.rbegin(), seq.rend());
//        if (seq != i && re_seq != i) {
//          cout << "solution disappear!" << endl;
//          for (int j : i) {
//            cout << j << " ";
//          }
//          if_wrong = true;
//          cout << endl;
//        }
//      }
//    }
//  }
  {
    unordered_map<yzzLong, int> idx;
    int NumCol = cvrp->NumCol;
    auto colpool = cvrp->ColPool4Pricing;
    for (int i = 0; i < NumCol; ++i) {
      yzzLong mem = 0;
      for (auto j = node->IdxCols[i] + 1;; ++j) {
        int cur_node = colpool[j];
        if (!cur_node) break;
        mem.set(cur_node);
      }
      idx[mem] = i;
    }

    int size_enu = node->SizeEnuColPool;
    for (int i = 0; i < size_enu; ++i) {
      yzzLong mem = 0;
      for (auto j = node->IdxColsInEnuColPool[i] + 1;; ++j) {
        int cur_node = colpool[j];
        if (!cur_node) break;
        mem.set(cur_node);
      }
      idx[mem] = i + NumCol;
    }

    if_wrong = false;
    for (auto &i : sol_) {
      yzzLong mem = 0;
      for (int j : i) {
        mem.set(j);
      }
      if (idx.find(mem) == idx.end()) {
        cout << "solution disappear!" << endl;
        for (int j : i) {
          cout << j << " ";
        }
        if_wrong = true;
        cout << endl;
      } else {
        data.emplace_back(idx[mem]);
        cout << "data is index meaning!" << endl;
        cout << idx[mem] << " ";
      }
    }
  }
  if (if_wrong) {
//    exit(0);
    cout << "This is the Mark!" << endl;
    exit(0);
  }
}
#endif

void CVRP::solveMIP(BBNODE *const node, bool if_inEnu) {
  //write MIP
//scheme one
#ifdef DEBUG_MIP_FILE
  //  writeMIP(node->solver, "LP");
  //  writeMIP(node->solver, "MPS");
  //  safe_solver(node->solver.SOLVERreoptimize())
  //  double objval;
  //  safe_solver(node->solver.SOLVERgetObjVal(&objval))
  //  cout << "objval= " << objval << endl;
#endif
  auto beg = high_resolution_clock::now();

  char *xtype = new char[NumCol];
#ifdef SolveMIPNotFixAll
  cout << "SolveMIPNotFixAll" << endl;
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetRC(0, NumCol, RC))
  int ccnt = 0;
  double half_Cap = calculateOptGap(node) / 2;
  for (int i = 0; i < NumCol; ++i) {
    if (RC[i] > half_Cap) {
      xtype[i] = SOLVER_CONTINUOUS;
      solver_ind[ccnt++] = i;
    } else {
      xtype[i] = SOLVER_BINARY;
    }
  }
  fill(solver_val, solver_val + ccnt, 1);
  safe_solver(node->solver.SOLVERaddconstr(ccnt, solver_ind, solver_val, SOLVER_LESS_EQUAL, 1, nullptr))
  ++NumRow;
  cout << "num of non fixed= " << ccnt << endl;
#else
  fill_n(xtype, NumCol, SOLVER_BINARY);
#endif
  safe_solver(node->solver.SOLVERsetenvCutoff(UB + roundUpTolerance))
  safe_solver(node->solver.SOLVERsetVTypearray(0, NumCol, xtype))

//  vector<double> obj(NumCol);
//  safe_solver(node->solver.SOLVERgetObj(0, NumCol, obj.data()))
//  iota(solver_ind, solver_ind + NumCol, 0);
//  GRBaddconstr(node->solver.model, NumCol, solver_ind, obj.data(), SOLVER_GREATER_EQUAL, UB - 1, nullptr);

  //write MIP
//  string str_ub = to_string(UB + roundUpTolerance) + ".mps";
//  safe_solver(node->solver.SOLVERwrite(str_ub.c_str()))

//  safe_solver(node->solver.SOLVERsetenvOutputFlag(1, false))
//  int bad_row = NumRow - Dim;
//  vector<int> cind(bad_row);
//  std::iota(cind.begin(), cind.begin() + bad_row, Dim);
//  safe_solver(node->solver.SOLVERdelconstrs(bad_row, cind.data()))
///delete the first column

  if (if_inEnu) {
    if_MIP_enumeration_suc = true;
    if (NumCol > CONFIG::MinNumRoute4MIP) {
      safe_solver(node->solver.SOLVERsetenvTimeLimit(CONFIG::MIPInEnumerationTimeLimit))
    }
    //reduce the rows in the enu
    //collect the coefficients of the rows
    if (CONFIG::MIPKeepPercentageRowsByDensity < 1 - TOLERANCE) {
      int all_size = int(node->RCCs.size() + node->R1Cs.size() + node->R1Cs_multi.size());
      unordered_map<int, int> row2Density(all_size);
      int MainCstr = Dim;
      int left = NumRow - MainCstr;
      if (left != all_size) {
        cout << "error in the number of rows!" << endl;
        exit(0);
      }
      size_t numnzP;
      safe_solver(node->solver.SOLVERXgetconstrs(&numnzP,
                                                 solver_beg,
                                                 solver_ind,
                                                 solver_val,
                                                 MainCstr,
                                                 left))
      int l_left = left - 1;
      for (int i = 0; i < l_left; ++i) {
        auto density = solver_beg[i + 1] - solver_beg[i];
        row2Density[i + MainCstr] = (int) density;
      }
      auto density = numnzP - solver_beg[l_left];
      row2Density[l_left + MainCstr] = (int) density;

      vector<pair<int, int>> row2DensityVec(row2Density.begin(), row2Density.end());
      int n = int(int(row2DensityVec.size()) * CONFIG::MIPKeepPercentageRowsByDensity);
      nth_element(row2DensityVec.begin(), row2DensityVec.begin() + n, row2DensityVec.end(),
                  [](const pair<int, int> &a, const pair<int, int> &b) {
                    return a.second < b.second;
                  });
//delete the rows
      vector<int> cstr_idx(NumRow);
      vector<int> deleted_cstrs(row2DensityVec.size() - n);
      int cnt = 0;
      std::iota(cstr_idx.begin(), cstr_idx.end(), 0);
      for (int i = n; i < row2DensityVec.size(); ++i) {
        cstr_idx[row2DensityVec[i].first] = -1;
        solver_ind[cnt++] = row2DensityVec[i].first;
        deleted_cstrs[i - n] = row2DensityVec[i].first;
      }
      std::sort(deleted_cstrs.begin(), deleted_cstrs.end());
      int delta = 0;
      auto stop_sign = deleted_cstrs.end() - 1;
      for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
        ++delta;
        for (int j = *i + 1; j < *(i + 1); ++j) cstr_idx[j] = j - delta;
      }
      ++delta;
      for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_idx[j] = j - delta;

      safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
      NumRow -= cnt;
      cout << "NumRow after deleting: " << NumRow << endl;
      //rcc
      for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
        if (cstr_idx[i->IdxRCC] == -1) {
          i = node->RCCs.erase(i);
        } else {
          i->IdxRCC = cstr_idx[i->IdxRCC];
          ++i;
        }
      }
      //r1c
      for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
        if (cstr_idx[i->IdxR1C] == -1) {
          i = node->R1Cs.erase(i);
        } else {
          i->IdxR1C = cstr_idx[i->IdxR1C];
          ++i;
        }
      }
      //r1c_multi
      for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
        if (cstr_idx[i->IdxR1C] == -1) {
          i = node->R1Cs_multi.erase(i);
        } else {
          i->IdxR1C = cstr_idx[i->IdxR1C];
          ++i;
        }
      }
    }
  }

#ifdef find_missed_solution
    vector<double> vec_rhs(NumRow);
    safe_solver(node->solver.SOLVERgetRHS(0, NumRow, vec_rhs.data()))
    vector<int> data;
    bool if_opt;
    findWhySolutionDisappear(node, this, data, if_opt);
    if (node->TreeLevel == 4 && node->BrCs.back().Edge == make_pair(10, 126) && node->BrCs.back().BrDir == false) {
      cout << "please delete this in feature use!" << endl;
      safe_solver(node->solver.SOLVERsetenvOutputFlag(1, false))
      safe_solver(GRBsetdblparam(GRBgetenv(node->solver.model), GRB_DBL_PAR_OPTIMALITYTOL, 1e-5))
      safe_solver(node->solver.SOLVERMIPoptimize())
      safe_solver(node->solver.SOLVERsetenvOutputFlag(0, false))
      safe_solver(node->solver.SOLVERwrite("unused.mps"))
      unordered_set<yzzLong> tmp_set;
      for (auto &s : sol_) {
        yzzLong tmp = 0;
        for (auto &i : s) tmp.set(i);
        tmp_set.emplace(tmp);
      }
      // get all variables of the node->solver
      size_t numnzP;
      safe_solver(node->solver.SOLVERXgetvars(&numnzP,
                                              solver_beg,
                                              solver_ind,
                                              solver_val,
                                              0,
                                              NumCol))
      //map the variable into yzzLong
      solver_beg[NumCol] = numnzP;
      vector<int> deleted_cols;
      for (int i = 0; i < NumCol; ++i) {
        yzzLong tmp = 0;
        for (auto j = solver_beg[i]; j < solver_beg[i + 1]; ++j) {
          auto idx = solver_ind[j];
          if (idx < RealDim) {
            tmp.set(idx + 1);
          }
        }
        if (tmp_set.find(tmp) == tmp_set.end()) {
          deleted_cols.emplace_back(i);
        }
      }
      if (!deleted_cols.empty()) {
        safe_solver(node->solver.SOLVERdelvars(deleted_cols.size(), deleted_cols.data()))
        NumCol -= (int) deleted_cols.size();
        cout << "NumCol after deleting: " << NumCol << endl;
      }
      safe_solver(node->solver.SOLVERwrite("optimal.mps"))
      exit(0);
    }
    if (if_opt) {
      safe_solver(node->solver.SOLVERsetenvOutputFlag(1, false))
      //    safe_solver(GRBsetdblparam(node->solver.env, "MIPGapAbs", 0))
      //    safe_solver(GRBsetdblparam(node->solver.env, "MIPGap", 0))
      vector<int> cind(NumRow - RealDim);
      std::iota(cind.begin(), cind.end(), RealDim);
      safe_solver(node->solver.SOLVERdelconstrs(NumRow - RealDim, cind.data()))
      safe_solver(node->solver.SOLVERupdatemodel())
      vector<double> obj_coeff(NumCol);
      safe_solver(GRBgetdblattrarray(node->solver.model, "Obj", 0, NumCol, obj_coeff.data()))
      //    int cnt = 0;
      for (int i = 1; i < NumCol; ++i) {
        double cost = 0;
        int past_node = 0;
        yzzLong mem = 0;
        for (auto j = node->IdxCols[i] + 1;; ++j) {
          int cur_node = ColPool4Pricing[j];
          cost += CostMat4Vertex[past_node][cur_node];
          if (!cur_node) break;
          mem.set(cur_node);
          past_node = cur_node;
        }
        if (cost != obj_coeff[i]) {
          cout << "--------------------" << endl;
          cout << "cost: " << cost << endl;
          cout << "obj_coeff[i]: " << obj_coeff[i] << endl;
        }
        if (std::find(data.begin(), data.end(), i) != data.end()) {
          int numnzP;
          cout << "------------------------" << endl;
          cout << "Real Seq: ";
          for (auto j = node->IdxCols[i] + 1;; ++j) {
            int cur_node = ColPool4Pricing[j];
            if (!cur_node) break;
            cout << cur_node << " ";
          }
          cout << endl;
          cout << "Fake Seq: ";
          GRBgetvars(node->solver.model, &numnzP, solver_ind2, solver_ind, solver_val, i, 1);
          for (int j = 0; j < numnzP; ++j) {
            if (!mem[solver_ind[j] + 1]) {
              cout << " " << solver_ind[j] + 1;
            }
          }
          cout << endl;
          //        solver_ind[cnt++] = i;
        }
      }
      //    cout << "we delete " << cnt << " columns!" << endl;
      //    safe_solver(node->solver.SOLVERdelvars(cnt, solver_ind))
      safe_solver(node->solver.SOLVERsetenvOutputFlag(0, false))
      safe_solver(node->solver.SOLVERwrite("test.mps"))
      cout << "here we check cuts!" << endl;
      cout << "we check rank-1 cuts!" << endl;
      for (auto &r1c : node->R1Cs) {
        if (int(r1c.InfoR1C.size() / 2) != r1c.RHS || vec_rhs[r1c.IdxR1C] != r1c.RHS) {
          cout << "---------------------" << endl;
          cout << "RHS: " << r1c.RHS << endl;
          cout << "InfoR1C: ";
          for (auto &i : r1c.InfoR1C) {
            cout << i << " ";
          }
          cout << endl;
          cout << "vec_rhs: " << vec_rhs[r1c.IdxR1C] << endl;
        }
      }
      cout << "we check rank-mul_cuts!" << endl;
      for (auto &r1c : node->R1Cs_multi) {
        const auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
        auto &multi = get<0>(plan);
        int denominator = get<1>(plan);
        double rhs = int(accumulate(multi.begin(), multi.end(), 0) / (denominator + 0.));
        if (rhs != r1c.RHS || vec_rhs[r1c.IdxR1C] != r1c.RHS) {
          cout << "---------------------" << endl;
          cout << "RHS: " << r1c.RHS << endl;
          cout << "InfoR1C: ";
          for (auto &i : r1c.InfoR1C.first) {
            cout << i << " ";
          }
          cout << endl;
          cout << "multi: ";
          for (auto &i : multi) {
            cout << i << " ";
          }
          cout << endl;
          cout << "vec_rhs: " << vec_rhs[r1c.IdxR1C] << endl;
        }
      }
      cout << "we check vehcle cuts!" << endl;
      if (vec_rhs[RealDim] != sol_.size()) {
        cout << "the optimal vehicle is " << sol_.size() << endl;
        cout << "the vehicle rhs= " << vec_rhs[RealDim] << endl;
      }
    }
#endif

  safe_solver(node->solver.SOLVERMIPoptimize())

//  if (if_opt) {
//    exit(0);
//  }
  ///
//  safe_solver(node->solver.SOLVERsetenvOutputFlag(0, false))
  ///
#ifdef DEBUG_MIP_FILE
  writeMIP(node->solver, "IP");
//  exit(0);
#endif

  int status;
  safe_solver(node->solver.SOLVERgetStatus(&status))
  if (status == SOLVER_INFEASIBLE || status == SOLVER_INF_OR_UNBD) {
    cout << "Model is infeasible after pre_solving!" << endl;
    cout << "status: " << status << endl;
    LPVal = 1e100;
  } else {
    safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  }
  auto end = high_resolution_clock::now();
  auto eps = duration_cast<milliseconds>(end - beg);
  cout << "MIP used " << (double) eps.count() * 1e-3 << " s!" << endl;

  safe_solver(node->solver.SOLVERsetenvCutoff(MAXCUTOFF))

  if (if_inEnu) {
    if (status == SOLVER_TIME_LIMIT) {
      if_MIP_enumeration_suc = false;
      cout << "MIP time limit reached!" << endl;
      MaxNumRoute4MIP = max(int(CONFIG::MIPInEnumerationTimeLimitChgFactor * NumCol), CONFIG::MinNumRoute4MIP);
      cout << "MaxNumRoute4MIP = " << MaxNumRoute4MIP << endl;
    }
    safe_solver(node->solver.SOLVERsetenvTimeLimit(MAXTIMELIMT4MIP))
    if (UB > ceil_transformed_number_related(LPVal - TOLERANCE) + TOLERANCE) {
      safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
      UB = ceil_transformed_number_related(LPVal - TOLERANCE);
      cout << "solve MIP get " << UB << endl;
      IPOptSol.clear();
      for (int i = 0; i < NumCol; ++i) {
        if (abs(X[i] - 1) < MIP_TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(MaxLengthEleRoute);
          tmp.emplace_back(0);
          for (auto j = node->IdxCols[i] + 1;; ++j) {
            if (!ColPool4Pricing[j]) break;
            tmp.emplace_back(ColPool4Pricing[j]);
          }
          tmp.emplace_back(0);
          IPOptSol.emplace_back(std::move(tmp));
        }
      }
    }

#ifdef MASTER_VALVE_ML
    updateUB_EdgeSol();
#endif

    if (status == SOLVER_TIME_LIMIT) {
      fill_n(xtype, NumCol, SOLVER_CONTINUOUS);
      safe_solver(node->solver.SOLVERsetVTypearray(0, NumCol, xtype))
      safe_solver(node->solver.SOLVERreoptimize())
    }
  } else {
    if (UB > ceil_transformed_number_related(LPVal - TOLERANCE) + TOLERANCE) {
      safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
      UB = ceil_transformed_number_related(LPVal - TOLERANCE);
      cout << "solve MIP get " << UB << endl;
      IPOptSol.clear();
      for (int i = 0; i < node->NumParentCols; ++i) {
        if (abs(X[i] - 1) < MIP_TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(MaxLengthEleRoute);
          tmp.emplace_back(0);
          for (auto j = node->IdxCols[i] + 1;; ++j) {
            if (!ColPool4Mem[j]) break;
            tmp.emplace_back(ColPool4Mem[j]);
          }
          tmp.emplace_back(0);
          IPOptSol.emplace_back(std::move(tmp));
        }
      }
      for (int i = node->NumParentCols; i < NumCol; ++i) {
        if (abs(X[i] - 1) < MIP_TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(MaxLengthEleRoute);
          tmp.emplace_back(0);
          for (auto j = node->IdxCols[i] + 1;; ++j) {
            if (!ColPool4Pricing[j]) break;
            tmp.emplace_back(ColPool4Pricing[j]);
          }
          tmp.emplace_back(0);
          IPOptSol.emplace_back(std::move(tmp));
        }
      }
    }

#ifdef MASTER_VALVE_ML
    updateUB_EdgeSol();
#endif

    fill_n(xtype, NumCol, SOLVER_CONTINUOUS);
    safe_solver(node->solver.SOLVERsetVTypearray(0, NumCol, xtype))
    safe_solver(node->solver.SOLVERreoptimize())
  }
  delete[]xtype;
}

void CVRP::enumerateMIP(BBNODE *&node) {
  tellIfEnumeration(node);
  if (!final_decision_4_enumeration) return;
  final_decision_4_enumeration = false;

  cout << BIG_PHASE_SEPARATION;
  cout << "try EnumerateRoutes..." << endl;

  if (abs(MeetPointResourceInBiDirEnu) < TOLERANCE) MeetPointResourceInBiDirEnu = MeetPointResourceInBiDir;
  //to use the info in arc_elimination, we cannot update the rc
  bool if_succeed;
  double time_labeling;

  PtrAllR1Cs ptrAllR1Cs(node, this);

  auto beg = high_resolution_clock::now();

  Rollback = 0;
  if_succeed = enumerateRoutes(node, ptrAllR1Cs);

  auto end = high_resolution_clock::now();
  auto eps = duration_cast<milliseconds>(end - beg);
  time_labeling = (double) eps.count() * 1e-3;

  if (if_succeed) {
    cout << "enumeration time= " << time_labeling << " and succeed!" << endl;
  } else {
    cout << "RollBack= " << Rollback << endl;
    cout << "enumeration time= " << time_labeling << " but failed!" << endl;
  }

  double gap = (UB - node->Val) / UB;
  if (if_succeed) {
    if (enumerationMode) {
      Count4Tolerance4tryEnumerationWhenArcEliminationFails = 0;
    }
    if (gap > GapTolerance4ArcEliminationNEnumeration) {
      GapTolerance4ArcEliminationNEnumeration = gap;
    }
    LastEnumerationFailGap /= CONFIG::EnumerationFailFactor;
#ifdef SYMMETRY_PROHIBIT
    double dif = abs(NumForwardLabelsInEnu - NumBackwardLabelsInEnu);
    double over = dif / min(NumForwardLabelsInEnu, NumBackwardLabelsInEnu);
#ifdef DETAILED_EXACT_PRINT_INFO
    cout << "over= " << over << endl;
#endif
    if (over > CONFIG::NumberOfOverLabelsInMeetPoint) {
#ifdef DETAILED_EXACT_PRINT_INFO
      cout << "we adjust the meetpoint!" << endl;
#endif
      if (NumForwardLabelsInEnu > NumBackwardLabelsInEnu) {
        MeetPointResourceInBiDirEnu *= (1 - CONFIG::MeetPointFactor);
      } else {
        MeetPointResourceInBiDirEnu *= (1 + CONFIG::MeetPointFactor);
      }
      cout << "MeetPointResourceInBiDirEnu= " << MeetPointResourceInBiDirEnu << endl;
#ifdef DETAILED_EXACT_PRINT_INFO
      cout << "MeetPointResourceInBiDirEnu= " << MeetPointResourceInBiDirEnu << endl;
#endif
    }
#endif
    if_enumeration_suc = true;
#ifdef BranchFashion_MemSaving
    terminateNodeMemSave(node);
#else
    terminateNode(node);
#endif
  } else {
    LastEnumerationFailGap = gap * CONFIG::EnumerationFailFactor;//set a smaller gap!
    if_enumeration_suc = false;
  }
  cout << "LastEnumerationFailGap= " << LastEnumerationFailGap << endl;
  cout << BIG_PHASE_SEPARATION;
}

bool CVRP::enumerateRoutes(BBNODE *const node, const PtrAllR1Cs &ptrAllR1Cs) {
  //in the enumeration phase, because the dominance rule is cost-related.
  //so we sort the labels in the bucket arcs by cost
  //but in concatenate phase, we sort the labels by rc
  int Max_labels = CONFIG::MaxNumLabelInEnumeration;
  int Max_routes_all = CONFIG::MaxNumRouteInEnumeration;
  int Max_routes_phase1 = (int) sqrt_self((float) Max_routes_all);
  int num_routes_now = 0;
  auto &cost_m = node->Cost4ColsInEnuColPool;
  auto &ptr = node->IdxColsInEnuColPool;
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;

  int index;
  OptGap = calculateOptGap(node);
  unordered_map<yzzLong, tuple<LABEL *, LABEL *, double>> Tags;
  Tags.reserve(Max_routes_all);

  //we still need the labels from arc elimination to accelerate
  //we copy all the labels into another graph
  auto copy_Forward_bucket = new vector<LABEL *> *[Dim];
#ifdef SYMMETRY_PROHIBIT
  auto copy_Backward_bucket = new vector<LABEL *> *[Dim];//these memory will be freed inside the forward function
#endif
  if (enumerationMode) {//do without arc elimination
    cout << "enumeration Mode : 1 !" << endl;
    auto tmp_label = AllLabel + IdxGlo - 1;
    tmp_label->RC = -numeric_limits<float>::max();
    tmp_label->Rank1CutMem = 0;
    tmp_label->numValidRank1Cut = 0;
    memset(tmp_label->Rank1CutMem_multi, 0, sizeof(int) * NumValidR1C_multi_InCG);
    tmp_label->numValidRank1Cut_multi = 0;
    for (int i = 0; i < Dim; ++i) {
      copy_Forward_bucket[i] = new vector<LABEL *>[NumBucketsPerVertex];
#ifdef SYMMETRY_PROHIBIT
      copy_Backward_bucket[i] = new vector<LABEL *>[NumBucketsPerVertex];
#endif
      for (int b = 0; b < NumBucketsPerVertex; ++b) {
        copy_Forward_bucket[i][b] = {tmp_label};
        RC2TillThisBinInForwardSense[i][b] = -numeric_limits<float>::max();
        RC2BinInForwardSense[i][b] = -numeric_limits<float>::max();
#ifdef SYMMETRY_PROHIBIT
        copy_Backward_bucket[i][b] = {tmp_label};
        RC2TillThisBinInBackwardSense[i][b] = -numeric_limits<float>::max();
        RC2BinInBackwardSense[i][b] = -numeric_limits<float>::max();
#endif
      }
    }
  } else {
    for (int i = 0; i < Dim; ++i) {
      copy_Forward_bucket[i] = new vector<LABEL *>[NumBucketsPerVertex];
#ifdef SYMMETRY_PROHIBIT
      copy_Backward_bucket[i] = new vector<LABEL *>[NumBucketsPerVertex];
#endif
      for (int b = 0; b < NumBucketsPerVertex; ++b) {
        copy_Forward_bucket[i][b].assign(LabelArrayInForwardSense[i][b].first.begin(),
                                         LabelArrayInForwardSense[i][b].first.begin()
                                             + LabelArrayInForwardSense[i][b].second);
#ifdef SYMMETRY_PROHIBIT
        copy_Backward_bucket[i][b].assign(LabelArrayInBackwardSense[i][b].first.begin(),
                                          LabelArrayInBackwardSense[i][b].first.begin()
                                              + LabelArrayInBackwardSense[i][b].second);
#endif
      }
    }
  }

  PriorPoolBeg4Pricing = PoolBeg4Pricing;
  if (checkPricingPool()) reallocatePricingPool();

  //half_extension
  auto beg = high_resolution_clock::now();
  auto end = high_resolution_clock::now();
  auto eps = duration_cast<milliseconds>(end - beg);

  int status =
#ifdef SYMMETRY_PROHIBIT
      enumerateHalfwardRoutes<true, false>(node, r1c_to_pi, r1c_multi_to_pi, Tags,
                                           copy_Backward_bucket, num_routes_now);
#else
      enumerateHalfwardRoutes<true, true>(node, r1c_to_pi, r1c_multi_to_pi, Tags, copy_Forward_bucket, num_routes_now);
#endif

  end = high_resolution_clock::now();
  eps = duration_cast<milliseconds>(end - beg);
  cout << "Half Forward time= " << (double) eps.count() * 1e-3 << endl;

  if (status || Rollback) {
    if (status == 1)
      cout << "the number of labels in Forward reached its limit!" << endl;
    return false;
  }

#ifdef SYMMETRY_PROHIBIT
  beg = high_resolution_clock::now();

  status = enumerateHalfwardRoutes<false, false>(node,
                                                 r1c_to_pi,
                                                 r1c_multi_to_pi,
                                                 Tags,
                                                 copy_Forward_bucket,
                                                 num_routes_now);

  end = high_resolution_clock::now();
  eps = duration_cast<milliseconds>(end - beg);
  cout << "Half Backward time= " << (double) eps.count() * 1e-3 << endl;

  if (status || Rollback) {
    if (status == 1)
      cout << "the number of labels in Backward reached its limit!" << endl;
    return false;
  }
#endif


  //begin to concat cols
  beg = high_resolution_clock::now();

  status = concatenateRoutes_prior_forward_InEnumeration(node, r1c_to_pi, r1c_multi_to_pi, Tags, num_routes_now);

  end = high_resolution_clock::now();
  eps = duration_cast<milliseconds>(end - beg);
  cout << "Concatenate time= " << (double) eps.count() * 1e-3 << endl;

  if (status) {
    cout << "the number of routes reached its limit!" << endl;
    return false;
  }

  //recover the rc of edges
  priceLabeling(node);
  //organize the cols
  organizeColsInMem2Pricing(node);//reset colPool for pricing! be careful!
  //assign mem for pool
  //no need for setting 0 in this case
  cost_m.resize(num_routes_now);
  ptr.resize(num_routes_now);
  node->SizeEnuColPool = num_routes_now;
  node->validSize = num_routes_now;

  auto PricingWarning = (size_t) (0.9 * (double) Mem4Pricing);
  //populate IdxColsInEnuColPool
  index = 0;
  LABEL *k, *l;
  for (auto &tag : Tags) {
    k = get<0>(tag.second);
    l = get<1>(tag.second);

    cost_m[index] = get<2>(tag.second);
    ptr[index++] = PoolBeg4Pricing;

    for (int m = 0; m <= k->IdxEndSeq; ++m) {
      ColPool4Pricing[PoolBeg4Pricing++] = *(k->Seq + m);
    }

    if (l) {
      for (int m = l->IdxEndSeq; m >= 0; --m) {
        ColPool4Pricing[PoolBeg4Pricing++] = *(l->Seq + m);
      }
    } else {
      ColPool4Pricing[PoolBeg4Pricing++] = 0;
    }
    if (PoolBeg4Pricing >= PricingWarning) {
      cout << SMALL_PHASE_SEPARATION;
      cout << "Warning: the pricing pool is almost full!" << endl;
      cout << "PoolBeg4Pricing=" << PoolBeg4Pricing << endl;
      cout << "Mem4Pricing=" << Mem4Pricing << endl;
      cout << "we reallocate the pricing pool!" << endl;
      reallocatePricingPool();
      PricingWarning = (size_t) (0.9 * (double) Mem4Pricing);
      cout << "the new Mem4Pricing=" << Mem4Pricing << endl;
      cout << SMALL_PHASE_SEPARATION;
    }
  }

  if (index != node->SizeEnuColPool) {
    cerr << "Wrong in enumerateRoutes! The index is not equal to the size_pool!" << endl;
    exit(0);
  }
  //clean the node's lp except the first one
  //clean node's lp
  cleanColsNonEle(node);
  deleteBrCsNR1C1s(node);
  recoverR1CsInEnu(node);

  //generateMatrix and delete bad br-cols
  generateVertex2IdxCols_N_Edge2IdxCols(node);
  return true;
}

//int CVRP::enumerateHalfwardRoutes(BBNODE *node,
//                                  const double *r1c_to_pi,
//                                  const double *r1c_multi_to_pi,
//                                  unordered_map<yzzLong, tuple<LABEL *, LABEL *, double>> &Tags,
//                                  vector<LABEL *> **copy_bucket,
//                                  int &num_routes_now) {
//  //no sort in this case
//  int status = 0;
//  int edgemap;
//  NumForwardLabelsInEnu = 0;
//  int Max_routes_phase1 = CONFIG::MaxNumRouteInEnumeration_half;
//  bool if_keep, if_break;
//  double path_rc, path_cost;
//  initializeLabels(node, 1, false, {true, 1, true});
//
//  auto beg = high_resolution_clock::now();
//  auto end = beg;
//  auto b4_end = beg;
//  auto af_end = beg;
//  double eps;
//  double eps2;
//  double left_time = CONFIG::HardTimeThresholdInAllEnumeration;
//
//  for (int b = 0; b < NumBucketsPerVertex; ++b) {
//    int i = 1;
//    STILL_EXIST:
//    for (; i < Dim; ++i) {
//      end = high_resolution_clock::now();
//      eps = duration<double>(end - b4_end).count();
//      if (eps > left_time) {
//        status = 2;
//        goto outside;
//      }
//      auto &valid_num = IfExistExtraLabelsInForwardSense[i][b].second;
//      if (!valid_num) continue;
//      auto &label_array = IfExistExtraLabelsInForwardSense[i][b].first;
//      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
//        auto &ki = label_array[vec_index];
//        if (ki->if_extended) continue;
//        ki->if_extended = true;
//        for (int j : node->AllForwardBuckets[i][b].BucketArcs) {
//          if (ki->PI[j]) continue;
//          auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
//          if (!increaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
//          if (tmp_mainResource > MeetPointResourceInBiDirEnu) {
//            concatenateLabelsInForwardCG[{i, j}].emplace_back(ki, tmp_mainResource);
//            continue;
//          }
//          auto &tmp_rc = AllLabel[IdxGlo].RC;
//          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
//          if_keep = false;
//#ifdef SYMMETRY_PROHIBIT
//          int arr_bj = int(tmp_mainResource / StepSize);
//          //first test
//          if (tmp_rc + RC2TillThisBinInBackwardSense[j][arr_bj] > OptGap) continue;
//          //second test
//          if (tmp_rc + RC2BinInBackwardSense[j][arr_bj] < OptGap) {
//            for (auto &kkj : copy_bucket[j][arr_bj]) {
//              if (tmp_mainResource > kkj->Sum_MainResource) continue;
//#else
//          int arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//          //first test
//          if (tmp_rc + RC2TillThisBinInForwardSense[j][arr_bj] > OptGap) continue;
//          //second test
//          if (tmp_rc + RC2BinInForwardSense[j][arr_bj] < OptGap) {
//            for (auto &kkj : copy_bucket[j][arr_bj]) {
//              if (tmp_mainResource + kkj->Sum_MainResource > MaxMainResource) continue;
//#endif
//              if ((ki->PI & kkj->PI).any()) continue;
//              path_rc = tmp_rc + kkj->RC;
//              if (path_rc > OptGap) break;
//              if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
//                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//                  if (kkj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//                }
//              } else {
//                for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
//                  if (ki->Rank1CutMem[kkj->validRank1Cut[l]])path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
//                }
//              }
//
//              if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
//                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = ki->validRank1Cut_multi[l];
//                  if (kkj->Rank1CutMem_multi[tmp_cut] +
//                      ki->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              } else {
//                for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = kkj->validRank1Cut_multi[l];
//                  if (ki->Rank1CutMem_multi[tmp_cut] +
//                      kkj->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              }
//
//              if (path_rc < OptGap) {
//                if_keep = true;
//                goto outside1;
//              }
//            }
//          }
//#ifdef SYMMETRY_PROHIBIT
//          //real test
//          for (++arr_bj; arr_bj < NumBucketsPerVertex; ++arr_bj) {
//            //first test
//            if (tmp_rc + RC2TillThisBinInBackwardSense[j][arr_bj] > OptGap) break;
//            //second test
//            if (tmp_rc + RC2BinInBackwardSense[j][arr_bj] > OptGap) continue;
//#else
//          //real test
//          for (--arr_bj; arr_bj >= 0; --arr_bj) {
//            //first test
//            if (tmp_rc + RC2TillThisBinInForwardSense[j][arr_bj] > OptGap) break;
//            //second test
//            if (tmp_rc + RC2BinInForwardSense[j][arr_bj] > OptGap) continue;
//#endif
//            //real test
//            for (auto &kkj : copy_bucket[j][arr_bj]) {
//              if ((ki->PI & kkj->PI).any()) continue;
//              path_rc = tmp_rc + kkj->RC;
//              if (path_rc > OptGap) break;
//              if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
//                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//                  if (kkj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//                }
//              } else {
//                for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
//                  if (ki->Rank1CutMem[kkj->validRank1Cut[l]])path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
//                }
//              }
//
//              if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
//                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = ki->validRank1Cut_multi[l];
//                  if (kkj->Rank1CutMem_multi[tmp_cut] +
//                      ki->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              } else {
//                for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
//                  int tmp_cut = kkj->validRank1Cut_multi[l];
//                  if (ki->Rank1CutMem_multi[tmp_cut] +
//                      kkj->Rank1CutMem_multi[tmp_cut]
//                      >= R1C_multi_denominator_InCG[tmp_cut]
//                      )
//                    path_rc -= r1c_multi_to_pi[tmp_cut];
//                }
//              }
//
//              if (path_rc < OptGap) {
//                if_keep = true;
//                goto outside1;
//              }
//            }
//          }
//          outside1:
//          if (!if_keep) continue;
//          int bj = int(tmp_mainResource / StepSize);
//          auto &labelList_j = LabelArrayInForwardSense[j][bj].first;
//          auto &valid_num_j = LabelArrayInForwardSense[j][bj].second;
//          auto &tmp_PI = AllLabel[IdxGlo].PI;
//          auto &tmp_Cost = AllLabel[IdxGlo].Cost;
//          auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
//          auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
//          auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
//          tmp_PI = ki->PI;
//          tmp_PI.set(j);
//          tmp_Cost = ki->Cost + CostMat4Vertex[i][j];
//          tmp_Rank1CutMem = ki->Rank1CutMem;
//          //tmp_num_valid_rank1_cut do not have to copy
//          for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
//            if (tmp_Rank1CutMem[l]) {
//              tmp_Rank1CutMem[l] = false;
//              tmp_rc -= r1c_to_pi[l];
//            } else tmp_Rank1CutMem[l] = true;
//          }
//          tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);
//
//          auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
//          auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
//          auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;
//          copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);
//          for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
//            int tmp_cut = get<0>(l);
//            tmp_Rank1CutMem_multi[tmp_cut] += get<1>(l);
//            if (tmp_Rank1CutMem_multi[tmp_cut] >= get<2>(l)) {
//              tmp_rc -= r1c_multi_to_pi[tmp_cut];
//              tmp_Rank1CutMem_multi[tmp_cut] -= get<2>(l);
//            }
//          }
//          for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;
//
//          if_break = false;
//
//          for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
//            auto &kj = labelList_j[vec_index_j];
//#ifdef CAPACITY_AS_MAIN_RESOURCE
//            if (abs(kj->Sum_MainResource - tmp_mainResource) > TOLERANCE) {
//              ++vec_index_j;
//              continue;
//            }
//            if ((kj->PI ^ tmp_PI).none()) {
//              if (kj->Cost > tmp_Cost) {
//                kj->if_extended = true;
//                kj = labelList_j[--valid_num_j];
//                --NumForwardLabelsInEnu;
//              } else {
//                if_break = true;
//                break;
//              }
//            } else ++vec_index_j;
//#else
//            if (kj->Cost > tmp_Cost) {
//              if (kj->Sum_MainResource > tmp_mainResource) {
//                if ((kj->PI ^ tmp_PI).none()) {
//                  kj->if_extended = true;
//                  kj = labelList_j[--valid_num_j];
//                  --NumForwardLabelsInEnu;
//                } else ++vec_index_j;
//              } else ++vec_index_j;
//            } else {
//              if (kj->Sum_MainResource < tmp_mainResource) {
//                if ((kj->PI ^ tmp_PI).none()) {
//                  if_break = true;
//                  break;
//                } else ++vec_index_j;
//              } else ++vec_index_j;
//            }
//#endif
//          }
//
//          if (if_break) continue;
//
//          tmp_num_valid_rank1_cut = 0;
//          for (auto l : get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
//            if (tmp_Rank1CutMem[l]) {
//              tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
//            }
//          }
//
//          tmp_num_valid_rank1_cut_multi = 0;
//          for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
//            if (tmp_Rank1CutMem_multi[l]) {
//              tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
//            }
//          }
//
//          labelList_j[valid_num_j++] = AllLabel + IdxGlo;
//          if (valid_num_j == labelList_j.size()) {
//            labelList_j.resize(labelList_j.size() * 2);
//          }
//
//          //Seq
//          AllLabel[IdxGlo].Seq = AllSeq + SeqBeg;
//          copy(ki->Seq, ki->Seq + ki->IdxEndSeq + 1, AllLabel[IdxGlo].Seq);
//          AllLabel[IdxGlo].IdxEndSeq = ki->IdxEndSeq + 1;
//          *(AllLabel[IdxGlo].Seq + AllLabel[IdxGlo].IdxEndSeq) = j;
//          SeqBeg += AllLabel[IdxGlo].IdxEndSeq + 1;
//          //EndVertex
//          AllLabel[IdxGlo].EndVertex = j;
//          //if_extended
//          AllLabel[IdxGlo].if_extended = false;
//          //bucket
//          auto &bucket = IfExistExtraLabelsInForwardSense[j][bj];
//          bucket.first[bucket.second++] = AllLabel + IdxGlo;
//          if (bucket.second == bucket.first.size()) {
//            bucket.first.resize(bucket.first.size() * 2);
//          }
//
//          if (tmp_rc + ChgCostMat4Vertex[j][0] < OptGap) {
//            path_cost = tmp_Cost + CostMat4Vertex[j][0];
//
//            if (Tags.find(tmp_PI) == Tags.end()) {
//              Tags[tmp_PI] = {AllLabel + IdxGlo, nullptr, path_cost};
//              ++num_routes_now;
//              if (num_routes_now > Max_routes_phase1) {
//                status = 2;//routes limit
//                goto outside;
//              }
//            } else if (get<2>(Tags[tmp_PI]) > path_cost) {
//              Tags[tmp_PI] = {AllLabel + IdxGlo, nullptr, path_cost};
//            }
//          }
//          if (++NumForwardLabelsInEnu > CONFIG::MaxNumLabelInEnumeration) {
//            status = 3;//all labels limit
//            goto outside;
//          }
//          ++IdxGlo;//can be put here, because once go to outside, the function will end
//          if (IdxGlo == LabelAssign) {
//            Rollback = 2;
//            goto outside;
//          }
//        }
//      }
//      valid_num = 0;
//    }
////test if all labels are extended
//    for (i = 1; i < Dim; ++i) {
//      if (IfExistExtraLabelsInForwardSense[i][b].second)
//        goto STILL_EXIST;
//    }
//    af_end = high_resolution_clock::now();
//    eps2 = duration<double>(af_end - b4_end).count();
//    eps = duration<double>(af_end - beg).count();
//    left_time = (CONFIG::HardTimeThresholdInAllEnumeration - eps) / (NumBucketsPerVertex - b);
//    if (eps2 > left_time) {
//      status = 2;
//      goto outside;
//    }
//    b4_end = af_end;
//  }
//  outside:
//  for (int i = 0; i < Dim; ++i) {
//    delete[]copy_bucket[i];
//  }
//  delete[] copy_bucket;
//  cout << "Half Forward labeling: num_labels= " << NumForwardLabelsInEnu << " num_routes= " << num_routes_now <<
//       endl;
//  if (status)return status;
////there is one situation that rc of one col could be smaller than optGap
////but could be dominated by cols whose rc is larger than optGap
////this could lead to a discrepancy between the number of routes
//
////for reasons above, we revise the number of routes and sort the routes
//  for (int i = 1; i < Dim; ++i) {
//    for (int b = 0; b < NumBucketsPerVertex; ++b) {
//      std::stable_sort(LabelArrayInForwardSense[i][b].first.begin(),
//                       LabelArrayInForwardSense[i][b].first.begin()
//                           + LabelArrayInForwardSense[i][b].second,
//                       CmpLabelRCLess);
//    }
//  }
////  populateRC2TillThisBinNRC2Bin(node, 1);//never use this function in phase1
//  return 0;
//}

//int CVRP::concatenateRoutes_prior_forward_InEnumeration(BBNODE *node,
//                                                        const double *r1c_to_pi,
//                                                        const double *r1c_multi_to_pi,
//                                                        unordered_map<yzzLong, tuple<LABEL *, LABEL *, double>> &Tags,
//                                                        int &num_routes_now) {
//  int status = 0;
//  double path_rc, path_cost;
//#ifdef SYMMETRY_PROHIBIT
//  populateRC2TillThisBinNRC2Bin<false>(node);
//#else
//  populateRC2TillThisBinNRC2Bin<true>(node);//use function here!
//#endif
//  for (auto &label_list : concatenateLabelsInForwardCG) {
//    int i = label_list.first.first;
//    int j = label_list.first.second;
//    auto &label_vec = label_list.second;
//    for (auto &pr : label_vec) {
//      auto &ki = pr.first;
//      double tmp_mainResource = pr.second;
//      double tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
//#ifdef SYMMETRY_PROHIBIT
//      int arr_bj = int((tmp_mainResource) / StepSize);
//      //most_negative_rc_till_this_bin
//      if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)
//        continue;
//
//      //most_negative_rc_in_this_bin
//      if (RC2BinInBackwardSense[j][arr_bj] + tmp_rc < OptGap) {
//        //add one more condition for testing capacity
//        auto &label_arr = LabelArrayInBackwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInBackwardSense[j][arr_bj].second;
//#else
//      int arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
//      //most_negative_rc_till_this_bin
//      if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)
//        continue;
//
//      //most_negative_rc_in_this_bin
//      if (RC2BinInForwardSense[j][arr_bj] + tmp_rc < OptGap) {
//        //add one more condition for testing capacity
//        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
//#endif
//        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//          auto &kj = label_arr[vec_index];
//          path_rc = kj->RC + tmp_rc;
//          if (path_rc > OptGap) break;
//#ifdef SYMMETRY_PROHIBIT
//          if (tmp_mainResource > kj->Sum_MainResource) continue;
//#else
//          if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
//#endif
//          if ((ki->PI & kj->PI).any()) continue;
//
//          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//              if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//            }
//          }
//
//          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = ki->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut] +
//                  ki->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (ki->Rank1CutMem_multi[tmp_cut] +
//                  kj->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          }
//
//          if (path_rc < OptGap) {
//            //begin to test
//            path_cost = ki->Cost + CostMat4Vertex[i][j] + kj->Cost;
//            auto tmp_PI = ki->PI | kj->PI;
//            if (Tags.find(tmp_PI) == Tags.end()) {
//              Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
//              ++num_routes_now;
//              if (num_routes_now > CONFIG::MaxNumRouteInEnumeration) {
//                status = 1;
//                goto QUIT;
//              }
//            } else if (get<2>(Tags[tmp_PI]) > path_cost) {
//              Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
//            }
//          }
//        }
//      }
//#ifdef SYMMETRY_PROHIBIT
//      //bj-1
//      for (++arr_bj; arr_bj < NumBucketsPerVertex; ++arr_bj) {
//        if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)
//          break;
//
//        if (RC2BinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)
//          continue;
//
//        auto &label_arr = LabelArrayInBackwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInBackwardSense[j][arr_bj].second;
//#else
//      for (--arr_bj; arr_bj >= 0; --arr_bj) {
//        if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)
//          break;
//
//        if (RC2BinInForwardSense[j][arr_bj] + tmp_rc > OptGap)
//          continue;
//
//        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
//        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
//#endif
//        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
//          auto &kj = label_arr[vec_index];
//          path_rc = kj->RC + tmp_rc;
//          if (path_rc > OptGap) break;
//
//          if ((ki->PI & kj->PI).any()) continue;
//
//          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
//            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
//              if (kj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
//              if (ki->Rank1CutMem[kj->validRank1Cut[l]])path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
//            }
//          }
//
//          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
//            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = ki->validRank1Cut_multi[l];
//              if (kj->Rank1CutMem_multi[tmp_cut] +
//                  ki->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          } else {
//            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
//              int tmp_cut = kj->validRank1Cut_multi[l];
//              if (ki->Rank1CutMem_multi[tmp_cut] +
//                  kj->Rank1CutMem_multi[tmp_cut]
//                  >= R1C_multi_denominator_InCG[tmp_cut]
//                  )
//                path_rc -= r1c_multi_to_pi[tmp_cut];
//            }
//          }
//
//          if (path_rc < OptGap) {
//            //begin to test
//            path_cost = ki->Cost + CostMat4Vertex[i][j] + kj->Cost;
//            auto tmp_PI = ki->PI | kj->PI;
//            if (Tags.find(tmp_PI) == Tags.end()) {
//              Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
//              ++num_routes_now;
//              if (num_routes_now > CONFIG::MaxNumRouteInEnumeration) {
//                status = 1;
//                goto QUIT;
//              }
//            } else if (get<2>(Tags[tmp_PI]) > path_cost) {
//              Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
//            }
//          }
//        }
//      }
//    }
//  }
//  QUIT:
//  cout << "num_routes_now= " << num_routes_now << endl;
//  return status;
//}

int CVRP::concatenateRoutes_prior_forward_InEnumeration(BBNODE *node,
                                                        const double *r1c_to_pi,
                                                        const double *r1c_multi_to_pi,
                                                        unordered_map<yzzLong, tuple<LABEL *, LABEL *, double>> &Tags,
                                                        int &num_routes_now) {
  int status = 0;
  double path_rc, path_cost;
  yzzLong tmp_PI;
#ifdef SYMMETRY_PROHIBIT
  populateRC2TillThisBinNRC2Bin<false>(node);
#else
  populateRC2TillThisBinNRC2Bin<true>(node);//use function here!
#endif
  for (auto &label_list : concatenateLabelsInForwardCG) {
    int i = label_list.first.first;
    int j = label_list.first.second;
    auto &label_vec = label_list.second;
    for (auto &pr : label_vec) {
      auto &ki = pr.first;
      double tmp_mainResource = pr.second;
      double tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];
#ifdef SYMMETRY_PROHIBIT
      int arr_bj = int((tmp_mainResource) / StepSize);
      //most_negative_rc_till_this_bin
      if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)
        continue;

      //most_negative_rc_in_this_bin
      if (RC2BinInBackwardSense[j][arr_bj] + tmp_rc < OptGap &&
          (std::find(node->AllBackwardBuckets[j][arr_bj].BucketArcs.begin(),
                     node->AllBackwardBuckets[j][arr_bj].BucketArcs.end(), i)
              != node->AllBackwardBuckets[j][arr_bj].BucketArcs.end())) {
        //add one more condition for testing capacity
        auto &label_arr = LabelArrayInBackwardSense[j][arr_bj].first;
        auto &label_valid_num = LabelArrayInBackwardSense[j][arr_bj].second;
#else
      int arr_bj = int((MaxMainResource - tmp_mainResource) / StepSize);
      //most_negative_rc_till_this_bin
      if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)
        continue;

      //most_negative_rc_in_this_bin
      if (RC2BinInForwardSense[j][arr_bj] + tmp_rc < OptGap
        //      &&
        //          (std::find(node->AllForwardBuckets[j][arr_bj].BucketArcs.begin(),
        //                     node->AllForwardBuckets[j][arr_bj].BucketArcs.end(), i)
        //              != node->AllForwardBuckets[j][arr_bj].BucketArcs.end())
          ) {
        //add one more condition for testing capacity
        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
#endif
        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
          auto &kj = label_arr[vec_index];
          path_rc = kj->RC + tmp_rc;
          if (path_rc > OptGap) break;
#ifdef SYMMETRY_PROHIBIT
          if (tmp_mainResource > kj->Sum_MainResource) continue;
#else
          if (tmp_mainResource + kj->Sum_MainResource > MaxMainResource) continue;
#endif
          if ((ki->PI & kj->PI).any()) continue;

          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
              if (kj->Rank1CutMem[ki->validRank1Cut[l]]) {
                path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
                if (path_rc > OptGap) goto here;
              }
            }
          } else {
            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
              if (ki->Rank1CutMem[kj->validRank1Cut[l]]) {
                path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
                if (path_rc > OptGap) goto here;
              }
            }
          }

          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
              int tmp_cut = ki->validRank1Cut_multi[l];
              if (kj->Rank1CutMem_multi[tmp_cut] +
                  ki->Rank1CutMem_multi[tmp_cut]
                  >= R1C_multi_denominator_InCG[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > OptGap) goto here;
              }
            }
          } else {
            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
              int tmp_cut = kj->validRank1Cut_multi[l];
              if (ki->Rank1CutMem_multi[tmp_cut] +
                  kj->Rank1CutMem_multi[tmp_cut]
                  >= R1C_multi_denominator_InCG[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > OptGap) goto here;
              }
            }
          }
          //begin to test
          path_cost = ki->Cost + CostMat4Vertex[i][j] + kj->Cost;
          tmp_PI = ki->PI | kj->PI;
          if (Tags.find(tmp_PI) == Tags.end()) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
            ++num_routes_now;
            if (num_routes_now > CONFIG::MaxNumRouteInEnumeration) {
              status = 1;
              goto QUIT;
            }
          } else if (get<2>(Tags[tmp_PI]) > path_cost) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
          }
          here:;
        }
      }
#ifdef SYMMETRY_PROHIBIT
      //bj-1
      for (++arr_bj; arr_bj < NumBucketsPerVertex; ++arr_bj) {
        if (RC2TillThisBinInBackwardSense[j][arr_bj] + tmp_rc > OptGap)
          break;

        if (RC2BinInBackwardSense[j][arr_bj] + tmp_rc > OptGap ||
            (std::find(node->AllBackwardBuckets[j][arr_bj].BucketArcs.begin(),
                       node->AllBackwardBuckets[j][arr_bj].BucketArcs.end(), i)
                == node->AllBackwardBuckets[j][arr_bj].BucketArcs.end()))
          continue;

        auto &label_arr = LabelArrayInBackwardSense[j][arr_bj].first;
        auto &label_valid_num = LabelArrayInBackwardSense[j][arr_bj].second;
#else
      for (--arr_bj; arr_bj >= 0; --arr_bj) {
        if (RC2TillThisBinInForwardSense[j][arr_bj] + tmp_rc > OptGap)
          break;

        if (RC2BinInForwardSense[j][arr_bj] + tmp_rc > OptGap
          //        ||
          //            (std::find(node->AllForwardBuckets[j][arr_bj].BucketArcs.begin(),
          //                       node->AllForwardBuckets[j][arr_bj].BucketArcs.end(), i)
          //                == node->AllForwardBuckets[j][arr_bj].BucketArcs.end())
            )
          continue;

        auto &label_arr = LabelArrayInForwardSense[j][arr_bj].first;
        auto &label_valid_num = LabelArrayInForwardSense[j][arr_bj].second;
#endif
        for (int vec_index = 0; vec_index < label_valid_num; ++vec_index) {
          auto &kj = label_arr[vec_index];
          path_rc = kj->RC + tmp_rc;
          if (path_rc > OptGap) break;

          if ((ki->PI & kj->PI).any()) continue;

          if (ki->numValidRank1Cut < kj->numValidRank1Cut) {
            for (int l = 0; l < ki->numValidRank1Cut; ++l) {
              if (kj->Rank1CutMem[ki->validRank1Cut[l]]) {
                path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
                if (path_rc > OptGap) goto here2;
              }
            }
          } else {
            for (int l = 0; l < kj->numValidRank1Cut; ++l) {
              if (ki->Rank1CutMem[kj->validRank1Cut[l]]) {
                path_rc -= r1c_to_pi[kj->validRank1Cut[l]];
                if (path_rc > OptGap) goto here2;
              }
            }
          }

          if (ki->numValidRank1Cut_multi < kj->numValidRank1Cut_multi) {
            for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
              int tmp_cut = ki->validRank1Cut_multi[l];
              if (kj->Rank1CutMem_multi[tmp_cut] +
                  ki->Rank1CutMem_multi[tmp_cut]
                  >= R1C_multi_denominator_InCG[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > OptGap) goto here2;
              }
            }
          } else {
            for (int l = 0; l < kj->numValidRank1Cut_multi; ++l) {
              int tmp_cut = kj->validRank1Cut_multi[l];
              if (ki->Rank1CutMem_multi[tmp_cut] +
                  kj->Rank1CutMem_multi[tmp_cut]
                  >= R1C_multi_denominator_InCG[tmp_cut]
                  ) {
                path_rc -= r1c_multi_to_pi[tmp_cut];
                if (path_rc > OptGap) goto here2;
              }
            }
          }
          //begin to test
          path_cost = ki->Cost + CostMat4Vertex[i][j] + kj->Cost;
          tmp_PI = ki->PI | kj->PI;
          if (Tags.find(tmp_PI) == Tags.end()) {
            Tags[tmp_PI] = make_tuple(ki, kj, path_cost);
            ++num_routes_now;
            if (num_routes_now > CONFIG::MaxNumRouteInEnumeration) {
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
  cout << "num_routes_now= " << num_routes_now << endl;
  return status;
}

#ifdef SYMMETRY_PROHIBIT
int CVRP::enumerateHalfBackwardRoutes(BBNODE *const node,
                                      const double *r1c_to_pi,
                                      const double *r1c_multi_to_pi,
                                      vector<LABEL *> **copy_bucket) {
  //no sort in this case
  int status = 0;
  int edgemap;
  NumBackwardLabelsInEnu = 0;
  bool if_keep, if_break;
  double path_rc, path_cost;
  initializeLabels(node, 2, false, {true, 2, false});

  auto beg = high_resolution_clock::now();
  auto end = beg;
  auto b4_end = beg;
  auto af_end = beg;
  double eps;
  double eps2;
  double left_time = CONFIG::HardTimeThresholdInAllEnumeration;

  for (int b = NumBucketsPerVertex - 1; b >= 0; --b) {
    int i = 1;
    STILL_EXIST:
    for (; i < Dim; ++i) {
      end = high_resolution_clock::now();
      eps = duration<double>(end - b4_end).count();
      if (eps > left_time) {
        status = 2;
        goto outside;
      }
      auto &valid_num = IfExistExtraLabelsInBackwardSense[i][b].second;
      if (!valid_num) continue;
      auto &label_array = IfExistExtraLabelsInBackwardSense[i][b].first;
      for (int vec_index = 0; vec_index < valid_num; ++vec_index) {
        auto &ki = label_array[vec_index];
        if (ki->if_extended) continue;
        ki->if_extended = true;
        for (int j : node->AllBackwardBuckets[i][b].BucketArcs) {
          if (ki->PI[j]) continue;
          auto &tmp_mainResource = AllLabel[IdxGlo].Sum_MainResource;
          if (!decreaseMainResourceConsumption(ki->Sum_MainResource, tmp_mainResource, i, j)) continue;
          if (tmp_mainResource < MeetPointResourceInBiDirEnu) {
            continue;
          }
          auto &tmp_rc = AllLabel[IdxGlo].RC;
          tmp_rc = ki->RC + ChgCostMat4Vertex[i][j];//real rc
          if_keep = false;

          int arr_bj = int(tmp_mainResource / StepSize);
          //first test
          if (tmp_rc + RC2TillThisBinInForwardSense[j][arr_bj] > OptGap) continue;
          //second test
          if (tmp_rc + RC2BinInForwardSense[j][arr_bj] < OptGap) {
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if (tmp_mainResource < kkj->Sum_MainResource) continue;
              if ((ki->PI & kkj->PI).any()) continue;
              path_rc = tmp_rc + kkj->RC;
              if (path_rc > OptGap) break;
              if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
                  if (kkj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
                  if (ki->Rank1CutMem[kkj->validRank1Cut[l]])path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
                }
              }
              if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = ki->validRank1Cut_multi[l];
                  if (kkj->Rank1CutMem_multi[tmp_cut] +
                      ki->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = kkj->validRank1Cut_multi[l];
                  if (ki->Rank1CutMem_multi[tmp_cut] +
                      kkj->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              }
              if (path_rc < OptGap) {
                if_keep = true;
                goto outside1;
              }
            }
          }
          //real test
          for (--arr_bj; arr_bj >= 0; --arr_bj) {
            //first test
            if (tmp_rc + RC2TillThisBinInForwardSense[j][arr_bj] > OptGap) break;
            //second test
            if (tmp_rc + RC2BinInForwardSense[j][arr_bj] > OptGap) continue;
            //real test
            for (auto &kkj : copy_bucket[j][arr_bj]) {
              if ((ki->PI & kkj->PI).any()) continue;
              path_rc = tmp_rc + kkj->RC;
              if (path_rc > OptGap) break;
              if (ki->numValidRank1Cut < kkj->numValidRank1Cut) {
                for (int l = 0; l < ki->numValidRank1Cut; ++l) {
                  if (kkj->Rank1CutMem[ki->validRank1Cut[l]])path_rc -= r1c_to_pi[ki->validRank1Cut[l]];
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut; ++l) {
                  if (ki->Rank1CutMem[kkj->validRank1Cut[l]])path_rc -= r1c_to_pi[kkj->validRank1Cut[l]];
                }
              }
              if (ki->numValidRank1Cut_multi < kkj->numValidRank1Cut_multi) {
                for (int l = 0; l < ki->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = ki->validRank1Cut_multi[l];
                  if (kkj->Rank1CutMem_multi[tmp_cut] +
                      ki->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              } else {
                for (int l = 0; l < kkj->numValidRank1Cut_multi; ++l) {
                  int tmp_cut = kkj->validRank1Cut_multi[l];
                  if (ki->Rank1CutMem_multi[tmp_cut] +
                      kkj->Rank1CutMem_multi[tmp_cut]
                      >= R1C_multi_denominator_InCG[tmp_cut]
                      )
                    path_rc -= r1c_multi_to_pi[tmp_cut];
                }
              }

              if (path_rc < OptGap) {
                if_keep = true;
                goto outside1;
              }
            }
          }
          outside1:
          if (!if_keep) continue;
          int bj = int(tmp_mainResource / StepSize);
          auto &labelList_j = LabelArrayInBackwardSense[j][bj].first;
          auto &valid_num_j = LabelArrayInBackwardSense[j][bj].second;
          auto &tmp_PI = AllLabel[IdxGlo].PI;
          auto &tmp_Cost = AllLabel[IdxGlo].Cost;
          auto &tmp_Rank1CutMem = AllLabel[IdxGlo].Rank1CutMem;
          auto &tmp_num_valid_rank1_cut = AllLabel[IdxGlo].numValidRank1Cut;
          auto &tmp_valid_rank1_cut = AllLabel[IdxGlo].validRank1Cut;
          tmp_PI = ki->PI;
          tmp_PI.set(j);
          tmp_Cost = ki->Cost + CostMat4Vertex[i][j];
          tmp_Rank1CutMem = ki->Rank1CutMem;
          //tmp_num_valid_rank1_cut do not have to copy
          for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_Rank1CutMem[l] = false;
              tmp_rc -= r1c_to_pi[l];
            } else tmp_Rank1CutMem[l] = true;
          }
          tmp_Rank1CutMem &= get<1>(Vertex2ActiveInOnePricingR1Cs[j]);

          auto &tmp_Rank1CutMem_multi = AllLabel[IdxGlo].Rank1CutMem_multi;
          auto &tmp_num_valid_rank1_cut_multi = AllLabel[IdxGlo].numValidRank1Cut_multi;
          auto &tmp_valid_rank1_cut_multi = AllLabel[IdxGlo].validRank1Cut_multi;
          copy(ki->Rank1CutMem_multi, ki->Rank1CutMem_multi + NumValidR1C_multi_InCG, tmp_Rank1CutMem_multi);
          for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[j])) {
            int tmp_cut = get<0>(l);
            tmp_Rank1CutMem_multi[tmp_cut] += get<1>(l);
            if (tmp_Rank1CutMem_multi[tmp_cut] >= get<2>(l)) {
              tmp_rc -= r1c_multi_to_pi[tmp_cut];
              tmp_Rank1CutMem_multi[tmp_cut] -= get<2>(l);
            }
          }
          for (auto l : get<1>(Vertex2ActiveInOnePricingR1C_multi[j])) tmp_Rank1CutMem_multi[l] = 0;

          if_break = false;

          for (int vec_index_j = 0; vec_index_j < valid_num_j;) {
            auto &kj = labelList_j[vec_index_j];
#ifdef CAPACITY_AS_MAIN_RESOURCE
            if (abs(kj->Sum_MainResource - tmp_mainResource) > TOLERANCE) {
              ++vec_index_j;
              continue;
            }
            if ((kj->PI ^ tmp_PI).none()) {
              if (kj->Cost > tmp_Cost) {
                kj->if_extended = true;
                kj = labelList_j[--valid_num_j];
                --NumBackwardLabelsInEnu;
              } else {
                if_break = true;
                break;
              }
            } else ++vec_index_j;
#else
            if (kj->Cost > tmp_Cost) {
              if (kj->Sum_MainResource < tmp_mainResource) {
                if ((kj->PI ^ tmp_PI).none()) {
                  kj->if_extended = true;
                  kj = labelList_j[--valid_num_j];
                  --NumBackwardLabelsInEnu;
                } else ++vec_index_j;
              } else ++vec_index_j;
            } else {
              if (kj->Sum_MainResource > tmp_mainResource) {
                if ((kj->PI ^ tmp_PI).none()) {
                  if_break = true;
                  break;
                } else ++vec_index_j;
              } else ++vec_index_j;
            }
#endif
          }
          if (if_break) continue;

          tmp_num_valid_rank1_cut = 0;
          for (auto l : get<2>(Vertex2ActiveInOnePricingR1Cs[j])) {
            if (tmp_Rank1CutMem[l]) {
              tmp_valid_rank1_cut[tmp_num_valid_rank1_cut++] = l;
            }
          }

          tmp_num_valid_rank1_cut_multi = 0;
          for (auto l : get<2>(Vertex2ActiveInOnePricingR1C_multi[j])) {
            if (tmp_Rank1CutMem_multi[l]) {
              tmp_valid_rank1_cut_multi[tmp_num_valid_rank1_cut_multi++] = l;
            }
          }

          labelList_j[valid_num_j++] = AllLabel + IdxGlo;
          if (valid_num_j == labelList_j.size()) {
            labelList_j.resize(labelList_j.size() * 2);
          }

          //Seq
          AllLabel[IdxGlo].Seq = AllSeq + SeqBeg;
          copy(ki->Seq, ki->Seq + ki->IdxEndSeq + 1, AllLabel[IdxGlo].Seq);
          AllLabel[IdxGlo].IdxEndSeq = ki->IdxEndSeq + 1;
          *(AllLabel[IdxGlo].Seq + AllLabel[IdxGlo].IdxEndSeq) = j;
          SeqBeg += AllLabel[IdxGlo].IdxEndSeq + 1;
          //EndVertex
          AllLabel[IdxGlo].EndVertex = j;
          //if_extended
          AllLabel[IdxGlo].if_extended = false;
          //bucket
          auto &bucket = IfExistExtraLabelsInBackwardSense[j][bj];
          bucket.first[bucket.second++] = AllLabel + IdxGlo;
          if (bucket.second == bucket.first.size()) {
            bucket.first.resize(bucket.first.size() * 2);
          }

          if (++NumBackwardLabelsInEnu > CONFIG::MaxNumLabelInEnumeration) {
            status = 3;//all labels limit
            goto outside;
          }
          ++IdxGlo;//can be put here, because once go outside, the function will end
          if (IdxGlo == LabelAssign) {
            Rollback = 2;
            goto outside;
          }
        }
      }
      valid_num = 0;
    }
//test if all labels are extended
    for (i = 1; i < Dim; ++i) {
      if (IfExistExtraLabelsInBackwardSense[i][b].second)
        goto STILL_EXIST;
    }
    af_end = high_resolution_clock::now();
    eps2 = duration<double>(af_end - b4_end).count();
    eps = duration<double>(af_end - beg).count();
    left_time = (CONFIG::HardTimeThresholdInAllEnumeration - eps) / (b + 1);
    if (eps2 > left_time) {
      status = 2;
      goto outside;
    }
    b4_end = af_end;
  }
  outside:
  for (int i = 0; i < Dim; ++i) {
    delete[]copy_bucket[i];
  }
  delete[] copy_bucket;
  cout << "Half Backward labeling: num_labels= " << NumBackwardLabelsInEnu << endl;
  if (status)return status;

  for (int i = 1; i < Dim; ++i) {
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      std::stable_sort(LabelArrayInBackwardSense[i][b].first.begin(),
                       LabelArrayInBackwardSense[i][b].first.begin()
                           + LabelArrayInBackwardSense[i][b].second,
                       CmpLabelRCLess);
    }
  }
//  populateRC2TillThisBinNRC2Bin(node, 2);//never use this function in phase 1
  return 0;
}
#endif

int CVRP::generateColsByInspection(BBNODE *node,
                                   bool if_only_need_value) {
  if (node->SizeEnuColPool == 0) return 0;
  if (if_only_need_value) {
    cout << "we only record value without doing anything to the col pool!" << endl;
  }
  int size_pool = node->SizeEnuColPool;
  auto &Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
  auto &mat = node->MatInEnu;

  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))

  RowVectorXd rc = node->Cost4ColsInEnuColPool;

  int num = 0;
  for (auto &it : mat) {
    RowVectorXd dual(it.rows());
    for (int i = 0; i < it.rows(); ++i)dual(i) = Pi[num++];
    rc -= dual * it;
  }

  vector<int> col_to_be_added;
  col_to_be_added.reserve(size_pool);

  int cnt = 0;
  if (if_only_need_value) {
    for (int i = 0; i < size_pool; ++i) {
      if (Deleted_ColsInEnuPool[i]) continue;
      if (rc(i) < RC_TOLERANCE) {
        col_to_be_added.emplace_back(i);
        //not delete cols.
        ++cnt;
        if (cnt == CONFIG::MaxNumRoutesInExact) break;
      }
    }
  } else {
    for (int i = 0; i < size_pool; ++i) {
      if (Deleted_ColsInEnuPool[i]) continue;
      if (rc(i) < RC_TOLERANCE) {
        col_to_be_added.emplace_back(i);
        Deleted_ColsInEnuPool[i] = true;
        ++cnt;
        if (cnt == CONFIG::MaxNumRoutesInExact) break;
      }
    }
  }

  if (col_to_be_added.empty()) {
    if (!if_only_need_value) {
      node->Val = LPVal;
      OptGap = calculateOptGap(node);
      for (int i = 0; i < size_pool; ++i) {
        if (rc(i) > OptGap) {
          Deleted_ColsInEnuPool[i] = true;
        }
      }
      //delete cols in lp
      cleanColsRCLargerThanOptGapInLP(node);
      //delete cols from pool
//#ifdef Use_heuristic_UB
////      if (if_use_heur_enumeration) {
////        //force to delete cols
////        refineArcsHeur(node);
//////          cleanColsInPool(node, nullptr, true);
//////          cout << "tes5!" << endl;
////        cleanColsInPool(node, nullptr);
////      } else {
////        cleanColsInPool(node, nullptr);
////      }
////        cout << "tes6!" << endl;
////        regenerateEnuMat(node, nullptr);
//#else
//#endif
      regenerateEnuMat(node, nullptr);
    }
  } else {
    addColsByInspection(node, col_to_be_added);
  }
  return cnt;
}

void CVRP::generateVertex2IdxCols_N_Edge2IdxCols(BBNODE *node) {
  int size_pool = node->SizeEnuColPool, curr_node;
  if (size_pool == 0) return;
  auto &ptr = node->IdxColsInEnuColPool;
  sparseRowMatrixXd mat(NumRow, size_pool);

  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(NumRow) * double(size_pool) * 0.1));

  auto &Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
  Deleted_ColsInEnuPool = new bool[size_pool]();

  auto &map = node->map_col_pool;
  map.clear();
  for (int i = 0; i < size_pool; ++i) {
    int past_node = 0;
    for (auto j = ptr[i] + 1;; ++j) {
      curr_node = ColPool4Pricing[j];
      if (past_node < curr_node)
        map[{past_node, curr_node}].emplace_back(i);
      else map[{curr_node, past_node}].emplace_back(i);
      if (!curr_node) break;
      triplets.emplace_back(curr_node - 1, i, 1);
      past_node = curr_node;
    }
  }

  sparseRowMatrixXd tmpMat(RealDim, size_pool);
  tmpMat.setFromTriplets(triplets.begin(), triplets.end());

  //vehicle
  for (int i = 0; i < size_pool; ++i) triplets.emplace_back(RealDim, i, 1);

  sparseRowMatrixXd sum(1, size_pool);
  unordered_map<int, double> tmp;
  tmp.reserve(size_pool);

  //for rcc
  for (auto &rcc : node->RCCs) {
    tmp.clear();
    if (rcc.FormRCC) {
      auto &info = rcc.InfoRCCCustomer;
      for (auto iter = info.begin(); iter != info.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != info.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
    } else {
      auto &infoRccCustomer = rcc.InfoRCCCustomer;
      auto &infoRccOutsideCustomer = rcc.InfoRCCOutsideCustomer;
      for (auto iter = infoRccOutsideCustomer.begin(); iter != infoRccOutsideCustomer.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != infoRccOutsideCustomer.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
      for (auto customer_it : infoRccOutsideCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] += 0.5;
      }
      for (auto customer_it : infoRccCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] -= 0.5;
      }
    }
    int row = rcc.IdxRCC;
    for (auto &it : tmp) {
      if (abs(it.second) > TOLERANCE) {
        triplets.emplace_back(row, it.first, it.second);
      }
    }
  }

  for (auto &r1c : node->R1Cs) {
    sum.setZero();
    auto &info = r1c.InfoR1C;
    for (auto j : info) {
      sum += tmpMat.row(j - 1);
    }
    sum /= 2;
    int row = r1c.IdxR1C;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }

  for (auto &r1c : node->R1Cs_multi) {
    sum.setZero();
    auto &info = r1c.InfoR1C;
    const auto &plan = map_rank1_multiplier[(int) info.first.size()][info.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto &j : info.first) {
      sum += tmpMat.row(j - 1) * multi[count++];
    }
    sum /= denominator;
    int row = r1c.IdxR1C;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }

  mat.setFromTriplets(triplets.begin(), triplets.end());
  node->MatInEnu.push_back((mat));

  //clean cols violated by branch constraints
  vector<bool> bad_cols(size_pool);
  for (auto &brc : node->BrCs) {
    if (brc.BrDir) {//must use: 1. use and only use one edge .2 use two but not next to each other
      //revise ColPool_related_info
      int ai = brc.Edge.first, aj = brc.Edge.second;
      if (ai) {
        sum = mat.row(ai - 1) + mat.row(aj - 1);
        std::fill(bad_cols.begin(), bad_cols.end(), true);
        for (int it : map[{ai, aj}]) bad_cols[it] = false;
        for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
          if (it.value() > 0.5) {
            if (it.value() < 1.5) {//==1
              Deleted_ColsInEnuPool[it.col()] = true;
            } else {//==2 further test
              if (bad_cols[it.col()]) Deleted_ColsInEnuPool[it.col()] = true;
            }
          }
        }
      } else {
        std::fill(bad_cols.begin(), bad_cols.end(), true);
        for (int it : map[{0, aj}]) bad_cols[it] = false;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, aj - 1); it; ++it) {
          if (it.value() > 0.5 && bad_cols[it.col()]) {
            Deleted_ColsInEnuPool[it.col()] = true;
          }
        }
      }
    }
  }
  //give to Deleted_ColsInEnuPool
  regenerateEnuMat(node, nullptr, false);
  if (size_pool != node->SizeEnuColPool)
    cout << "now size= " << node->SizeEnuColPool << endl;
}

void CVRP::addColsByInspection(BBNODE *const node, const vector<int> &Col_added) {
  if (Col_added.empty()) return;
#ifdef DEBUG_LP_FILE
  cout << "before add cols, NumCol=" << NumCol << endl;
  writeLP_N_tell_if_LP_corrected(node->solver);
#endif
  auto &ptr = node->IdxColsInEnuColPool;
  auto &ptr_cost = node->Cost4ColsInEnuColPool;
  auto &mat = node->MatInEnu;
  auto col_idx = node->IdxCols;
  size_t nzcnt = 0;
  int ccnt = 0;
  int index = NumCol;

  vector<int> if_col_added(node->SizeEnuColPool, -1);
  int cnt = 0;
  for (auto col : Col_added) {
    if_col_added[col] = cnt++;
  }

  double nonZeros = 0;
  SparseMatrix<double, ColMajor> tmp_mat(NumRow, (int) Col_added.size());
  vector<Triplet<double>> triplets;
  triplets.reserve(Col_added.size() * NumRow * 0.1);
  size_t num = 0;
  for (auto &it : mat) {
    nonZeros += (double) it.nonZeros();
    for (int i = 0; i < it.rows(); ++i) {
      for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator inner_it(it, i); inner_it; ++inner_it) {
        if (if_col_added[inner_it.col()] != -1) {
          triplets.emplace_back(num, if_col_added[inner_it.col()], inner_it.value());
        }
      }
      ++num;
    }
  }

  tmp_mat.setFromTriplets(triplets.begin(), triplets.end());

  num = size_t(min(nonZeros / node->SizeEnuColPool * (int) Col_added.size() * 2,
                   nonZeros));
  if (checkSolverPtr(num)) {
    reallocateSolverPtr(num);
  }

//  auto beg = high_resolution_clock::now();
//  for (int row = 0; row < NumRow; ++row) {
//    int enu_row = LPRow2MatRow[row];
//    for (auto col : Col_added) {
//      if (mat.coeff(enu_row, col) != 0) {
//        idx_map[col].emplace_back(row, mat.coeff(enu_row, col));
//      }
//    }
//  }
//  auto end = high_resolution_clock::now();
//  cout << "time add cols= " << duration<double>(end - beg).count() << endl;
//
//  idx_map.clear();
//  vector<bool> if_col_added(mat.cols(), false);
//  for (auto col : Col_added)if_col_added[col] = true;
//  beg = high_resolution_clock::now();
//  for (int row = 0; row < NumRow; ++row) {
//    int enu_row = LPRow2MatRow[row];
//    for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(mat, enu_row); it; ++it) {
//      if (if_col_added[it.col()]) {
//        idx_map[it.col()].emplace_back(row, it.value());
//      }
//    }
//  }
//  end = high_resolution_clock::now();
//  cout << "time2 add cols= " << duration<double>(end - beg).count() << endl;


  for (auto col : Col_added) {
    col_idx[index++] = ptr[col];
    solver_obj[ccnt] = ptr_cost[col];
    solver_beg[ccnt] = nzcnt;
    for (Eigen::SparseMatrix<double, ColMajor>::InnerIterator it(tmp_mat, ccnt); it; ++it) {
      solver_ind[nzcnt] = (int) it.row();
      solver_val[nzcnt++] = it.value();
    }
    ++ccnt;
  }
  solver_beg[ccnt] = nzcnt;

  if (nzcnt >= MaxNonZeroEntry) throw runtime_error("ccnt>=MaxNonZeroEntry");

  safe_solver(node->solver.SOLVERXaddvars(ccnt,
                                          nzcnt,
                                          solver_beg,
                                          solver_ind,
                                          solver_val,
                                          solver_obj,
                                          nullptr,
                                          nullptr,
                                          nullptr,
                                          nullptr))
  safe_solver(node->solver.SOLVERupdatemodel())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
#ifdef DEBUG_LP_FILE
  cout << "after add cols, NumCol=" << NumCol << endl;
  writeLP_N_tell_if_LP_corrected(node->solver);
#endif
}

void CVRP::regenerateEnuMat(BBNODE *node, BBNODE *node2, bool if_force) {
  bool *Deleted_ColsInEnuPool;
  int size_pool = (int) node->SizeEnuColPool;
  if (!size_pool) {
    delete[] node->Deleted_ColsInEnuPool;
    node->Deleted_ColsInEnuPool = nullptr;
    return;
  }
  BBNODE *out_node;
  int del_size = 0;
  if (!node2) {
    Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
    out_node = node;
    //count the number of true
    for (int i = 0; i < size_pool; ++i) {
      if (Deleted_ColsInEnuPool[i]) ++del_size;
    }
    node->validSize = size_pool - del_size;
    if (!if_force) {
      auto left = double(node->SizeEnuColPool - del_size) / node->SizeEnuColPool;
      if (left > CONFIG::LeftThresholdRCFixing4EnumerationPool) {
        cout << "stash_size= " << del_size << endl;
        return;
      } else {
        cout << "del_size= " << del_size << endl;
      }
    }
  } else {
    Deleted_ColsInEnuPool = node2->Deleted_ColsInEnuPool;
    out_node = node2;
    //count the number of true
    for (int i = 0; i < size_pool; ++i) {
      if (Deleted_ColsInEnuPool[i]) ++del_size;
    }
    node2->validSize = size_pool - del_size;
  }

  int colIndex = 0;
  vector<int> new_col_map(size_pool, -1);
  for (int i = 0; i < size_pool; ++i) {
    if (!Deleted_ColsInEnuPool[i]) {
      new_col_map[i] = colIndex++;
    }
  }

  if (CstrIndex.empty()) {
    safe_solver(node->solver.SOLVERoptimize())
    findNonactiveCuts(node);
    if (CstrIndex.empty()) {
      CstrIndex.resize(NumRow);
      iota(CstrIndex.begin(), CstrIndex.end(), 0);
    }
  }

  //take cols from tmp - cnt and take the rows from LPRow2MatRow - NumRow from node->MatInEnu
  int new_size_pool = size_pool - del_size;
  sparseRowMatrixXd tmpMatrix(NumRow, new_size_pool);

  int rowIndex = 0;
  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(NumRow) * double(new_size_pool) * 0.1));

  for (auto &it : node->MatInEnu) {
    for (int i = 0; i < it.rows(); ++i) {
      if (CstrIndex[rowIndex] != -1) {
        int j = CstrIndex[rowIndex];
        for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator inner_it(it, i); inner_it; ++inner_it) {
          if (new_col_map[inner_it.col()] != -1) {
            triplets.emplace_back(j, new_col_map[inner_it.col()], inner_it.value());
          }
        }
      }
      ++rowIndex;
    }
  }

  tmpMatrix.setFromTriplets(triplets.begin(), triplets.end());

  out_node->MatInEnu.clear();
  out_node->MatInEnu.push_back(std::move(tmpMatrix));

  colIndex = 0;
  for (int i = 0; i < out_node->SizeEnuColPool; ++i) {
    if (!out_node->Deleted_ColsInEnuPool[i]) {
      out_node->Cost4ColsInEnuColPool[colIndex] = out_node->Cost4ColsInEnuColPool[i];
      out_node->IdxColsInEnuColPool[colIndex++] = out_node->IdxColsInEnuColPool[i];
    }
  }
  out_node->SizeEnuColPool = new_size_pool;
  out_node->Cost4ColsInEnuColPool.conservativeResize(out_node->SizeEnuColPool);
  out_node->IdxColsInEnuColPool.conservativeResize(out_node->SizeEnuColPool);
  delete[]Deleted_ColsInEnuPool;
  out_node->Deleted_ColsInEnuPool = new bool[out_node->SizeEnuColPool]();


  //update
  auto &map_col_pool = out_node->map_col_pool;
  map_col_pool.clear();
  auto ptr_col = out_node->IdxColsInEnuColPool;
  for (int i = 0; i < out_node->SizeEnuColPool; ++i) {
    int past_node = 0;
    for (auto j = ptr_col[i] + 1;; ++j) {
      int curr_node = ColPool4Pricing[j];
      if (past_node < curr_node) map_col_pool[{past_node, curr_node}].emplace_back(i);
      else map_col_pool[{curr_node, past_node}].emplace_back(i);
      if (!curr_node) break;
      past_node = curr_node;
    }
  }
  CstrIndex.clear();
}

void CVRP::solveLPByInspection(BBNODE *const node, bool if_only_need_value,
                               bool if_heuristic, bool if_record_sol) {
  node->if_Int = false;
  int ccnt = 0, iter_exact = 0, old_ncol = NumCol;

  time_point<high_resolution_clock, duration<long long, ratio<1L, 1000000000LL>>> mt_beg,
      mt_end, spt_beg, spt_end;
  duration<long long int, ratio<1LL, 1000000000LL>> mt_elap{0}, spt_elap{0};
  spt_beg = high_resolution_clock::now();

  optimizeLP4OneIterInEnu(node);

  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.SOLVERgetenvMethod(&env_method))
  if (env_method != SOLVER_PRIMAL_SIMPLEX) {
    safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_PRIMAL_SIMPLEX))
    if_changed = true;
  }

  //Exact Phase
  //loop : cg->opt_lp->delete column->reOpt_lp
  if (if_heuristic)
    cout << "Heuristic phase begin...\n";
  else
    cout << "Exact phase begin...\n";

//  safe_solver(node->solver.SOLVERsetenvOutputFlag(1, false))
  //print cuts
//  cout << "r1c: " << endl;
//  for (auto &r1c : node->R1Cs) {
//    for (auto &i : r1c.InfoR1C) {
//      cout << i << " ";
//    }
//    cout << endl;
//  }
//  cout << "r2c: " << endl;
//  for (auto &r2c : node->R1Cs_multi) {
//    for (auto &i : r2c.InfoR1C.first) {
//      cout << i << " ";
//    }
//    cout << endl;
//  }
//  cout << endl;

//get coefficient of the first column and compare it with the rhs


  while (true) {

    ++iter_exact;

    spt_beg = high_resolution_clock::now();

    ccnt = generateColsByInspection(node, if_only_need_value);

    if (!ccnt) {
      --iter_exact;
      break;
    }

    spt_end = high_resolution_clock::now();
    spt_elap += duration_cast<milliseconds>(spt_end - spt_beg);
    mt_beg = high_resolution_clock::now();

    optimizeLP4OneIterInEnu(node);

    mt_end = high_resolution_clock::now();
    mt_elap += duration_cast<milliseconds>(mt_end - mt_beg);

    if (!(iter_exact % PRINT_LABELING_STEP_SIZE)) {
      GloEnd = high_resolution_clock::now();
      GloEps = duration<double>(GloEnd - GloBeg).count();
      printInfoLabeling(iter_exact, NumCol - old_ncol, NumCol, NumRow, double(mt_elap.count()) * 1e-9,
                        double(spt_elap.count()) * 1e-9, GloEps,
                        LPVal, LB, UB);
      if (!if_only_need_value) {
        cout << "col_pool= " << node->SizeEnuColPool << "  remain "
             << (MaxNumEnuColPool ? (double(node->SizeEnuColPool) / MaxNumEnuColPool * 100) : 0) << "%\n";
      }
      spt_elap = duration_values<milliseconds>::zero();
      mt_elap = duration_values<milliseconds>::zero();
      old_ncol = NumCol;
    }
  }

//  if (!if_heuristic && !ml.EnumEdgeDeleteCols.empty()) {
//    node->Val = LPVal;
//    OptGap = calculateOptGap(node);
//    int size_pool = (int) node->SizeEnuColPool;
//    auto &edge = node->BrCs.back().Edge;
//    vector<int> *vec_ptr{};
//    if (node->BrCs.back().BrDir) {
//      vec_ptr = &ml.EnumEdgeDeleteCols[edge].Idx4RightBrCol;
//    } else {
//      vec_ptr = &ml.EnumEdgeDeleteCols[edge].Idx4LeftBrCol;
//    }
//    for (int i = 0; i < size_pool; ++i) {
//      if (rc(i) > OptGap) {
//        (*vec_ptr)[i] = 1;
//      }
//    }
////    int stepsize = 10;
////    int length = int(OptGap / stepsize);
////    vector<double> vec(length + 1, 0);
////
////    for (int i = 0; i < size_pool; ++i) {
////      int seq = int(rc(i) / stepsize);
////      if (seq < length)
////        vec[seq] += 1;
////      else {
////        if (rc(i) < OptGap) {
////          cerr << "rc(i) < OptGap";
////        }
////        vec[length] += 1;
////      }
////    }
////    //divide by size_pool
////    transform(vec.begin(), vec.end(), vec.begin(), [size_pool](double x) { return x / size_pool; });
////    for (int i = 0; i < length; ++i) {
////      cout << "rc= " << (i + 1) * stepsize << "  num= " << vec[i] << endl;
////    }
////    cout << "Eliminate num= " << vec[length] << endl;
//  }

//  if (!if_heuristic && if_populate_DumpIdx4BrCol) {
//    node->Val = LPVal;
//    OptGap = calculateOptGap(node);
//    int size_pool = (int) node->SizeEnuColPool;
//    auto &edge = node->BrCs.back().Edge;
//    vector<size_t> *vec_ptr{};
//    if (node->BrCs.back().BrDir) {
//      vec_ptr = &(DumpIdx4RightBrCol);
//    } else {
//      vec_ptr = &(DumpIdx4LeftBrCol);
//    }
//    safe_solver(node->solver.SOLVERgetRC(0, NumCol, RC))
//    int keep = 1;
//    safe_solver(node->solver.SOLVERgetObj(0, NumCol, solver_obj))
//    for (int i = keep; i < NumCol; ++i) {
//      if (abs(solver_obj[i] - UB) < 1e-6) continue;//it is ok, even if the UB is updated!
//      if (RC[i] > OptGap) {
//        (*vec_ptr).emplace_back(node->IdxCols[i]);
//      }
//    }
//    for (int i = 0; i < size_pool; ++i) {
//      if (node->Deleted_ColsInEnuPool[i]) continue;
//      if (rc(i) > OptGap) {
//        (*vec_ptr).emplace_back(node->IdxColsInEnuColPool[i]);
//      }
//    }
//    cout << "DumpIdx4BrCol size= " << (*vec_ptr).size() << endl;
//  }

//#ifdef ENV_ENU_ML
//  if (!if_heuristic) ml.rc = rc;
//#endif

  //update the node->Val
  if (if_record_sol) {
    recordOptCol(node);

//    cout << "Test here!" << endl;
//
//    unordered_map<int, int> map1;
//    {
//      node->Val = LPVal;
//      OptGap = calculateOptGap(node);
//      int size_pool = (int) node->SizeEnuColPool;
//      double half_gap = OptGap / 2;
//      int tmp_all = 0;
//      for (int i = 0; i < size_pool; ++i) {
//        if (node->Deleted_ColsInEnuPool[i]) continue;
//        ++tmp_all;
//        if (rc(i) > half_gap) {
//          ++map1[2];
//        }
//      }
//      cout << "half_larger= " << map1[2] << endl;
//      cout << "smaller= " << tmp_all - map1[2] << endl;
//    }
  }

  if (if_changed) safe_solver(node->solver.SOLVERsetenvMethod(env_method))
  if (iter_exact % PRINT_LABELING_STEP_SIZE || !iter_exact) {
    GloEnd = high_resolution_clock::now();
    GloEps = duration<double>(GloEnd - GloBeg).count();
    printInfoLabeling(iter_exact, NumCol - old_ncol, NumCol, NumRow, double(mt_elap.count()) * 1e-9,
                      double(spt_elap.count()) * 1e-9, GloEps,
                      LPVal, LB, UB);
    if (!if_only_need_value) {
      cout << "col_pool= " << node->SizeEnuColPool << "  remain "
           << (MaxNumEnuColPool ? (double(node->SizeEnuColPool) / MaxNumEnuColPool * 100) : 0) << "%\n";
    }
  }
}

void CVRP::cleanColsNonEle(BBNODE *const node) {
  //used in enumerated state and only once and after the cols are enumerated in pricing pool
  int len = 0, keep = 0, current_node;
  bool if_break;
  yzzLong PI;
  for (int i = 0; i < NumCol; ++i) {
    if_break = false;
    PI = 0;
    for (size_t j = node->IdxCols[i] + 1;; ++j) {
      current_node = ColPool4Pricing[j];
      if (!current_node) break;
      if (PI[current_node]) {
        solver_ind[len++] = i;
        if_break = true;
        break;
      }
      PI.set(current_node);
    }
    if (!if_break) {
      node->IdxCols[keep++] = node->IdxCols[i];
    }
  }

  if (len) {
    safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  }

  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  cout << "after clean non-ele routes lpval= " << LPVal << endl;
}

void CVRP::deleteBrCsNR1C1s(BBNODE *const node) {
  //never touch the first column
  //remember always clean enu_col_pool as well, because the rc <= optGap not just 0
  if (node->BrCs.empty() && node->R1Cs.empty()) return;
  if (!node->BrCs.empty()) {
    set<int> delete_col;
    int *ai_col = new int[NumCol];
    int *aj_col = new int[NumCol];
    int tmp;
    size_t numnzP;
    int ai, aj;
    for (auto &brc : node->BrCs) {
      ai = brc.Edge.first;
      aj = brc.Edge.second;
      tmp = ai * Dim + aj;
      //check if coefficient is zero, if i and j has one 1, delete them directly
      if (brc.BrDir) {
        //in lp
        for (int i = 0; i < NumCol; ++i) {
          ai_col[i] = 0;
          aj_col[i] = 0;
        }
        if (ai) {
          safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, ai - 1, 1))
          for (size_t i = 0; i < numnzP; ++i)ai_col[solver_ind[i]] = 1;
        }
        safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, aj - 1, 1))
        for (size_t i = 0; i < numnzP; ++i)aj_col[solver_ind[i]] = 1;
        safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, brc.IdxBrC, 1))
        for (int j = 0; j < solver_ind[0]; ++j)if (aj_col[j] || ai_col[j]) delete_col.insert(j);
        for (size_t i = 1; i < numnzP; ++i)
          for (int j = solver_ind[i - 1] + 1; j < solver_ind[i]; ++j)
            if (aj_col[j] || ai_col[j])delete_col.insert(j);
        for (int j = solver_ind[numnzP - 1] + 1; j < NumCol; ++j)if (aj_col[j] || ai_col[j])delete_col.insert(j);
      } else {//delete all cstrs that coefficient == 1
        //in lp
        safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, brc.IdxBrC, 1))
        for (size_t i = 0; i < numnzP; ++i) delete_col.insert(solver_ind[i]);
      }
    }

    //begin to delete all cstr
    //remove all cols whose coefficients are all 1
    //in lp
    for (int i = 0; i < NumCol; ++i) ai_col[i] = 0;
    for (auto col : delete_col) ai_col[col] = 1;
    //keep the first col
    ai_col[0] = 0;

    int len = 0, keep = 0;
    for (int i = keep; i < NumCol; ++i) {
      if (ai_col[i]) solver_ind[len++] = i;
      else node->IdxCols[keep++] = node->IdxCols[i];
    }

    if (len) {
      safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
    }
    delete[]ai_col;
    delete[]aj_col;
    cout << "delete " << len << " cols due to replacing the old columns" << endl;
  }

  //delete cstrs
  int len = 0, keep;
  auto cstr_index = new int[NumRow];
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);
  iota(cstr_index, cstr_index + NumRow, 0);
  for (auto &brc : node->BrCs) {
    if (brc.IdxBrC == -1) continue;
    keep = brc.IdxBrC;
    solver_ind[len++] = keep;
    cstr_index[keep] = -1;
    deleted_cstrs.emplace_back(keep);
  }
  for (auto &r1c : node->R1Cs) {
    if (r1c.InfoR1C.size() == 1) {
      keep = r1c.IdxR1C;
      solver_ind[len++] = keep;
      cstr_index[keep] = -1;
      deleted_cstrs.emplace_back(keep);
    }
  }
  if (deleted_cstrs.empty()) {
    delete[]cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  //recover RCCs
  for (auto &rcc : node->RCCs) rcc.IdxRCC = cstr_index[rcc.IdxRCC];

  //delete r1c
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //recover r1c_multi
  for (auto &r1c : node->R1Cs_multi) r1c.IdxR1C = cstr_index[r1c.IdxR1C];

  //recover brc
  for (auto &brc : node->BrCs) brc.IdxBrC = -1;

  safe_solver(node->solver.SOLVERdelconstrs(len, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  cout << "after clean unsatisfied-br-routes and rank1-1 cuts lpval= " << LPVal << endl;
  cout << "NumCol= " << NumCol << ", NumRow= " << NumRow << endl;
  delete[]cstr_index;
}

void CVRP::cleanColsRCLargerThanOptGapInLP(BBNODE *const node) {
  //the first artificial cols will be saved
  OptGap = calculateOptGap(node);
  safe_solver(node->solver.SOLVERgetRC(0, NumCol, RC))
  int len = 0, keep = 1;
  for (int i = keep; i < NumCol; ++i) {
    if (RC[i] > OptGap) {
      solver_ind[len++] = i;
    } else {
      node->IdxCols[keep++] = node->IdxCols[i];
    }
  }
  safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
}

void CVRP::organizeColsInMem2Pricing(BBNODE *const node) {

#ifdef Use_heuristic_UB
  if (if_use_heur_enumeration) {
    cout << "this is a heuristic enumeration!" << endl;
    if (copyColPool4Pricing) throw runtime_error("copyColPool4Pricing is not null!");
    copyColPool4Pricing = new int[PoolBeg4Pricing];
    memcpy(copyColPool4Pricing, ColPool4Pricing, sizeof(int) * PoolBeg4Pricing);
    PoolBeg4copyPricing = PoolBeg4Pricing;
  }
#endif
  int num_parent_cols = node->NumParentCols;
  size_t start;
  PoolBeg4Pricing = 0;
  //first we rearrange the cols in pricing to pricing
  for (int i = node->NumParentCols; i < NumCol; ++i) {
    start = PoolBeg4Pricing;
    ColPool4Pricing[PoolBeg4Pricing] = 0;
    ++PoolBeg4Pricing;
    for (size_t j = node->IdxCols[i] + 1;; ++j) {
      if (!(ColPool4Pricing[j]))break;
      ColPool4Pricing[PoolBeg4Pricing] = ColPool4Pricing[j];
      ++PoolBeg4Pricing;
    }
    ColPool4Pricing[PoolBeg4Pricing] = 0;
    ++PoolBeg4Pricing;
    node->IdxCols[i] = start;
  }
  //second we rearrange the cols in mem to pricing
  for (int i = 0; i < num_parent_cols; ++i) {
    start = PoolBeg4Pricing;
    ColPool4Pricing[PoolBeg4Pricing] = 0;
    ++PoolBeg4Pricing;
    for (size_t j = node->IdxCols[i] + 1;; ++j) {
      if (!(ColPool4Mem[j]))break;
      ColPool4Pricing[PoolBeg4Pricing] = ColPool4Mem[j];
      ++PoolBeg4Pricing;
    }
    ColPool4Pricing[PoolBeg4Pricing] = 0;
    ++PoolBeg4Pricing;
    node->IdxCols[i] = start;
  }
}

void CVRP::recoverR1CsInEnu(BBNODE *const node) {
  //just traverse the previous sequence (and then modify the corresponding coefficients)
  if (node->R1Cs.empty() && node->R1Cs_multi.empty()) return;
  int *ai_col = new int[NumCol];
  auto *sum_col = new double[NumCol];
  vector<double> val_(NumCol);
  size_t nzcnt = 0;
  int index;
  size_t numnzP;
  //r1c
  for (auto &r1c : node->R1Cs) {
    index = r1c.IdxR1C;
    memset(sum_col, 0, sizeof(double) * NumCol);
    for (auto i : r1c.InfoR1C) {
      safe_solver(node->solver.SOLVERXgetconstrs(
          &numnzP, solver_beg, ai_col, val_.data(), i - 1, 1))
      for (size_t j = 0; j < numnzP; ++j) ++sum_col[ai_col[j]];
    }
    for (int i = 0; i < NumCol; ++i) {
      int val = int(sum_col[i] / 2 + TOLERANCE);
      if (val) {
        solver_ind[nzcnt] = index;
        solver_ind2[nzcnt] = i;
        solver_val[nzcnt++] = val;
      }
    }
  }

  //r1c_multi
  for (auto &r1c : node->R1Cs_multi) {
    index = r1c.IdxR1C;
    memset(sum_col, 0, sizeof(double) * NumCol);
    const auto &plan = map_rank1_multiplier[(int) r1c.InfoR1C.first.size()][r1c.InfoR1C.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto &i : r1c.InfoR1C.first) {
      safe_solver(node->solver.SOLVERXgetconstrs(
          &numnzP, solver_beg, ai_col, val_.data(), i - 1, 1))
      for (size_t j = 0; j < numnzP; ++j) sum_col[ai_col[j]] += multi[count];
      ++count;
    }
    for (int i = 0; i < NumCol; ++i) {
      int val = int(sum_col[i] / denominator + TOLERANCE);
      if (val) {
        solver_ind[nzcnt] = index;
        solver_ind2[nzcnt] = i;
        solver_val[nzcnt++] = val;
      }
    }
  }

  safe_solver(node->solver.SOLVERXchgcoeffs(nzcnt, solver_ind, solver_ind2, solver_val))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  cout << "after recover rank1c lpval= " << LPVal << endl;
  delete[]ai_col;
  delete[]sum_col;
#ifdef check_if_cuts_are_legal
  checkIfCutsLegal(node);
#endif
}

void CVRP::recordOptColInEnu(BBNODE *const node) {

  safe_solver(node->solver.SOLVERupdatemodel())
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  node->Val = LPVal;
  if (ceil_transformed_number_related(node->Val - TOLERANCE) + TOLERANCE >= UB) {
    node->if_terminated = true;
    return;
  }

  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))

  vector<pair<size_t, double>> tmp_solindex;

  vector<pair<int, double>> tmp_first_dim;

//  auto &tmp_opt= node->optSolInIdxCol;
//  tmp_opt.clear();

  for (int i = 0; i < Dim; ++i) {
    if (X[i] > TOLERANCE) {
      tmp_first_dim.emplace_back(i, X[i]);
      tmp_solindex.emplace_back(node->IdxCols[i], X[i]);
//      tmp_opt.emplace_back(i, X[i]);
    }
  }
  for (int i = Dim; i < NumCol; ++i) {
    if (X[i] > TOLERANCE) {
      tmp_solindex.emplace_back(node->IdxCols[i], X[i]);
//      tmp_opt.emplace_back(i, X[i]);
    }
  }

  node->Idx4LPSolsInColPool = tmp_solindex;

  //update the edge info
  for (int i = 0; i < Dim; ++i) for (int j = 0; j < Dim; ++j) ArcGraph[i][j] = 0;

  for (auto &i : tmp_solindex) {
    for (size_t j = i.first;;) {
      ArcGraph[ColPool4Pricing[j]][ColPool4Pricing[j + 1]] += i.second;
      if (!ColPool4Pricing[++j]) break;
    }
  }

  for (int i = 0; i < Dim; ++i)
    for (int j = i + 1; j < Dim; ++j) {
      ArcGraph[i][j] += ArcGraph[j][i];
      ArcGraph_revised[i][j] = ArcGraph[i][j];
    }

  //revise the arti_cols from 1 to fx_Dimension, we always keep the first col
  for (auto &i : tmp_first_dim) {
    if (i.first) {
      int valid_length = 0, curr_node;
      for (size_t j = node->IdxCols[i.first] + 1;; ++j) {
        curr_node = ColPool4Pricing[j];
        if (!curr_node) break;
        ++valid_length;
      }
      curr_node = ColPool4Pricing[node->IdxCols[i.first] + 1];
      if (valid_length == 1)ArcGraph_revised[0][curr_node] -= i.second;
      else break;
    }
  }

  bool if_reformulate = true;
  for (int i = 0; i < Dim; ++i) {
    for (int j = i + 1; j < Dim; ++j) {
      if (ArcGraph_revised[i][j] > TOLERANCE && abs(ArcGraph_revised[i][j] - 1) > TOLERANCE) {
        if_reformulate = false;
        goto outside;
      }
    }
  }

  outside:
  if (if_reformulate) {
    reformulateIPSol(node);
  }

  NumEdge = 0;

  for (int i = 0; i < Dim; ++i) for (int j = i + 1; j < Dim; ++j) if (ArcGraph[i][j] > TOLERANCE) ++NumEdge;

  if (NumEdge > MaxNumEdge) MaxNumEdge *= 2;
}

void CVRP::optimizeLP4OneIterInEnu(BBNODE *const node) {

  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))

  node->if_Int = true;
  for (int i = 0; i < NumCol; ++i)
    if (X[i] > TOLERANCE && abs(X[i] - 1) > TOLERANCE) {
      node->if_Int = false;
      break;
    }

  if (node->if_Int) {
    if (ceil_transformed_number_related(LPVal - TOLERANCE) + TOLERANCE < UB) {
      UB = ceil_transformed_number_related(LPVal - TOLERANCE);

      IPOptSol.clear();

      for (int i = 0; i < NumCol; ++i)
        if (X[i] > TOLERANCE) {
          vector<int> tmp;
          tmp.reserve(MaxLengthEleRoute);
          tmp.emplace_back(0);
          for (auto j = node->IdxCols[i] + 1;; ++j) {
            if (!ColPool4Pricing[j])
              break;
            tmp.emplace_back(ColPool4Pricing[j]);
          }
          tmp.emplace_back(0);
          IPOptSol.emplace_back(std::move(tmp));
        }

#ifdef MASTER_VALVE_ML
      updateUB_EdgeSol();
#endif

      if (ceil_transformed_number_related(LB_transformed - TOLERANCE) + TOLERANCE >= UB) {
        node->Val = LPVal;
        node->if_terminated = true;
        cout << TERMINATED_MESSAGE_PROMISING_UPDATE_UB;
        return;
      }
    }
  }
}

void CVRP::terminateByMIP(BBNODE *node) {
  cout << "b4_col= " << NumCol << endl;
  regenerateEnuMat(node, nullptr, true);
  auto &Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
  vector<int> added_cols;
  added_cols.reserve(node->SizeEnuColPool);
  for (int i = 0; i < node->SizeEnuColPool; ++i) {
    if (Deleted_ColsInEnuPool[i])continue;
    added_cols.emplace_back(i);
  }
  addColsByInspection(node, added_cols);
//    auto name = "UB_" + to_string(UB) + "_col_" + to_string(NumCol) + ".lp";
//    safe_solver(node->solver.SOLVERwrite(name.c_str()))
  cout << "after_cols= " << NumCol << endl;
//    exit(-1);
  //second change lp to mip
  solveMIP(node, true);
  if (if_MIP_enumeration_suc) {
    node->if_terminated = true;
  } else {
    cout << "MIP enumeration failed, seek for branch!\n";
    node->if_terminated = false;
    node->map_col_pool.clear();
    node->SizeEnuColPool = 0;
    node->MatInEnu.clear();
    node->Cost4ColsInEnuColPool.resize(0);
    node->IdxColsInEnuColPool.resize(0);
  }
}

void CVRP::terminateNode(BBNODE *&root_node) {
  If_in_Enu_State = true;
  int root_index = root_node->Idx;
  time_point<high_resolution_clock> beg, end;
  double eps;

  cout << "To terminate the current node, " <<
       "chg node selection strategy by exploring the children nodes of the node first!\n";

  //second we construct a BBT tree locally...
  subBBT.push(root_node);
  --NumExploredNodes;
  for (auto &r1c : root_node->R1Cs) {
    r1c.Mem.clear();
  }
  for (auto &r1c : root_node->R1Cs_multi) {
    r1c.Mem.clear();
  }

#ifdef useM_dynamically
  bool if_reset{false};
#endif

  while (!subBBT.empty()
      ) {

    auto node = subBBT.top();
    subBBT.pop();

    ++NumExploredNodes;

    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))

    if (ceil_transformed_number_related(node->Val - TOLERANCE) + TOLERANCE >= UB) {
      cout << "The subtree rooted by " << root_index <<
           " has been terminated for tree_lb= " << node->Val << " but ub= " << UB << endl;
      delete node;
      goto DELETE;
    }

    cout << BIG_PHASE_SEPARATION;
    GloEnd = std::chrono::high_resolution_clock::now();
    GloEps = duration<double>(GloEnd - GloBeg).count();
    cout << "nd_ind= " << node->Idx << "  nd_col= " << NumCol << "  nd_val= " << node->Val << "  nd_dep= "
         << node->TreeLevel << "  et= " << GloEps << "  lb= " << LB << "  ub= "
         << UB << "  sub_nd_rmn= " << subBBT.size() << "  ColPool= " << node->SizeEnuColPool << endl;
    for (auto &brc : node->BrCs) {
      cout << (brc.BrDir ? "true" : "false") << "(" << brc.Edge.first << "," << brc.Edge.second << ")" << endl;
    }

#ifdef find_missed_solution
    vector<int> data;
    bool if_suc;
    findWhySolutionDisappear(node, this, data, if_suc);
    if (!if_suc) {
      cout << "This is not the optimal node!" << endl;
      delete node;
      continue;
    }
#endif

#ifdef check_enumeration_pool_unsatisfied_cols_by_br
    unordered_map<pair<int, int>, vector<int>, PairHasher> unsatisfied_cols;
    for (int i = 0; i < node->SizeEnuColPool; ++i) {
      auto idx = node->IdxColsInEnuColPool[i];
      int past_node = 0;
      for (auto j = idx + 1;; ++j) {
        int cur_node = ColPool4Pricing[j];
        if (past_node < cur_node)
          unsatisfied_cols[{past_node, cur_node}].emplace_back(i);
        else
          unsatisfied_cols[{cur_node, past_node}].emplace_back(i);
        if (!cur_node) {
          break;
        }
        past_node = cur_node;
      }
    }
    for (auto &brc : node->BrCs) {
      if (brc.BrDir) {
        int ai = brc.Edge.first, aj = brc.Edge.second;
        //only check aj is enough
        set<int> un_cols;
        for (int i = 0; i < aj; ++i) {
          if (i != ai) {
            for (auto &col : unsatisfied_cols[{i, aj}]) {
              if (un_cols.find(col) != un_cols.end()) {
                cout << "The edge (" << brc.Edge.first << "," << brc.Edge.second
                     << ") dir: true! is not satisfied by the enumeration pool!" << endl;
                cout << "col= " << col << endl;
                safe_Hyperparameter(1)
              }
              un_cols.emplace(col);
            }
          }
        }
        for (int i = aj + 1; i < Dim; ++i) {
          for (auto &col : unsatisfied_cols[{aj, i}]) {
            if (un_cols.find(col) != un_cols.end()) {
              cout << "The edge (" << brc.Edge.first << "," << brc.Edge.second
                   << ") dir: true! is not satisfied by the enumeration pool!" << endl;
              cout << "col= " << col << endl;
              safe_Hyperparameter(1)
            }
            un_cols.emplace(col);
          }
        }
      } else {
        if (unsatisfied_cols.find(brc.Edge) != unsatisfied_cols.end()) {
          cout << "The edge (" << brc.Edge.first << "," << brc.Edge.second
               << ") dir: false! is not satisfied by the enumeration pool!" << endl;
          cout << "The unsatisfied cols size is: " << unsatisfied_cols[brc.Edge].size() << endl;
          safe_Hyperparameter(1)
        }
      }
    }
#endif

    auto old_val = node->Val;

#ifdef TrainingDataTreeLevel
    if (node->TreeLevel > TrainingDataTreeLevel) {
      delete node;
      continue;
    }
#endif

#ifdef check_if_cuts_are_legal
    checkIfCutsLegal(node);
//    if (node->Idx == 14) {
//      cout << "we write lp!" << endl;
//      auto &Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
//      vector<int> added_cols;
//      added_cols.reserve(node->SizeEnuColPool);
//      for (int i = 0; i < node->SizeEnuColPool; ++i) {
//        if (Deleted_ColsInEnuPool[i])continue;
//        added_cols.emplace_back(i);
//      }
//      addColsByInspection(node, added_cols);
//      safe_solver(node->solver.SOLVERwrite("check.lp"))
//      exit(0);
//    }
#endif

#ifdef useM_dynamically
    auto beg_c = std::chrono::high_resolution_clock::now();
#endif

    MaxNumEnuColPool = node->SizeEnuColPool;

    beg = high_resolution_clock::now();

    solveLPByInspection(node, false, false, true);

    end = high_resolution_clock::now();
    eps = duration<double>(end - beg).count();

    cout << "time for inspection is " << eps << endl;

#ifdef check_lower_bound
    if (node->Val < old_value - TOLERANCE) {
      cerr << "The lower bound is not correct!" << endl;
      cerr << "node->Val= " << node->Val << " old_value= " << old_value << endl;
      safe_Hyperparameter(1);
    }
#endif

    if (node->Idx) {
      auto &brc = node->BrCs.back();
      if (brc.BrDir) {
        RealImprovement_up[brc.Edge].first += node->Val - old_val;
        ++RealImprovement_up[brc.Edge].second;
      } else {
        RealImprovement_down[brc.Edge].first += node->Val - old_val;
        ++RealImprovement_down[brc.Edge].second;
      }
    }

    if (node->if_terminated) {
      delete node;
      continue;
    } else if (NumCol + node->validSize <= MaxNumRoute4MIP) {
      terminateByMIP(node);
      if (node->if_terminated) {
        delete node;
        continue;
      }
    }

#ifndef CutAtNonRootNotAllowed
    //    sepHybridCutsInEnu(node);
    sepHybridCuts(node);
//    terminateByMIP(node);
//    if (node->if_terminated) {
//      delete node;
//      continue;
//    }
    if (node->if_terminated) {
      delete node;
      continue;
    } else if (NumCol + node->validSize <= MaxNumRoute4MIP) {
      terminateByMIP(node);
      if (node->if_terminated) {
        delete node;
        continue;
      }
    }
#else
    //    findNonactiveCuts(node); no need here, transfor to add BrC
    //    printCutsInfo(node);
#endif

#ifdef BranchInEnuNotAllowed
    cout << "In this version, branching in enumeration is not allowed!" << endl;
#ifndef Resolve_Ins_with_Optimal_Val
    exit(0);
#endif
#endif

#ifdef DELUXING_APPLIED
    cout << "Trying to apply RCF..." << endl;
    applyRCF(node, DELUXING_ROUND, true);
    solveLPByInspection(node, false, false, true);
    if (NumCol + node->validSize <= MaxNumRoute4MIP) {
      terminateByMIP(node);
      if (node->if_terminated) {
        delete node;
        continue;
      }
    }
#endif

#ifdef writeEnumerationTrees
    writeEnuTree(node);
    delete node;
    continue;
#endif

    BRANCH:
#ifdef useM_dynamically
    if (node->Idx) {
      auto end_c = std::chrono::high_resolution_clock::now();
      auto eps = duration<double>(end_c - beg_c).count();
      //reset c in this case! (since the pricing for enumeration can be different!)
      if (if_reset)
        node->updateState(eps, node->c, node->BrCs.size() - 1);
      else {
        node->c = eps;
        if_reset = true;
      }
      double new_r;
      node->calculateR_star(node->Val - old_val, new_r, this);
      node->updateState(new_r, node->geo_r_star, node->BrCs.size() - 1);
      cout << "node->c= " << node->c << " node->geo_r_star= " << node->geo_r_star << endl;
    }
#endif

    cout << "Begin branching...\n";

    regenerateEnuMat(node, nullptr, true);

    recordOptCol(node);

    getInfoEdge(node, true);

    writeMap_Edge_ColIdx_in_Enu(node);

    pair<int, int> info;

    do_SB(node, info);

    addBrCut2UnsolvedInEnu(node, info);

    ++NumBr;
  }

  DELETE:
  while (!subBBT.empty()) {
    BBNODE *extra_node = subBBT.top();
    subBBT.pop();
    delete extra_node;
  }

  root_node = nullptr;
  If_in_Enu_State = false;
}

void
CVRP::addBrCut2UnsolvedInEnu(BBNODE *const node, const std::pair<int, int> &info) {
  if (node->if_terminated) return;
  BrC bf;
  bf.Edge = info;
  ++node->TreeLevel;

  bf.BrDir = true;
  bf.IdxBrC = -1;

  CstrIndex.resize(NumRow);
  iota(CstrIndex.begin(), CstrIndex.end(), 0);

  auto node2 = new BBNODE(node, NumCol, IdxNode + 2, bf);
  reviseEnuColInfoByBrC(node, node2, bf);
  //if we use the last step the num_col will be changed
  safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  //保留左node
  bf.BrDir = false;
  node->BrCs.emplace_back(bf);

  CstrIndex.resize(NumRow);
  iota(CstrIndex.begin(), CstrIndex.end(), 0);
  reviseEnuColInfoByBrC(node, node, bf);
  node->Idx = ++IdxNode;
  ++IdxNode;

  subBBT.push(node2);
  subBBT.push(node);
}

void CVRP::reviseEnuColInfoByBrC(BBNODE *node, BBNODE *out_node, const BrC &bf) {
  int ai = bf.Edge.first, aj = bf.Edge.second;
  if (ai >= aj) cerr << "Wrong in reviseEnuColInfoByBrC!" << endl;
  int col_idx;
  int size_pool = int(node->SizeEnuColPool);
  auto Deleted_ColsInEnuPool = out_node->Deleted_ColsInEnuPool;
  //revise lp info
  cleanColsInLPByNewBrC(out_node, bf);
  if (!node->SizeEnuColPool) return;
  auto &mat0 = *node->MatInEnu.begin();

  if (bf.BrDir) {//must use: 1. use and only use one edge .2 use two but not next to each other
    //revise ColPool_related_info
    vector<bool> bad_cols(size_pool, true);
    for (int it : node->map_col_pool[{ai, aj}]) bad_cols[it] = false;
    if (ai) {
      sparseRowMatrixXd tmp = mat0.row(ai - 1) + mat0.row(aj - 1);
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, 0); it; ++it) {
        if (it.value() > 0.5) {
          if (it.value() < 1.5) {//==1
            Deleted_ColsInEnuPool[it.col()] = true;
          } else {//==2 further test
            if (bad_cols[it.col()]) Deleted_ColsInEnuPool[it.col()] = true;
          }
        }
      }
    } else {
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat0, aj - 1); it;
           ++it) {
        if (it.value() > 0.5 && bad_cols[it.col()]) Deleted_ColsInEnuPool[it.col()] = true;
      }
    }
  } else {//use two and two adjacent
    for (int it : node->map_col_pool[{ai, aj}]) Deleted_ColsInEnuPool[it] = true;
  }
  regenerateEnuMat(node, out_node);
}

void CVRP::cleanColsInLPByNewBrC(BBNODE *const node, const BrC &bf) {
  //always keep the first col
  vector<int> delete_col;
  int *sum_col = new int[NumCol]();
  int col, curr_node;
  size_t numnzP;
  bool if_keep;
  int ai = bf.Edge.first, aj = bf.Edge.second;
  if (ai >= aj) cerr << "Wrong in cleanColsInLPByNewBrC\n";
  if (!ai) {//=0
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, aj - 1, 1))
    if (bf.BrDir) {//一定要使用 // 找到有aj的col，然后判断其sequence，如果不是在开头或者结尾，就去除
      for (size_t i = 0; i < numnzP; ++i) {
        col = solver_ind[i];
        auto j = node->IdxCols[col] + 1;
        curr_node = ColPool4Pricing[j];
        if (curr_node == aj)continue;
        ++j;
        if_keep = true;
        for (;; ++j) {
          curr_node = ColPool4Pricing[j];
          if (!curr_node) {
            if (ColPool4Pricing[j - 1] == aj) {
              if_keep = false;
            }
            break;
          }
        }
        if (!if_keep) continue;
        delete_col.emplace_back(col);
      }
    } else {//find col that traverse aj, check its sequence, if it's in the beginning or end, remove it
      for (size_t i = 0; i < numnzP; ++i) {
        col = solver_ind[i];
        auto j = node->IdxCols[col] + 1;
        curr_node = ColPool4Pricing[j];
        if (curr_node == aj) {
          delete_col.emplace_back(col);
          continue;
        }
        ++j;
        if_keep = false;
        for (;; ++j) {
          curr_node = ColPool4Pricing[j];
          if (!curr_node) {
            if (ColPool4Pricing[j - 1] == aj) {
              if_keep = true;
            }
            break;
          }
        }
        if (!if_keep) continue;
        delete_col.emplace_back(col);
      }
    }
  } else {
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, ai - 1, 1))
    for (size_t i = 0; i < numnzP; ++i)++sum_col[solver_ind[i]];
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, aj - 1, 1))
    for (size_t i = 0; i < numnzP; ++i)++sum_col[solver_ind[i]];
    if (bf.BrDir) {//if 1 exists, remove it directly, if 2 exists, check if connect together
      for (int i = 1; i < NumCol; ++i) {
        if (sum_col[i]) {
          if (sum_col[i] == 1) {
            delete_col.emplace_back(i);
          } else {//==2
            //double check
            if_keep = true;
            for (auto j = node->IdxCols[i] + 1;; ++j) {
              curr_node = ColPool4Pricing[j];
              if (!curr_node) break;
              if (curr_node == ai) {
                if (ColPool4Pricing[j + 1] == aj) {
                  if_keep = false;
                  break;
                }
              } else if (curr_node == aj) {
                if (ColPool4Pricing[j + 1] == ai) {
                  if_keep = false;
                  break;
                }
              }
            }
            if (!if_keep) continue;
            delete_col.emplace_back(i);
          }
        }
      }
    } else {//cstr only check if 2 exists, check if connect together
      for (int i = 1; i < NumCol; ++i) {
        //double check
        if (sum_col[i] == 2) {
          //double check
          if_keep = false;
          for (auto j = node->IdxCols[i] + 1;; ++j) {
            curr_node = ColPool4Pricing[j];
            if (!curr_node) break;
            if (curr_node == ai) {
              if (ColPool4Pricing[j + 1] == aj) {
                if_keep = true;
                break;
              }
            } else if (curr_node == aj) {
              if (ColPool4Pricing[j + 1] == ai) {
                if_keep = true;
                break;
              }
            }
          }
          if (if_keep)delete_col.emplace_back(i);
        }
      }
    }
  }

  //remove cols whose coefficients are 1 but keep the first col
  for (int i = 0; i < NumCol; ++i) sum_col[i] = 0;
  for (auto c : delete_col) sum_col[c] = 1;
  sum_col[0] = 0;//keep the first col

  int len = 0, keep = 1;
  for (int i = 1; i < NumCol; ++i) {
    if (sum_col[i]) solver_ind[len++] = i;
    else node->IdxCols[keep++] = node->IdxCols[i];
  }

  if (len) {
    safe_solver(node->solver.SOLVERdelvars(len, solver_ind))
    safe_solver(node->solver.SOLVERreoptimize())
    safe_solver(node->solver.SOLVERgetNumCol(&NumCol))
  }

  delete[]sum_col;
}

void CVRP::generateRCCsInEnu(BBNODE *const node) {
  int numnz, cnt = 0;
  auto &ptr_col = node->IdxColsInEnuColPool;
  int curr_node;
  set<int> tmp_cutinfo;

  CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

  CMGR_CreateCMgr(&MyCutsCMP, 100);
  CMGR_CreateCMgr(&MyOldCutsCMP, 100);

  getInfoEdge(node, false);
  //Because the separation functor considers n+1 as the depot, we do the corresponding converse
  for (int i = 1; i <= node->NumEdges; ++i) if (!node->EdgeTail[i]) node->EdgeTail[i] = Dim; else break;

  CAPSEP_SeparateCapCuts(RealDim, Demand, Cap, node->NumEdges, node->EdgeTail,
                         node->EdgeHead, node->EdgeVal, MyOldCutsCMP,
                         MAXNOOFCUTS, TOLERANCE, TOLERANCE,
                         &if_IntNFeasible, &MaxVio, MyCutsCMP);
  //Because the separation functor considers n+1 as the depot, we do the corresponding converse
  for (int i = 1; i <= node->NumEdges; ++i) if (node->EdgeTail[i] == Dim) node->EdgeTail[i] = 0; else break;
  unordered_map<pair<int, int>, vector<int>, PairHasher> map_node_lp;
  int old_num_rcc = (int) node->RCCs.size();
  if (!MyCutsCMP->Size) goto QUIT;

  for (int col = 0; col < NumCol; ++col) {
    int past_node = 0;
    for (auto j = node->IdxCols[col] + 1;; ++j) {
      curr_node = ColPool4Pricing[j];
      if (past_node < curr_node) {
        map_node_lp[{past_node, curr_node}].emplace_back(col);
      } else {
        map_node_lp[{curr_node, past_node}].emplace_back(col);
      }
      if (!curr_node) break;
      past_node = curr_node;
    }
  }

  //Read cuts and add them to lp
  for (int i = 0; i < MyCutsCMP->Size; ++i) {
    RCC rcc;
    auto &tmp_customerInfo = rcc.InfoRCCCustomer;
    for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; ++j) {
      tmp_customerInfo.emplace_back(MyCutsCMP->CPL[i]->IntList[j]);
    }
    if (tmp_customerInfo.size() <= Dim / 2) {
      rcc.FormRCC = true;
      rcc.RHS = MyCutsCMP->CPL[i]->RHS;

      if (std::find(node->RCCs.begin(), node->RCCs.end(), rcc) != node->RCCs.end()) {
        continue;
      }

      ++cnt;

      vector<double> col_val(NumCol, 0);
      for (auto iter = tmp_customerInfo.begin(); iter != tmp_customerInfo.end(); ++iter) {
        int ai = *iter;
        auto inner_iter = iter;
        ++inner_iter;
        for (; inner_iter != tmp_customerInfo.end(); ++inner_iter) {
          int aj = *inner_iter;
          if (ai < aj) {
            for (auto it : map_node_lp[{ai, aj}]) {
              ++col_val[it];
            }
          } else {
            for (auto it : map_node_lp[{aj, ai}]) {
              ++col_val[it];
            }
          }
        }
      }
      col_val[0] = rcc.RHS;
      numnz = 0;
      for (int j = 0; j < NumCol; ++j) {
        if (col_val[j] != 0) {
          solver_ind[numnz] = j;
          solver_val[numnz++] = col_val[j];
        }
      }
    } else {
      rcc.FormRCC = false;
      rcc.RHS = RealDim - double(2 * tmp_customerInfo.size()) + MyCutsCMP->CPL[i]->RHS;

      auto &tmp_NoCustomerInfo = rcc.InfoRCCOutsideCustomer;
      vector<bool> tmp_supp(Dim, false);
      for (int j : tmp_customerInfo) {
        tmp_supp[j] = true;
      }
      for (int j = 1; j < Dim; ++j) {
        if (!tmp_supp[j]) {
          tmp_NoCustomerInfo.emplace_back(j);
        }
      }

      if (std::find(node->RCCs.begin(), node->RCCs.end(), rcc) != node->RCCs.end()) {
        continue;
      }

      ++cnt;

      vector<double> col_val(NumCol, 0);
      for (auto iter = tmp_NoCustomerInfo.begin(); iter != tmp_NoCustomerInfo.end(); ++iter) {
        int ai = *iter;
        auto inner_iter = iter;
        ++inner_iter;
        for (; inner_iter != tmp_NoCustomerInfo.end(); ++inner_iter) {
          int aj = *inner_iter;
          if (ai < aj) {
            for (auto it : map_node_lp[{ai, aj}]) {
              ++col_val[it];
            }
          } else {
            for (auto it : map_node_lp[{aj, ai}]) {
              ++col_val[it];
            }
          }
        }
      }
      for (auto customer_it : tmp_NoCustomerInfo) {
        for (auto it : map_node_lp[{0, customer_it}]) col_val[it] += 0.5;
      }
      for (auto customer_it : tmp_customerInfo) {
        for (auto it : map_node_lp[{0, customer_it}]) col_val[it] -= 0.5;
      }
      col_val[0] = rcc.RHS;
      numnz = 0;
      for (int j = 0; j < NumCol; ++j) {
        if (col_val[j] != 0) {
          solver_ind[numnz] = j;
          solver_val[numnz++] = col_val[j];
        }
      }
    }
    rcc.IdxRCC = NumRow;
    node->RCCs.emplace_back(rcc);
    //diff with generateRCCs in terms of without using getCoeffRCC instead we traverse the sequence
    safe_solver(node->solver.SOLVERaddconstr(numnz, solver_ind, solver_val, SOLVER_LESS_EQUAL, rcc.RHS, nullptr))
    safe_solver(node->solver.SOLVERupdatemodel())
    safe_solver(node->solver.SOLVERgetNumRow(&NumRow))
  }

  QUIT:
  for (int i = 0; i < MyCutsCMP->Size; ++i) CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
  MyCutsCMP->Size = 0;

  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  CMGR_FreeMemCMgr(&MyCutsCMP);
}

void CVRP::deleteNonActiveCutsSafely(BBNODE *const node, int old_num) {
  int cnt = 0;
  safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack))
  vector<int> cstr_index(NumRow);
  iota(cstr_index.begin(), cstr_index.end(), 0);
  int delta = 0;
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);
  decltype(deleted_cstrs.begin()) stop_sign;
  for (int i = old_num; i < NumRow; ++i) {
    if (abs(Slack[i]) > TOLERANCE) {
      solver_ind[cnt++] = i;
      cstr_index[i] = -1;
      deleted_cstrs.emplace_back(i);
    }
  }
  if (deleted_cstrs.empty()) {
    goto QUIT;
  }
  std::sort(deleted_cstrs.begin(), deleted_cstrs.end());
  stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))//cannot be deleted

  //start the delete operation and the assignment operation synchronously,
  //but BrC is not needed because we only add cuts here, old cuts are protected!
  //RCC
  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //R1Cs
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }
  //R1Cs_multi
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  QUIT:
  if (If_in_Enu_State) {
    if (node->SizeEnuColPool) {
      chgEnuMatByCuts(node);
    }
  }
}

void CVRP::generateAllR1CsInEnu(BBNODE *node, bool if_use_MIP) {

  int old_num = NumRow;
  vector<pair<vector<int>, double>> cut_info;
  vector<tuple<vector<int>, int, double>> multi_cut_info;
  vector<vector<int>> routes;
  vector<double> frac_routes;
  vector<int> route;
  route.reserve(Dim);

  for (auto &i : node->Idx4LPSolsInColPool) {
    route.clear();
    for (auto j = i.first + 1;; ++j) {
      int curr_node = ColPool4Pricing[j];
      if (!curr_node)break;
      route.emplace_back(curr_node);
    }
    routes.emplace_back(route);
    frac_routes.emplace_back(i.second);
  }

  generateR1C3s(routes, {}, frac_routes, {}, cut_info);

  if (if_use_MIP) {
    vector<int> in_cut_type;
    for (int i = 5; i <= CONFIG::MaxRowRank1; i += 2) {
      in_cut_type.emplace_back(i);
    }
    generateHighDimR1Cs(routes, frac_routes, in_cut_type, cut_info);
  } else {
    findR1C_multi(routes, frac_routes, routes, cut_info, multi_cut_info);
  }

  std::vector<std::tuple<int, std::set<int>, double, int, int>> cut_info_set(cut_info.size());
  for (int i = 0; i < cut_info.size(); ++i) {
    cut_info_set[i] = make_tuple(i, set < int > {}, cut_info[i].second, 0, 0);
  }
  addSelectedR1C_N_multiCuts<vector<pair<vector<int>, double>>, true>(node, cut_info_set, cut_info);
  cut_info_set.resize(multi_cut_info.size());
  for (int i = 0; i < multi_cut_info.size(); ++i) {
    cut_info_set[i] = make_tuple(i, set < int > {}, get<2>(multi_cut_info[i]), 0, 0);
  }
  addSelectedR1C_N_multiCuts<vector<tuple<vector<int>, int, double>>, true>(node, cut_info_set, multi_cut_info);
  //
//  for (auto &cut : cut_info) {
//    //add cuts into lp
//    addR1CInEnu(node, cut.first);
//  }
//
//  for (auto &cut : multi_cut_info) {
//    //add cuts into lp
//    auto new_cut = make_pair(get<0>(cut), get<1>(cut));
//    addR1C_multiInEnu(node, new_cut);
//  }

  safe_Hyperparameter(checkMaxNum_R1Cs((int) node->R1Cs.size()))
  safe_Hyperparameter((int) node->R1Cs_multi.size() > MaxNum_R1C_multi)
  safe_Hyperparameter(NumRow > CST_LIMIT)
}

void CVRP::deleteNewAddedNonActiveCutsBySlackInEnu(BBNODE *node, int oldNum, bool if_delete_rcc) {
  int cnt = 0;

  safe_solver(node->solver.SOLVERgetSlack(0, NumRow, Slack))
  int *cstr_index = new int[NumRow];
  iota(cstr_index, cstr_index + NumRow, 0);
  vector<int> deleted_cstrs;
  deleted_cstrs.reserve(NumRow);

  if (if_delete_rcc) {
    for (auto &cut : node->RCCs) {
      int i = cut.IdxRCC;
      if (i < oldNum)continue;
      if (abs(Slack[i]) > TOLERANCE) {
        solver_ind[cnt++] = i;
        cstr_index[i] = -1;
        deleted_cstrs.emplace_back(i);
      }
    }
  }

  //we don't delete rcc
  for (auto &cut : node->R1Cs) {
    int i = cut.IdxR1C;
    if (i < oldNum)continue;
    if (abs(Slack[i]) > TOLERANCE) {
      solver_ind[cnt++] = i;
      cstr_index[i] = -1;
      deleted_cstrs.emplace_back(i);
    }
  }

  for (auto &cut : node->R1Cs_multi) {
    int i = cut.IdxR1C;
    if (i < oldNum)continue;
    if (abs(Slack[i]) > TOLERANCE) {
      solver_ind[cnt++] = i;
      cstr_index[i] = -1;
      deleted_cstrs.emplace_back(i);
    }
  }

  if (deleted_cstrs.empty()) {
    delete[]cstr_index;
    return;
  }
  std::stable_sort(deleted_cstrs.begin(), deleted_cstrs.end());
  int delta = 0;
  auto stop_sign = deleted_cstrs.end() - 1;
  for (auto i = deleted_cstrs.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;

  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  safe_solver(node->solver.SOLVERreoptimize())
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))

  //rcc must be here!
  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //r1c
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //r1c_multi
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  delete[]cstr_index;
}

void CVRP::addR1C_multiInEnu(BBNODE *node, const pair<vector<int>, int> &cut) {
  size_t num_nz = 0;
  vector<int> ai_col(NumCol, 0);
  size_t numnzP;
  const auto &plan = map_rank1_multiplier[(int) cut.first.size()][cut.second];
  const auto &multi = get<0>(plan);
  int denominator = get<1>(plan);
  int count = 0;
  for (auto i : cut.first) {
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, i - 1, 1))
    for (size_t j = 0; j < numnzP; ++j) {
      ai_col[solver_ind[j]] += multi[count];
    }
    ++count;
  }
  int rhs = get<2>(plan);
  if (rhs > 0.1) {// is not 0
    solver_ind[num_nz] = 0;
    solver_val[num_nz++] = rhs;
  }
  for (int i = 1; i < NumCol; ++i) {
    int tmp = int(ai_col[i] / denominator);
    if (tmp) {
      solver_ind[num_nz] = i;
      solver_val[num_nz++] = tmp;
    }
  }

  safe_solver(node->solver.SOLVERaddconstr((int) num_nz,
                                           solver_ind,
                                           solver_val,
                                           SOLVER_LESS_EQUAL,
                                           rhs,
                                           nullptr))
  R1C_multi r1c;
  r1c.InfoR1C = cut;
  r1c.IdxR1C = NumRow++;
  r1c.RHS = rhs;
  node->R1Cs_multi.emplace_back(r1c);
}

void CVRP::addR1CInEnu(BBNODE *node, const vector<int> &cut) {
  size_t num_nz = 0;
  vector<int> ai_col(NumCol, 0);
  size_t numnzP;
  for (auto i : cut) {
    safe_solver(node->solver.SOLVERXgetconstrs(&numnzP, solver_beg, solver_ind, solver_val, i - 1, 1))
    for (size_t j = 0; j < numnzP; ++j)++ai_col[solver_ind[j]];
  }
  int rhs = int(cut.size() / 2);
  if (rhs > 0.1) {// is not 0
    solver_ind[num_nz] = 0;
    solver_val[num_nz++] = rhs;
  }
  for (int i = 1; i < NumCol; ++i) {
    int tmp = int(ai_col[i] / 2);
    if (tmp) {
      solver_ind[num_nz] = i;
      solver_val[num_nz++] = tmp;
    }
  }

  safe_solver(node->solver.SOLVERaddconstr((int) num_nz,
                                           solver_ind,
                                           solver_val,
                                           SOLVER_LESS_EQUAL,
                                           rhs,
                                           nullptr))
  R1C r1c3;
  r1c3.InfoR1C = cut;
  r1c3.IdxR1C = NumRow++;
  r1c3.RHS = rhs;
  node->R1Cs.emplace_back(r1c3);
}

void CVRP::chgEnuMatByCuts(BBNODE *node) {
  if (!node->SizeEnuColPool) return;

  auto &mat0 = node->MatInEnu.front();
  auto &map = node->map_col_pool;
  int oldNum = 0;
  for (auto &it : node->MatInEnu) oldNum += it.rows();

  sparseRowMatrixXd mat(NumRow - oldNum, node->SizeEnuColPool);
  //only add those whose col is not marked as deleted! they will be totally zero in the matrix
  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(NumRow) * double(node->SizeEnuColPool) * 0.1));
  auto &deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;

  sparseRowMatrixXd sum(1, mat.cols());
  unordered_map<int, double> tmp;
  tmp.reserve(mat.cols());

  //rcc
  for (auto &rcc : node->RCCs) {
    if (rcc.IdxRCC < oldNum) continue;
    tmp.clear();
    if (rcc.FormRCC) {
      auto &info = rcc.InfoRCCCustomer;
      for (auto iter = info.begin(); iter != info.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != info.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
    } else {
      auto &infoRccCustomer = rcc.InfoRCCCustomer;
      auto &infoRccOutsideCustomer = rcc.InfoRCCOutsideCustomer;
      for (auto iter = infoRccOutsideCustomer.begin(); iter != infoRccOutsideCustomer.end(); ++iter) {
        auto inn_iter = iter;
        ++inn_iter;
        int ai = *iter;
        for (; inn_iter != infoRccOutsideCustomer.end(); ++inn_iter) {
          int aj = *inn_iter;
          if (ai < aj) {
            for (auto it : map[{ai, aj}]) {
              ++tmp[it];
            }
          } else {
            for (auto it : map[{aj, ai}]) {
              ++tmp[it];
            }
          }
        }
      }
      for (auto customer_it : infoRccOutsideCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] += 0.5;
      }
      for (auto customer_it : infoRccCustomer) {
        for (auto it : map[{0, customer_it}]) tmp[it] -= 0.5;
      }
    }
    int row = rcc.IdxRCC - oldNum;
    for (auto &it : tmp) {
      if (abs(it.second) > TOLERANCE && !deleted_ColsInEnuPool[it.first]) {
        triplets.emplace_back(row, it.first, it.second);
      }
    }
  }

  for (auto &r1c : node->R1Cs) {
    if (r1c.IdxR1C < oldNum)continue;
    sum.setZero();
    auto &info = r1c.InfoR1C;
    for (auto j : info) {
      sum += mat0.row(j - 1);
    }
    sum /= 2;
    int row = r1c.IdxR1C - oldNum;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_ && !deleted_ColsInEnuPool[it.col()]) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }

  for (auto &r1c : node->R1Cs_multi) {
    if (r1c.IdxR1C < oldNum)continue;
    sum.setZero();
    auto &info = r1c.InfoR1C;
    const auto &plan = map_rank1_multiplier[(int) info.first.size()][info.second];
    const auto &multi = get<0>(plan);
    int denominator = get<1>(plan);
    int count = 0;
    for (auto &j : info.first) {
      sum += mat0.row(j - 1) * multi[count++];
    }
    sum /= denominator;
    int row = r1c.IdxR1C - oldNum;
    for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
      int val_ = int(it.value() + TOLERANCE);
      if (val_ && !deleted_ColsInEnuPool[it.col()]) {
        triplets.emplace_back(row, it.col(), val_);
      }
    }
  }
  mat.setFromTriplets(triplets.begin(), triplets.end());
  node->MatInEnu.push_back(std::move(mat));
}


//void CVRP::chgEnuMatByCuts(BBNODE *node) {
//  auto &mat = node->MatInEnu;
//  if (!mat.cols()) return;
//
//  auto &LPRow2MatRow = node->LPRow2MatRow;
//  auto &MatRow2LPRow = node->MatRow2LPRow;
//  auto &MaxUsedRowInMat = node->MaxUsedRowInMat;
//  auto &map = node->map_col_pool;
//  int oldNum = (int) MatRow2LPRow.size();
//  if (NumRow >= mat.rows()) {
//    mat.conservativeResize(2 * NumRow, mat.cols());
//  }
//  LPRow2MatRow.resize(NumRow);
//  vector<tuple<int, int, int>> vec(NumRow - oldNum);
//  int cnt = 0;
//  for (int i = 0; i < node->RCCs.size(); ++i) {
//    if (node->RCCs[i].IdxRCC < oldNum)continue;
//    vec[cnt++] = make_tuple(node->RCCs[i].IdxRCC, i, 0);
//  }
//  for (int i = 0; i < node->R1Cs.size(); ++i) {
//    if (node->R1Cs[i].IdxR1C < oldNum)continue;
//    vec[cnt++] = make_tuple(node->R1Cs[i].IdxR1C, i, 1);
//  }
//  for (int i = 0; i < node->R1Cs_multi.size(); ++i) {
//    if (node->R1Cs_multi[i].IdxR1C < oldNum)continue;
//    vec[cnt++] = make_tuple(node->R1Cs_multi[i].IdxR1C, i, 2);
//  }
//
//  std::sort(vec.begin(), vec.end(), [](const tuple<int, int, int> &a, const tuple<int, int, int> &b) {
//    return get<0>(a) < get<0>(b);
//  });
//
//  std::vector<Eigen::Triplet<double>> triplets;
//  triplets.reserve(1000000);
//
//  sparseRowMatrixXd sum(1, mat.cols());
//  unordered_map<int, double> tmp;
//  tmp.reserve(mat.cols());
//
//  auto beg = high_resolution_clock::now();
//  for (auto &i : vec) {
//    LPRow2MatRow[get<0>(i)] = MaxUsedRowInMat;
//    MatRow2LPRow[MaxUsedRowInMat] = get<0>(i);
//    if (get<2>(i) == 0) {
//      tmp.clear();
//      auto &rcc = node->RCCs[get<1>(i)];
//      if (rcc.FormRCC) {
//        auto &info = rcc.InfoRCCCustomer;
//        for (auto iter = info.begin(); iter != info.end(); ++iter) {
//          auto inn_iter = iter;
//          ++inn_iter;
//          int ai = *iter;
//          for (; inn_iter != info.end(); ++inn_iter) {
//            int aj = *inn_iter;
//            if (ai < aj) {
//              for (auto it : map[{ai, aj}]) {
//                ++tmp[it];
//              }
//            } else {
//              for (auto it : map[{aj, ai}]) {
//                ++tmp[it];
//              }
//            }
//          }
//        }
//      } else {
//        auto &infoRccCustomer = rcc.InfoRCCCustomer;
//        auto &infoRccOutsideCustomer = rcc.InfoRCCOutsideCustomer;
//        for (auto iter = infoRccOutsideCustomer.begin(); iter != infoRccOutsideCustomer.end(); ++iter) {
//          auto inn_iter = iter;
//          ++inn_iter;
//          int ai = *iter;
//          for (; inn_iter != infoRccOutsideCustomer.end(); ++inn_iter) {
//            int aj = *inn_iter;
//            if (ai < aj) {
//              for (auto it : map[{ai, aj}]) {
//                ++tmp[it];
//              }
//            } else {
//              for (auto it : map[{aj, ai}]) {
//                ++tmp[it];
//              }
//            }
//          }
//        }
//        for (auto customer_it : infoRccOutsideCustomer) {
//          for (auto it : map[{0, customer_it}]) tmp[it] += 0.5;
//        }
//        for (auto customer_it : infoRccCustomer) {
//          for (auto it : map[{0, customer_it}]) tmp[it] -= 0.5;
//        }
//      }
//      vector<pair<int, double>> vec_tmp(tmp.begin(), tmp.end());
//      std::sort(vec_tmp.begin(), vec_tmp.end(), [](const pair<int, double> &a, const pair<int, double> &b) {
//        return a.first < b.first;
//      });
//      for (auto &it : vec_tmp) {
//        if (abs(it.second) > TOLERANCE) {
//          triplets.emplace_back(MaxUsedRowInMat, it.first, it.second);
////          mat.insert(MaxUsedRowInMat, it.first) = it.second;
//        }
//      }
//    } else {
//      sum.setZero();
//      if (get<2>(i) == 1) {
//        auto &info = node->R1Cs[get<1>(i)].InfoR1C;
//        for (auto j : info) {
//          sum += mat.row(j - 1);
//        }
//        sum /= 2;
//      } else {
//        auto &info = node->R1Cs_multi[get<1>(i)].InfoR1C;
//        const auto &plan = map_rank1_multiplier[(int) info.first.size()][info.second];
//        const auto &multi = get<0>(plan);
//        int denominator = get<1>(plan);
//        int count = 0;
//        for (auto &j : info.first) {
//          sum += mat.row(j - 1) * multi[count++];
//        }
//        sum /= denominator;
//      }
//      for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
//        int val_ = int(it.value() + TOLERANCE);
//        if (val_) {
////          mat.insert(MaxUsedRowInMat, it.col()) = val_;
//          triplets.emplace_back(MaxUsedRowInMat, it.col(), val_);
//        }
//      }
//    }
//    ++MaxUsedRowInMat;
//  }
//  mat.setFromTriplets(triplets.begin(), triplets.end());
//}

void CVRP::refineArcsHeur(BBNODE *node) const {
  unordered_set < pair<int, int>, PairHasher > arcs;
  auto Deleted_ColsInEnuPool = node->Deleted_ColsInEnuPool;
  arcs.reserve(Dim * Dim / 2);
  // get lp solution here
  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
  for (int i = 1; i < NumCol; ++i) {
//    if (X[i] > TOLERANCE) {
    // traverse the sol
    int past_node = 0;
    for (auto j = node->IdxCols[i] + 1;; ++j) {
      int curr_node = ColPool4Pricing[j];
      if (past_node < curr_node) arcs.emplace(past_node, curr_node);
      else arcs.emplace(curr_node, past_node);//one side
      if (!curr_node) break;
      past_node = curr_node;
    }
//    }
  }
  //examine the enumeration col pool
  for (auto &arc_pr : node->map_col_pool) {
    if (arcs.find(arc_pr.first) == arcs.end()) {
      for (auto &it : arc_pr.second) {
        Deleted_ColsInEnuPool[it] = true;
      }
    }
  }
}

//void CVRP::changeModel4DeletingSlackDual_zeroCuts(BBNODE *node) {
//  if (!If_in_Enu_State) throw runtime_error("changeModel4DeletingSlackDual_zeroCuts not in enumeration state");
//  safe_solver(node->solver.SOLVERchgobj(0, 1, &LPVal))
//  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//  cout << "old_dual= " << endl;
//  for (auto &i : target_cuts) {
//    cout << Pi[i] << ' ';
//  }
//  cout << endl;
//  vector<double> old_rhs(NumRow);
//  safe_solver(node->solver.SOLVERgetRHS(0, NumRow, old_rhs.data()))
//  vector<double> rhs(NumRow, 0);
//  for (auto &cut : target_cuts) {
//    rhs[cut] = -1;
//  }
//  safe_solver(node->solver.SOLVERchgrhs(0, NumRow, rhs.data()))
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_UB, 0, 0);
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_LB, 0, -GRB_INFINITY);
//  safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_DUAL_SIMPLEX))
//  int status;
//  bool if_new_dual = true;
//  safe_solver(node->solver.SOLVERgetStatus(&status))
//  if (status == SOLVER_UNBOUNDED || status == SOLVER_INF_OR_UNBD) {
//    if_new_dual = false;
//  }
//  if (if_new_dual) {
//    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//    cout << "new dual" << endl;
//    for (auto &i : target_cuts) {
//      cout << Pi[i] << ' ';
//    }
//    cout << endl;
//  }
//  safe_solver(node->solver.SOLVERchgobj(0, 1, &UB))
//  safe_solver(node->solver.SOLVERchgrhs(0, NumRow, old_rhs.data()))
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_LB, 0, 0);
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_UB, 0, GRB_INFINITY);
//  safe_solver(node->solver.SOLVERreoptimize())
//  if (!if_new_dual) {
//    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//  }
//  safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_PRIMAL_SIMPLEX))
//}

void CVRP::findNonactiveCuts(BBNODE *node) {
  safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
  vector<int> nonactive_cuts;
  nonactive_cuts.reserve(NumRow);
  for (auto &rcc : node->RCCs) {
    int idx = rcc.IdxRCC;
    if (abs(Pi[idx]) < TOLERANCE) {
      nonactive_cuts.emplace_back(idx);
    }
  }
  for (auto &r1c : node->R1Cs) {
    int idx = r1c.IdxR1C;
    if (abs(Pi[idx]) < TOLERANCE) {
      nonactive_cuts.emplace_back(idx);
    }
  }
  for (auto &r1c_mul : node->R1Cs_multi) {
    int idx = r1c_mul.IdxR1C;
    if (abs(Pi[idx]) < TOLERANCE) {
      nonactive_cuts.emplace_back(idx);
    }
  }
  deleteNonactiveCuts(node, nonactive_cuts);

//  safe_solver(node->solver.SOLVERgetObjVal(&LPVal))
//  safe_solver(node->solver.SOLVERchgobj(0, 1, &LPVal))
//
//  vector<double> old_rhs(NumRow);
//  safe_solver(node->solver.SOLVERgetRHS(0, NumRow, old_rhs.data()))
//  vector<double> rhs(NumRow, 0);
//  for (auto &cut : nonactive_cuts) {
//    rhs[cut] = -1;
//  }
//  safe_solver(node->solver.SOLVERchgrhs(0, NumRow, rhs.data()))
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_UB, 0, 0);
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_LB, 0, -GRB_INFINITY);
//
//  int env_method;
//  bool if_changed = false;
//  safe_solver(node->solver.SOLVERgetenvMethod(&env_method))
//  if (env_method != SOLVER_DUAL_SIMPLEX) {
//    safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_DUAL_SIMPLEX))
//    if_changed = true;
//  }
//  safe_solver(node->solver.SOLVERreoptimize())
//  int status;
//  bool if_new_dual = true;
//  safe_solver(node->solver.SOLVERgetStatus(&status))
//  if (status == SOLVER_UNBOUNDED || status == SOLVER_INF_OR_UNBD) {
//    if_new_dual = false;
//  }
//  if (if_new_dual) {
//    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//    for (auto beg = nonactive_cuts.begin(); beg != nonactive_cuts.end();) {
//      if (abs(Pi[*beg]) > TOLERANCE) {
//        beg = nonactive_cuts.erase(beg);
//        cout << "?????" << endl;
//      } else {
//        beg++;
//      }
//    }
//  }
//  safe_solver(node->solver.SOLVERchgobj(0, 1, &Obj4FirstCol))
//  safe_solver(node->solver.SOLVERchgrhs(0, NumRow, old_rhs.data()))
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_LB, 0, 0);
//  GRBsetdblattrelement(node->solver.model, GRB_DBL_ATTR_UB, 0, GRB_INFINITY);
//  safe_solver(node->solver.SOLVERreoptimize())
//  if (!if_new_dual) {
//    safe_solver(node->solver.SOLVERgetDual(0, NumRow, Pi))
//    nonactive_cuts.clear();
//  }
//  deleteNonactiveCuts(node, nonactive_cuts);
//  auto end = high_resolution_clock::now();
//  cout << "findNonactiveCuts time= " << duration<double>(end - beg).count() << endl;
//  if (if_changed) {
//    safe_solver(node->solver.SOLVERsetenvMethod(env_method))
//  }
}

void CVRP::deleteNonactiveCuts(BBNODE *node, std::vector<int> &nonactive_cuts) {
  if (nonactive_cuts.empty()) return;
  int old_num_row = NumRow;
  sort(nonactive_cuts.begin(), nonactive_cuts.end());
  vector<int> cstr_index(NumRow);
  iota(cstr_index.begin(), cstr_index.end(), 0);
  int cnt = 0;
  for (auto &cut : nonactive_cuts) {
    solver_ind[cnt++] = cut;
    cstr_index[cut] = -1;
  }
  int delta = 0;
  auto stop_sign = nonactive_cuts.end() - 1;
  for (auto i = nonactive_cuts.begin(); i < stop_sign; ++i) {
    ++delta;
    for (int j = *i + 1; j < *(i + 1); ++j) cstr_index[j] = j - delta;
  }
  ++delta;
  for (int j = *stop_sign + 1; j < NumRow; ++j) cstr_index[j] = j - delta;
  safe_solver(node->solver.SOLVERdelconstrs(cnt, solver_ind))
  int env_method;
  bool if_changed = false;
  safe_solver(node->solver.SOLVERgetenvMethod(&env_method))
  if (env_method != SOLVER_DUAL_SIMPLEX) {
    safe_solver(node->solver.SOLVERsetenvMethod(SOLVER_DUAL_SIMPLEX))
    if_changed = true;
  }
  safe_solver(node->solver.SOLVERreoptimize())
  if (if_changed) {
    safe_solver(node->solver.SOLVERsetenvMethod(env_method))
  }
  safe_solver(node->solver.SOLVERgetNumRow(&NumRow))

  //start the delete operation and the assignment operation synchronously,
  //but BrC does not need it, BrC only needs to perform the reassignment operation.
  //RCCs are the first

  for (auto i = node->RCCs.begin(); i < node->RCCs.end();) {
    if (cstr_index[i->IdxRCC] == -1) {
      i = node->RCCs.erase(i);
    } else {
      i->IdxRCC = cstr_index[i->IdxRCC];
      ++i;
    }
  }

  //R1C3s next
  for (auto i = node->R1Cs.begin(); i < node->R1Cs.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  //R1C3s_multi next
  for (auto i = node->R1Cs_multi.begin(); i < node->R1Cs_multi.end();) {
    if (cstr_index[i->IdxR1C] == -1) {
      i = node->R1Cs_multi.erase(i);
    } else {
      i->IdxR1C = cstr_index[i->IdxR1C];
      ++i;
    }
  }

  if (If_in_Enu_State) {
    CstrIndex = std::move(cstr_index);
  } else {
    for (auto &i : node->BrCs) i.IdxBrC = cstr_index[i.IdxBrC];
  }
//  if (If_in_Enu_State) {
//#ifdef NominalBranchingInEnu
//    for (auto &i : node->NBrCs) {
//      if (i.IdxNBrC != -2)i.IdxNBrC = cstr_index[i.IdxNBrC];
//    }
//#endif
//    regenerateEnuMat(node, nullptr, false);
//  } else {
//    //BrC finally, BrC does not need to delete, only need to reassign
//
//  }
}
