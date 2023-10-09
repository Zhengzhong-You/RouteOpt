//
// Created by Zhengzhong You on 7/14/22.
//

#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace Eigen;
using namespace chrono;

void CVRP::addCols(BBNODE *node, int &ccnt) {

  int ccnt_cnt;
  int curr_node, past_node;
  double cost_sum;

  vector<int> r1c_eff(node->R1Cs.size(), 0);
  vector<int> r1c_multi_eff(node->R1Cs_multi.size(), 0);
  vector<int> r1c_multi_state(node->R1Cs_multi.size(), 0);
  unordered_map<int, vector<int>> map_node_lp;
  map_node_lp.reserve(Dim * Dim);
  R1CMem Rank1CutMem = 0;
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(NumRow, ccnt);
  Eigen::RowVectorXd cost(ccnt);

  ccnt_cnt = 0;
  past_node = 0;
  cost_sum = 0;
  for (size_t start = PriorPoolBeg4Pricing + 1; start < PoolBeg4Pricing; ++start) {
	curr_node = ColPool4Pricing[start];
	cost_sum += CostMat4Vertex[past_node][curr_node];
	//rcc & br
	map_node_lp[past_node * Dim + curr_node].emplace_back(ccnt_cnt);

//    int cnt = 0;
//    bool if_open = false;
//    if (ccnt_cnt == 0 && NumRow == 470) {
//      for (auto &r1c : node->R1Cs) {
//        if (r1c.IdxR1C == 310) {
//          for (auto &i : r1c.InfoR1C) {
//            cout << i << " ";
//          }
//          cout << endl;
//          break;
//        }
//        ++cnt;
//      }
//      cout << "cnt= " << cnt << endl;
//      cout << "l= ";
//      if_open = true;
//    }

	//r1c

	for (auto l : Vertex2AllInOneLPR1Cs[curr_node].first) {
	  if (Rank1CutMem[l]) {
		Rank1CutMem[l] = false;
		++r1c_eff[l];
	  } else Rank1CutMem[l] = true;
	}
//    if (if_open) {
//      for (auto l : Vertex2AllInOneLPR1Cs[curr_node].first) {
//        cout << l << " ";
//      }
//      cout << endl;
//    }
//    if (if_open) {
//      cout << "curr_node: " << curr_node << " r1c(ccnt)= " << Rank1CutMem[cnt] << " r1c_eff= " << r1c_eff[cnt] << endl;
//    }
	Rank1CutMem &= Vertex2AllInOneLPR1Cs[curr_node].second;
//    if (if_open) {
//      cout << "curr_node: " << curr_node << " r1c(ccnt)= " << Rank1CutMem[cnt] << " r1c_eff= " << r1c_eff[cnt] << endl;
//    }

	//r1c_multi
	for (auto &l : Vertex2AllInOneLPR1C_multi[curr_node].first) {
	  int tmp_cut = l.first;
	  r1c_multi_state[tmp_cut] += l.second;
	  if (r1c_multi_state[tmp_cut] >= R1C_multi_denominator_InLP[tmp_cut]) {
		r1c_multi_state[tmp_cut] -= R1C_multi_denominator_InLP[tmp_cut];
		++r1c_multi_eff[tmp_cut];
	  }
	}

	for (auto l : Vertex2AllInOneLPR1C_multi[curr_node].second) r1c_multi_state[l] = 0;

	if (!curr_node) {
	  cost(ccnt_cnt) = cost_sum;
	  //vehicle
	  mat(RealDim, ccnt_cnt) = 1;
	  //r1c
	  for (int j = 0; j < r1c_eff.size(); ++j) {
		if (r1c_eff[j]) {
		  mat(node->R1Cs[j].IdxR1C, ccnt_cnt) = r1c_eff[j];
		}
	  }
	  //r1c_multi
	  for (int j = 0; j < r1c_multi_eff.size(); ++j) {
		if (r1c_multi_eff[j]) {
		  mat(node->R1Cs_multi[j].IdxR1C, ccnt_cnt) = r1c_multi_eff[j];
		}
	  }
	  ++start;
	  cost_sum = 0;
	  memset(r1c_eff.data(), 0, sizeof(int) * r1c_eff.size());
	  memset(r1c_multi_eff.data(), 0, sizeof(int) * r1c_multi_eff.size());
	  memset(r1c_multi_state.data(), 0, sizeof(int) * r1c_multi_state.size());
	  Rank1CutMem = 0;
	  ++ccnt_cnt;
	} else {
	  ++mat(curr_node - 1, ccnt_cnt);
	}
	past_node = curr_node;
  }

  //rcc
  for (auto &rcc : node->RCCs) {
	if (rcc.FormRCC) {
	  auto &info = rcc.InfoRCCCustomer;
	  int idx = rcc.IdxRCC;
	  for (auto it = info.begin(); it != info.end(); ++it) {
		int ai = *it;
		auto it_inner = it;
		++it_inner;
		for (; it_inner != info.end(); ++it_inner) {
		  int aj = *it_inner;
		  for (auto it_map : map_node_lp[ai * Dim + aj]) ++mat(idx, it_map);
		  for (auto it_map : map_node_lp[aj * Dim + ai]) ++mat(idx, it_map);
		}
	  }
	} else {
	  auto &customer_info = rcc.InfoRCCCustomer;
	  auto &outside_customer_info = rcc.InfoRCCOutsideCustomer;
	  int idx = rcc.IdxRCC;
	  for (auto it = outside_customer_info.begin(); it != outside_customer_info.end(); ++it) {
		int ai = *it;
		auto it_inner = it;
		++it_inner;
		for (; it_inner != outside_customer_info.end(); ++it_inner) {
		  int aj = *it_inner;
		  for (auto it_map : map_node_lp[ai * Dim + aj]) ++mat(idx, it_map);
		  for (auto it_map : map_node_lp[aj * Dim + ai]) ++mat(idx, it_map);
		}
	  }
	  for (int aj : outside_customer_info) {
		for (auto it_map : map_node_lp[aj]) mat(idx, it_map) += 0.5;
		for (auto it_map : map_node_lp[aj * Dim]) mat(idx, it_map) += 0.5;
	  }
	  for (int aj : customer_info) {
		for (auto it_map : map_node_lp[aj]) mat(idx, it_map) -= 0.5;
		for (auto it_map : map_node_lp[aj * Dim]) mat(idx, it_map) -= 0.5;
	  }
	}
  }

  //br
  for (auto &br : node->BrCs) {
	int ai = br.Edge.first;
	int aj = br.Edge.second;
	int idx = br.IdxBrC;
	for (auto it_map : map_node_lp[ai * Dim + aj]) ++mat(idx, it_map);
	for (auto it_map : map_node_lp[aj * Dim + ai]) ++mat(idx, it_map);
  }

  if (ccnt != ccnt_cnt) {
	cerr << "Wrong in exactLabeling 1122" << endl;
	exit(0);
  }

  Map<RowVectorXd> pi(Pi4Labeling.data(), NumRow);
  RowVectorXd rc = cost - pi * mat;

  ccnt_cnt = 0;
  size_t nzcnt = 0;
  for (int i = 0; i < ccnt; ++i) {
	if (rc(i) < RC_TOLERANCE) {
	  solver_beg[ccnt_cnt] = nzcnt;
	  solver_obj[ccnt_cnt] = cost(i);
	  ++ccnt_cnt;
	  for (int row = 0; row < NumRow; ++row) {
		if (mat(row, i) != 0) {
		  solver_ind[nzcnt] = row;
		  solver_val[nzcnt++] = mat(row, i);
		}
	  }
	} else {
	  cout << "col " << i << " is not allowed! The rc= " << rc(i) << endl;
	}
  }

  solver_beg[ccnt_cnt] = nzcnt;

  ccnt = ccnt_cnt;
  if (!ccnt) return;

  safe_solver(node->solver.SOLVERXaddvars(ccnt_cnt,
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
  writeLP_N_tell_if_LP_corrected(node->solver);
#endif
#ifdef checkR1CTotally
  //  checkR1Ctotally(node);
#endif
}

void CVRP::writeColsInPricingPool(BBNODE *const node, int &index) {
  LABEL *p;
  if (checkPricingPool()) reallocatePricingPool();
  auto seq = new int[int(MaxMainResource) + 3];
  int cnt;
  for (auto &i : NegativeRCLabelTuple) {
	auto &ki = get<0>(i);
	auto &kj = get<1>(i);
	node->IdxCols[++index] = PoolBeg4Pricing;
	cnt = 0;
	p = ki;
	while (p) {
	  seq[cnt++] = p->EndVertex;
	  p = p->PLabel;
	}
	for (int k = 0; k < cnt; ++k) {
	  ColPool4Pricing[PoolBeg4Pricing++] = seq[cnt - 1 - k];
	}
	if (kj) {
	  p = kj;
	  while (p) {
		ColPool4Pricing[PoolBeg4Pricing++] = p->EndVertex;
		p = p->PLabel;
	  }
	} else {
	  ColPool4Pricing[PoolBeg4Pricing++] = 0;
	}
  }
  delete[]seq;
}

void CVRP::initializeLabels(BBNODE *const node,
							int mode,
							bool if_resetLabelPoint,
							tuple<bool, int, bool> control_cleanAllPtr) {

  if (if_resetLabelPoint) {
	//if you want to insert artificial variables, you need to pay attention to this
#ifdef SYMMETRY_PROHIBIT
	IdxGlo = 2 * Dim - 1;
#else
	IdxGlo = Dim;
#endif
	SeqBeg = SeqSizeArtiVars;
	rc_std = RC_TOLERANCE;
	NumDominanceChecks = 0;
	NegativeRCLabelTuple.clear();
	Map4EachNegativeRCRoute.clear();
  }

  if (get<0>(control_cleanAllPtr)) {
	// clean AllPtr pool
	cleanAllPtr(node, get<1>(control_cleanAllPtr), get<2>(control_cleanAllPtr));
  }

  if (mode == 1) {
	unordered_set<int> depot_set;
	for (auto j : node->AllForwardBuckets[0][0].BucketArcs) depot_set.emplace(j);
	for (int i = 1; i < Dim; ++i) {
	  if (depot_set.find(i) == depot_set.end()) continue;
	  AllLabel[i].RC = ChgCostMat4Vertex[0][i];
	  AllLabel[i].if_extended = false;
	  //update r1c
	  AllLabel[i].Rank1CutMem = 0;
	  AllLabel[i].numValidRank1Cut = 0;
	  for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[i])) {
		AllLabel[i].Rank1CutMem.set(l);
		AllLabel[i].validRank1Cut[AllLabel[i].numValidRank1Cut++] = l;
	  }
	  //update r1c_multi
	  memset(AllLabel[i].Rank1CutMem_multi, 0, sizeof(int) * NumValidR1C_multi_InCG);
	  AllLabel[i].numValidRank1Cut_multi = 0;
	  for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[i])) {
		int tmp_cut = get<0>(l);
		AllLabel[i].Rank1CutMem_multi[tmp_cut] = get<1>(l);
		AllLabel[i].validRank1Cut_multi[AllLabel[i].numValidRank1Cut_multi++] = tmp_cut;
	  }
	  //reset initial reduced cost
	  int bin = int(AllLabel[i].Sum_MainResource / StepSize);
	  auto &bucket = LabelArrayInForwardSense[i][bin];
	  bucket.first[bucket.second++] = AllLabel + i;
	  auto &bucket2 = IfExistExtraLabelsInForwardSense[i][bin];
	  bucket2.first[bucket2.second++] = AllLabel + i;
	}
  }
#ifdef SYMMETRY_PROHIBIT
  else if (mode == 2) {
	 unordered_set<int> depot_set;
	for (auto j : node->AllBackwardBuckets[0][0].BucketArcs) depot_set.emplace(j);
	int max_num = 2 * Dim - 1;
	for (int i = Dim; i < max_num; ++i) {
	  int point = i - Dim + 1;
	   if (depot_set.find(point) == depot_set.end()) continue;
	  AllLabel[i].RC = ChgCostMat4Vertex[0][point];
	  AllLabel[i].if_extended = false;
	  //update r1c
	  AllLabel[i].Rank1CutMem = 0;
	  AllLabel[i].numValidRank1Cut = 0;
	  for (auto l : get<0>(Vertex2ActiveInOnePricingR1Cs[point])) {
		AllLabel[i].Rank1CutMem.set(l);
		AllLabel[i].validRank1Cut[AllLabel[i].numValidRank1Cut++] = l;
	  }
	  //update r1c_multi
	  memset(AllLabel[i].Rank1CutMem_multi, 0, sizeof(int) * NumValidR1C_multi_InCG);
	  AllLabel[i].numValidRank1Cut_multi = 0;
	  for (auto &l : get<0>(Vertex2ActiveInOnePricingR1C_multi[point])) {
		int tmp_cut = get<0>(l);
		AllLabel[i].Rank1CutMem_multi[tmp_cut] = get<1>(l);
		AllLabel[i].validRank1Cut_multi[AllLabel[i].numValidRank1Cut_multi++] = tmp_cut;
	  }
	  int bin = int(AllLabel[i].Sum_MainResource / StepSize);
	  auto &bucket = LabelArrayInBackwardSense[point][bin];
	  bucket.first[bucket.second++] = AllLabel + i;
	  auto &bucket2 = IfExistExtraLabelsInBackwardSense[point][bin];
	  bucket2.first[bucket2.second++] = AllLabel + i;
	}
  }
#endif
}

void CVRP::addPathByRC(double path_rc, LABEL *ki, LABEL *kj, int num) {
  if (path_rc < rc_std) {
	auto p = ki;
	yzzLong tmp_seq = 0;
	while (p) {
	  tmp_seq.set(p->EndVertex);
	  p = p->PLabel;
	}
	p = kj;
	while (p) {
	  tmp_seq.set(p->EndVertex);
	  p = p->PLabel;
	}
	if (Map4EachNegativeRCRoute.find(tmp_seq) == Map4EachNegativeRCRoute.end()) {
	  Map4EachNegativeRCRoute[tmp_seq] = {ki, kj, path_rc};
	  auto it = std::lower_bound(NegativeRCLabelTuple.begin(), NegativeRCLabelTuple.end(), path_rc,
								 [](const std::tuple<LABEL *, LABEL *, double> &a, double b) {
								   return get<2>(a) < b;
								 });
	  NegativeRCLabelTuple.insert(it, {ki, kj, path_rc});
	  if (NegativeRCLabelTuple.size() > num) NegativeRCLabelTuple.resize(num);
	  rc_std = get<2>(NegativeRCLabelTuple.back());
	} else {
	  auto &tmp_route = Map4EachNegativeRCRoute[tmp_seq];
	  if (path_rc < get<2>(tmp_route)) {
		auto tmp_old = std::find(NegativeRCLabelTuple.begin(), NegativeRCLabelTuple.end(), tmp_route);
		if (tmp_old != NegativeRCLabelTuple.end()) {
		  auto &_old = *tmp_old;
		  get<0>(_old) = ki;
		  get<1>(_old) = kj;
		  get<2>(_old) = path_rc;
		  std::sort(NegativeRCLabelTuple.begin(), NegativeRCLabelTuple.end(),
					[](const std::tuple<LABEL *, LABEL *, double> &a, const std::tuple<LABEL *, LABEL *, double> &b) {
					  return get<2>(a) < get<2>(b);
					});
		} else {
		  auto it = std::lower_bound(NegativeRCLabelTuple.begin(), NegativeRCLabelTuple.end(), path_rc,
									 [](const std::tuple<LABEL *, LABEL *, double> &a, double b) {
									   return get<2>(a) < b;
									 });
		  NegativeRCLabelTuple.insert(it, {ki, kj, path_rc});
		  NegativeRCLabelTuple.pop_back();
		}
		get<0>(tmp_route) = ki;
		get<1>(tmp_route) = kj;
		get<2>(tmp_route) = path_rc;
		//sort
		rc_std = get<2>(NegativeRCLabelTuple.back());
	  }
	}
  }
}

