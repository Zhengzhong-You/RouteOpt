
#include "CVRP.hpp"
#include "templateFunctors.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::eliminateArcs(BbNode *const node) {
  determineIfArcElimination(node);
  if (!final_decision_4_arc_elimination) {
    if_arc_elimination_succeed = false;
    return;
  }
  final_decision_4_arc_elimination = false;
  if_arc_elimination_tried_but_failed = false;

  double time_labeling, time_elimination, time_obtainjump;

  cout << BIG_PHASE_SEPARATION;
  cout << "run ArcElimination..." << endl;

  populateTellWhichBin4ArcElimination<true>();
#ifdef SYMMETRY_PROHIBIT
  populateTellWhichBin4ArcElimination<false>();
#endif

  PtrAllR1CS ptrAllR1Cs(node, this);

  auto beg = chrono::high_resolution_clock::now();
  auto end = chrono::high_resolution_clock::now();
  auto eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  rollback = 0;

  runLabelingForArcElimination(node, ptrAllR1Cs);

  if (rollback == 2) {
    cout << "rollback: " << rollback << endl;
    if_arc_elimination_succeed = false;
    if_arc_elimination_tried_but_failed = true;
    goto QUIT;
  } else if (rollback == 1) {
    cout << "rollback: " << rollback << endl;
    if_arc_elimination_succeed = false;
    if_arc_elimination_tried_but_failed = true;
    cout << "Reach the hard limit!" << endl;
    goto QUIT;
  }
  safe_Hyperparameter(rollback)

  end = chrono::high_resolution_clock::now();
  eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  time_labeling = (double) eps.count() * 1e-3;

  beg = chrono::high_resolution_clock::now();
  eliminateBucketArcs(node, ptrAllR1Cs);

  end = chrono::high_resolution_clock::now();
  eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  time_elimination = (double) eps.count() * 1e-3;

  beg = chrono::high_resolution_clock::now();
  obtainJumpArcs(node);
  end = chrono::high_resolution_clock::now();
  eps = chrono::duration_cast<chrono::milliseconds>(end - beg);

  time_obtainjump = (double) eps.count() * 1e-3;

  cout << "time summary= " << endl;
  cout << "runLabelingForArcElimination= " << time_labeling << endl;
  cout << "eliminateBucketArcs= " << time_elimination << endl;
  cout << "obtainJumpArcs= " << time_obtainjump << endl;
  if_arc_elimination_succeed = true;
  QUIT:
  cout << BIG_PHASE_SEPARATION;
}

void CVRP::runLabelingForArcElimination(BbNode *const node, const PtrAllR1CS &ptrAllR1Cs) {
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  double dif;
  bool if_break, if_find, if_keep;
  double path_rc;
  int edgemap;
  opt_gap = calculateOptimalGap(node);

  auto time_beg = chrono::high_resolution_clock::now();
  auto time_end = time_beg;
  double time_eps;

#ifdef SYMMETRY_PROHIBIT
  concatenatePhaseInArcElimination<true, false>(r1c_to_pi, r1c_multi_to_pi);
  if (rollback == 1 || rollback == 2) goto QUIT;
  runLabeling<true, true, false>(node, ptrAllR1Cs);
  if (rollback == 1 || rollback == 2) goto QUIT;
  concatenatePhaseInArcElimination<false, false>(r1c_to_pi, r1c_multi_to_pi);
  if (rollback == 1 || rollback == 2) goto QUIT;
  runLabeling<false, true, false>(node, ptrAllR1Cs);
  if (rollback == 1 || rollback == 2) goto QUIT;
#else
  concatenatePhaseInArcElimination<true, true>(r1c_to_pi, r1c_multi_to_pi);
  if (rollback == 1 || rollback == 2) goto QUIT;
  runLabeling<true, true, true>(node, ptrAllR1Cs);
  if (rollback == 1 || rollback == 2) goto QUIT;
#endif

  time_end = chrono::high_resolution_clock::now();
  time_eps = chrono::duration<double>(time_end - time_beg).count();
  cout << "Last half forward labeling= " << time_eps << endl;

  QUIT:
  return;
}


void CVRP::eliminateBucketArcs(BbNode *const node, const PtrAllR1CS &ptrAllR1Cs) {
  auto r1c_to_pi = ptrAllR1Cs.r1c_to_pi;
  auto r1c_multi_to_pi = ptrAllR1Cs.r1c_multi_to_pi;
  int dim_sq = dim * dim;
  auto stateBetween2Buckets = new bool[dim_sq * num_buckets_per_vertex];
  auto latest_bucket = new int[dim_sq];
  opt_gap = calculateOptimalGap(node);

#ifdef SYMMETRY_PROHIBIT
  eliminateBucketArcs<true, false>(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);
  eliminateBucketArcs<false, false>(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);
#else
  eliminateBucketArcs<true, true>(node, r1c_to_pi, r1c_multi_to_pi, dim_sq, stateBetween2Buckets, latest_bucket);
#endif

  delete[]stateBetween2Buckets;
  delete[]latest_bucket;
}


void CVRP::obtainJumpArcs(BbNode *const node) const {
  auto tmp = new bitset<2> *[dim];
  for (int i = 0; i < dim; ++i) tmp[i] = new bitset<2>[num_buckets_per_vertex];
  obtainjumpArcs<true>(node, tmp);
#ifdef SYMMETRY_PROHIBIT
  obtainjumpArcs<false>(node, tmp);
#endif
  for (int i = 0; i < dim; ++i) {
    delete[]tmp[i];
  }
  delete[]tmp;
}

void CVRP::checkBucketIJ(BbNode *node, int i, int j) const {
  cout << "from i = " << i << " to j = " << j << endl;
  int num_buckets = 0, num_jump = 0;
  for (int bin = 0; bin < num_buckets_per_vertex; ++bin) {
    auto if_find = std::find(node->all_forward_buckets[i][bin].bucket_arcs.begin(),
                             node->all_forward_buckets[i][bin].bucket_arcs.end(), j);
    if (if_find == node->all_forward_buckets[i][bin].bucket_arcs.end()) {
      auto iff = std::find_if(node->all_forward_buckets[i][bin].jump_arcs.begin(),
                              node->all_forward_buckets[i][bin].jump_arcs.end(),
                              [&](const pair<double, int> &p) { return p.second == j; });
      if (iff == node->all_forward_buckets[i][bin].jump_arcs.end()) {
        continue;
        cout << "(" << bin;
        cout << " n) ";
      } else {
        cout << "(" << bin;
        cout << " j) ";
        ++num_jump;
      }
    } else {
      cout << "(" << bin;
      cout << " B) ";
      ++num_buckets;
    }
  }
  cout << endl;
}
