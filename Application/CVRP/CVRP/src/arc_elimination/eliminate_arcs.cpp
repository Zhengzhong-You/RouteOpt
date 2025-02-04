#include "cvrp.hpp"
#include "template_functors.hpp"

using namespace std;
using namespace std::chrono;

void CVRP::eliminateArcs(BbNode *const node) {
    bool if_ban_arc_elimination = false;
    setting_2_call(banArcElimination(node->index == 0, if_ban_arc_elimination))
    if (if_ban_arc_elimination) return;
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

    priceLabeling(node, optimal_dual_vector);
    getRank1DualsInCG(node, optimal_dual_vector);

    chrono::time_point<chrono::high_resolution_clock> beg, end;
    beg = chrono::high_resolution_clock::now();

    runLabelingForArcElimination(node);

    if (if_short_memory) {
        if_arc_elimination_succeed = false;
        if_arc_elimination_tried_but_failed = true;
        goto QUIT;
    } else if (if_roll_back) {
        if_arc_elimination_succeed = false;
        if_arc_elimination_tried_but_failed = true;
        cout << "Reach the hard limit!" << endl;
        goto QUIT;
    }

    end = chrono::high_resolution_clock::now();
    time_labeling = chrono::duration<double>(end - beg).count();

    beg = chrono::high_resolution_clock::now();
    eliminateBucketArcs(node);
    end = chrono::high_resolution_clock::now();
    time_elimination = chrono::duration<double>(end - beg).count();

    beg = chrono::high_resolution_clock::now();
    obtainJumpArcs(node);
    end = chrono::high_resolution_clock::now();
    time_obtainjump = chrono::duration<double>(end - beg).count();

#if VERBOSE_MODE == 1
    cout << "time summary: " << endl;
    cout << "runLabelingForArcElimination= " << time_labeling << endl;
    cout << "eliminateBucketArcs= " << time_elimination << endl;
    cout << "obtainJumpArcs= " << time_obtainjump << endl;
#endif
    if_arc_elimination_succeed = true;
QUIT:
    cout << BIG_PHASE_SEPARATION;
}

void CVRP::runLabelingForArcElimination(BbNode *const node) {
    double dif;
    bool if_break, if_find, if_keep;
    double path_rc;
    int edgemap;
    if_roll_back = false;
    if_short_memory = false;
    opt_gap = calculateOptimalGap(node);

    auto time_beg = chrono::high_resolution_clock::now();
    auto time_end = time_beg;
    double time_eps;

#ifdef SYMMETRY_PROHIBIT
  auto beg=chrono::high_resolution_clock::now();
  concatenatePhaseInArcElimination<true, false>(node);
  auto end=chrono::high_resolution_clock::now();
  cout<<"concatenatePhaseInArcElimination= "<<chrono::duration<double>(end-beg).count()<<endl;
  if (if_roll_back || if_short_memory ) goto QUIT;
  runLabeling<true, true, false, false, 0>(node);
  if (if_roll_back || if_short_memory ) goto QUIT;
  concatenatePhaseInArcElimination<false, false>(node);
  if (if_roll_back || if_short_memory ) goto QUIT;
  runLabeling<false, true,  false, false, 0>(node);
  if (if_roll_back || if_short_memory ) goto QUIT;
#else
    concatenatePhaseInArcElimination<true, true>(node);
    if (if_roll_back || if_short_memory) goto QUIT;
    runLabeling<true, true, false, true, 0>(node);
    if (if_roll_back || if_short_memory) goto QUIT;
#endif

    time_end = chrono::high_resolution_clock::now();
    time_eps = chrono::duration<double>(time_end - time_beg).count();
#if VERBOSE_MODE == 1
    cout << "Last half forward labeling= " << time_eps << endl;
#endif

QUIT:
    if (if_roll_back) {
        gap_improved_4_arc_elimination_n_enumeration += Config::GapBarIncreased4ArcEliminationNEnumeration;
        arc_elimination_time += Config::ArcEliminationTimeIncreased;
    }
}

void CVRP::eliminateBucketArcs(BbNode *const node) {
    int dim_sq = dim * dim;
    auto stateBetween2Buckets = new bool[dim_sq * num_buckets_per_vertex];
    auto latest_bucket = new int[dim_sq];
    opt_gap = calculateOptimalGap(node);

#ifdef SYMMETRY_PROHIBIT
  eliminateBucketArcs<true, false>(node, dim_sq, stateBetween2Buckets, latest_bucket);
  eliminateBucketArcs<false, false>(node, dim_sq, stateBetween2Buckets, latest_bucket);
#else
    eliminateBucketArcs<true, true>(node, dim_sq, stateBetween2Buckets, latest_bucket);
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
