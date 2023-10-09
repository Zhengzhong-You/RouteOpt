//
// Created by Zhengzhong You on 5/31/22.
//

#include "BBNODE.hpp"
#include "CVRP.hpp"

using namespace std;

BBNODE::BBNODE(int num, int p_col, CVRP *cvrp) {
  allocateMem(num);
  IdxCols = new size_t[CONFIG::MaxNumCols];
  SizeAllocatedMem = num;
  NumParentCols = p_col;

  NumRowInBucketGraph = cvrp->Dim;

  AllForwardBuckets = new Bucket *[NumRowInBucketGraph];
  for (int i = 0; i < NumRowInBucketGraph; ++i) {
    AllForwardBuckets[i] = new Bucket[cvrp->NumBucketsPerVertex];
  }
//  VBasis.resize(CONFIG::MaxNumCols);
//  CBasis.resize(CST_LIMIT);
#ifdef SYMMETRY_PROHIBIT
  AllBackwardBuckets = new Bucket *[NumRowInBucketGraph];
  for (int i = 0; i < NumRowInBucketGraph; ++i) {
    AllBackwardBuckets[i] = new Bucket[cvrp->NumBucketsPerVertex];
  }
#endif
}

#ifdef readEnumerationTrees
BBNODE::BBNODE(int num) {
  allocateMem(num);
  IdxCols = new size_t[CONFIG::MaxNumCols];
  SizeAllocatedMem = num;
}
#endif

BBNODE::BBNODE(BBNODE *node, BidirLinkList *ptr, int p_col, int idx, const BrC &bf, int NumBucketsPerVertex) {
  //copy all data
  Ptr = ptr;
  NumParentCols = p_col;
  Idx = idx;

  SizeAllocatedMem = node->SizeAllocatedMem;
  allocateMem(SizeAllocatedMem);
  IdxCols = new size_t[CONFIG::MaxNumCols];
  for (int i = 0; i < NumParentCols; ++i) {
    IdxCols[i] = node->IdxCols[i];
  }

  solver.SOLVERgetsolver(&node->solver);
  TreeLevel = node->TreeLevel;
  RCCs = node->RCCs;
  R1Cs = node->R1Cs;
  R1Cs_multi = node->R1Cs_multi;
  BrCs = node->BrCs;
  BrCs.emplace_back(bf);
  Val = node->Val;
  NumRowInBucketGraph = node->NumRowInBucketGraph;
  NumForwardBucketArcs = node->NumForwardBucketArcs;
  NumForwardJumpArcs = node->NumForwardJumpArcs;
#ifdef SYMMETRY_PROHIBIT
  NumBackwardBucketArcs = node->NumBackwardBucketArcs;
  NumBackwardJumpArcs = node->NumBackwardJumpArcs;
#endif
  LastGap = node->LastGap;

  AllForwardBuckets = new Bucket *[NumRowInBucketGraph];
  for (int i = 0; i < NumRowInBucketGraph; ++i) {
    AllForwardBuckets[i] = new Bucket[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      AllForwardBuckets[i][j] = node->AllForwardBuckets[i][j];
    }
  }
#ifdef useM_dynamically
  t4oneLP = node->t4oneLP;
  geo_r_star = node->geo_r_star;
  c = node->c;
  obj_change = node->obj_change;
  l_r_ratio = node->l_r_ratio;
#endif
//  VBasis.resize(CONFIG::MaxNumCols);
//  CBasis.resize(CST_LIMIT);
#ifdef SYMMETRY_PROHIBIT
  AllBackwardBuckets = new Bucket *[NumRowInBucketGraph];
  for (int i = 0; i < NumRowInBucketGraph; ++i) {
    AllBackwardBuckets[i] = new Bucket[NumBucketsPerVertex];
    for (int j = 0; j < NumBucketsPerVertex; ++j) {
      AllBackwardBuckets[i][j] = node->AllBackwardBuckets[i][j];
    }
  }
#endif
}

BBNODE::BBNODE(BBNODE *node, int p_col, int idx, const BrC &bf) {
  //copy all data
  Idx = idx;
  //do not delete this, since it will be used later, for safety, require for every function
  NumParentCols = p_col;
  SizeAllocatedMem = node->SizeAllocatedMem;
  allocateMem(SizeAllocatedMem);
  TreeLevel = node->TreeLevel;
  RCCs = node->RCCs;
  R1Cs = node->R1Cs;
  R1Cs_multi = node->R1Cs_multi;
  BrCs = node->BrCs;
  LastGap = node->LastGap;
  BrCs.emplace_back(bf);


  //according to bf, we revise col_pool info Vertex2IdxCols and Edge2IdxCols and lp.
  solver.SOLVERgetsolver(&node->solver);
  SizeEnuColPool = node->SizeEnuColPool;
  validSize = node->validSize;
  IdxCols = new size_t[NumParentCols + SizeEnuColPool];
  IdxColsInEnuColPool = node->IdxColsInEnuColPool;
  Cost4ColsInEnuColPool = node->Cost4ColsInEnuColPool;
//  MatInEnu = node->MatInEnu;
  Deleted_ColsInEnuPool = new bool[SizeEnuColPool];
  copy(node->Deleted_ColsInEnuPool, node->Deleted_ColsInEnuPool + SizeEnuColPool, Deleted_ColsInEnuPool);
  copy(node->IdxCols, node->IdxCols + NumParentCols, IdxCols);
  Val = node->Val;
  NumRowInBucketGraph = node->NumRowInBucketGraph;
#ifdef useM_dynamically
  t4oneLP = node->t4oneLP;
  geo_r_star = node->geo_r_star;
  c = node->c;
  obj_change = node->obj_change;
  l_r_ratio = node->l_r_ratio;
#endif
  //no need for graph info
//  VBasis.resize(CONFIG::MaxNumCols);
//  CBasis.resize(CST_LIMIT);
}

BBNODE::~BBNODE() {
  freeMem();
  delete[]IdxCols;
  delete[]Deleted_ColsInEnuPool;

  if (AllForwardBuckets) {
    for (int i = 0; i < NumRowInBucketGraph; ++i) {
      delete[]AllForwardBuckets[i];
    }
    delete[]AllForwardBuckets;
#ifdef SYMMETRY_PROHIBIT
    for (int i = 0; i < NumRowInBucketGraph; ++i) {
      delete[]AllBackwardBuckets[i];
    }
    delete[]AllBackwardBuckets;
#endif
    if (Ptr)
      this->Ptr->deleteSelf();
  }

  solver.SOLVERfreemodel();
}

void BBNODE::allocateMem(int num) {
  if (EdgeTail) {
    freeMem();
  }

  EdgeTail = new int[num];
  EdgeHead = new int[num];
  EdgeVal = new double[num];
}

void BBNODE::freeMem() const {
  delete[] EdgeTail;
  delete[] EdgeHead;
  delete[] EdgeVal;
}

void BidirLinkList::becomeParent(BBNODE *node, CVRP *cvrp) {
  //At this time, all the cols are in the mem
  //mapping rule is i*fx_dimension+j (where i < j)
  auto tmp_p = node->IdxCols;
  int *pool = cvrp->ColPool4Mem;
  int fx_Dimension = cvrp->Dim;
  int num_col = cvrp->NumCol;

  if (!Edge2Cols.empty()) {
    cout << "Edge2Cols is not empty" << endl;
    exit(0);
  }
  Edge2Cols.clear();

  //note that the 1 - NumArtiVars columns need special treatment
  for (int i = 1; i < fx_Dimension; ++i) {
    Edge2Cols[i].emplace_back(i);
    Edge2Cols[i].emplace_back(i);
  }

  //normal cols
  for (int i = fx_Dimension; i < num_col; ++i) {
    for (size_t j = tmp_p[i];;) {
      int ai = pool[j], aj = pool[++j];
      if (ai > aj) {
        Edge2Cols[aj * fx_Dimension + ai].emplace_back(i);
      } else {
        Edge2Cols[ai * fx_Dimension + aj].emplace_back(i);
      }
      if (!aj) break;
    }
  }

  //delete self
  //the rest of the unclean is completely processed in the destructor of node
  if (this == PNode->LNode) {
    PNode->LNode = nullptr;
  } else {
    PNode->RNode = nullptr;
  }
  PNode->deleteSelf();
  //become parent self
  PNode = nullptr;
  node->NumParentCols = num_col;
}

void BidirLinkList::giveBirth(BidirLinkList *&lnode, BidirLinkList *&rnode) {
  lnode = new BidirLinkList(this);
  rnode = new BidirLinkList(this);
  LNode = lnode;
  RNode = rnode;
}

void BidirLinkList::deleteSelf() {
  //the delete function will judge by itself whether it can delete the current self
  //You must first judge whether you can delete the parent.
  //If the parent can't delete, just disconnect each other directly.
  //Iteratively delete the parent
  if ((!LNode) && (!RNode)) {
    //indicate that the current node can be deleted
    if (PNode) {
      if (this == PNode->LNode) {
        PNode->LNode = nullptr;
      } else {
        PNode->RNode = nullptr;
      }
      PNode->deleteSelf();
    }
    delete this;
  }
}

PtrAllR1Cs::PtrAllR1Cs(BBNODE *node, CVRP *cvrp) {
  int size_r1c = int(node->R1Cs.size());
  r1c_to_pi = new double[size_r1c];
  num_r1c_nzero = 0;

  auto &pi = cvrp->Pi4Labeling;

  for (int i = 0; i < size_r1c; ++i) {
    if (abs(pi[node->R1Cs[i].IdxR1C]) < DUAL_TOLERANCE) continue;
    r1c_to_pi[num_r1c_nzero++] = pi[node->R1Cs[i].IdxR1C];
  }

  int size_r1c_multi = int(node->R1Cs_multi.size());
  r1c_multi_to_pi = new double[size_r1c_multi];
  num_r1c_multi_nzero = 0;

  for (int i = 0; i < size_r1c_multi; ++i) {
    if (abs(pi[node->R1Cs_multi[i].IdxR1C]) < DUAL_TOLERANCE) continue;
    r1c_multi_to_pi[num_r1c_multi_nzero++] = pi[node->R1Cs_multi[i].IdxR1C];
  }
  cvrp->convertVertex2R1CsInPricing(node);
}

PtrAllR1Cs::~PtrAllR1Cs() {
  delete[]r1c_to_pi;
  delete[]r1c_multi_to_pi;
}

BBNODE::BBNODE(BBNODE *node, int NumBucketsPerVertex, int NumCol, const bool *if_use_arc) {
  //copy all data
  NumParentCols = node->NumParentCols;
  Idx = node->Idx;

  SizeAllocatedMem = node->SizeAllocatedMem;
  allocateMem(SizeAllocatedMem);
  IdxCols = new size_t[CONFIG::MaxNumCols];
  for (int i = 0; i < NumCol; ++i) {
    IdxCols[i] = node->IdxCols[i];
  }

  solver.SOLVERgetsolver(&node->solver);
  safe_solver(solver.SOLVERreoptimize())
  TreeLevel = node->TreeLevel;
  RCCs = node->RCCs;
  R1Cs = node->R1Cs;
  R1Cs_multi = node->R1Cs_multi;
  BrCs = node->BrCs;
  Val = node->Val;
//  VBasis = node->VBasis;
//  CBasis = node->CBasis;
  NumRowInBucketGraph = node->NumRowInBucketGraph;
  NumForwardBucketArcs = node->NumForwardBucketArcs;
  NumForwardJumpArcs = node->NumForwardJumpArcs;
#ifdef SYMMETRY_PROHIBIT
  NumBackwardBucketArcs = node->NumBackwardBucketArcs;
  NumBackwardJumpArcs = node->NumBackwardJumpArcs;
#endif
  LastGap = node->LastGap;

  AllForwardBuckets = new Bucket *[NumRowInBucketGraph];
  for (int i = 0; i < NumRowInBucketGraph; ++i) {
    AllForwardBuckets[i] = new Bucket[NumBucketsPerVertex];
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      auto &bucket = AllForwardBuckets[i][b];
      for (auto &j : node->AllForwardBuckets[i][b].BucketArcs) {
        if (if_use_arc[i * NumRowInBucketGraph + j]) {
          bucket.BucketArcs.emplace_back(j);
        }
      }
      for (auto &arc : node->AllForwardBuckets[i][b].JumpArcs) {
        if (if_use_arc[i * NumRowInBucketGraph + arc.second]) {
          bucket.JumpArcs.emplace_back(arc);
        }
      }
    }
  }
#ifdef useM_dynamically
  t4oneLP = node->t4oneLP;
  geo_r_star = node->geo_r_star;
  c = node->c;
  obj_change = node->obj_change;
  l_r_ratio = node->l_r_ratio;
#endif
#ifdef SYMMETRY_PROHIBIT
  AllBackwardBuckets = new Bucket *[NumRowInBucketGraph];
  for (int i = 0; i < NumRowInBucketGraph; ++i) {
    AllBackwardBuckets[i] = new Bucket[NumBucketsPerVertex];
    for (int b = 0; b < NumBucketsPerVertex; ++b) {
      auto &bucket = AllBackwardBuckets[i][b];
      for (auto &j : node->AllBackwardBuckets[i][b].BucketArcs) {
        if (if_use_arc[i * NumRowInBucketGraph + j]) {
          bucket.BucketArcs.emplace_back(j);
        }
      }
      for (auto &arc : node->AllBackwardBuckets[i][b].JumpArcs) {
        if (if_use_arc[i * NumRowInBucketGraph + arc.second]) {
          bucket.JumpArcs.emplace_back(arc);
        }
      }
    }
  }
#endif
}

#ifdef NominalBranchingInEnu
BBNODE::BBNODE(BBNODE *node, int p_col, const std::vector<int> &BranchingColSet, int idx, const NBrC &nbf) {
  //copy basic all data
  Idx = idx;
  NumParentCols = p_col;
  SizeAllocatedMem = node->SizeAllocatedMem;
  allocateMem(SizeAllocatedMem);
  TreeLevel = node->TreeLevel;
  RCCs = node->RCCs;
  R1Cs = node->R1Cs;
  R1Cs_multi = node->R1Cs_multi;
  BrCs = node->BrCs;
  NBrCs = node->NBrCs;
  NBrCs.emplace_back(nbf);
  LastGap = node->LastGap;

  solver.SOLVERgetsolver(&node->solver);
  SizeEnuColPool = node->SizeEnuColPool;
  IdxCols = new size_t[NumParentCols + SizeEnuColPool];
  Deleted_ColsInEnuPool = new bool[SizeEnuColPool];
  copy(node->Deleted_ColsInEnuPool, node->Deleted_ColsInEnuPool + SizeEnuColPool, Deleted_ColsInEnuPool);
  for (auto &i : BranchingColSet) {
    Deleted_ColsInEnuPool[i] = true;
  }
  copy(node->IdxCols, node->IdxCols + NumParentCols, IdxCols);
  Val = node->Val;
  NumRowInBucketGraph = node->NumRowInBucketGraph;
}
#endif

#ifdef useM_dynamically
void BBNODE::updateState(double new_value, double &old_value, int n) {
  //use geometric mean
//  G^{\prime}=G \times(a)^{1 /(n+1)} \times\left(G^{-1 / n}\right)
  if (n == 0) {
    old_value = new_value;
  } else {
    old_value = old_value * pow(new_value, 1.0 / (n + 1)) * pow(1.0 / old_value, 1.0 / (n + 1));
  }
}

void BBNODE::calculateR_star(double lift, double &new_r_star, CVRP *cvrp) {
  if (Idx == 0) throw runtime_error("calculateR_star: Idx==0");
  auto &edge = BrCs.back().Edge;
  auto dir = BrCs.back().BrDir;
  auto &info = obj_change[edge];

  if (dir) {
    if (get<0>(info) == 0) {
      get<0>(info) = max(lift * l_r_ratio, TOLERANCE);
    } else {
      if (get<1>(info) < TOLERANCE) {
        get<0>(info) = TOLERANCE;
      } else {
        get<0>(info) *= max(lift / get<1>(info), TOLERANCE);
      }
    }
    get<1>(info) = max(lift, TOLERANCE);
  } else {
    get<0>(info) = lift;
    if (get<1>(info) == 0) {
      get<1>(info) = max(lift / l_r_ratio, TOLERANCE);
    } else {
      if (get<0>(info) < TOLERANCE) {
        get<1>(info) = TOLERANCE;
      } else {
        get<1>(info) *= max(lift / get<0>(info), TOLERANCE);
      }
    }
  }

  cvrp->if_use_full_k = false;
  if (get<0>(info) == TOLERANCE || get<1>(info) == TOLERANCE) {
    cvrp->if_use_full_k = true;
    new_r_star = geo_r_star;
    return;
  }

  auto l = get<0>(info);
  auto r = get<1>(info);

  auto bar_r = sqrt(l * r);
  int m = cvrp->m;
  int n = cvrp->ml.giveInitialScreeningNumCandidates();
  auto k = min(get<2>(info), m);// generalize the root node branching decision
  new_r_star = bar_r / ((m + 1.0) / n * k / (k + 1.0) + 1.0 - m / n);
  cout << "new_r_star = " << new_r_star << endl;
  cout << "l_r_ratio1 = " << l_r_ratio << endl;
  auto new_l_r = get<0>(info) / get<1>(info);
  updateState(new_l_r, l_r_ratio, BrCs.size() - 1);
  cout << "l_r_ratio2 = " << l_r_ratio << endl;
}
#endif
