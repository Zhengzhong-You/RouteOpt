#include "branching_node.hpp"
#include "cvrp.hpp"

using namespace std;

BbNode::BbNode(int num, CVRP *cvrp) {
    allocateMem(num);
    num_rows_in_bucket_graph = cvrp->dim;
    all_forward_buckets = new Bucket *[num_rows_in_bucket_graph];
    for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
        all_forward_buckets[i] = new Bucket[cvrp->num_buckets_per_vertex];
    }
#ifdef SYMMETRY_PROHIBIT
  all_backward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_backward_buckets[i] = new Bucket[cvrp->num_buckets_per_vertex];
  }
#endif
}

BbNode::BbNode(int num) {
    allocateMem(num);
}

BbNode::BbNode(BbNode *node,
               int idx,
               const Brc &bf,
               int num_buckets_per_vertex) {
    index = idx;
    cols = node->getCols();
    solver.getSolver(&node->getSolver());
    tree_level = node->getTreeLevel();
    rccs = node->rccs;
    r1cs = node->r1cs;
    brcs = node->getBrCs();
    brcs.emplace_back(bf);
    value = node->getCurrentNodeVal();
    num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
    num_forward_bucket_arcs = node->num_forward_bucket_arcs;
    num_forward_jump_arcs = node->num_forward_jump_arcs;
    topological_order_forward = node->topological_order_forward;
#ifdef SYMMETRY_PROHIBIT
  topological_order_backward = node->topological_order_backward;
  num_backward_bucket_arcs = node->num_backward_bucket_arcs;
  num_backward_jump_arcs = node->num_backward_jump_arcs;
#endif
    last_gap = node->last_gap;

    all_forward_buckets = new Bucket *[num_rows_in_bucket_graph];
    for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
        all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
        for (int j = 0; j < num_buckets_per_vertex; ++j) {
            all_forward_buckets[i][j] = node->all_forward_buckets[i][j];
        }
    }
    dynamic_call(Dynamics::copyDynamicData4EachNode(node->dynamics, this->dynamics))
    symmetry_prohibit_call(all_backward_buckets = new Bucket *[num_rows_in_bucket_graph];
        for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
        all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
        for (int j = 0; j < num_buckets_per_vertex; ++j) {
        all_backward_buckets[i][j] = node->all_backward_buckets[i][j];
        }
        })
}

BbNode::BbNode(BbNode *node, int idx, const Brc &bf) {
    index = idx;
    allocateMem((int) node->getEdgeHead().size());
    tree_level = node->getTreeLevel();
    rccs = node->rccs;
    r1cs = node->r1cs;
    brcs = node->getBrCs();
    last_gap = node->last_gap;
    brcs.emplace_back(bf);

    solver.getSolver(&node->getSolver());
    size_enumeration_col_pool = node->size_enumeration_col_pool;
    valid_size = node->valid_size;

    index_columns_in_enumeration_column_pool = node->index_columns_in_enumeration_column_pool;
    cost_for_columns_in_enumeration_column_pool = node->cost_for_columns_in_enumeration_column_pool;
    deleted_columns_in_enumeration_pool.assign(node->deleted_columns_in_enumeration_pool.begin(),
                                               node->deleted_columns_in_enumeration_pool.begin()
                                               + size_enumeration_col_pool);
    cols = node->getCols();
    value = node->getCurrentNodeVal();
    num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
    dynamic_call(Dynamics::copyDynamicData4EachNode(node->dynamics, this->dynamics))
}

BbNode::~BbNode() {
    if (all_forward_buckets) {
        for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
            delete[]all_forward_buckets[i];
        }
        delete[]all_forward_buckets;
        symmetry_prohibit_call(for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
            delete[]all_backward_buckets[i];
            }
            delete[]all_backward_buckets;)
    }
    solver.freeModel();
}

void BbNode::allocateMem(int num) {
    edge_head.resize(num);
    edge_tail.resize(num);
    edge_value.resize(num);
}

BbNode::BbNode(BbNode *node, int num_buckets_per_vertex, const bool *if_use_arc) {
    index = node->index;

    cols = node->getCols();
    allocateMem((int) node->getEdgeHead().size());
    solver.getSolver(&node->getSolver());
    safe_solver(solver.reoptimize())
    tree_level = node->getTreeLevel();
    rccs = node->rccs;
    r1cs = node->r1cs;
    brcs = node->getBrCs();
    value = node->getCurrentNodeVal();
    num_rows_in_bucket_graph = node->num_rows_in_bucket_graph;
    num_forward_bucket_arcs = node->num_forward_bucket_arcs;
    num_forward_jump_arcs = node->num_forward_jump_arcs;
    topological_order_forward = node->topological_order_forward;
    symmetry_prohibit_call(
        topological_order_backward = node->topological_order_backward;
        num_backward_bucket_arcs = node->num_backward_bucket_arcs;
        num_backward_jump_arcs = node->num_backward_jump_arcs;)
    last_gap = node->last_gap;

    all_forward_buckets = new Bucket *[num_rows_in_bucket_graph];
    for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
        all_forward_buckets[i] = new Bucket[num_buckets_per_vertex];
        for (int b = 0; b < num_buckets_per_vertex; ++b) {
            auto &bucket = all_forward_buckets[i][b];
            for (auto &j: node->all_forward_buckets[i][b].bucket_arcs) {
                if (if_use_arc[i * num_rows_in_bucket_graph + j]) {
                    bucket.bucket_arcs.emplace_back(j);
                }
            }
            for (auto &arc: node->all_forward_buckets[i][b].jump_arcs) {
                if (if_use_arc[i * num_rows_in_bucket_graph + arc.second]) {
                    bucket.jump_arcs.emplace_back(arc);
                }
            }
        }
    }
    dynamic_call(Dynamics::copyDynamicData4EachNode(node->dynamics, this->dynamics))
#ifdef SYMMETRY_PROHIBIT
  all_backward_buckets = new Bucket *[num_rows_in_bucket_graph];
  for (int i = 0; i < num_rows_in_bucket_graph; ++i) {
	all_backward_buckets[i] = new Bucket[num_buckets_per_vertex];
	for (int b = 0; b < num_buckets_per_vertex; ++b) {
	  auto &bucket = all_backward_buckets[i][b];
	  for (auto &j : node->all_backward_buckets[i][b].bucket_arcs) {
		if (if_use_arc[i * num_rows_in_bucket_graph + j]) {
		  bucket.bucket_arcs.emplace_back(j);
		}
	  }
	  for (auto &arc : node->all_backward_buckets[i][b].jump_arcs) {
		if (if_use_arc[i * num_rows_in_bucket_graph + arc.second]) {
		  bucket.jump_arcs.emplace_back(arc);
		}
	  }
	}
  }
#endif
}


int &BbNode::getTreeLevel() {
    return tree_level;
}


double BbNode::obtainParentNodeValue() {
    if (tree_level == 0 && !if_just_enter_enu) {
        return 0;
    }
    return value;
}

double &BbNode::getCurrentNodeVal() {
    return value;
}

int &BbNode::getNumEdges() {
    return num_edges;
}

std::vector<double> &BbNode::getEdgeValue() {
    return edge_value;
}

std::vector<int> &BbNode::getEdgeTail() {
    return edge_tail;
}

std::vector<int> &BbNode::getEdgeHead() {
    return edge_head;
}


int &BbNode::getNodeIdx() {
    return index;
}

void BbNode::printBrCInfo() {
    if (tree_level) {
        for (const auto &brc: brcs) {
            cout << (brc.br_dir ? "true" : "false") << "(" << brc.edge.first << "," << brc.edge.second << ")" <<
                    " ";
        }
        cout << endl;
    }
}

Solver &BbNode::getSolver() {
    return solver;
}


bool &BbNode::getIfJustEnterEnu() {
    return if_just_enter_enu;
}

double &BbNode::getNodeBrValueImproved() {
    return br_value_improved;
}

std::vector<SequenceInfo> &BbNode::getCols() {
    return cols;
}


std::vector<Brc> &BbNode::getBrCs() {
    return brcs;
}

bool &BbNode::getIfTerminated() {
    return is_terminated;
}



int &BbNode::getDynamicKNode() {
    return dynamics.getKNode();
}

