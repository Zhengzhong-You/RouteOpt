
#include "CVRP.hpp"

using namespace std;
using namespace Eigen;

void CVRP::generateRCMatrix4Label(std::vector<Eigen::Triplet<double> > &triplets,
								  Label *label,
								  RowVectorXd &cost,
								  std::unordered_map<std::pair<int, int>,
													 std::vector<int>,
													 PairHasher> &edge_map) {
  auto &aux_label = label->aux_label;
  int index = aux_label->index4_rc_matrix;
  int past_node = 0, curr_node;
  double cost_sum = 0;

  std::unordered_map<int, int> cnt;
  cnt.reserve(lp_r1c_denominator.size());
  auto &state = aux_label->states;
  auto &sparse_rep = aux_label->sparse_lp_states;
  auto &valid_sparse_num = aux_label->sparse_num;

  vector<int> col;
  col.reserve(dim);
  auto p = label;
  while (p) {
	col.emplace_back(p->end_vertex);
	p = p->p_label;
  }

  vector<double> vis(dim, 0);
  for (auto it = col.rbegin() + 1; it != col.rend(); ++it) {
	curr_node = *it;
	cost_sum += cost_mat4_vertex[past_node][curr_node];
	++vis[curr_node];
	auto pr = past_node < curr_node ? make_pair(past_node, curr_node) : make_pair(curr_node, past_node);
	edge_map[pr].emplace_back(index);

	getCoefficientExtendR1C(state, sparse_rep, cnt, valid_sparse_num, past_node, curr_node);

	past_node = curr_node;
  }
  vis[curr_node] -= 0.5;
  cost(index) = cost_sum;

  int old_size = (int)triplets.size();
  triplets.resize(old_size + cnt.size() + dim);
  for (auto &pr : cnt) {
	triplets[old_size++] = {lp_r1c_map.at(pr.first), index, (double)pr.second};
  }
  for (int i = 1; i < dim; ++i) {
	if (vis[i] > 0.4) {
	  triplets[old_size++] = {i - 1, index, vis[i]};
	}
  }
  triplets.resize(old_size);
}




