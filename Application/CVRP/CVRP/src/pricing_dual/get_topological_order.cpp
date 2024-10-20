#include "cvrp.hpp"

using namespace std;

void fillOrder(const unordered_map<int, vector<int>> &graph,
			   int node,
			   unordered_set<int> &visited,
			   stack<int> &finishStack) {
  vector<bool> onStack(graph.size() + 1, false);//since we start from 1
  stack<int> dfsStack;
  dfsStack.push(node);

  while (!dfsStack.empty()) {
	int curr = dfsStack.top();
	if (!visited.count(curr)) {
	  visited.insert(curr);
	  onStack[curr] = true;
	  if (graph.find(curr) != graph.end()) {
		for (int neighbor : graph.at(curr)) {
		  if (!visited.count(neighbor) && !onStack[neighbor]) {
			dfsStack.push(neighbor);
		  }
		}
	  }
	} else {
	  dfsStack.pop();
	  if (onStack[curr]) {
		finishStack.push(curr);
		onStack[curr] = false;
	  }
	}
  }
}

void getTranspose(const unordered_map<int, vector<int>> &graph,
				  unordered_map<int, vector<int>> &transposeGraph) {
  transposeGraph.clear();
  for (auto &pair : graph) {
	for (int neighbor : pair.second) {
	  transposeGraph[neighbor].push_back(pair.first);
	}
  }
}

void dfsTransposed(const unordered_map<int, vector<int>> &graph,
				   int node,
				   unordered_set<int> &visited,
				   vector<int> &component) {
  stack<int> dfsStack;
  dfsStack.push(node);

  while (!dfsStack.empty()) {
	int curr = dfsStack.top();
	dfsStack.pop();
	if (!visited.count(curr)) {
	  visited.insert(curr);
	  component.push_back(curr);
	  if (graph.count(curr)) {
		for (int neighbor : graph.at(curr)) {
		  if (!visited.count(neighbor)) {
			dfsStack.push(neighbor);
		  }
		}
	  }
	}
  }
}

void kosaraju(const unordered_map<int, vector<int>> &graph,
			  vector<vector<int>> &scc) {

  stack<int> orderStack;
  unordered_set<int> visited;
  scc.clear();

  for (auto &pair : graph) {
	if (!visited.count(pair.first)) {
	  fillOrder(graph, pair.first, visited, orderStack);
	}
  }

  unordered_map<int, vector<int>> transposedGraph;
  getTranspose(graph, transposedGraph);

  visited.clear();
  while (!orderStack.empty()) {
	int curr = orderStack.top();
	orderStack.pop();
	if (!visited.count(curr)) {
	  vector<int> component;
	  dfsTransposed(transposedGraph, curr, visited, component);
	  scc.emplace_back(component);
	}
  }
}

void buildCondensedGraph(const vector<vector<int>> &scc,
						 const unordered_map<int, vector<int>> &graph,
						 unordered_map<int, vector<int>> &condensedGraph) {

  unordered_map<int, int> nodeToComponent;
  int index = 0;
  for (const auto &component : scc) {
	for (int node : component) {
	  nodeToComponent[node] = index;
	}
	index++;
  }

  condensedGraph.clear();
  for (int i = 0; i < scc.size(); ++i) {
	unordered_set<int> seen;
	for (int node : scc[i]) {
	  if (graph.find(node) != graph.end()) {
		for (int neighbor : graph.at(node)) {
		  int neighborComponent = nodeToComponent[neighbor];
		  if (neighborComponent != i && seen.find(neighborComponent) == seen.end()) {
			condensedGraph[i].push_back(neighborComponent);
			seen.insert(neighborComponent);
		  }
		}
	  }
	}
  }
}

void topologicalSort(const unordered_map<int, vector<int>> &graph, vector<int> &order) {
  unordered_map<int, int> inDegree;
  for (const auto &pair : graph) {
	for (int neighbor : pair.second) {
	  inDegree[neighbor]++;
	}
  }

  queue<int> queue;
  for (const auto &pair : graph) {
	if (inDegree[pair.first] == 0) {
	  queue.push(pair.first);
	}
  }

  order.clear();
  while (!queue.empty()) {
	int node = queue.front();
	queue.pop();
	order.push_back(node);
	if (graph.find(node) == graph.end())continue;
	for (int neighbor : graph.at(node)) {
	  if (--inDegree[neighbor] == 0) {
		queue.push(neighbor);
	  }
	}
  }
}

template<bool dir>
void CVRP::getTopologicalOrder4OneBin(BbNode *node, int b) {
  auto &order = dir ? node->topological_order_forward[b] : node->topological_order_backward[b];
  order.clear();
  ResTuple first_res;
  dir ? first_res.first_res = {b * step_size} : first_res.first_res = {(b + 1) * step_size - 1};
  res_int max_res;
  dir ? max_res = (b + 1) * step_size : max_res = b * step_size - 1;
  ResTuple tmp_res;
  unordered_map<int, vector<int>> graph;
  for (int i = 1; i < dim; ++i) {
	yzzLong connect = 0;
	for (auto &j : dir ? node->all_forward_buckets[i][b].bucket_arcs :
				   node->all_backward_buckets[i][b].bucket_arcs
		)
	  connect.set(j);
	graph[i] = vector<int>();
	for (int j = 1; j < dim; ++j) {
	  if (!connect.test(j)) continue;
	  if (i == j) continue;
	  dir ? increaseMainResourceConsumption(first_res, tmp_res, i, j) :
	  decreaseMainResourceConsumption(first_res, tmp_res, i, j);
	  if (dir ? tmp_res.first_res >= max_res : tmp_res.first_res <= max_res) continue;
	  graph[i].emplace_back(j);
	}
  }
  vector<vector<int>> scc;
  kosaraju(graph, scc);
  unordered_map<int, vector<int>> condensedGraph;
  buildCondensedGraph(scc, graph, condensedGraph);
  vector<int> topoOrder;
  topologicalSort(condensedGraph, topoOrder);
  unordered_set<int> tmpSet;
  tmpSet.reserve(topoOrder.size());
  order.resize(scc.size());
  for (int i = 0; i < topoOrder.size(); ++i) {
	int index = topoOrder[i];
	order[i] = scc[index];
	tmpSet.emplace(index);
  }
  vector<vector<int>> tmp;
  tmp.reserve(scc.size());
  for (int i = 0; i < scc.size(); ++i) {
	if (tmpSet.find(i) == tmpSet.end()) {
	  tmp.emplace_back(scc[i]);
	}
  }
  sort(tmp.begin(), tmp.end(), [](const vector<int> &a, const vector<int> &b) {
	return a.size() > b.size();
  });

  auto index = topoOrder.size();
  for (auto &vec : tmp)order[index++] = vec;
}

void CVRP::getTopologicalOrder(BbNode *node) {
  node->topological_order_forward.resize(num_buckets_per_vertex);
  for (int b = 0; b < num_buckets_per_vertex; ++b) {
	getTopologicalOrder4OneBin<true>(node, b);
  }
  symmetry_prohibit_call(node->topological_order_backward.resize(num_buckets_per_vertex);
							 for (int b = 0; b < num_buckets_per_vertex; ++b) {
							   getTopologicalOrder4OneBin<false>(node, b);
							 })
}
