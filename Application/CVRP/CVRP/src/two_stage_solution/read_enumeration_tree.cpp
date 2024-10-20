

#include "read_enumeration_tree.hpp"
#include "branching.hpp"
#if SOLUTION_TYPE == 1
#include "best_bound_first_branching.hpp"
#elif SOLUTION_TYPE == 2
#include "depth_first_branching.hpp"
#endif
using namespace std;
using namespace chrono;
using namespace Eigen;

CVRP *ReadEnumerationTree::cvrp{};
BbNode *ReadEnumerationTree::node{};
string ReadEnumerationTree::tree_path{};
string ReadEnumerationTree::col_pool_path{};

void ReadEnumerationTree::getPath() {
  tree_path = TREE_FOLDER + "/" + Config::tree_path;
  col_pool_path = COL_POOL_FOLDER + "/" + Config::col_pool_path;
  auto pos = tree_path.find_last_of('/') + 1;
  auto end = tree_path.find_last_of('.');
  cout << tree_path.substr(pos, end) << endl;
  cvrp->file_name = tree_path.substr(pos, end);
}

void parseRCC(Rcc &rcc, const std::string &line);
void parseR1C(R1c &r1c, const std::string &line);
void parseBrC(Brc &brc, const std::string &line);

void ReadEnumerationTree::createNodeModel() {
  cvrp->setSolverEnv();
  node->getSolver().getEnv(&cvrp->solver);

  cout << SMALL_PHASE_SEPARATION;
  cout << "<Instance  " << cvrp->file_name << "  Capacity  " << cvrp->cap << ">" << endl;

  int cnt = cvrp->dim;
  vector<double> rhs(cnt, 1);
  vector<char> sense(cnt, SOLVER_EQUAL);
  rhs[cvrp->real_dim] = cvrp->num_vehicle;
  sense[cvrp->real_dim] = SOLVER_GREATER_EQUAL;
  for (auto &rcc : node->rccs) {
	rcc.idx_rcc = cnt++;
	double k;
	if (rcc.form_rcc) {
	  k = (int)rcc.info_rcc_customer.size() - rcc.rhs;
	} else {
	  k = (int)rcc.info_rcc_outside_customer.size() - rcc.rhs;
	}
	rhs.emplace_back(k);
	sense.emplace_back(SOLVER_GREATER_EQUAL);
  }
  for (auto &r1c : node->r1cs) {
	r1c.idx_r1c = cnt++;
	rhs.emplace_back(r1c.rhs);
	sense.emplace_back(SOLVER_LESS_EQUAL);
  }
  if (cnt != cvrp->getNumRow()) throw runtime_error("cnt!=cvrp->getNumRow()");
  cvrp->getVCutMapLP(node);

  const char *model_name = "CVRP.lp";
  safe_solver(node->getSolver().newModel(model_name, 0, nullptr, nullptr, nullptr, nullptr, nullptr))
  safe_solver(node->getSolver().addConstraints(cnt, 0, nullptr, nullptr, nullptr, sense.data(), rhs.data(), nullptr))
  safe_solver(node->getSolver().updateModel())
}

void ReadEnumerationTree::restoreModel() {

  node = new BbNode(500);

  vector<pair<size_t, double>> col_idx;

  readEnumTree(col_idx);
  readColumnPool(col_idx);

  node->getTreeLevel() = 0;
  node->index = 0;

#if SOLUTION_TYPE == 1
  BestBoundFirstBranching::sub_bbt.push(node);
#elif SOLUTION_TYPE == 2
  DepthFirstBranching::addNodeIn(DepthFirstBranching::sub_bbt, node);
#endif
  BaseBranching::lb = 0;
  BaseBranching::lb_transformed = 0;
}

void parseRCC(Rcc &rcc, const std::string &line) {
  std::stringstream ss(line);
  std::string token;

  int num;
  ss >> rcc.form_rcc;
  ss >> token;
  ss >> rcc.rhs;
  ss >> token;
  ss >> rcc.idx_rcc;
  ss >> token;
  while (ss >> token) {
	if (token == "|") {
	  break;
	}
	rcc.info_rcc_customer.emplace_back(stoi(token));
  }
  while (ss >> num) {
	rcc.info_rcc_outside_customer.emplace_back(num);
  }
}

void parseR1C(R1c &r1c, const std::string &line) {
  std::stringstream ss(line);
  std::string token;
  int num;
  ss >> r1c.rhs;
  ss >> token;
  ss >> r1c.idx_r1c;
  ss >> token;
  while (ss >> num) {
	r1c.info_r1c.first.emplace_back(num);
	ss >> std::ws;
	if (ss.peek() == '|') {
	  ss >> token;
	  break;
	}
  }

  ss >> r1c.info_r1c.second;
}

void parseBrC(Brc &brc, const std::string &line) {
  std::stringstream ss(line);
  std::string token;
  int num;

  ss >> brc.edge.first;
  ss >> token;
  ss >> brc.edge.second;
  ss >> token;
  ss >> brc.br_dir;
  brc.idx_br_c = -1;
}

void ReadEnumerationTree::init(CVRP *pr_cvrp) {
  ReadEnumerationTree::cvrp = pr_cvrp;
}

void ReadEnumerationTree::readEnumTree(std::vector<pair<size_t, double>> &col_idx) {
  std::ifstream file(tree_path);
  std::string line;

  if (!file.is_open()) throw runtime_error("cannot open file: " + tree_path);

  while (std::getline(file, line)) {
	if (line.find("nd_ind=") != std::string::npos) {
	  std::istringstream ss(line);
	  std::string key;
	  while (ss >> key) {
		if (key == "nd_col=") {
		  ss >> cvrp->getNumCol();
		} else if (key == "nd_cstr=") {
		  ss >> cvrp->getNumRow();
		} else if (key == "ub=") {
		  ss >> BaseBranching::ub;
		} else if (key == "ColPool=") {
		  int all_col;
		  ss >> all_col;
		  all_col += cvrp->getNumCol();
		  col_idx.resize(all_col);
		}
	  }
	} else if (line.find("ColPool:") != std::string::npos) {
	  int cnt = 0;
	  std::istringstream ss(line);
	  std::string key;
	  ss >> key;//skip ColPool:
	  size_t index, next_index;
	  double value;
	  ss >> key;
	  index = std::stoull(key);
	  while (ss >> key) {
		if (key.find(',') != std::string::npos) {
		  size_t comma_pos = key.find(',');
		  value = std::stod(key.substr(0, comma_pos));
		  if (comma_pos + 1 < key.size()) {
			next_index = std::stoull(key.substr(comma_pos + 1));
		  }
		  col_idx[cnt++] = {index, value};
		  index = next_index;
		}
	  }
	  col_idx.resize(cnt);
	} else if (line.find("rccs") != std::string::npos) {
	  while (std::getline(file, line)) {
		if (line.find("r1cs") != std::string::npos) {
		  break;
		}
		Rcc rcc;
		parseRCC(rcc, line);
		node->rccs.emplace_back(std::move(rcc));
	  }
	  while (std::getline(file, line)) {
		if (line.find("Brc") != std::string::npos) {
		  break;
		}
		R1c r1c;
		parseR1C(r1c, line);
		node->r1cs.emplace_back(std::move(r1c));
	  }
	  while (std::getline(file, line)) {
		Brc brc;
		parseBrC(brc, line);
		node->getBrCs().emplace_back(std::move(brc));
	  }
	}
  }


  file.close();
  cout << "read enumeration tree successfully!" << endl;
}

void ReadEnumerationTree::readColumnPool(
	const std::vector<std::pair<size_t, double>> &col_idx) {
  auto &col_pool4_pricing = cvrp->col_pool4_pricing;
  auto &pool_beg4_pricing = cvrp->pool_beg4_pricing;
  auto &mem4_pricing = cvrp->mem4_pricing;
  unordered_map<size_t, double> col_set;
  col_set.reserve(col_idx.size());
  for (const auto &i : col_idx) {
	if (col_set.find(i.first) == col_set.end()) {
	  col_set[i.first] = i.second;
	} else {
	  if (col_set[i.first] > i.second) {
		col_set[i.first] = i.second;
	  }
	}
  }
  node->size_enumeration_col_pool = (int)col_set.size();
  node->cost_for_columns_in_enumeration_column_pool.resize(node->size_enumeration_col_pool);
  node->index_columns_in_enumeration_column_pool.resize(node->size_enumeration_col_pool);
  std::ifstream file(col_pool_path);
  std::string line;
  size_t lineNumber = 0;
  int cnt = 0;
  int numCol = 0, numPool = 0;
  if (!file.is_open()) throw runtime_error("cannot open file: " + col_pool_path);

  vector<size_t> tmp_col;
  tmp_col.reserve(node->size_enumeration_col_pool);
  double value;
  int num_col_skips = 0;
  int num;

  double k = ceil(accumulate(cvrp->demand + 1, cvrp->demand + cvrp->dim, 0.0) / cvrp->cap);
  double aver_len = ceil(cvrp->dim / k) + 2;
  auto allocate_size = size_t(node->size_enumeration_col_pool * aver_len);
  cvrp->reallocatePricingPool(allocate_size);
  auto PricingWarning = (size_t)(0.9 * (double)mem4_pricing) - cvrp->dim;

  while (std::getline(file, line)) {
	if (col_set.find(lineNumber) != col_set.end()) {
	  std::stringstream ss(line);
	  char delimiter;
	  ss >> value;
	  if (value < col_set[lineNumber] - TOLERANCE) {
		++num_col_skips;
		goto HERE;
	  }
	  node->cost_for_columns_in_enumeration_column_pool[cnt] = value;
	  node->index_columns_in_enumeration_column_pool[cnt] = pool_beg4_pricing;
	  ss >> delimiter;
	  col_pool4_pricing[pool_beg4_pricing++] = 0;
	  while (ss >> num) {
		col_pool4_pricing[pool_beg4_pricing++] = num;
	  }
	  col_pool4_pricing[pool_beg4_pricing++] = 0;
	  ++cnt;
	}
HERE:
	++lineNumber;
	cvrp->resizePoolWarning(PricingWarning);
  }
  file.close();
  node->getCols().clear();

  node->size_enumeration_col_pool = cnt;
  node->valid_size = cnt;
  node->cost_for_columns_in_enumeration_column_pool.conservativeResize(node->size_enumeration_col_pool);
  node->index_columns_in_enumeration_column_pool.conservativeResize(node->size_enumeration_col_pool);
  cvrp->getIfInEnuState() = true;

  createNodeModel();

  cvrp->generateVertex2IndexColsAndEdge2IndexCols(node);

  vector<int> idx(cvrp->getNumCol());
  iota(idx.begin(), idx.end(), 0);
  cvrp->getNumCol() = 0;
  cvrp->addColumnsByInspection(node, idx);

  cout << "num_col_skips= " << num_col_skips << endl;
  cout << "cnt= " << cnt << endl;
  cout << "node->size_enumeration_col_pool= " << node->size_enumeration_col_pool << endl;
  cout << "read column pool successfully!" << endl;
}
