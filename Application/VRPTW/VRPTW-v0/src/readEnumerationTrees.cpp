

#include "CVRP.hpp"

#ifdef READ_ENUMERATION_TREES
using namespace std;
using namespace chrono;
using namespace Eigen;


void CVRP::solveEnumTree() {

  cout << "\n<solveEnumTree>\n" << endl;

  cout << "turn on only read enumeration mode!" << endl;

  lateProcessing();
  glo_beg = std::chrono::high_resolution_clock::now();

  BbNode *node = nullptr;
  restoreModel(node);

  node->value = 0;//reset
  cout << "reset node->value= " << node->value << endl;
  cout<<"old brc= ";
  for (auto &brc : node->brcs) {
	cout << (brc.br_dir ? "true" : "false") << "(" << brc.edge.first << "," << brc.edge.second << ")" << " ";
  }
  cout << endl;
  node->brcs.clear();
  ++num_explored_nodes;
  terminateNode(node);

  solver.freeEnv();

  global_gap = (ub - lb_transformed) / ub;

#ifdef PrintOptimalSolution
  printOptIntSol();
#endif
}

void addEnuCol(BbNode *node, CVRP *cvrp);

void CVRP::restoreModel(BbNode *&node) {

  node = new BbNode(500);

  vector<pair<size_t, double>> col_idx;
  vector<double> obj;

  readEnumTree(node, col_idx);
  readColumnPool(node, col_idx, obj);
  obj4_first_col = obj[0];

  safe_solver(solver.loadEnv(nullptr))
  safe_solver(solver.setEnvThreads(NUM_THREADS_LP, true))
  safe_solver(solver.setEnvOutputFlag(0, true))
  safe_solver(solver.setEnvInfUnbdInfo(1, true))
  safe_solver(solver.setEnvMIPGap(MIP_GAP_TOLERANCE, true))
  node->solver.getEnv(&solver);

  cout << SMALL_PHASE_SEPARATION;
  cout << "<Instance  " << file_name << "  Capacity  " << cap << ">" << endl;

  int cnt = dim;
  for (auto &rcc : node->rccs) rcc.idx_rcc = cnt++;
  for (auto &r1c : node->r1cs) r1c.idx_r1c = cnt++;
  for (auto &r1c_mul : node->r1cs_multi) r1c_mul.idx_r1c = cnt++;

  const char *model_name = "CVRP.lp";
  safe_solver(node->solver.newModel(model_name, num_col, obj.data(), nullptr, nullptr, nullptr, nullptr))
  safe_solver(node->solver.optimize())
  vector<double> X(num_col);
  safe_solver(node->solver.getX(0, num_col, X.data()))
  addEnuCol(node, this);
  generateVertex2IndexColsAndEdge2IndexCols(node);

  idx_node = 0;
  node->tree_level = 0;
  node->index = idx_node;
}

void parseRCC(Rcc &rcc, const std::string &line);
void parseR1C(R1c &r1c, const std::string &line);
void parseR1C_multi(R1cMulti &r1c, const std::string &line);
void parseBrC(Brc &brc, const std::string &line);

void CVRP::readEnumTree(BbNode *node, std::vector<pair<size_t, double>> &col_idx) {
  std::ifstream file(tree_path);
  std::string line;

  if (!file.is_open()) {
	cout << "cannot open file: " << tree_path << endl;
	exit(1);
  }

  while (std::getline(file, line)) {
	if (line.find("nd_ind=") != std::string::npos) {
	  std::istringstream ss(line);
	  std::string key;
	  while (ss >> key) {
		if (key == "nd_col=") {
		  ss >> num_col;
		} else if (key == "nd_cstr=") {
		  ss >> num_row;
		} else if (key == "ub=") {
		  ss >> ub;
		} else if (key == "ColPool=") {
		  int all_col;
		  ss >> all_col;
		  all_col += num_col;
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
		if (line.find("R1cMulti") != std::string::npos) {
		  break;
		}
		R1c r1c;
		parseR1C(r1c, line);
		node->r1cs.emplace_back(std::move(r1c));
	  }
	  while (std::getline(file, line)) {
		if (line.find("Brc") != std::string::npos) {
		  break;
		}
		R1cMulti r1c_multi;
		parseR1C_multi(r1c_multi, line);
		node->r1cs_multi.emplace_back(std::move(r1c_multi));
	  }
	  while (std::getline(file, line)) {
		Brc brc;
		parseBrC(brc, line);
		node->brcs.emplace_back(std::move(brc));
	  }
	}
  }

  file.close();
  cout << "read enumeration tree successfully!" << endl;
}

void CVRP::readColumnPool(BbNode *node, const std::vector<std::pair<size_t, double>> &col_idx, std::vector<double> &objs) {
  unordered_map<size_t, double> col_set;
  col_set.reserve(col_idx.size());
  for (int i = 0; i < col_idx.size(); ++i) {
	if(col_set.find(col_idx[i].first)==col_set.end()) {
	  col_set[col_idx[i].first] = col_idx[i].second;
	}else{
	  if(col_set[col_idx[i].first]>col_idx[i].second){
		col_set[col_idx[i].first] = col_idx[i].second;
	  }
	}
  }
  node->size_enumeration_col_pool = col_set.size() - num_col;
  cout << "size of enumeration column pool: " << node->size_enumeration_col_pool << endl;
  node->cost_for_columns_in_enumeration_column_pool.resize(node->size_enumeration_col_pool);
  node->index_columns_in_enumeration_column_pool.resize(node->size_enumeration_col_pool);
  std::ifstream file(col_pool_path);
  std::string line;
  size_t lineNumber = 0;
  int cnt = 0;
  int numCol=0, numPool=0;
  if (!file.is_open()) {
	cout << "cannot open file: " << col_pool_path << endl;
	exit(1);
  }

  vector<size_t> tmp_col;
  tmp_col.reserve(Config::MaxNumRouteInEnumeration);
  objs.resize(num_col);
  double value;
  int num_col_skips = 0;
  int num;
  while (std::getline(file, line)) {
	if (col_set.find(lineNumber) != col_set.end()) {
	  std::stringstream ss(line);
	  char delimiter;
	  ss >> value;
	  if (value < col_set[lineNumber] - TOLERANCE) {
		++num_col_skips;
		goto HERE;
	  }
	  if (cnt >= num_col) {
		node->cost_for_columns_in_enumeration_column_pool[cnt - num_col] = value;
		node->index_columns_in_enumeration_column_pool[cnt - num_col] = pool_beg4_pricing;
	  } else {
		objs[cnt] = value;
		tmp_col.emplace_back(pool_beg4_pricing);
	  }
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
  }
  if (cnt>num_col) {
	numPool=cnt-num_col;
	numCol=num_col;
	node->size_enumeration_col_pool=numPool;
  }else{
	numPool=0;
	numCol=cnt;
	node->size_enumeration_col_pool=0;
  }
  node->index_columns.assign(tmp_col.begin(), tmp_col.end());
  num_col=numCol;
  objs.resize(numCol);
  node->cost_for_columns_in_enumeration_column_pool.conservativeResize(numPool);
  node->index_columns_in_enumeration_column_pool.conservativeResize(numPool);

  cout<< "num_col_skips= "<<num_col_skips<<endl;
  cout << "cnt= " << cnt << endl;
  cout << "node->size_enumeration_col_pool= " << node->size_enumeration_col_pool << endl;

  if (pool_beg4_pricing >= mem4_pricing) {
	throw runtime_error("the pricing pool is full!");
  }

  file.close();
  cout << "read column pool successfully!" << endl;
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
	r1c.info_r1c.emplace_back(num);
  }
}

void parseR1C_multi(R1cMulti &r1c_multi, const std::string &line) {
  std::stringstream ss(line);
  std::string token;
  int num;

  ss >> r1c_multi.rhs;
  ss >> token;
  ss >> r1c_multi.idx_r1c;
  ss >> token;

  auto &vec = r1c_multi.info_r1c.first;
  while (ss >> token) {
	if (token == "|") {
	  break;
	}
	vec.emplace_back(stoi(token));
  }
  ss >> r1c_multi.info_r1c.second;
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

void addEnuCol(BbNode *node, CVRP *cvrp) {
  int size_pool = cvrp->num_col, curr_node, num_row = cvrp->num_row, real_dim = cvrp->real_dim;
  auto &map_rank1_multiplier = cvrp->map_rank1_multiplier;
  int *col_pool4_pricing = cvrp->col_pool4_pricing;
  auto &ptr = node->index_columns;
  sparseRowMatrixXd mat(num_row, size_pool);

  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(num_row) * double(size_pool) * 0.1));

  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> map;
  for (int i = 0; i < size_pool; ++i) {
	int past_node = 0;
	for (auto j = ptr[i] + 1;; ++j) {
	  curr_node = col_pool4_pricing[j];
	  if (past_node < curr_node)
		map[{past_node, curr_node}].emplace_back(i);
	  else map[{curr_node, past_node}].emplace_back(i);
	  if (!curr_node) break;
	  triplets.emplace_back(curr_node - 1, i, 1);
	  past_node = curr_node;
	}
  }

  sparseRowMatrixXd tmpMat(real_dim, size_pool);
  tmpMat.setFromTriplets(triplets.begin(), triplets.end());

  for (int i = 0; i < size_pool; ++i) triplets.emplace_back(real_dim, i, 1);

  sparseRowMatrixXd sum(1, size_pool);
  unordered_map<int, double> tmp;
  tmp.reserve(size_pool);

  for (auto &rcc : node->rccs) {
	tmp.clear();
	if (rcc.form_rcc) {
	  auto &info = rcc.info_rcc_customer;
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
	  auto &infoRccCustomer = rcc.info_rcc_customer;
	  auto &infoRccOutsideCustomer = rcc.info_rcc_outside_customer;
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
	int row = rcc.idx_rcc;
	for (auto &it : tmp) {
	  if (abs(it.second) > TOLERANCE) {
		triplets.emplace_back(row, it.first, it.second);
	  }
	}
  }
  for (auto &r1c : node->r1cs) {
	sum.setZero();
	auto &info = r1c.info_r1c;
	for (auto j : info) {
	  sum += tmpMat.row(j - 1);
	}
	sum /= 2;
	int row = r1c.idx_r1c;
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
	  int val_ = int(it.value() + TOLERANCE);
	  if (val_) {
		triplets.emplace_back(row, it.col(), val_);
	  }
	}
  }
  for (auto &r1c : node->r1cs_multi) {
	sum.setZero();
	auto &info = r1c.info_r1c;
	const auto &plan = map_rank1_multiplier[(int)info.first.size()][info.second];
	const auto &multi = get<0>(plan);
	int denominator = get<1>(plan);
	int count = 0;
	for (auto &j : info.first) {
	  sum += tmpMat.row(j - 1) * multi[count++];
	}
	sum /= denominator;
	int row = r1c.idx_r1c;
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sum, 0); it; ++it) {
	  int val_ = int(it.value() + TOLERANCE);
	  if (val_) {
		triplets.emplace_back(row, it.col(), val_);
	  }
	}
  }
  mat.setFromTriplets(triplets.begin(), triplets.end());

  vector<char> sense(num_row, SOLVER_LESS_EQUAL);
  vector<double> rhs(num_row);

  vector<size_t> solver_beg(num_row + 1);
  vector<int> solver_ind;
  vector<double> solver_val;
  size_t numRe=size_pool*cvrp->aver_route_length;
  solver_ind.reserve(numRe);
  solver_val.reserve(numRe);

  for (int i = 0; i < num_row; ++i) {
	solver_beg[i] = solver_ind.size();
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, i); it; ++it) {
	  solver_ind.emplace_back(it.col());
	  solver_val.emplace_back(it.value());
	}
  }

  int i = 0;
  for (; i < real_dim; ++i) {
	sense[i] = SOLVER_EQUAL;
	rhs[i] = 1;
  }
  sense[i] = SOLVER_GREATER_EQUAL;
  rhs[i] = cvrp->num_vehicle;

  for (auto &rcc : node->rccs) rhs[rcc.idx_rcc] = rcc.rhs;
  for (auto &r1c : node->r1cs) rhs[r1c.idx_r1c] = r1c.rhs;
  for (auto &r1c_mul : node->r1cs_multi) rhs[r1c_mul.idx_r1c] = r1c_mul.rhs;

  safe_solver(node->solver.XaddContraints(num_row,
											 solver_ind.size(),
											 solver_beg.data(),
											 solver_ind.data(),
											 solver_val.data(),
											 sense.data(),
											 rhs.data(),
											 nullptr))
  vector<int> cind(num_row), vind(num_row, 0);
  iota(cind.begin(), cind.end(), 0);
  safe_solver(node->solver.changeCoeffs(num_row, cind.data(), vind.data(), rhs.data()))
  safe_solver(node->solver.optimize())
  safe_solver(node->solver.getObjVal(&node->value))
}

#endif
