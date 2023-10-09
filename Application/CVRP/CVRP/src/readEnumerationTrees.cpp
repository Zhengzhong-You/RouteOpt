//
// Created by You, Zhengzhong on 8/15/23.
//

#include "CVRP.hpp"

using namespace std;
using namespace chrono;
using namespace Eigen;

#ifdef readEnumerationTrees

void CVRP::solveEnuTree() {

  cout << "\n<solveEnuTree>\n" << endl;

  cout << "turn on only read enumeration mode!" << endl;

  if_only_read_enumerationTree = true;

  lateProcessing();
  GloBeg = std::chrono::high_resolution_clock::now();
//  buildModel();
  BBNODE *node = nullptr;
  restoreModel(node);

  node->Val = 0;//reset
  cout << "reset node->Val= " << node->Val << endl;
  ++NumExploredNodes;
  terminateNode(node);

  Solver.SOLVERfreeenv();

  GlobalGap = (UB - LB_transformed) / UB;

#ifdef PrintOptimalSolution
  printOptIntSol();
#endif
}

void addEnuCol(BBNODE *node, CVRP *cvrp);

void CVRP::restoreModel(BBNODE *&node) {
  //initialize a binary tree p0

  node = new BBNODE(MaxNumEdge);

  vector<size_t> col_idx;
  vector<double> obj;

  readEnuTree(node, col_idx);
  readColPool(node, col_idx, obj);
  Obj4FirstCol = obj[0];

  //build the initial model
  safe_solver(Solver.SOLVERloadenv(nullptr))
  safe_solver(Solver.SOLVERsetenvThreads(NUM_THREAD_LP, true))
  safe_solver(Solver.SOLVERsetenvOutputFlag(0, true))
  safe_solver(Solver.SOLVERsetenvInfUnbdInfo(1, true))
  safe_solver(Solver.SOLVERsetenvMIPGap(MIPGap, true))
  node->solver.SOLVERgetenv(&Solver);

  cout << SMALL_PHASE_SEPARATION;
  cout << "<Instance  " << FileName << "  Capacity  " << Cap << ">" << endl;

  //rearrange the rows
  int cnt = Dim;
  for (auto &rcc : node->RCCs) rcc.IdxRCC = cnt++;
  for (auto &r1c : node->R1Cs) r1c.IdxR1C = cnt++;
  for (auto &r1c_mul : node->R1Cs_multi) r1c_mul.IdxR1C = cnt++;

  const char *model_name = "CVRP.lp";
  safe_solver(node->solver.SOLVERnewmodel(model_name, NumCol, obj.data(), nullptr, nullptr, nullptr, nullptr))
  safe_solver(node->solver.SOLVERoptimize())
  safe_solver(node->solver.SOLVERgetX(0, NumCol, X))
  addEnuCol(node, this);
  deleteBrCs_spec4readEnuTree(node);
  generateVertex2IdxCols_N_Edge2IdxCols(node);

  IdxNode = 0;
  node->TreeLevel = 0;
  node->Idx = IdxNode;
}

void CVRP::deleteBrCs_spec4readEnuTree(BBNODE *const node) {
  if (node->BrCs.empty()) return;
  set<int> delete_col;
  unordered_map<pair<int, int>, vector<int>, PairHasher> col_map;//record the index;
  unordered_map<int, vector<int>> i_map;
  for (int i = 1; i < NumCol; ++i) {
	int past = 0;
	for (auto j = node->IdxCols[i] + 1;; ++j) {
	  int curr_node = ColPool4Pricing[j];
	  i_map[curr_node].emplace_back(i);
	  if (!curr_node) break;
	  if (past < curr_node) {
		col_map[{past, curr_node}].emplace_back(i);
	  } else col_map[{curr_node, past}].emplace_back(i);
	  past = curr_node;
	}
  }
  for (auto it = col_map.begin(); it != col_map.end(); ++it) {
	auto &vec = it->second;
	std::sort(vec.begin(), vec.end());
  }
  vector<int> ai_col(NumCol), aj_col(NumCol);
  for (auto &brc : node->BrCs) {
	int ai = brc.Edge.first, aj = brc.Edge.second;
	if (brc.BrDir) {
	  //in lp
	  std::fill(ai_col.begin(), ai_col.end(), 0);
	  std::fill(aj_col.begin(), aj_col.end(), 0);
	  if (ai) {
		for (auto i : i_map[ai])ai_col[i] = 1;
	  }
	  for (auto j : i_map[aj])aj_col[j] = 1;
	  if (col_map.find({ai, aj}) != col_map.end()) {
		auto &cols = col_map[{ai, aj}];
		for (int j = 0; j < cols[0]; ++j) if (aj_col[j] || ai_col[j]) delete_col.insert(j);
		for (size_t i = 1; i < cols.size(); ++i)
		  for (int j = cols[i - 1] + 1; j < cols[i]; ++j)
			if (aj_col[j] || ai_col[j])delete_col.insert(j);
		for (int j = cols.back() + 1; j < NumCol; ++j)if (aj_col[j] || ai_col[j])delete_col.insert(j);
	  } else {
		for (int j = 0; j < NumCol; ++j) if (aj_col[j] || ai_col[j]) delete_col.insert(j);
	  }
	} else {//delete all cstrs that coefficient == 1
	  if (col_map.find({ai, aj}) != col_map.end())
		for (auto i : col_map[{ai, aj}])delete_col.insert(i);
	}
  }

  //in lp
  std::fill(ai_col.begin(), ai_col.end(), 0);
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
  cout << "delete " << len << " cols due to replacing the old columns" << endl;
}

void parseRCC(RCC &rcc, const std::string &line);
void parseR1C(R1C &r1c, const std::string &line);
void parseR1C_multi(R1C_multi &r1c, const std::string &line);
void parseBrC(BrC &brc, const std::string &line);

void CVRP::readEnuTree(BBNODE *node, std::vector<size_t> &col_idx) {
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
		  ss >> NumCol;
		} else if (key == "nd_cstr=") {
		  ss >> NumRow;
		} else if (key == "ub=") {
		  ss >> UB;
		} else if (key == "ColPool=") {
		  int all_col;
		  ss >> all_col;
		  all_col += NumCol;
		  col_idx.resize(all_col);
		}
	  }
	} else if (line.find("ColPool:") != std::string::npos) {
	  std::istringstream ss(line);
	  int value;
	  int cnt = 0;
	  string key;
	  ss >> key;
	  while (ss >> value) {
		col_idx[cnt++] = value;
	  }
	} else if (line.find("RCCs") != std::string::npos) {
	  while (std::getline(file, line)) {
		if (line.find("R1Cs") != std::string::npos) {
		  break;
		}
		RCC rcc;
		parseRCC(rcc, line);
		node->RCCs.emplace_back(std::move(rcc));
	  }
	  while (std::getline(file, line)) {
		if (line.find("R1C_multi") != std::string::npos) {
		  break;
		}
		R1C r1c;
		parseR1C(r1c, line);
		node->R1Cs.emplace_back(std::move(r1c));
	  }
	  while (std::getline(file, line)) {
		if (line.find("BrC") != std::string::npos) {
		  break;
		}
		R1C_multi r1c_multi;
		parseR1C_multi(r1c_multi, line);
		node->R1Cs_multi.emplace_back(std::move(r1c_multi));
	  }
	  while (std::getline(file, line)) {
		BrC brc;
		parseBrC(brc, line);
		node->BrCs.emplace_back(std::move(brc));
	  }
	}
  }

  file.close();
  cout << "read enumeration tree successfully!" << endl;
}

void CVRP::readColPool(BBNODE *node, const std::vector<size_t> &col_idx, std::vector<double> &objs) {
  unordered_set<int> col_set(col_idx.begin(), col_idx.end());

  node->SizeEnuColPool = col_set.size() - NumCol;
  cout << "size of enumeration column pool: " << node->SizeEnuColPool << endl;
  node->Cost4ColsInEnuColPool.resize(node->SizeEnuColPool);
  node->IdxColsInEnuColPool.resize(node->SizeEnuColPool);

  std::ifstream file(colPool_path);
  std::string line;
  size_t lineNumber = 0;
  int cnt = 0;

  if (!file.is_open()) {
	cout << "cannot open file: " << colPool_path << endl;
	exit(1);
  }
  objs.resize(NumCol);
  while (std::getline(file, line)) {
	if (col_set.find(lineNumber) != col_set.end()) {
	  std::stringstream ss(line);
	  char delimiter;
	  if (cnt >= NumCol) {
		ss >> node->Cost4ColsInEnuColPool[cnt - NumCol];
		node->IdxColsInEnuColPool[cnt - NumCol] = PoolBeg4Pricing;
	  } else {
		ss >> objs[cnt];
		node->IdxCols[cnt] = PoolBeg4Pricing;
	  }
	  ss >> delimiter;
	  int num;
	  ColPool4Pricing[PoolBeg4Pricing++] = 0;
	  while (ss >> num) {
		ColPool4Pricing[PoolBeg4Pricing++] = num;
	  }
	  ColPool4Pricing[PoolBeg4Pricing++] = 0;
	  ++cnt;
	}
	++lineNumber;
  }

  cout << "cnt= " << cnt << endl;
  cout << "node->SizeEnuColPool= " << node->SizeEnuColPool << endl;

  if (PoolBeg4Pricing >= Mem4Pricing) {
	throw runtime_error("the pricing pool is full!");
  }

  file.close();
  cout << "read column pool successfully!" << endl;
}

void parseRCC(RCC &rcc, const std::string &line) {
  std::stringstream ss(line);
  std::string token;

  int num;
  ss >> rcc.FormRCC;
  ss >> token;
  ss >> rcc.RHS;
  ss >> token;
  ss >> rcc.IdxRCC;
  ss >> token;
  while (ss >> token) {
	if (token == "|") {
	  break;
	}
	rcc.InfoRCCCustomer.emplace_back(stoi(token));
  }
  while (ss >> num) {
	rcc.InfoRCCOutsideCustomer.emplace_back(num);
  }
}

void parseR1C(R1C &r1c, const std::string &line) {
  std::stringstream ss(line);
  std::string token;
  int num;
  ss >> r1c.RHS;
  ss >> token;
  ss >> r1c.IdxR1C;
  ss >> token;
  while (ss >> num) {
	r1c.InfoR1C.emplace_back(num);
  }
}

void parseR1C_multi(R1C_multi &r1c_multi, const std::string &line) {
  std::stringstream ss(line);
  std::string token;
  int num;

  ss >> r1c_multi.RHS;
  ss >> token;
  ss >> r1c_multi.IdxR1C;
  ss >> token;

  auto &vec = r1c_multi.InfoR1C.first;
  while (ss >> token) {
	if (token == "|") {
	  break;
	}
	vec.emplace_back(stoi(token));
  }
  ss >> r1c_multi.InfoR1C.second;
}

void parseBrC(BrC &brc, const std::string &line) {
  std::stringstream ss(line);
  std::string token;
  int num;

  ss >> brc.Edge.first;
  ss >> token;
  ss >> brc.Edge.second;
  ss >> token;
  ss >> brc.BrDir;
  brc.IdxBrC = -1;
}

void addEnuCol(BBNODE *node, CVRP *cvrp) {
  int size_pool = cvrp->NumCol, curr_node, NumRow = cvrp->NumRow, RealDim = cvrp->RealDim;
  auto &map_rank1_multiplier = cvrp->map_rank1_multiplier;
  int *ColPool4Pricing = cvrp->ColPool4Pricing;
  auto &ptr = node->IdxCols;
  sparseRowMatrixXd mat(NumRow, size_pool);

  vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(int(double(NumRow) * double(size_pool) * 0.1));

  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> map;
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
	const auto &plan = map_rank1_multiplier[(int)info.first.size()][info.second];
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

  vector<char> sense(NumRow, SOLVER_LESS_EQUAL);
  vector<double> rhs(NumRow);

  auto solver_beg = cvrp->solver_beg;
  auto solver_ind = cvrp->solver_ind;
  auto solver_val = cvrp->solver_val;

  int nzcnt = 0;
  for (int i = 0; i < NumRow; ++i) {
	solver_beg[i] = nzcnt;
	for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, i); it; ++it) {
	  solver_ind[nzcnt] = it.col();
	  solver_val[nzcnt] = it.value();
	  ++nzcnt;
	}
  }

  int i = 0;
  for (; i < RealDim; ++i) {
	sense[i] = SOLVER_EQUAL;
	rhs[i] = 1;
  }
  sense[i] = SOLVER_GREATER_EQUAL;
  rhs[i] = cvrp->K;

  for (auto &rcc : node->RCCs) rhs[rcc.IdxRCC] = rcc.RHS;
  for (auto &r1c : node->R1Cs) rhs[r1c.IdxR1C] = r1c.RHS;
  for (auto &r1c_mul : node->R1Cs_multi) rhs[r1c_mul.IdxR1C] = r1c_mul.RHS;

  safe_solver(node->solver.SOLVERXaddconstrs(NumRow,
											 nzcnt,
											 solver_beg,
											 solver_ind,
											 solver_val,
											 sense.data(),
											 rhs.data(),
											 nullptr))
//change the coefficient of the first column
  vector<int> cind(NumRow), vind(NumRow, 0);
  iota(cind.begin(), cind.end(), 0);
  safe_solver(node->solver.SOLVERchgcoeffs(NumRow, cind.data(), vind.data(), rhs.data()))
  safe_solver(node->solver.SOLVERoptimize())
  safe_solver(node->solver.SOLVERgetObjVal(&node->Val))
}

#endif
