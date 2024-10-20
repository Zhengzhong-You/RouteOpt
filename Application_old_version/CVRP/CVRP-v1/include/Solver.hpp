
#ifndef CVRP_SOLVER_HPP
#define CVRP_SOLVER_HPP

//#define SOLVER_STATISTICS

#define SOLVERTYPE 0

#if SOLVERTYPE == 0

#include "gurobi_c.h"
#include "gurobi_c++.h"

using SOLVERmodel = GRBmodel;
using SOLVERenv = GRBenv;
#define SOLVER_EQUAL GRB_EQUAL
#define SOLVER_GREATER_EQUAL GRB_GREATER_EQUAL
#define SOLVER_LESS_EQUAL GRB_LESS_EQUAL
#define SOLVER_BINARY GRB_BINARY
#define SOLVER_INTEGER GRB_INTEGER
#define SOLVER_CONTINUOUS GRB_CONTINUOUS
#define SOLVER_INFEASIBLE 3
#define SOLVER_INF_OR_UNBD 4
#define SOLVER_UNBOUNDED 5
#define SOLVER_PRIMAL_SIMPLEX 0
#define SOLVER_DUAL_SIMPLEX 1
#define SOLVER_BARRIER 2
#define SOLVER_CROSSOVER_CLOSE 0
#define SOLVER_MAX_SENSE (-1)
#define SOLVER_TIME_LIMIT 9
#define SOLVER_INFINITY GRB_INFINITY
#elif SOLVERTYPE == 1
# SOLVER_INFINITY 1e100
#elif SOLVERTYPE == 2

#include "Milp\LinearModel\LinearModel.hpp"

#define SOLVER_EQUAL 0
#define SOLVER_GREATER_EQUAL 2
#define SOLVER_LESS_EQUAL 1
#define SOLVER_BINARY 2
#define SOLVER_INTEGER 1
#define SOLVER_CONTINUOUS 0
#define VARTYPE Mindo::VariableType
#define SOLVER_INFEASIBLE 1
#define SOLVER_UNBOUNDED 2
#define SOLVER_INF_OR_UNBD 3
#define SOLVER_PRIMAL_SIMPLEX 0
#define SOLVER_DUAL_SIMPLEX 1
#define SOLVER_MAX_SENSE (-1)
#define SOLVER_TIME_LIMIT 9

#endif

class Solver {
 public://for test here

#if SOLVERTYPE == 0
  SOLVERmodel *model;
  SOLVERenv *env;
#elif SOLVERTYPE == 2

  Mindo::LinearModel model;

#endif
 public:

  int getSolFromPool(int i, int first, int len, double *values) const;
  int getObjFromPool(int i, double *valueP) const;
  int getSolCount(int *valueP) const;

  int setVBasis(int first, int len, int *vbasis);
  int setCBasis(int first, int len, int *cbasis);

  int getVBasis(int first, int len, int *vbasis) const;

  int getCBasis(int first, int len, int *cbasis) const;

  int getSlack(int first, int len, double *values) const;

  int getDual(int first, int len, double *values) const;

  int getNumRow(int *valueP) const;

  int getNumCol(int *valueP) const;

  int delConstraints(int len, int *cind);

  [[nodiscard]] int reoptimize();

  [[nodiscard]] int optimize();

  [[nodiscard]] int mipOptimize();

  int setVTypeArray(int first, int len, char *newvalues);

#if SOLVERTYPE == 0
  [[nodiscard]] SOLVERmodel *copyModel() const;
  [[nodiscard]] SOLVERenv *passEnv() const;
#elif SOLVERTYPE == 2
  [[nodiscard]] Solver* copyModel() const;
#endif

  int freeModel();

  int newModel(const char *Pname, int numvars,
			   double *obj, double *lb, double *ub, char *vtype,
			   char **varnames);

  int addConstraints(int numconstrs, int numnz,
					 int *cbeg, int *cind, double *cval,
					 char *sense, double *rhs, char **constrnames);

  [[nodiscard]] int updateModel();

  int write(const char *filename) const;

  int delVars(int len, int *ind);

  void freeEnv() const;

  int getObjVal(double *valueP) const;

  int getX(int first, int len, double *values) const;

  int getObj(int first, int len, double *values) const;

  int getRhs(int first, int len, double *values) const;

  void getSolver(Solver *solver);

  void getEnv(Solver *solver);

  int getRC(int first, int len, double *values) const;

  int changeObj(int first, int len, double *values);

  int changeCoeffs(int cnt, int *cind, int *vind, double *val);

  int XchangeCoeffs(size_t cnt, int *cind, int *vind, double *val);

  int setRhs(int first, int len, char *sense, double *newvalues);

  int setRhs(int first, int len, double *newvalues);

  int setColLower(int col, double value);

  int setColUpper(int col, double value);

  int removeColLower(int col);

  int removeColUpper(int col);

  int addConstraint(int numnz, int *cind, double *cval,
					char sense, double rhs, const char *constrname);

  int addVars(int numvars, int numnz,
			  int *vbeg, int *vind, double *vval,
			  double *obj, double *lb, double *ub, char *vtype,
			  char **varnames);

  int addVar(int numnz, int *vind, double *vval,
			 double obj, double lb, double ub, char vtype,
			 const char *varname);

  int XaddVars(int numvars, size_t numnz,
			   size_t *vbeg, int *vind, double *vval,
			   double *obj, double *lb, double *ub, char *vtype,
			   char **varnames);

  int XaddContraints(int numconstrs, size_t numnz,
					 size_t *cbeg, int *cind, double *cval,
					 char *sense, double *rhs, char **constrnames);

  int XgetConstraints(size_t *numnzP, size_t *cbeg,
					  int *cind, double *cval, int start, int len) const;

  int XgetVars(size_t *numnzP, size_t *vbeg,
			   int *vind, double *vval, int start, int len) const;

  int getConstraints(int *numnzP, int *cbeg,
					 int *cind, double *cval, int start, int len) const;

  int getVars(int *numnzP, int *vbeg, int *vind, double *vval, int start, int len) const;

  int getCoeff(int row, int col, double &value) const;

  int getStatus(int *valueP) const;

  [[nodiscard]] int setEnvMethod(int value);

  [[nodiscard]] int setEnvCrossOver(int value);

  [[nodiscard]] int getEnvMethod(int *value) const;

  int loadEnv(const char *logfilename);

  [[nodiscard]] int setEnvThreads(int value, bool if_model_free);

  [[nodiscard]] int setEnvOutputFlag(int value, bool if_model_free);

  [[nodiscard]] int setEnvInfUnbdInfo(int value, bool if_model_free);

  [[nodiscard]] int setEnvMIPGap(double value, bool if_model_free);

  [[nodiscard]] int setEnvCutoff(double value);

  [[nodiscard]] int setEnvTimeLimit(double value);

  [[nodiscard]] int setModelSense(int value);

  [[nodiscard]] int getSense(int first, int len, char *sense);
};

#ifdef SOLVER_STATISTICS
class STATISTICS {
 public:
  static int calls_SOLVERgetSlack;
  static int calls_SOLVERgetDual;
  static int calls_SOLVERgetNumRow;
  static int calls_SOLVERgetNumCol;
  static int calls_SOLVERdelconstrs;
  static int calls_SOLVERreoptimize;
  static int calls_SOLVERoptimize;
  static int calls_SOLVERsetVTypearray;
  static int calls_SOLVERsetSensearray;
  static int calls_SOLVERcopymodel;
  static int calls_SOLVERpassenv;
  static int calls_SOLVERfreemodel;
  static int calls_SOLVERnewmodel;
  static int calls_SOLVERaddconstrs;
  static int calls_SOLVERupdatemodel;
  static int calls_SOLVERwrite;
  static int calls_SOLVERdelvars;
  static int calls_SOLVERfreeenv;
  static int calls_SOLVERgetObjVal;
  static int calls_SOLVERgetX;
  static int calls_SOLVERgetsolver;
  static int calls_SOLVERgetenv;
  static int calls_SOLVERgetRC;
  static int calls_SOLVERchgcoeffs;
  static int calls_SOLVERXchgcoeffs;
  static int calls_SOLVERsetRHS;
  static int calls_SOLVERaddconstr;
  static int calls_SOLVERaddvars;
  static int calls_SOLVERXaddvars;
  static int calls_SOLVERXaddconstrs;
  static int calls_SOLVERXgetconstrs;
  static int calls_SOLVERgetconstrs;
  static int calls_SOLVERgetStatus;
  static int calls_SOLVERsetenvMethod;
  static int calls_SOLVERloadenv;
  static int calls_SOLVERsetenvThreads;
  static int calls_SOLVERsetenvOutputFlag;
  static int calls_SOLVERsetenvInfUnbdInfo;
  static int calls_SOLVERsetenvCutoff;
  static int calls_SOLVERsetenvTimeLimit;
  static int calls_SOLVERsetModelSense;
  static int calls_SOLVERsetVBasis;
  static int calls_SOLVERsetCBasis;
  static int calls_SOLVERgetVBasis;
  static int calls_SOLVERgetCBasis;
  static int calls_SOLVERXgetvars;
  static int calls_SOLVERgetObj;
  static int calls_SOLVERgetRhs;
  static int calls_SOLVERMIPoptimize;
  static int calls_SOLVERsetenvMIPGap;
  static int calls_SOLVERgetSense;
  static int calls_SOLVERchgObj;
  static int calls_SOLVERchgRhs;
  static int calls_SOLVERgetenvMethod;
  static double time_SOLVERgetSlack;
  static double time_SOLVERgetDual;
  static double time_SOLVERgetNumRow;
  static double time_SOLVERgetNumCol;
  static double time_SOLVERdelconstrs;
  static double time_SOLVERreoptimize;
  static double time_SOLVERoptimize;
  static double time_SOLVERsetVTypearray;
  static double time_SOLVERsetSensearray;
  static double time_SOLVERcopymodel;
  static double time_SOLVERpassenv;
  static double time_SOLVERfreemodel;
  static double time_SOLVERnewmodel;
  static double time_SOLVERaddconstrs;
  static double time_SOLVERupdatemodel;
  static double time_SOLVERwrite;
  static double time_SOLVERdelvars;
  static double time_SOLVERfreeenv;
  static double time_SOLVERgetObjVal;
  static double time_SOLVERgetX;
  static double time_SOLVERgetsolver;
  static double time_SOLVERgetenv;
  static double time_SOLVERgetRC;
  static double time_SOLVERchgcoeffs;
  static double time_SOLVERXchgcoeffs;
  static double time_SOLVERsetRHS;
  static double time_SOLVERaddconstr;
  static double time_SOLVERaddvars;
  static double time_SOLVERXaddvars;
  static double time_SOLVERXaddconstrs;
  static double time_SOLVERXgetconstrs;
  static double time_SOLVERgetconstrs;
  static double time_SOLVERgetStatus;
  static double time_SOLVERsetenvMethod;
  static double time_SOLVERloadenv;
  static double time_SOLVERsetenvThreads;
  static double time_SOLVERsetenvOutputFlag;
  static double time_SOLVERsetenvInfUnbdInfo;
  static double time_SOLVERsetenvCutoff;
  static double time_SOLVERsetenvTimeLimit;
  static double time_SOLVERsetModelSense;
  static double time_SOLVERsetVBasis;
  static double time_SOLVERsetCBasis;
  static double time_SOLVERgetVBasis;
  static double time_SOLVERgetCBasis;
  static double time_SOLVERXgetvars;
  static double time_SOLVERgetObj;
  static double time_SOLVERgetRhs;
  static double time_SOLVERMIPoptimize;
  static double time_SOLVERsetenvMIPGap;
  static double time_SOLVERgetSense;
  static double time_SOLVERchgObj;
  static double time_SOLVERchgRhs;
  static double time_SOLVERgetenvMethod;
  static void giveStatistics(const std::string &StatisticsName);
};
#endif

#endif //CVRP_SOLVER_HPP
