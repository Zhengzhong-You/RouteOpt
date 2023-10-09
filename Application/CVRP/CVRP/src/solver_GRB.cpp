//
// Created by Zhengzhong You on 10/21/22.
//
#include "solver.hpp"
#include <chrono>
using namespace std::chrono;

#if SOLVERTYPE == 0

int SOLVER::SOLVERgetSlack(int first, int len, double *const values) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetSlack;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattrarray(model, "Slack", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetSlack += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetDual(int first, int len, double *const values) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetDual;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattrarray(model, "Pi", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetDual += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetNumRow(int *const valueP) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetNumRow;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetintattr(model, "NumConstrs", valueP);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetNumRow += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetNumCol(int *const valueP) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetNumCol;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetintattr(model, "NumVars", valueP);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetNumCol += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetStatus(int *const valueP) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetStatus;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetintattr(model, "Status", valueP);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetStatus += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERdelconstrs(int len, int *const cind) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERdelconstrs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBdelconstrs(model, len, cind);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERdelconstrs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERreoptimize() {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERreoptimize;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBoptimize(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERreoptimize += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERoptimize() {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERoptimize;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBoptimize(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERoptimize += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERMIPoptimize() {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERMIPoptimize;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBoptimize(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERMIPoptimize += (double) duration.count() * 1e-3;
#endif
  return error;
}

SOLVERmodel *SOLVER::SOLVERcopymodel() const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERcopymodel;
  auto beg = high_resolution_clock::now();
#endif
  auto newModel = GRBcopymodel(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERcopymodel += (double) duration.count() * 1e-3;
#endif
  return newModel;
}

int SOLVER::SOLVERfreemodel() const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERfreemodel;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBfreemodel(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERfreemodel += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERloadenv(const char *const logfilename) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERloadenv;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBloadenv(&env, logfilename);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERloadenv += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetenvThreads(int value, bool if_model_free) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvThreads;
  auto beg = high_resolution_clock::now();
#endif
  int error;
  if (if_model_free)
    error = GRBsetintparam(env, "Threads", value);
  else
    error = GRBsetintparam(GRBgetenv(model), "Threads", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvThreads += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetenvOutputFlag(int value, bool if_model_free) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvOutputFlag;
  auto beg = high_resolution_clock::now();
#endif
  int error;
  if (if_model_free)
    error = GRBsetintparam(env, "OutputFlag", value);
  else
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvOutputFlag += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetenvInfUnbdInfo(int value, bool if_model_free) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvInfUnbdInfo;
  auto beg = high_resolution_clock::now();
#endif
  int error;
  if (if_model_free)
    error = GRBsetintparam(env, "InfUnbdInfo", value);
  else
    error = GRBsetintparam(GRBgetenv(model), "InfUnbdInfo", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvInfUnbdInfo += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetenvMIPGap(double value, bool if_model_free) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvMIPGap;
  auto beg = high_resolution_clock::now();
#endif
  int error;
  if (if_model_free)
    error = GRBsetdblparam(env, "MIPGap", value);
  else
    error = GRBsetdblparam(GRBgetenv(model), "MIPGap", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvMIPGap += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERnewmodel(const char *const Pname, int numvars,
                           double *const obj, double *const lb, double *const ub, char *const vtype,
                           char **const varnames) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERnewmodel;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBnewmodel(env, &model, Pname, numvars, obj, lb, ub, vtype, varnames);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERnewmodel += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERaddconstrs(int numconstrs, int numnz,
                             int *const cbeg, int *const cind, double *const cval,
                             char *const sense, double *const rhs, char **const constrnames) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERaddconstrs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBaddconstrs(model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERaddconstrs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERupdatemodel() {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERupdatemodel;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBupdatemodel(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERupdatemodel += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERwrite(const char *filename) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERwrite;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBwrite(model, filename);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERwrite += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERdelvars(int len, int *const ind) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERdelvars;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBdelvars(model, len, ind);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERdelvars += (double) duration.count() * 1e-3;
#endif
  return error;
}

void SOLVER::SOLVERfreeenv() const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERfreeenv;
  auto beg = high_resolution_clock::now();
#endif
  GRBfreeenv(env);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERfreeenv += (double) duration.count() * 1e-3;
#endif
}

int SOLVER::SOLVERgetObjVal(double *const valueP) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetObjVal;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattr(model, "ObjVal", valueP);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetObjVal += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetObj(int first, int len, double *const values) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetObj;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattrarray(model, "Obj", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetObj += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetRhs(int first, int len, double *const values) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetRhs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattrarray(model, "RHS", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetRhs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetX(int first, int len, double *const values) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetX;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattrarray(model, "X", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetX += (double) duration.count() * 1e-3;
#endif
  return error;
}

void SOLVER::SOLVERgetsolver(SOLVER *const solver) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetsolver;
  auto beg = high_resolution_clock::now();
#endif
  model = solver->SOLVERcopymodel();
  env = GRBgetenv(model);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetsolver += (double) duration.count() * 1e-3;
#endif
}

void SOLVER::SOLVERgetenv(SOLVER *const solver) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetenv;
  auto beg = high_resolution_clock::now();
#endif
  env = solver->SOLVERpassenv();
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetenv += (double) duration.count() * 1e-3;
#endif
}

SOLVERenv *SOLVER::SOLVERpassenv() const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERpassenv;
  auto beg = high_resolution_clock::now();
#endif
  auto error = env;
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERpassenv += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetRC(int first, int len, double *const values) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetRC;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetdblattrarray(model, "RC", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetRC += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERchgObj(int first, int len, double *const values) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERchgObj;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetdblattrarray(model, "Obj", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERchgObj += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetRhs(int first, int len, char *sense, double *const values) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetRHS;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetdblattrarray(model, "RHS", first, len, values);
  error += GRBsetcharattrarray(model, "Sense", first, len, sense);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetRHS += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetRhs(int first, int len, double *const values) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetRHS;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetdblattrarray(model, "RHS", first, len, values);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetRHS += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERchgcoeffs(int cnt, int *const cind, int *const vind, double *const val) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERchgcoeffs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBchgcoeffs(model, cnt, cind, vind, val);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERchgcoeffs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERXchgcoeffs(size_t cnt, int *const cind, int *const vind, double *const val) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERXchgcoeffs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBXchgcoeffs(model, cnt, cind, vind, val);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERXchgcoeffs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERaddconstr(int numnz, int *const cind, double *const cval, char sense, double rhs,
                            const char *const constrname) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERaddconstr;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBaddconstr(model, numnz, cind, cval, sense, rhs, constrname);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERaddconstr += (double) duration.count() * 1e-3;
#endif
  return error;
}

int
SOLVER::SOLVERaddvars(int numvars, int numnz, int *const vbeg, int *const vind, double *const vval, double *const obj,
                      double *const lb, double *const ub, char *const vtype, char **const varnames) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERaddvars;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBaddvars(model, numvars, numnz, vbeg, vind, vval, obj, lb, ub, vtype, varnames);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERaddvars += (double) duration.count() * 1e-3;
#endif
  return error;
}

int
SOLVER::SOLVERXaddvars(int numvars,
                       size_t numnz,
                       size_t *const vbeg,
                       int *const vind,
                       double *const vval,
                       double *const obj,
                       double *const lb,
                       double *const ub,
                       char *const vtype,
                       char **const varnames) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERXaddvars;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBXaddvars(model, numvars, numnz, vbeg, vind, vval, obj, lb, ub, vtype, varnames);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERXaddvars += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERXaddconstrs(int numconstrs, size_t numnz,
                              size_t *const cbeg, int *const cind, double *const cval,
                              char *const sense, double *const rhs, char **const constrnames) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERXaddconstrs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBXaddconstrs(model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERXaddconstrs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int
SOLVER::SOLVERgetconstrs(int *const numnzP,
                         int *const cbeg,
                         int *const cind,
                         double *const cval,
                         int start,
                         int len) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetconstrs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetconstrs(model, numnzP, cbeg, cind, cval, start, len);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetconstrs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERXgetconstrs(size_t *const numnzP, size_t *const cbeg,
                              int *const cind, double *const cval, int start, int len) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERXgetconstrs;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBXgetconstrs(model, numnzP, cbeg, cind, cval, start, len);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERXgetconstrs += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERXgetvars(size_t *const numnzP, size_t *const vbeg, int *const vind, double *const vval,
                           int start, int len) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERXgetvars;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBXgetvars(model, numnzP, vbeg, vind, vval, start, len);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERXgetvars += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetenvCutoff(double value) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvCutoff;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetdblparam(GRBgetenv(model), "Cutoff", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvCutoff += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetVTypearray(int first, int len, char *const newvalues) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetVTypearray;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetcharattrarray(model, "VType", first, len, newvalues);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetVTypearray += (double) duration.count() * 1e-3;
#endif
  return error;
}

//int SOLVER::SOLVERsetSensearray(int first, int len, char *newvalues) const {
//#ifdef SOLVER_STATISTICS
//  ++STATISTICS::calls_SOLVERsetSensearray;
//  auto beg = high_resolution_clock::now();
//#endif
//  auto error = GRBsetcharattrarray(model, "Sense", first, len, newvalues);
//#ifdef SOLVER_STATISTICS
//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<milliseconds>(stop - beg);
//  STATISTICS::time_SOLVERsetSensearray += (double) duration.count() * 1e-3;
//#endif
//  return error;
//}

int SOLVER::SOLVERsetenvMethod(int value) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvMethod;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetintparam(GRBgetenv(model), "Method", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvMethod += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetenvMethod(int *value) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetenvMethod;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetintparam(GRBgetenv(model), "Method", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto elap = duration<double>(stop - beg).count();
  STATISTICS::time_SOLVERgetenvMethod += elap;
#endif
  return error;
}

int SOLVER::SOLVERsetenvTimeLimit(double value) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetenvTimeLimit;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetenvTimeLimit += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetModelSense(int value) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetModelSense;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetintattr(model, "ModelSense", value);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetModelSense += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetVBasis(int first, int len, int *vbasis) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetVBasis;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetintattrarray(model, "VBasis", first, len, vbasis);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetVBasis += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetCBasis(int first, int len, int *cbasis) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERsetCBasis;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBsetintattrarray(model, "CBasis", first, len, cbasis);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERsetCBasis += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetVBasis(int first, int len, int *vbasis) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetVBasis;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetintattrarray(model, "VBasis", first, len, vbasis);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetVBasis += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetCBasis(int first, int len, int *cbasis) const {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetCBasis;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetintattrarray(model, "CBasis", first, len, cbasis);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetCBasis += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERgetSense(int first, int len, char *sense) {
#ifdef SOLVER_STATISTICS
  ++STATISTICS::calls_SOLVERgetSense;
  auto beg = high_resolution_clock::now();
#endif
  auto error = GRBgetcharattrarray(model, "Sense", first, len, sense);
#ifdef SOLVER_STATISTICS
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - beg);
  STATISTICS::time_SOLVERgetSense += (double) duration.count() * 1e-3;
#endif
  return error;
}

int SOLVER::SOLVERsetColLower(int col, double value) {
  return GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, col, value);
}

int SOLVER::SOLVERsetColUpper(int col, double value) {
  return GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, col, value);
}

int SOLVER::SOLVERremoveColLower(int col) {
  return GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, col, -GRB_INFINITY);
}

int SOLVER::SOLVERremoveColUpper(int col) {
  return GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, col, GRB_INFINITY);
}

int SOLVER::SOLVERgetcoeff(int row, int col, double &value) const {
  return GRBgetcoeff(model, row, col, &value);
}

int SOLVER::SOLVERgetSolCount(int *valueP) const {
  return GRBgetintattr(model, "SolCount", valueP);
}

int SOLVER::SOLVERgetObjFromPool(int sol, double *valueP) const {
  auto error = GRBsetintparam(GRBgetenv(model), "SolutionNumber", sol);
  error += GRBgetdblattr(model, "PoolObjVal", valueP);
  return error;
}

int SOLVER::SOLVERgetSolFromPool(int i, int first, int len, double *values) const {
  return GRBgetdblattrarray(model, "Xn", first, len, values);
}

#endif




