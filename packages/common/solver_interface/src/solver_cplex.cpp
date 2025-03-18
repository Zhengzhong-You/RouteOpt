/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "solver.hpp"
#include <vector>


namespace RouteOpt {
#if SOLVER_TYPE == 1

    int Solver::getSlack(int first, int len, double *const values) const {
        return CPXgetslack(env, model, values, first, first + len - 1);
    }

    int Solver::getDual(int first, int len, double *const values) const {
        return CPXgetpi(env, model, values, first, first + len - 1);
    }

    int Solver::getCrossOver(int *const valueP) const {
        return CPXgetintparam(env, CPXPARAM_SolutionType, valueP);
    }

    int Solver::getNumRow(int *const valueP) const {
        *valueP = CPXgetnumrows(env, model);
        return 0;
    }

    int Solver::getNumCol(int *const valueP) const {
        *valueP = CPXgetnumcols(env, model);
        return 0;
    }

    int Solver::getStatus(int *const valueP) const {
        *valueP = CPXgetstat(env, model);
        return 0;
    }

    int Solver::delConstraints(int len, int *const cind) {
        int numRows = CPXgetnumrows(env, model);

        std::vector<int> delstat(numRows, 0);

        for (int i = 0; i < len; i++) delstat[cind[i]] = 1;

        int status = CPXdelsetrows(env, model, delstat.data());
        return status;
    }

    int Solver::reoptimize(int method) const {
        int error;
        int now_method;
        CPXgetintparam(env, CPX_PARAM_LPMETHOD, &now_method);
        if (now_method != method) {
            error = CPXsetintparam(env, CPX_PARAM_LPMETHOD, method);
            if (error) return error;
        }

        error = CPXlpopt(env, model);

        if (error) return error;

        if (now_method != method) {
            error = CPXsetintparam(env, CPX_PARAM_LPMETHOD, now_method);
            if (error) return error;
        }
        return 0;
    }

    int Solver::optimize() {
        return CPXlpopt(env, model);
    }

    int Solver::mipOptimize() {
        return CPXmipopt(env, model);
    }

    SOLVERmodel Solver::copyModel() const {
        int status_p;
        CPXLPptr newModel = CPXcloneprob(env, model, &status_p);
        if (status_p)
          throw std::runtime_error("Failed to std::copy model: " + std::to_string(status_p));
        return newModel;
    }

    int Solver::freeModel() {
        int error = CPXfreeprob(env, &model);
        model = nullptr;
        return error;
    }

    int Solver::loadEnv(const char *const logfilename) {
        int status;

        env = CPXopenCPLEX(&status);
        if (status) return status;

        status = CPXsetlogfilename(env, logfilename, "w");
        return status;
    }

    int Solver::setEnvThreads(int value, bool if_model_free) {
        return CPXsetintparam(env, CPXPARAM_Threads, value);
    }

    int Solver::getThreads(int *const valueP) const {
        return CPXgetintparam(env, CPXPARAM_Threads, valueP);
    }

    int Solver::setEnvOutputFlag(int value, bool if_model_free) {
        return CPXsetintparam(env, CPXPARAM_ScreenOutput, value);
    }

    int Solver::setEnvInfUnbdInfo(int value, bool if_model_free) {
        return 0; //disable this one in cplex
    }

    int Solver::setEnvMIPGap(double value, bool if_model_free) {
        return CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, value);
    }

    int Solver::setEnvFeasibilityTol(double value, bool if_model_free){
      return CPXsetdblparam(env, CPXPARAM_EpRHS, value);
      }

    int Solver::newModel(const char *const Pname, int numvars,
                         double *const obj, double *const lb, double *const ub, char *const vtype,
                         char **const varnames) {
        int status;
        model = CPXcreateprob(env, &status, Pname);
        if (status) return status;

        status = CPXnewcols(env, model, numvars, obj, lb, ub, vtype, varnames);
        return status;
    }

    int Solver::addConstraints(int numconstrs, int numnz,
                               int *const cbeg, int *const cind, double *const cval,
                               char *const sense, double *const rhs, char **const constrnames) {
        return CPXaddrows(env, model, 0, numconstrs, numnz, rhs, sense, cbeg, cind, cval, nullptr, constrnames);
    }

    int Solver::updateModel() {
        return 0; //do nothing
    }

    int Solver::write(const char *filename) const {
        return CPXwriteprob(env, model, filename, nullptr);
    }

    int Solver::delVars(int len, int *const ind) {
        int numCols = CPXgetnumcols(env, model);

        std::vector<int> delstat(numCols, 0);

        for (int i = 0; i < len; i++) delstat[ind[i]] = 1;

        int status = CPXdelsetcols(env, model, delstat.data());
        return status;
    }

    void Solver::freeEnv() {
        CPXcloseCPLEX(&env);
    }

    int Solver::getObjVal(double *const valueP) const {
        return CPXgetobjval(env, model, valueP);
    }

    int Solver::getObj(int first, int len, double *const values) const {
        return CPXgetobj(env, model, values, first, first + len - 1);
    }

    int Solver::getRhs(int first, int len, double *const values) const {
        return CPXgetrhs(env, model, values, first, first + len - 1);
    }

    int Solver::getX(int first, int len, double *const values) const {
        return CPXgetx(env, model, values, first, first + len - 1);
    }

    void Solver::getSolver(Solver *const solver) {
        model = solver->copyModel();
        env = solver->env;
    }

    void Solver::getEnv(Solver *const solver) {
        env = solver->env;
    }

    int Solver::getRC(int first, int len, double *const values) const {
        return CPXgetdj(env, model, values, first, first + len - 1);
    }

    int Solver::changeObj(int first, int len, double *const values) {
        std::vector<int> indices(len);
        for (int i = 0; i < len; i++) indices[i] = first + i;
        return CPXchgobj(env, model, len, indices.data(), values);
    }

    int Solver::setRhs(int first, int len, char *sense, double *const values) {
        std::vector<int> indices(len);
        for (int i = 0; i < len; i++) indices[i] = first + i;
        int error = CPXchgrhs(env, model, len, indices.data(), values);
        error += CPXchgsense(env, model, len, indices.data(), sense);
        return error;
    }

    int Solver::setRhs(int first, int len, double *const values) {
        std::vector<int> indices(len);
        for (int i = 0; i < len; i++) indices[i] = first + i;
        return CPXchgrhs(env, model, len, indices.data(), values);
    }

    int Solver::changeCoeffs(int cnt, int *const cind, int *const vind, double *const val) {
        return CPXchgcoeflist(env, model, cnt, cind, vind, val);
    }

    int Solver::XchangeCoeffs(size_t cnt, int *const cind, int *const vind, double *const val) {
        return CPXXchgcoeflist(env, model, cnt, cind, vind, val);
    }

    int Solver::addConstraint(int numnz, int *const cind, double *const cval, char sense, double rhs,
                              const char *const constrname) {
        std::vector<int> cmatbeg(1, 0);
        return CPXaddrows(env, model, 0, 1, numnz, &rhs, &sense, cmatbeg.data(), cind, cval, nullptr, nullptr);
    }

    int Solver::addVars(int numvars, int numnz, int *const vbeg, int *const vind, double *const vval, double *const obj,
                        double *const lb, double *const ub, char *const vtype, char **const varnames) {
        return CPXaddcols(env, model, numvars, numnz, obj, vbeg, vind, vval, lb, ub, varnames);
    }

    int Solver::addVar(int numnz, int *const vind, double *const vval, double obj, double lb, double ub, char vtype,
                       const char *const varname) {
        std::vector<int> vbeg(1, 0);
        return CPXaddcols(env, model, 1, numnz, &obj, vbeg.data(), vind, vval, &lb, &ub, nullptr);
    }

    int Solver::XaddVars(int numvars,
                         size_t numnz,
                         size_t *const vbeg,
                         int *const vind,
                         double *const vval,
                         double *const obj,
                         double *const lb,
                         double *const ub,
                         char *const vtype,
                         char **const varnames) {
        return CPXXaddcols(env, model, numvars, numnz, obj,
                           reinterpret_cast<const CPXNNZ *>(vbeg), vind, vval, lb, ub, varnames);
    }

    int Solver::XaddConstraints(int numconstrs, size_t numnz, size_t *const cbeg, int *const cind, double *const cval,
                               char *const sense, double *const rhs, char **const constrnames) {
        return CPXXaddrows(env, model, 0, numconstrs, (int) numnz, rhs, sense,
                           reinterpret_cast<const CPXNNZ *>(cbeg), cind, cval, nullptr, constrnames);
    }

    int Solver::getConstraints(int *const numnzP,
                               int *const cbeg,
                               int *const cind,
                               double *const cval,
                               int start,
                               int len) const {
        int matspace, surplus_p;
        if (cbeg == nullptr) {
            matspace = 0;
            int error = CPXgetrows(env,
                                   model,
                                   numnzP,
                                   cbeg,
                                   cind,
                                   cval,
                                   matspace,
                                   &surplus_p,
                                   start,
                                   start + len - 1);
            if (error == CPXERR_NEGATIVE_SURPLUS) {
                *numnzP = -surplus_p;
                return 0;
            } else {
                return error;
            }
        } else {
            matspace = *numnzP;
            int error = CPXgetrows(env,
                                   model,
                                   numnzP,
                                   cbeg,
                                   cind,
                                   cval,
                                   matspace,
                                   &surplus_p,
                                   start,
                                   start + len - 1);
            return error;
        }
    }

    int Solver::getVars(int *const numnzP, int *const vbeg, int *const vind, double *const vval, int start,
                        int len) const {
        int matspace, surplus_p;

        if (vbeg == nullptr) {
            matspace = 0;
            int error = CPXgetcols(env,
                                   model,
                                   numnzP,
                                   vbeg,
                                   vind,
                                   vval,
                                   matspace,
                                   &surplus_p,
                                   start,
                                   start + len - 1);
            if (error == CPXERR_NEGATIVE_SURPLUS) {
                *numnzP = -surplus_p;
                return 0;
            } else {
                return error;
            }
        } else {
            matspace = *numnzP;
            int error = CPXgetcols(env,
                                   model,
                                   numnzP,
                                   vbeg,
                                   vind,
                                   vval,
                                   matspace,
                                   &surplus_p,
                                   start,
                                   start + len - 1);
            return error;
        }
    }

    int Solver::XgetConstraints(size_t *const numnzP,
                                size_t *const cbeg,
                                int *const cind,
                                double *const cval,
                                int start,
                                int len) const {
        CPXNNZ matspace, surplus_p;

        if (cbeg == nullptr) {
            matspace = 0;
            int error = CPXXgetrows(env,
                                    model,
                                    (CPXNNZ *) numnzP,
                                    reinterpret_cast<CPXNNZ *>(cbeg),
                                    cind,
                                    cval,
                                    matspace,
                                    &surplus_p,
                                    start,
                                    start + len - 1);
            if (error == CPXERR_NEGATIVE_SURPLUS) {
                *numnzP = -surplus_p;
                return 0;
            } else {
                return error;
            }
        } else {
            matspace = *numnzP;
            int error = CPXXgetrows(env,
                                    model,
                                    (CPXNNZ *) numnzP,
                                    reinterpret_cast<CPXNNZ *>(cbeg),
                                    cind,
                                    cval,
                                    matspace,
                                    &surplus_p,
                                    start,
                                    start + len - 1);
            return error;
        }
    }

    int Solver::XgetVars(size_t *const numnzP,
                         size_t *const vbeg,
                         int *const vind,
                         double *const vval,
                         int start,
                         int len) const {
        CPXNNZ matspace, surplus_p;

        if (vbeg == nullptr) {
            matspace = 0;
            int error = CPXXgetcols(env,
                                    model,
                                    (CPXNNZ *) numnzP,
                                    reinterpret_cast<CPXNNZ *>(vbeg),
                                    vind,
                                    vval,
                                    matspace,
                                    &surplus_p,
                                    start,
                                    start + len - 1);
            if (error == CPXERR_NEGATIVE_SURPLUS) {
                *numnzP = -surplus_p;
                return 0;
            } else {
                return error;
            }
        } else {
            matspace = *numnzP;
            int error = CPXXgetcols(env,
                                    model,
                                    (CPXNNZ *) numnzP,
                                    reinterpret_cast<CPXNNZ *>(vbeg),
                                    vind,
                                    vval,
                                    matspace,
                                    &surplus_p,
                                    start,
                                    start + len - 1);
            return error;
        }
    }

    int Solver::setEnvCrossOver(int value) const {
        if (value == SOLVER_CROSSOVER_DOWN) {
            return CPXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
        }
        CPXsetintparam(env, CPXPARAM_SolutionType, CPX_BASIC_SOLN);
        return CPXsetintparam(env, CPXPARAM_Barrier_Crossover, value);
    }

    int Solver::setEnvCutoff(double value) {
        return CPXsetdblparam(env, CPX_PARAM_CUTUP, value);
    }

    int Solver::setVTypeArray(int first, int len, char *const newvalues) {
        std::vector<int> indices(len);
        for (int i = 0; i < len; i++) indices[i] = first + i;
        return CPXchgctype(env, model, len, indices.data(), newvalues);
    }

    int Solver::setEnvMethod(int value) const {
        return CPXsetintparam(env, CPXPARAM_LPMethod, value);
    }

    int Solver::getEnvMethod(int *value) const {
        return CPXgetintparam(env, CPXPARAM_LPMethod, value);
    }

    int Solver::setEnvTimeLimit(double value) {
        return CPXsetdblparam(env, CPXPARAM_TimeLimit, value);
    }

    int Solver::setModelSense(int value) {
        return CPXchgobjsen(env, model, value);
    }

    int Solver::setVBasis(int first, int len, int *vbasis) {
        return CPXcopybase(env, model, vbasis, nullptr);
    }

    int Solver::setCBasis(int first, int len, int *cbasis) {
        return CPXcopybase(env, model, nullptr, cbasis);
    }

    int Solver::getVBasis(int first, int len, int *vbasis) const {
        return CPXgetbase(env, model, vbasis, nullptr);
    }

    int Solver::getCBasis(int first, int len, int *cbasis) const {
        return CPXgetbase(env, model, nullptr, cbasis);
    }

    int Solver::getSense(int first, int len, char *sense) {
        return CPXgetsense(env, model, sense, first, first + len - 1);
    }

    int Solver::setColLower(int col, double value) {
        return CPXchgbds(env, model, 1, &col, "L", &value);
    }

    int Solver::setColUpper(int col, double value) {
        return CPXchgbds(env, model, 1, &col, "U", &value);
    }

    int Solver::removeColLower(int col) {
        double value = -CPX_INFBOUND;
        return CPXchgbds(env, model, 1, &col, "L", &value);
    }

    int Solver::removeColUpper(int col) {
        double value = CPX_INFBOUND;
        return CPXchgbds(env, model, 1, &col, "U", &value);
    }

    int Solver::getCoeff(int row, int col, double &value) const {
        return CPXgetcoef(env, model, row, col, &value);
    }

    int Solver::getSolCount(int *valueP) const {
        *valueP = CPXgetsolnpoolnumsolns(env, model);
        return 0;
    }

    int Solver::getObjFromPool(int sol, double *valueP) const {
        return CPXgetsolnpoolobjval(env, model, sol, valueP);
    }

    int Solver::getSolFromPool(int i, int first, int len, double *values) const {
        return CPXgetsolnpoolx(env, model, i, values, first, first + len - 1);
    }

    int Solver::readModel(const char *filename) {
        int status;
        model = CPXcreateprob(env, &status, nullptr);
        if (status) return status;
        status = CPXreadcopyprob(env, model, filename, nullptr);
        return status;
    }
#endif
}
