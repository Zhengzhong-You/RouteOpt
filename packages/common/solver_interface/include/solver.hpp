/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file solver.hpp
 * @brief Solver interface for RouteOpt.
 *
 * This header defines the Solver class, which provides a unified interface for different optimization solvers.
 * It abstracts the interactions with various solver backends (Gurobi, CPLEX, MindOpt) based on the SOLVER_TYPE macro.
 */

#ifndef ROUTE_OPT_SOLVER_HPP
#define ROUTE_OPT_SOLVER_HPP

#include "solver_macro.hpp"
#if SOLVER_TYPE == 0
#include "solver_grb.hpp"
#elif SOLVER_TYPE == 1
#include "solver_cplex.hpp"
#elif SOLVER_TYPE == 2
#include "solver_mindopt.hpp"
#endif

namespace RouteOpt {
    /**
     * @brief The Solver class provides a unified interface for different optimization solvers.
     *
     * This class abstracts the interactions with various solver backends (Gurobi, CPLEX, MindOpt)
     * based on the SOLVER_TYPE macro. It offers methods for model creation, updating, optimization,
     * and querying solution information.
     */
    class Solver {
    public:
#if SOLVER_TYPE == 0
        SOLVERmodel *model; ///< Pointer to the solver model (Gurobi model).
        SOLVERenv *env; ///< Pointer to the solver environment.
#elif SOLVER_TYPE == 1
    SOLVERmodel model;   ///< Solver model instance (CPLEX model).
    SOLVERenv env;       ///< Solver environment instance.
#elif SOLVER_TYPE == 2
    Mindo::LinearModel model;  ///< MindOpt linear model.
#endif

        /**
         * @brief Retrieve solution values from the solution pool.
         * @param i Index of the solution in the pool.
         * @param first Starting index for retrieving solution variables.
         * @param len Number of variables to retrieve.
         * @param values Array to store the retrieved solution values.
         * @return Status code indicating success or failure.
         */
        int getSolFromPool(int i, int first, int len, double *values) const;

        /**
         * @brief Retrieve the objective value of a solution from the pool.
         * @param i Index of the solution.
         * @param valueP Pointer to store the objective value.
         * @return Status code.
         */
        int getObjFromPool(int i, double *valueP) const;

        /**
         * @brief Get the number of solutions available in the solution pool.
         * @param valueP Pointer to store the solution count.
         * @return Status code.
         */
        int getSolCount(int *valueP) const;

        /**
         * @brief Set the basis for decision variables.
         * @param first Starting index.
         * @param len Number of variables.
         * @param vbasis Array containing the new basis indices.
         * @return Status code.
         */
        int setVBasis(int first, int len, int *vbasis);

        /**
         * @brief Set the basis for constraints.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param cbasis Array containing the new basis indices.
         * @return Status code.
         */
        int setCBasis(int first, int len, int *cbasis);

        /**
         * @brief Retrieve the current basis for decision variables.
         * @param first Starting index.
         * @param len Number of variables.
         * @param vbasis Array to store the current basis.
         * @return Status code.
         */
        int getVBasis(int first, int len, int *vbasis) const;

        /**
         * @brief Retrieve the current basis for constraints.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param cbasis Array to store the current basis.
         * @return Status code.
         */
        int getCBasis(int first, int len, int *cbasis) const;

        /**
         * @brief Retrieve slack values for constraints.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param values Array to store slack values.
         * @return Status code.
         */
        int getSlack(int first, int len, double *values) const;

        /**
         * @brief Retrieve dual variable values.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param values Array to store dual values.
         * @return Status code.
         */
        int getDual(int first, int len, double *values) const;

        /**
         * @brief Retrieve crossover information.
         * @param valueP Pointer to store the crossover status.
         * @return Status code.
         */
        int getCrossOver(int *valueP) const;

        /**
         * @brief Get the number of rows (constraints) in the model.
         * @param valueP Pointer to store the row count.
         * @return Status code.
         */
        int getNumRow(int *valueP) const;

        /**
         * @brief Get the number of columns (variables) in the model.
         * @param valueP Pointer to store the column count.
         * @return Status code.
         */
        int getNumCol(int *valueP) const;

        /**
         * @brief Delete constraints from the model.
         * @param len Number of constraints to delete.
         * @param cind Array of constraint indices to delete.
         * @return Status code.
         */
        int delConstraints(int len, int *cind);

        /**
         * @brief Re-optimize the model using the specified method.
         * @param method Optimization method (default is SOLVER_DUAL_SIMPLEX).
         * @return Status code.
         */
        int reoptimize(int method = SOLVER_DUAL_SIMPLEX) const;

        /**
         * @brief Optimize the model.
         * @return Status code.
         */
        int optimize();

        /**
         * @brief Optimize a mixed-integer programming (MIP) model.
         * @return Status code.
         */
        int mipOptimize();

        /**
         * @brief Set the variable type array for the model.
         * @param first Starting index.
         * @param len Number of variables.
         * @param newvalues Array containing new variable types.
         * @return Status code.
         */
        int setVTypeArray(int first, int len, char *newvalues);

        /**
         * @brief Free the current model and release its resources.
         * @return Status code.
         */
        int freeModel();

        /**
         * @brief Create a new model.
         * @param Pname Model name.
         * @param numvars Number of variables.
         * @param obj Array of objective coefficients.
         * @param lb Array of lower bounds.
         * @param ub Array of upper bounds.
         * @param vtype Array of variable types.
         * @param varnames Array of variable names.
         * @return Status code.
         */
        int newModel(const char *Pname, int numvars,
                     double *obj, double *lb, double *ub, char *vtype,
                     char **varnames);

        /**
         * @brief Add multiple constraints to the model.
         * @param numconstrs Number of constraints.
         * @param numnz Total number of non-zero coefficients.
         * @param cbeg Array indicating the starting index of each constraint.
         * @param cind Array of variable indices corresponding to coefficients.
         * @param cval Array of coefficient values.
         * @param sense Array indicating the sense (e.g., '=', '<', '>') for each constraint.
         * @param rhs Array of right-hand side values.
         * @param constrnames Array of constraint names.
         * @return Status code.
         */
        int addConstraints(int numconstrs, int numnz,
                           int *cbeg, int *cind, double *cval,
                           char *sense, double *rhs, char **constrnames);

        /**
         * @brief Update the model after making changes.
         * @return Status code.
         */
        int updateModel();

        /**
         * @brief Write the model to a file.
         * @param filename Output file name.
         * @return Status code.
         */
        int write(const char *filename) const;

        /**
         * @brief Delete variables from the model.
         * @param len Number of variables to delete.
         * @param ind Array of variable indices to delete.
         * @return Status code.
         */
        int delVars(int len, int *ind);

        /**
         * @brief Free the solver environment resources.
         */
        void freeEnv();

        /**
         * @brief Retrieve the objective value of the optimized model.
         * @param valueP Pointer to store the objective value.
         * @return Status code.
         */
        int getObjVal(double *valueP) const;

        /**
         * @brief Retrieve solution values for decision variables.
         * @param first Starting index.
         * @param len Number of variables.
         * @param values Array to store the solution values.
         * @return Status code.
         */
        int getX(int first, int len, double *values) const;

        /**
         * @brief Retrieve objective coefficients for a subset of variables.
         * @param first Starting index.
         * @param len Number of variables.
         * @param values Array to store the objective coefficients.
         * @return Status code.
         */
        int getObj(int first, int len, double *values) const;

        /**
         * @brief Retrieve the right-hand side values for constraints.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param values Array to store the RHS values.
         * @return Status code.
         */
        int getRhs(int first, int len, double *values) const;

        /**
         * @brief Copy the current solver instance into another solver.
         * @param solver Pointer to the solver instance to be populated.
         */
        void getSolver(Solver *solver);

        /**
         * @brief Copy the current solver environment into another solver.
         * @param solver Pointer to the solver instance to be populated.
         */
        void getEnv(Solver *solver);

        /**
         * @brief Retrieve the reduced costs for decision variables.
         * @param first Starting index.
         * @param len Number of variables.
         * @param values Array to store reduced cost values.
         * @return Status code.
         */
        int getRC(int first, int len, double *values) const;

        /**
         * @brief Change the objective coefficients for specified variables.
         * @param first Starting index.
         * @param len Number of variables.
         * @param values Array of new objective coefficients.
         * @return Status code.
         */
        int changeObj(int first, int len, double *values);

        /**
         * @brief Change multiple coefficient values in the model.
         * @param cnt Number of coefficients to change.
         * @param cind Array of constraint indices.
         * @param vind Array of variable indices.
         * @param val Array of new coefficient values.
         * @return Status code.
         */
        int changeCoeffs(int cnt, int *cind, int *vind, double *val);

        /**
         * @brief Change coefficient values using an extended interface with size_t.
         * @param cnt Number of coefficients to change.
         * @param cind Array of constraint indices.
         * @param vind Array of variable indices.
         * @param val Array of new coefficient values.
         * @return Status code.
         */
        int XchangeCoeffs(size_t cnt, int *cind, int *vind, double *val);

        /**
         * @brief Set the right-hand side values with specified senses.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param sense Array of constraint senses.
         * @param newvalues Array of new RHS values.
         * @return Status code.
         */
        int setRhs(int first, int len, char *sense, double *newvalues);

        /**
         * @brief Set the right-hand side values for constraints.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param newvalues Array of new RHS values.
         * @return Status code.
         */
        int setRhs(int first, int len, double *newvalues);

        /**
         * @brief Set the lower bound for a specific variable.
         * @param col Variable (column) index.
         * @param value New lower bound value.
         * @return Status code.
         */
        int setColLower(int col, double value);

        /**
         * @brief Set the upper bound for a specific variable.
         * @param col Variable (column) index.
         * @param value New upper bound value.
         * @return Status code.
         */
        int setColUpper(int col, double value);

        /**
          * @brief Retrieve the lower bound of a specific variable (column).
          * @param col Index of the column whose lower bound is queried.
          * @param value Pointer to a double where the retrieved lower bound will be stored.
          * @return Integer status code (0 if successful, otherwise an error code).
          */
        int getColLower(int col, double *value) const;

        /**
         * @brief Retrieve the upper bound of a specific variable (column).
         * @param col Index of the column whose upper bound is queried.
         * @param value Pointer to a double where the retrieved upper bound will be stored.
         * @return Integer status code (0 if successful, otherwise an error code).
         */
        int getColUpper(int col, double *value) const;


        /**
         * @brief Remove the lower bound restriction from a variable.
         * @param col Variable (column) index.
         * @return Status code.
         */
        int removeColLower(int col);

        /**
         * @brief Remove the upper bound restriction from a variable.
         * @param col Variable (column) index.
         * @return Status code.
         */
        int removeColUpper(int col);

        /**
         * @brief Add a single constraint to the model.
         * @param numnz Number of non-zero coefficients.
         * @param cind Array of variable indices.
         * @param cval Array of coefficient values.
         * @param sense Constraint sense (e.g., '=', '<', '>').
         * @param rhs Right-hand side value.
         * @param constrname Name of the constraint.
         * @return Status code.
         */
        int addConstraint(int numnz, int *cind, double *cval,
                          char sense, double rhs, const char *constrname);

        /**
         * @brief Add multiple variables to the model.
         * @param numvars Number of variables.
         * @param numnz Total number of non-zero entries.
         * @param vbeg Array indicating the starting index for each variable's coefficients.
         * @param vind Array of row indices corresponding to non-zero entries.
         * @param vval Array of coefficient values.
         * @param obj Array of objective coefficients.
         * @param lb Array of lower bounds.
         * @param ub Array of upper bounds.
         * @param vtype Array of variable types.
         * @param varnames Array of variable names.
         * @return Status code.
         */
        int addVars(int numvars, int numnz,
                    int *vbeg, int *vind, double *vval,
                    double *obj, double *lb, double *ub, char *vtype,
                    char **varnames);

        /**
         * @brief Add a single variable to the model.
         * @param numnz Number of non-zero coefficients.
         * @param vind Array of constraint indices where the variable appears.
         * @param vval Array of coefficient values.
         * @param obj Objective coefficient.
         * @param lb Lower bound.
         * @param ub Upper bound.
         * @param vtype Variable type (e.g., continuous, integer).
         * @param varname Variable name.
         * @return Status code.
         */
        int addVar(int numnz, int *vind, double *vval,
                   double obj, double lb, double ub, char vtype,
                   const char *varname);

        /**
         * @brief Add multiple variables using an extended interface with size_t for indexing.
         * @param numvars Number of variables.
         * @param numnz Total number of non-zero entries.
         * @param vbeg Array indicating the starting index for each variable's coefficients.
         * @param vind Array of row indices corresponding to non-zero entries.
         * @param vval Array of coefficient values.
         * @param obj Array of objective coefficients.
         * @param lb Array of lower bounds.
         * @param ub Array of upper bounds.
         * @param vtype Array of variable types.
         * @param varnames Array of variable names.
         * @return Status code.
         */
        int XaddVars(int numvars, size_t numnz,
                     size_t *vbeg, int *vind, double *vval,
                     double *obj, double *lb, double *ub, char *vtype,
                     char **varnames);

        /**
         * @brief Add multiple constraints using an extended interface with size_t for indexing.
         * @param numconstrs Number of constraints.
         * @param numnz Total number of non-zero entries.
         * @param cbeg Array indicating the starting index for each constraint's coefficients.
         * @param cind Array of variable indices corresponding to non-zero entries.
         * @param cval Array of coefficient values.
         * @param sense Array of constraint senses.
         * @param rhs Array of right-hand side values.
         * @param constrnames Array of constraint names.
         * @return Status code.
         */
        int XaddConstraints(int numconstrs, size_t numnz,
                            size_t *cbeg, int *cind, double *cval,
                            char *sense, double *rhs, char **constrnames);

        /**
         * @brief Retrieve constraints using the extended interface.
         * @param numnzP Pointer to store the number of non-zero entries.
         * @param cbeg Array to store beginning indices.
         * @param cind Array to store variable indices.
         * @param cval Array to store coefficient values.
         * @param start Starting index for retrieval.
         * @param len Number of constraints to retrieve.
         * @return Status code.
         */
        int XgetConstraints(size_t *numnzP, size_t *cbeg,
                            int *cind, double *cval, int start, int len) const;

        /**
         * @brief Retrieve variables using the extended interface.
         * @param numnzP Pointer to store the number of non-zero entries.
         * @param vbeg Array to store beginning indices.
         * @param vind Array to store row indices.
         * @param vval Array to store coefficient values.
         * @param start Starting index for retrieval.
         * @param len Number of variables to retrieve.
         * @return Status code.
         */
        int XgetVars(size_t *numnzP, size_t *vbeg,
                     int *vind, double *vval, int start, int len) const;

        /**
         * @brief Retrieve constraints using the standard interface.
         * @param numnzP Pointer to store the number of non-zero entries.
         * @param cbeg Array to store beginning indices.
         * @param cind Array to store variable indices.
         * @param cval Array to store coefficient values.
         * @param start Starting index for retrieval.
         * @param len Number of constraints to retrieve.
         * @return Status code.
         */
        int getConstraints(int *numnzP, int *cbeg,
                           int *cind, double *cval, int start, int len) const;

        /**
         * @brief Retrieve variables using the standard interface.
         * @param numnzP Pointer to store the number of non-zero entries.
         * @param vbeg Array to store beginning indices.
         * @param vind Array to store row indices.
         * @param vval Array to store coefficient values.
         * @param start Starting index for retrieval.
         * @param len Number of variables to retrieve.
         * @return Status code.
         */
        int getVars(int *numnzP, int *vbeg, int *vind, double *vval, int start, int len) const;

        /**
         * @brief Retrieve a specific coefficient value from the model.
         * @param row Row index.
         * @param col Column index.
         * @param value Reference to store the coefficient value.
         * @return Status code.
         */
        int getCoeff(int row, int col, double &value) const;

        /**
         * @brief Get the current optimization status of the model.
         * @param valueP Pointer to store the status code.
         * @return Status code.
         */
        int getStatus(int *valueP) const;

        /**
         * @brief Set the optimization method for the solver environment.
         * @param value Method value.
         * @return Status code.
         */
        int setEnvMethod(int value) const;

        /**
         * @brief Set the crossover option for the solver environment.
         * @param value Crossover option value.
         * @return Status code.
         */
        int setEnvCrossOver(int value) const;

        /**
         * @brief Retrieve the current optimization method from the environment.
         * @param value Pointer to store the method value.
         * @return Status code.
         */
        int getEnvMethod(int *value) const;

        /**
         * @brief Load the solver environment with the specified log file.
         * @param logfilename Name of the log file.
         * @return Status code.
         */
        int loadEnv(const char *logfilename);

        /**
         * @brief Set the number of threads for the solver environment.
         * @param value Number of threads.
         * @param if_model_free Flag indicating if the model has not bound to the environment.
         * @return Status code.
         */
        int setEnvThreads(int value, bool if_model_free);

        /**
         * @brief Retrieve the number of threads used by the solver environment.
         * @param value Pointer to store the thread count.
         * @return Status code.
         */
        int getThreads(int *value) const;

        /**
         * @brief Set the output flag for the solver environment.
         * @param value Output flag value.
         * @param if_model_free Flag indicating if the model has not bound to the environment.
         * @return Status code.
         */
        int setEnvOutputFlag(int value, bool if_model_free);

        /**
         * @brief Set the infeasibility/unbounded information option in the environment.
         * @param value Option value.
         * @param if_model_free Flag indicating if the model has not bound to the environment.
         * @return Status code.
         */
        int setEnvInfUnbdInfo(int value, bool if_model_free);

        /**
         * @brief Set the MIP gap tolerance for the solver environment.
         * @param value Gap tolerance value.
         * @param if_model_free Flag indicating if the model has not bound to the environment.
         * @return Status code.
         */
        int setEnvMIPGap(double value, bool if_model_free);

        /**
         * @brief Set the feasibility tolerance for the solver environment.
         * @param value Tolerance value.
         * @param if_model_free Flag indicating if the model has not bound to the environment.
         * @return Status code.
         */
        int setEnvFeasibilityTol(double value, bool if_model_free);

        /**
         * @brief Set the cutoff value for the solver environment.
         * @param value Cutoff value.
         * @return Status code.
         */
        int setEnvCutoff(double value);

        /**
         * @brief Set the time limit for optimization.
         * @param value Time limit value.
         * @return Status code.
         */
        int setEnvTimeLimit(double value);

        /**
         * @brief Set the optimization sense for the model (minimize or maximize).
         * @param value Sense value (e.g., SOLVER_MAX_SENSE).
         * @return Status code.
         */
        int setModelSense(int value);

        /**
         * @brief Retrieve the constraint senses for a subset of constraints.
         * @param first Starting index.
         * @param len Number of constraints.
         * @param sense Array to store the constraint senses.
         * @return Status code.
         */
        int getSense(int first, int len, char *sense);

        /**
         * @brief Read a model from a file.
         * @param filename File name containing the model.
         * @return Status code.
         */
        int readModel(const char *filename);

#if SOLVER_TYPE == 0
        /**
         * @brief Create a copy of the current solver model.
         * @return Pointer to the copied solver model.
         */
        [[nodiscard]] SOLVERmodel *copyModel() const;

        /**
         * @brief Retrieve the current solver environment.
         * @return Pointer to the solver environment.
         */
        [[nodiscard]] SOLVERenv *passEnv() const;
#elif SOLVER_TYPE == 1
    /**
     * @brief Create a copy of the current solver model.
     * @return A copy of the solver model.
     */
    SOLVERmodel copyModel() const;
#elif SOLVER_TYPE == 2
    /**
     * @brief Create a copy of the current solver.
     * @return Pointer to the copied Solver instance.
     */
    Solver* copyModel() const;
#endif
    };
} // namespace RouteOpt

#endif // ROUTE_OPT_SOLVER_HPP
