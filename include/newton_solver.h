#ifndef newton_solver
#define newton_solver

#include <sundials/sundials_dense.h>
#include <nvector/nvector_serial.h> // DlsMat
#include <sundials/sundials_sparse.h> // SlsMat
#include <klu.h>

class NewtonSolverDense;
class NewtonSolverSparse;
class NewtonSolverIterative;
class UserData;
class ReturnData;
class TempData;
class Model;

/**
 * @brief The NewtonSolver class sets up the linear solver for the Newton method.
 */

class NewtonSolver {
    
public:
    NewtonSolver(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata);
    
    static NewtonSolver* getSolver(int linsolType, Model *model, ReturnData *rdata, UserData *udata, TempData *tdata, int *status);

    int getStep(int ntry, int nnewt, N_Vector delta);
    
    int getSensis(int it);
    
    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear solver
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     * @return stats integer flag indicating success of the method
     */
    virtual int prepareLinearSystem(int ntry, int nnewt) = 0;

    /**
     * Solves the linear system for the Newton step
     *
     * @param[in,out] rhs containing the RHS of the linear system, will be overwritten by solution to the linear system @type N_Vector
     * @return stats integer flag indicating success of the method
     */
    virtual int solveLinearSystem(N_Vector rhs) = 0;
    
    virtual ~NewtonSolver();
    
protected:
    /** pointer to the AMICI model object */
    Model *model;
    /** pointer to the return data object */
    ReturnData *rdata;
    /** pointer to the user data object */
    UserData *udata;
    /** pointer to the temporary data object */
    TempData *tdata;
};

/**
 * @brief The NewtonSolverDense provides access to the dense linear solver for the Newton method.
 */

class NewtonSolverDense : public NewtonSolver {
    
public:
    NewtonSolverDense(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata);
    int solveLinearSystem(N_Vector rhs);
    int prepareLinearSystem(int ntry, int nnewt);
    ~NewtonSolverDense();
    
private:
    /** temporary storage of pivot array */
    long int *pivots;
    /** temporary N_Vector storage  */
    N_Vector tmp1;
    /** temporary N_Vector storage  */
    N_Vector tmp2;
    /** temporary N_Vector storage  */
    N_Vector tmp3;
};

/**
 * @brief The NewtonSolverSparse provides access to the sparse linear solver for the Newton method.
 */

class NewtonSolverSparse : public NewtonSolver {
    
public:
    NewtonSolverSparse(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata);
    int solveLinearSystem(N_Vector rhs);
    int prepareLinearSystem(int ntry, int nnewt);
    ~NewtonSolverSparse();
    
private:
    /** temporary N_Vector storage  */
    N_Vector tmp1;
    /** temporary N_Vector storage  */
    N_Vector tmp2;
    /** temporary N_Vector storage  */
    N_Vector tmp3;
    /** klu common storage? */
    klu_common common;
    /** klu symbolic storage? */
    klu_symbolic *symbolic;
    /** klu numeric stoarge? */
    klu_numeric *numeric;
    /** klu status flag  */
    int klu_status = 0;
};

/**
 * @brief The NewtonSolverIterative provides access to the iterative linear solver for the Newton method.
 */

class NewtonSolverIterative : public NewtonSolver {
    
public:
    NewtonSolverIterative(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata);
    int solveLinearSystem(N_Vector rhs);
    int prepareLinearSystem(int ntry, int nnewt);
    ~NewtonSolverIterative();
    
private:
    /** number of tries  */
    int newton_try;
    /** number of iterations  */
    int i_newton;
};

#endif // NEWTON_SOLVER