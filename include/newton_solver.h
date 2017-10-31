#ifndef amici_newton_solver_h
#define amici_newton_solver_h

#include <klu.h>
#include <nvector/nvector_serial.h> // DlsMat
#include <sundials/sundials_dense.h>
#include <sundials/sundials_sparse.h> // SlsMat

namespace amici {

class UserData;
class ReturnData;
class TempData;
class Model;

/**
 * @brief The NewtonSolver class sets up the linear solver for the Newton
 * method.
 */

class NewtonSolver {

  public:
    NewtonSolver(Model *model, ReturnData *rdata, const UserData *udata,
                 TempData *tdata);

    static NewtonSolver *getSolver(int linsolType, Model *model,
                                   ReturnData *rdata, const UserData *udata,
                                   TempData *tdata);

    void getStep(int ntry, int nnewt, N_Vector delta);

    void getSensis(int it);

    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     * @return stats integer flag indicating success of the method
     */
    virtual void prepareLinearSystem(int ntry, int nnewt) = 0;

    /**
     * Solves the linear system for the Newton step
     *
     * @param[in,out] rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system @type N_Vector
     * @return stats integer flag indicating success of the method
     */
    virtual void solveLinearSystem(N_Vector rhs) = 0;

    ~NewtonSolver() {
    if(sx_ip)
        N_VDestroy_Serial(sx_ip);
    };

  protected:
    /** pointer to the AMICI model object */
    Model *model;
    /** pointer to the return data object */
    ReturnData *rdata;
    /** pointer to the user data object */
    const UserData *udata;
    /** pointer to the temporary data object */
    TempData *tdata;
    /** sensitivity N_Vector */
    N_Vector sx_ip = nullptr;
};

/**
 * @brief The NewtonSolverDense provides access to the dense linear solver for
 * the Newton method.
 */

class NewtonSolverDense : public NewtonSolver {

  public:
    NewtonSolverDense(Model *model, ReturnData *rdata, const UserData *udata,
                      TempData *tdata);
    void solveLinearSystem(N_Vector rhs);
    void prepareLinearSystem(int ntry, int nnewt);
    ~NewtonSolverDense();

  private:
    /** temporary storage of pivot array */
    long int *pivots;
    /** temporary N_Vector storage  */
    N_Vector tmp1 = nullptr;
    /** temporary N_Vector storage  */
    N_Vector tmp2 = nullptr;
    /** temporary N_Vector storage  */
    N_Vector tmp3 = nullptr;
};

/**
 * @brief The NewtonSolverSparse provides access to the sparse linear solver for
 * the Newton method.
 */

class NewtonSolverSparse : public NewtonSolver {

  public:
    NewtonSolverSparse(Model *model, ReturnData *rdata, const UserData *udata,
                       TempData *tdata);
    void solveLinearSystem(N_Vector rhs);
    void prepareLinearSystem(int ntry, int nnewt);
    ~NewtonSolverSparse();

  private:
    /** temporary N_Vector storage  */
    N_Vector tmp1 = nullptr;
    /** temporary N_Vector storage  */
    N_Vector tmp2 = nullptr;
    /** temporary N_Vector storage  */
    N_Vector tmp3 = nullptr;
    /** klu common storage? */
    klu_common common;
    /** klu symbolic storage? */
    klu_symbolic *symbolic = nullptr;
    /** klu numeric stoarge? */
    klu_numeric *numeric = nullptr;
    /** klu status flag  */
    int klu_status = 0;
};

/**
 * @brief The NewtonSolverIterative provides access to the iterative linear
 * solver for the Newton method.
 */

class NewtonSolverIterative : public NewtonSolver {

  public:
    NewtonSolverIterative(Model *model, ReturnData *rdata,
                          const UserData *udata, TempData *tdata);
    void solveLinearSystem(N_Vector rhs);
    void prepareLinearSystem(int ntry, int nnewt);
    void linsolveSPBCG(int ntry,int nnewt, N_Vector ns_delta);
    ~NewtonSolverIterative();

  private:
    /** number of tries  */
    int newton_try;
    /** number of iterations  */
    int i_newton;
    /** ???  */
    N_Vector ns_p = nullptr;
    /** ???  */
    N_Vector ns_h = nullptr;
    /** ???  */
    N_Vector ns_t = nullptr;
    /** ???  */
    N_Vector ns_s = nullptr;
    /** ???  */
    N_Vector ns_r = nullptr;
    /** ???  */
    N_Vector ns_rt = nullptr;
    /** ???  */
    N_Vector ns_v = nullptr;
    /** ???  */
    N_Vector ns_Jv = nullptr;
    /** ???  */
    N_Vector ns_tmp = nullptr;
    /** ???  */
    N_Vector ns_Jdiag = nullptr;
};


} // namespace amici

#endif // NEWTON_SOLVER
