#ifndef amici_newton_solver_h
#define amici_newton_solver_h

#include <klu.h>
#include <nvector/nvector_serial.h> // DlsMat
#include <sundials/sundials_dense.h>
#include <sundials/sundials_sparse.h> // SlsMat

namespace amici {

class UserData;
class ReturnData;
class Model;
class AmiVector;

/**
 * @brief The NewtonSolver class sets up the linear solver for the Newton
 * method.
 */

class NewtonSolver {

  public:
    NewtonSolver(realtype *t, AmiVector *x, Model *model, ReturnData *rdata, const UserData *udata);

    static NewtonSolver *getSolver(realtype *t, AmiVector *x, int linsolType, Model *model,
                                   ReturnData *rdata, const UserData *udata);

    void getStep(int ntry, int nnewt, AmiVector *delta);

    void getSensis(const int it, AmiVectorArray *sx);

    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    virtual void prepareLinearSystem(int ntry, int nnewt) = 0;

    /**
     * Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system 
     */
    virtual void solveLinearSystem(AmiVector *rhs) = 0;
    virtual ~NewtonSolver() {};

  protected:
    /** time variable */
    realtype *t;
    /** pointer to the AMICI model object */
    Model *model;
    /** pointer to the return data object */
    ReturnData *rdata;
    /** pointer to the user data object */
    const UserData *udata;
    /** right hand side AmiVector */
    AmiVector xdot;
    /** current state*/
    AmiVector *x;
    /** current state time derivative (DAE) */
    AmiVector dx;
};

/**
 * @brief The NewtonSolverDense provides access to the dense linear solver for
 * the Newton method.
 */

class NewtonSolverDense : public NewtonSolver {

  public:
    NewtonSolverDense(realtype *t, AmiVector *x, Model *model, ReturnData *rdata, const UserData *udata);
    void solveLinearSystem(AmiVector *rhs) override;
    void prepareLinearSystem(int ntry, int nnewt) override;
    virtual ~NewtonSolverDense();

  private:
    /** temporary storage of pivot array */
    long int *pivots = nullptr;
    /** temporary storage of Jacobian */
    DlsMat Jtmp = nullptr;
};

/**
 * @brief The NewtonSolverSparse provides access to the sparse linear solver for
 * the Newton method.
 */

class NewtonSolverSparse : public NewtonSolver {

  public:
    NewtonSolverSparse(realtype *t, AmiVector *x, Model *model, ReturnData *rdata, const UserData *udata);
    void solveLinearSystem(AmiVector *rhs) override;
    void prepareLinearSystem(int ntry, int nnewt) override;
    virtual ~NewtonSolverSparse();

  private:
    /** klu common storage? */
    klu_common common;
    /** klu symbolic storage? */
    klu_symbolic *symbolic = nullptr;
    /** klu numeric stoarge? */
    klu_numeric *numeric = nullptr;
    /** klu status flag  */
    int klu_status = 0;
    /** temporary storage of Jacobian */
    SlsMat Jtmp = nullptr;
};

/**
 * @brief The NewtonSolverIterative provides access to the iterative linear
 * solver for the Newton method.
 */

class NewtonSolverIterative : public NewtonSolver {

  public:
    NewtonSolverIterative(realtype *t, AmiVector *x, Model *model, ReturnData *rdata,
                          const UserData *udata);
    void solveLinearSystem(AmiVector *rhs);
    void prepareLinearSystem(int ntry, int nnewt);
    void linsolveSPBCG(int ntry, int nnewt, AmiVector *ns_delta);
    virtual ~NewtonSolverIterative();

  private:
    /** number of tries  */
    int newton_try;
    /** number of iterations  */
    int i_newton;
    /** ???  */
    AmiVector ns_p;
    /** ???  */
    AmiVector ns_h;
    /** ???  */
    AmiVector ns_t;
    /** ???  */
    AmiVector ns_s;
    /** ???  */
    AmiVector ns_r;
    /** ???  */
    AmiVector ns_rt;
    /** ???  */
    AmiVector ns_v;
    /** ???  */
    AmiVector ns_Jv;
    /** ???  */
    AmiVector ns_tmp;
    /** ???  */
    AmiVector ns_Jdiag;
};


} // namespace amici

#endif // NEWTON_SOLVER
