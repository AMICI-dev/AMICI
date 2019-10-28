#ifndef amici_newton_solver_h
#define amici_newton_solver_h

#include "amici/vector.h"
#include "amici/defines.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/sundials_linsol_wrapper.h"

#include <memory>



namespace amici {

class ReturnData;
class Model;
class AmiVector;

/**
 * @brief The NewtonSolver class sets up the linear solver for the Newton
 * method.
 */

class NewtonSolver {

  public:
    /**
     * Initializes all members with the provided objects
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */
    NewtonSolver(realtype *t, AmiVector *x, Model *model, ReturnData *rdata);

    /**
     * Factory method to create a NewtonSolver based on linsolType
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param linsolType integer indicating which linear solver to use
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     * @param maxlinsteps maximum number of allowed linear steps per Newton step for steady state computation
     * @param maxsteps maximum number of allowed Newton steps for steady state computation
     * @param atol absolute tolerance
     * @param rtol relative tolerance
     * @return solver NewtonSolver according to the specified linsolType
     */
    static std::unique_ptr<NewtonSolver> getSolver(realtype *t, AmiVector *x,
                                                   LinearSolver linsolType,
                                                   Model *model,
                                                   ReturnData *rdata,
                                                   int maxlinsteps,
                                                   int maxsteps,
                                                   double atol, double rtol);

    /**
     * Computes the solution of one Newton iteration
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     * @param delta containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void getStep(int ntry, int nnewt, AmiVector &delta);

    /**
     * Computes steady state sensitivities
     *
     * @param sx pointer to state variable sensitivities
     */
    void computeNewtonSensis(AmiVectorArray &sx);

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
    virtual void solveLinearSystem(AmiVector &rhs) = 0;

    virtual ~NewtonSolver() = default;

    /** maximum number of allowed linear steps per Newton step for steady state
     * computation */
    int maxlinsteps = 0;
    /** maximum number of allowed Newton steps for steady state computation */
    int maxsteps = 0;
    /** absolute tolerance */
    double atol = 1e-16;
    /** relative tolerance */
    double rtol = 1e-8;

  protected:
    /** time variable */
    realtype *t;
    /** pointer to the AMICI model object */
    Model *model;
    /** pointer to the return data object */
    ReturnData *rdata;
    /** right hand side AmiVector */
    AmiVector xdot;
    /** current state */
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
    /**
     * Constructor, initializes all members with the provided objects
     * and initializes temporary storage objects
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */

    NewtonSolverDense(realtype *t, AmiVector *x, Model *model, ReturnData *rdata);
    ~NewtonSolverDense() override;

    /**
     * Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystem(int ntry, int nnewt) override;

  private:
    /** temporary storage of Jacobian */
    SUNMatrixWrapper Jtmp;

    /** dense linear solver */
    SUNLinearSolver linsol = nullptr;
};

/**
 * @brief The NewtonSolverSparse provides access to the sparse linear solver for
 * the Newton method.
 */

class NewtonSolverSparse : public NewtonSolver {

  public:
    /**
     * Constructor, initializes all members with the provided objects,
     * initializes temporary storage objects and the klu solver
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */
    NewtonSolverSparse(realtype *t, AmiVector *x, Model *model, ReturnData *rdata);
    ~NewtonSolverSparse() override;

    /**
     * Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystem(int ntry, int nnewt) override;

  private:
    /** temporary storage of Jacobian */
    SUNMatrixWrapper Jtmp;

    /** sparse linear solver */
    SUNLinearSolver linsol = nullptr;
};

/**
 * @brief The NewtonSolverIterative provides access to the iterative linear
 * solver for the Newton method.
 */

class NewtonSolverIterative : public NewtonSolver {

  public:
    /**
     * Constructor, initializes all members with the provided objects
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */
    NewtonSolverIterative(realtype *t, AmiVector *x, Model *model, ReturnData *rdata);
    ~NewtonSolverIterative() override = default;

    /**
     * Solves the linear system for the Newton step by passing it to
     * linsolveSPBCG
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver.
     * Also wraps around getSensis for iterative linear solver.
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystem(int ntry, int nnewt) override;

    /**
     * Iterative linear solver created from SPILS BiCG-Stab.
     * Solves the linear system within each Newton step if iterative solver is
     * chosen.
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     * @param ns_delta Newton step
     */
    void linsolveSPBCG(int ntry, int nnewt, AmiVector &ns_delta);

  private:
    /** number of tries  */
    int newton_try = 0;
    /** number of iterations  */
    int i_newton = 0;
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
