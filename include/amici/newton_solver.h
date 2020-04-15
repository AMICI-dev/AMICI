#ifndef amici_newton_solver_h
#define amici_newton_solver_h

#include "amici/vector.h"
#include "amici/defines.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/sundials_linsol_wrapper.h"

#include <memory>

namespace amici {

class Model;
class Solver;
class AmiVector;

/**
 * @brief The NewtonSolver class sets up the linear solver for the Newton
 * method.
 */

class NewtonSolver {

  public:
    /**
     * @brief Initializes all members with the provided objects
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the model object
     */
    NewtonSolver(realtype *t, AmiVector *x, Model *model);

    /**
     * @brief Factory method to create a NewtonSolver based on linsolType
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param simulationSolver solver with settings
     * @param model pointer to the model object
     * @return solver NewtonSolver according to the specified linsolType
     */
    static std::unique_ptr<NewtonSolver> getSolver(
    realtype *t, AmiVector *x, Solver &simulationSolver, Model *model);

    /**
     * @brief Computes the solution of one Newton iteration
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     * @param delta containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void getStep(int ntry, int nnewt, AmiVector &delta);

    /**
     * @brief Computes steady state sensitivities
     *
     * @param sx pointer to state variable sensitivities
     */
    void computeNewtonSensis(AmiVectorArray &sx);
    
    /**
     * @brief Accessor for numlinsteps
     *
     * @return numlinsteps
     */
    const std::vector<int> &getNumLinSteps() const {
        return numlinsteps;
    }

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    virtual void prepareLinearSystem(int ntry, int nnewt) = 0;

    /**
     * @brief Solves the linear system for the Newton step
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
    /** damping factor flag */
    NewtonDampingFactorMode dampingFactorMode = NewtonDampingFactorMode::on;
    /** damping factor lower bound */
    double dampingFactorLowerBound = 1e-8;

  protected:
    /** time variable */
    realtype *t;
    /** pointer to the model object */
    Model *model;
    /** right hand side AmiVector */
    AmiVector xdot;
    /** current state */
    AmiVector *x;
    /** current state time derivative (DAE) */
    AmiVector dx;
    
    /** history of number of linear steps */
    std::vector<int> numlinsteps;
};

/**
 * @brief The NewtonSolverDense provides access to the dense linear solver for
 * the Newton method.
 */

class NewtonSolverDense : public NewtonSolver {

  public:
    /**
     * @brief Constructor, initializes all members with the provided objects
     * and initializes temporary storage objects
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the model object
     */

    NewtonSolverDense(realtype *t, AmiVector *x, Model *model);
    ~NewtonSolverDense() override;

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the linear
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
     * @brief Constructor, initializes all members with the provided objects,
     * initializes temporary storage objects and the klu solver
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the model object
     */
    NewtonSolverSparse(realtype *t, AmiVector *x, Model *model);
    ~NewtonSolverSparse() override;

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the linear
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
     * @brief Constructor, initializes all members with the provided objects
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the model object
     */
    NewtonSolverIterative(realtype *t, AmiVector *x, Model *model);
    ~NewtonSolverIterative() override = default;

    /**
     * @brief Solves the linear system for the Newton step by passing it to
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
