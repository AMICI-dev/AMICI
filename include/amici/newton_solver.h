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
    realtype *t, AmiVector *x, const Solver &simulationSolver, Model *model);

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
        return num_lin_steps_;
    }

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the
     * linear solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    virtual void prepareLinearSystem(int ntry, int nnewt) = 0;

    /**
     * Writes the Jacobian (JB) for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    virtual void prepareLinearSystemB(int ntry, int nnewt) = 0;

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    virtual void solveLinearSystem(AmiVector &rhs) = 0;
    
    void updateModel(Model *model) {model_ = model;}

    virtual ~NewtonSolver() = default;

    /** maximum number of allowed linear steps per Newton step for steady state
     * computation */
    int max_lin_steps_ {0};
    /** maximum number of allowed Newton steps for steady state computation */
    int max_steps {0};
    /** damping factor flag */
    NewtonDampingFactorMode damping_factor_mode_ {NewtonDampingFactorMode::on};
    /** damping factor lower bound */
    realtype damping_factor_lower_bound {1e-8};

  protected:
    /** time variable */
    realtype *t_;
    /** pointer to the model object */
    Model *model_;
    /** right hand side AmiVector */
    AmiVector xdot_;
    /** current state */
    AmiVector *x_;
    /** current state time derivative (DAE) */
    AmiVector dx_;
    /** history of number of linear steps */
    std::vector<int> num_lin_steps_;
    /** current adjoint state */
    AmiVector xB_;
    /** current adjoint state time derivative (DAE) */
    AmiVector dxB_;
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

    NewtonSolverDense(const NewtonSolverDense&) = delete;

    NewtonSolverDense& operator=(const NewtonSolverDense& other) = delete;

    ~NewtonSolverDense() override;

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the
     * linear solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystem(int ntry, int nnewt) override;

    /**
     * Writes the Jacobian (JB) for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystemB(int ntry, int nnewt) override;

  private:
    /** temporary storage of Jacobian */
    SUNMatrixWrapper Jtmp_;

    /** dense linear solver */
    SUNLinearSolver linsol_ {nullptr};
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

    NewtonSolverSparse(const NewtonSolverSparse&) = delete;

    NewtonSolverSparse& operator=(const NewtonSolverSparse& other) = delete;

    ~NewtonSolverSparse() override;

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector &rhs) override;

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the
     * linear solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystem(int ntry, int nnewt) override;

    /**
     * Writes the Jacobian (JB) for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */
    void prepareLinearSystemB(int ntry, int nnewt) override;

  private:
    /** temporary storage of Jacobian */
    SUNMatrixWrapper Jtmp_;

    /** sparse linear solver */
    SUNLinearSolver linsol_ {nullptr};
};


} // namespace amici

#endif // NEWTON_SOLVER
