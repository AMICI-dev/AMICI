#ifndef amici_newton_solver_h
#define amici_newton_solver_h

#include "amici/solver.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/vector.h"

#include <memory>

namespace amici {

class Model;
class Solver;
class AmiVector;
struct SimulationState;

/**
 * @brief The NewtonSolver class sets up the linear solver for the Newton
 * method.
 */

class NewtonSolver {

  public:
    /**
     * @brief Initializes solver according to the dimensions in the provided
     * model
     *
     * @param model pointer to the model object
     */
    explicit NewtonSolver(Model const& model);

    /**
     * @brief Factory method to create a NewtonSolver based on linsolType
     *
     * @param simulationSolver solver with settings
     * @param model pointer to the model instance
     * @return solver NewtonSolver according to the specified linsolType
     */
    static std::unique_ptr<NewtonSolver>
    getSolver(Solver const& simulationSolver, Model const& model);

    /**
     * @brief Computes the solution of one Newton iteration
     *
     * @param delta containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     * @param model pointer to the model instance
     * @param state current simulation state
     */
    void getStep(AmiVector& delta, Model& model, SimulationState const& state);

    /**
     * @brief Computes steady state sensitivities
     *
     * @param sx pointer to state variable sensitivities
     * @param model pointer to the model instance
     * @param state current simulation state
     */
    void computeNewtonSensis(
        AmiVectorArray& sx, Model& model, SimulationState const& state
    );

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the
     * linear solver
     *
     * @param model pointer to the model instance
     * @param state current simulation state
     */
    virtual void prepareLinearSystem(Model& model, SimulationState const& state)
        = 0;

    /**
     * Writes the Jacobian (JB) for the Newton iteration and passes it to the
     * linear solver
     *
     * @param model pointer to the model instance
     * @param state current simulation state
     */
    virtual void
    prepareLinearSystemB(Model& model, SimulationState const& state)
        = 0;

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    virtual void solveLinearSystem(AmiVector& rhs) = 0;

    /**
     * @brief Reinitialize the linear solver
     *
     */
    virtual void reinitialize() = 0;

    /**
     * @brief Checks whether linear system is singular
     *
     * @param model pointer to the model instance
     * @param state current simulation state
     * @return boolean indicating whether the linear system is singular
     * (condition number < 1/machine precision)
     */
    virtual bool is_singular(Model& model, SimulationState const& state) const
        = 0;

    virtual ~NewtonSolver() = default;

  protected:
    /** dummy rhs, used as dummy argument when computing J and JB */
    AmiVector xdot_;
    /** dummy state, attached to linear solver */
    AmiVector x_;
    /** dummy adjoint state, used as dummy argument when computing JB */
    AmiVector xB_;
    /** dummy differential adjoint state, used as dummy argument when computing
     * JB */
    AmiVector dxB_;
};

/**
 * @brief The NewtonSolverDense provides access to the dense linear solver for
 * the Newton method.
 */

class NewtonSolverDense : public NewtonSolver {

  public:
    /**
     * @brief constructor for sparse solver
     *
     * @param model model instance that provides problem dimensions
     */
    explicit NewtonSolverDense(Model const& model);

    NewtonSolverDense(NewtonSolverDense const&) = delete;

    NewtonSolverDense& operator=(NewtonSolverDense const& other) = delete;

    ~NewtonSolverDense() override;

    void solveLinearSystem(AmiVector& rhs) override;

    void
    prepareLinearSystem(Model& model, SimulationState const& state) override;

    void
    prepareLinearSystemB(Model& model, SimulationState const& state) override;

    void reinitialize() override;

    bool is_singular(Model& model, SimulationState const& state) const override;

  private:
    /** temporary storage of Jacobian */
    SUNMatrixWrapper Jtmp_;

    /** dense linear solver */
    SUNLinearSolver linsol_{nullptr};
};

/**
 * @brief The NewtonSolverSparse provides access to the sparse linear solver for
 * the Newton method.
 */

class NewtonSolverSparse : public NewtonSolver {

  public:
    /**
     * @brief constructor for dense solver
     *
     * @param model model instance that provides problem dimensions
     */
    explicit NewtonSolverSparse(Model const& model);

    NewtonSolverSparse(NewtonSolverSparse const&) = delete;

    NewtonSolverSparse& operator=(NewtonSolverSparse const& other) = delete;

    ~NewtonSolverSparse() override;

    void solveLinearSystem(AmiVector& rhs) override;

    void
    prepareLinearSystem(Model& model, SimulationState const& state) override;

    void
    prepareLinearSystemB(Model& model, SimulationState const& state) override;

    bool is_singular(Model& model, SimulationState const& state) const override;

    void reinitialize() override;

  private:
    /** temporary storage of Jacobian */
    SUNMatrixWrapper Jtmp_;

    /** sparse linear solver */
    SUNLinearSolver linsol_{nullptr};
};

} // namespace amici

#endif // NEWTON_SOLVER
