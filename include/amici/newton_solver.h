#ifndef amici_newton_solver_h
#define amici_newton_solver_h

#include "amici/defines.h"
#include "amici/model_state.h"
#include "amici/sundials_linsol_wrapper.h"
#include "amici/vector.h"

namespace amici {

class Model;
class Solver;
class AmiVector;

/**
 * @brief A non-owning container that groups references to the dynamical state
 * of differential equation system.
 *
 * It's intended to passing around the current state of a system while
 * avoiding repetetive argument passing in function signatures.
 */
struct DEStateView {
    /**
     * @brief Reference to the current simulation time \( t \).
     */
    realtype& t;

    /**
     * @brief Reference to the state vector \( x(t) \).
     */
    AmiVector& x;

    /**
     * @brief Reference to the time derivative of the state vector (DAE)
     */
    AmiVector& dx;

    /**
     * @brief Construct a DEStateView from references to state, derivative, and
     * time.
     *
     * @param t_ Reference to the simulation time.
     * @param x_ Reference to the state vector.
     * @param dx_ Reference to the state derivative.
     */
    DEStateView(realtype& t_, AmiVector& x_, AmiVector& dx_)
        : t(t_)
        , x(x_)
        , dx(dx_) {}

    /**
     * @brief Construct a DEStateView from a SolutionState reference.
     * @param sol
     */
    DEStateView(SolutionState& sol)
        : t(sol.t)
        , x(sol.x)
        , dx(sol.dx) {}
};

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
     * @param model the model object
     * @param linsol_type type of linear solver to use
     * @param sunctx SUNDIALS context
     */
    explicit NewtonSolver(
        Model const& model, LinearSolver linsol_type, SUNContext sunctx
    );

    NewtonSolver(NewtonSolver const&) = delete;

    NewtonSolver& operator=(NewtonSolver const& other) = delete;

    /**
     * @brief Computes the solution of one Newton iteration
     *
     * @param delta containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     * @param model the model instance
     * @param state current simulation state
     */
    void getStep(AmiVector& delta, Model& model, DEStateView const& state);

    /**
     * @brief Computes steady state sensitivities
     *
     * @param sx state variable sensitivities
     * @param model the model instance
     * @param state current simulation state
     */
    void computeNewtonSensis(
        AmiVectorArray& sx, Model& model, DEStateView const& state
    );

    /**
     * @brief Writes the Jacobian for the Newton iteration and passes it to the
     * linear solver
     *
     * @param model the model instance
     * @param state current simulation state
     */
    void prepareLinearSystem(Model& model, DEStateView const& state);

    /**
     * Writes the Jacobian (JB) for the Newton iteration and passes it to the
     * linear solver
     *
     * @param model the model instance
     * @param state current simulation state
     */
    void prepareLinearSystemB(Model& model, DEStateView const& state);

    /**
     * @brief Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */
    void solveLinearSystem(AmiVector& rhs);

    /**
     * @brief Reinitialize the linear solver
     *
     */
    void reinitialize();

    /**
     * @brief Checks whether the linear system is singular
     *
     * @param model the model instance
     * @param state current simulation state
     * @return boolean indicating whether the linear system is singular
     * (condition number < 1/machine precision)
     */
    bool is_singular(Model& model, DEStateView const& state) const;

  private:
    /** dummy rhs, used as dummy argument when computing J and JB */
    AmiVector xdot_;
    /** dummy state, attached to linear solver */
    AmiVector x_;
    /** dummy adjoint state, used as dummy argument when computing JB */
    AmiVector xB_;
    /** dummy differential adjoint state, used as dummy argument when computing
     * JB */
    AmiVector dxB_;

    /** linear solver */
    std::unique_ptr<SUNLinSolWrapper> linsol_;
};

} // namespace amici

#endif // NEWTON_SOLVER
