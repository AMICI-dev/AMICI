#ifndef AMICI_STEADYSTATEPROBLEM_H
#define AMICI_STEADYSTATEPROBLEM_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/solver_cvodes.h"
#include "amici/forwardproblem.h"
#include <amici/newton_solver.h>

#include <nvector/nvector_serial.h>

#include <functional>
#include <memory>

namespace amici {

class Solver;
class Model;

/**
 * @brief The SteadystateProblem class solves a steady-state problem using
 * Newton's method and falls back to integration on failure.
 */
class SteadystateProblem {
  public:
    /**
     * @brief constructor
     * @param solver Solver instance
     * @param model Model instance
     */
    explicit SteadystateProblem(const Solver &solver,
                                const Model &model);

    /**
     * Tries to determine the steady state of the ODE system by a Newton
     * solver, uses forward intergration, if the Newton solver fails,
     * restarts Newton solver, if integration fails.
     * Computes steady state sensitivities
     *
     * @param solver pointer to the solver object
     * @param model pointer to the model object
     * @param it integer with the index of the current time step
     */
    void workSteadyStateProblem(Solver *solver, Model *model, int it);

    /**
     * Computes the weighted root mean square of xdot
     * the weights are computed according to x:
     * w_i = 1 / ( rtol * x_i + atol )
     *
     * @param x current state
     * @param xdot current rhs
     * @param atol absolute tolerance
     * @param rtol relative tolerance
     * @return root-mean-square norm
     */
    realtype getWrmsNorm(AmiVector const &x,
                         AmiVector const &xdot,
                         realtype atol,
                         realtype rtol);

    /**
     * Checks convergence for state and respective sensitivities
     *
     * @param solver Solver instance
     * @param model instance
     * @return boolean indicating convergence
     */
    bool checkConvergence(const Solver *solver,
                          Model *model);

    /**
     * Runs the Newton solver iterations and checks for convergence to steady
     * state
     *
     * @param model pointer to the model object
     * @param newtonSolver pointer to the NewtonSolver object @type
     * NewtonSolver
     * @param steadystate_try start status of Newton solver
     */
    void applyNewtonsMethod(Model *model, NewtonSolver *newtonSolver,
                            NewtonStatus steadystate_try);

    /**
     * Forward simulation is launched, if Newton solver fails in first try
     *
     * @param solver pointer to the solver object
     * @param model pointer to the model object
     */
    void getSteadystateSimulation(Solver *solver, Model *model);

    /**
     * initialize CVodeSolver instance for preequilibration simulation
     *
     * @param solver pointer to the solver object
     * @param model pointer to the model object
     * @return solver instance
     */
    std::unique_ptr<Solver> createSteadystateSimSolver(const Solver *solver,
                                                       Model *model) const;

    /**
     * @brief store carbon copy of current simulation state variables as SimulationState
     * @param model model carrying the ModelState to be used
     * @param storesensi flag to enable storage of sensitivities
     */
    void storeSimulationState(Model *model, bool storesensi);

    /**
     * @brief returns the stored SimulationState
     * @return stored SimulationState
     */
    const SimulationState getSimulationState() const {
        return state;
    };

    /**
     * @brief returns state at steadystate
     * @return x
     */
    const AmiVector &getState() const {
        return x;
    };


    /**
     * @brief returns state sensitivity at steadystate
     * @return sx
     */
    const AmiVectorArray &getStateSensitivity() const {
        return sx;
    };

     /**
      * @brief Accessor for dJydx
      * @return dJydx
      */
    std::vector<realtype> const& getDJydx() const {
         return dJydx;
     }

    /**
     * @brief Accessor for run_time
     * @return run_time
     */
    double getCPUTime() const { return cpu_time; }

    /**
     * @brief Accessor for newton_status
     * @return newton_status
     */
    NewtonStatus getNewtonStatus() const { return newton_status; }

    /**
     * @brief Accessor for t
     * @return t
     */
    realtype getSteadyStateTime() const { return t; }

    /**
     * @brief Accessor for wrms
     * @return wrms
     */
    realtype getResidualNorm() const { return wrms; }

    /**
     * @brief Accessor for numsteps
     * @return numsteps
     */
    const std::vector<int> &getNumSteps() const { return numsteps; }

    /**
     * @brief Accessor for numlinsteps
     * @return numlinsteps
     */
    const std::vector<int> &getNumLinSteps() const { return numlinsteps; }

    /**
     * @brief computes adjoint updates dJydx according to provided model and expdata
     * @param model Model instance
     * @param edata experimental data
     */
    void getAdjointUpdates(Model &model,
                           const ExpData &edata);



  private:
    /** time variable for simulation steadystate finding */
    realtype t;
    /** newton step */
    AmiVector delta;
    /** error weights */
    AmiVector ewt;
    /** container for relative error calcuation? */
    AmiVector rel_x_newton;
    /** container for absolute error calcuation? */
    AmiVector x_newton;
    /** state vector */
    AmiVector x;
    /** old state vector */
    AmiVector x_old;
    /** differential state vector */
    AmiVector dx;
    /** time derivative state vector */
    AmiVector xdot;
    /** old time derivative state vector */
    AmiVector xdot_old;
    /** state sensitivities */
    AmiVectorArray sx;
    /** state differential sensitivities */
    AmiVectorArray sdx;

    /** weighted root-mean-square error */
    realtype wrms = NAN;

    /** state derivative of data likelihood
     * (dimension nJ x nx x nt, ordering =?) */
    std::vector<realtype> dJydx;

    SimulationState state;

    /** stores diagnostic information about employed number of steps */
    std::vector<int> numsteps;

    /** stores diagnostic information about employed number of linear steps */
    std::vector<int> numlinsteps;

    /** stores diagnostic information about runtime */
    double cpu_time;

    /** stores diagnostic information about execution success*/
    NewtonStatus newton_status = NewtonStatus::failed;
};

} // namespace amici
#endif // STEADYSTATEPROBLEM_H
