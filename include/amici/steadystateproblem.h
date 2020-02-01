#ifndef AMICI_STEADYSTATEPROBLEM_H
#define AMICI_STEADYSTATEPROBLEM_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/solver_cvodes.h"
#include <amici/newton_solver.h>

#include <nvector/nvector_serial.h>

#include <functional>
#include <memory>

namespace amici {

class ExpData;
class ReturnData;
class Solver;
class Model;

/**
 * @brief The SteadystateProblem class solves a steady-state problem using
 * Newton's method and falls back to integration on failure.
 */

class SteadystateProblem {
  public:
    void workSteadyStateProblem(ReturnData *rdata, Solver *solver,
                                      Model *model, int it);

    void workSteadyStateBackwardProblem(ReturnData *rdata, Solver *solver,
                                        const ExpData *edata, Model *model);

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
                         realtype rtol
                         );

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
     * @param rdata pointer to the return data object
     * @param model pointer to the AMICI model object
     * @param newtonSolver pointer to the NewtonSolver object @type
     * NewtonSolver
     * @param steadystate_try start status of Newton solver
     */
    void applyNewtonsMethod(ReturnData *rdata, Model *model,
                            NewtonSolver *newtonSolver,
                            NewtonStatus steadystate_try);
    /**
     * Stores output of workSteadyStateProblem in return data
     *
     * @param newton_status integer flag indicating when a steady state was
     * found
     * @param run_time double coputation time of the solver in milliseconds
     * @param rdata pointer to the return data instance
     * @param model pointer to the model instance
     * @param it current timepoint index, <0 indicates preequilibration
     */
    void writeNewtonOutput(ReturnData *rdata, const Model *model,
                           NewtonStatus newton_status, double run_time, int it);

    /**
     * Forward simulation is launched, if Newton solver fails in first try
     *
     * @param solver pointer to the AMICI solver object
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */
    void getSteadystateSimulation(ReturnData *rdata, Solver *solver,
                                  Model *model);

    /**
     * initialize CVodeSolver instance for preequilibration simulation
     *
     * @param solver pointer to the AMICI solver object
     * @param model pointer to the AMICI model object
     * @return solver instance
     */
    std::unique_ptr<Solver> createSteadystateSimSolver(const Solver *solver,
                                                       Model *model) const;

    /**
     * @brief constructor
     * @param solver pointer to Solver instance
     * @param x0 initial state
     */
    explicit SteadystateProblem(const Solver *solver, const AmiVector &x0);
    
    /**
     * @brief routine that writes solutions of steadystate problem to target
     vectors
     * @param t final timepoint
     * @param x steadystate state
     * @param sx steadystate state sensitivity
     */
    void writeSolution(realtype *t, AmiVector &x, AmiVectorArray &sx) const;

    /**
     * @brief routine that writes computed quadrature state vector
     * @param xQB quadrature state vector
     */
    void writeSolutionBackward(AmiVector &xQB) const;

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
    /** quadrature state vector */
    AmiVector xQB;

    /** weighted root-mean-square error */
    realtype wrms = NAN;

};

} // namespace amici
#endif // STEADYSTATEPROBLEM_H
