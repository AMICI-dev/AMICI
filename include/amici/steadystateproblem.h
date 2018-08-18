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

    /**
     * applyNewtonsMethod applies Newtons method to the current state x to
     * find the steady state
     */
    void applyNewtonsMethod(ReturnData *rdata, Model *model,
                                  NewtonSolver *newtonSolver, int newton_try);

    void getNewtonOutput(ReturnData *rdata, const Model *model,
                         int newton_status, double run_time, int it);

    void getSteadystateSimulation(ReturnData *rdata, Solver *solver,
                                  Model *model, int it);
    
    std::unique_ptr<CVodeSolver> createSteadystateSimSolver(Solver *solver, Model *model, realtype tstart);

    /** default constructor
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param sx pointer to state sensitivity variables
     */
    SteadystateProblem(realtype *t, AmiVector *x, AmiVectorArray *sx) :
    delta(x->getLength()),
    rel_x_newton(x->getLength()),
    x_newton(x->getLength()),
    x_old(x->getLength()),
    dx(x->getLength()),
    xdot(x->getLength()),
    xdot_old(x->getLength()),
    sdx(x->getLength(),sx->getLength())
    {
        this->t = t;
        this->x = x;
        this->sx = sx;
    }

  private:
    realtype *t;
    /** newton step? */
    AmiVector delta;
    /** container for relative error calcuation? */
    AmiVector rel_x_newton;
    /** container for absolute error calcuation? */
    AmiVector x_newton;
    /** state vector */
    AmiVector *x;
    /** old state vector */
    AmiVector x_old;
    /** differential state vector */
    AmiVector dx;
    /** time derivative state vector */
    AmiVector xdot;
    /** old time derivative state vector */
    AmiVector xdot_old;
    /** state sensitivities */
    AmiVectorArray *sx;
    /** state differential sensitivities */
    AmiVectorArray sdx;
    
};

} // namespace amici
#endif // STEADYSTATEPROBLEM_H
