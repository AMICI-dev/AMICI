#ifndef AMICI_STEADYSTATEPROBLEM_H
#define AMICI_STEADYSTATEPROBLEM_H

#include "include/amici_defines.h"
#include "include/amici_vector.h"
#include <nvector/nvector_serial.h>
#include <include/newton_solver.h>

namespace amici {

class UserData;
class TempData;
class ReturnData;
class ExpData;
class Solver;
class Model;
class NewtonSolver;
class NewtonSolverDense;
class NewtonSolverSparse;
class NewtonSolverIterative;

/**
 * @brief The SteadystateProblem class solves a steady-state problem using
 * Newton's method and falls back to integration on failure.
 */

class SteadystateProblem {
  public:
    void workSteadyStateProblem(const UserData *udata,
                                      ReturnData *rdata, Solver *solver,
                                      Model *model, int it);

    /**
     * applyNewtonsMethod applies Newtons method to the current state x to
     * find the steady state
     */
    void applyNewtonsMethod(const UserData *udata, ReturnData *rdata, Model *model,
                                  NewtonSolver *newtonSolver, int newton_try);

    void getNewtonOutput(ReturnData *rdata,
                                Model *model, int newton_status,
                                double run_time, int it);

    void getNewtonSimulation(const UserData *udata,
                                   ReturnData *rdata, Solver *solver,
                                   Model *model, int it);
    
    
    /** default constructor
     * @param[in] nx number of state variables
     */
    SteadystateProblem(const int nx, const int nplist) :
    delta(nx), rel_x_newton(nx), x_newton(nx),
    x(nx), dx(nx), xdot(nx),
    x_old(nx), xdot_old(nx), sx(nx,nplist), sdx(nx,nplist){};
    
    /** default destructor */
    ~SteadystateProblem(){};
  private:
    realtype t;
    /** newton step? */
    AmiVector delta;
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
    
};

} // namespace amici
#endif // STEADYSTATEPROBLEM_H
