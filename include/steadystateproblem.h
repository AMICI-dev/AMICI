#ifndef AMICI_STEADYSTATEPROBLEM_H
#define AMICI_STEADYSTATEPROBLEM_H

#include "include/amici_defines.h"
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
    void workSteadyStateProblem(const UserData *udata, TempData *tdata,
                                      ReturnData *rdata, Solver *solver,
                                      Model *model, int it);

    /**
     * applyNewtonsMethod applies Newtons method to the current state x to
     * find the steady state
     */
    void applyNewtonsMethod(const UserData *udata, ReturnData *rdata,
                                  TempData *tdata, Model *model,
                                  NewtonSolver *newtonSolver, int newton_try);

    void getNewtonOutput(TempData *tdata, ReturnData *rdata,
                                Model *model, int newton_status,
                                double run_time, int it);

    void getNewtonSimulation(const UserData *udata, TempData *tdata,
                                   ReturnData *rdata, Solver *solver,
                                   Model *model, int it);
    
    
    /** default constructor
     * @param[in] nx number of state variables
     */
    SteadystateProblem(const int nx) {
        delta = N_VNew_Serial(nx);
        rel_x_newton = N_VNew_Serial(nx);
        x_newton = N_VNew_Serial(nx);
    }
    /** default destructor */
    ~SteadystateProblem(){
        N_VDestroy_Serial(delta);
        N_VDestroy_Serial(rel_x_newton);
        N_VDestroy_Serial(x_newton);
    };
  private:
    /** newton step? */
    N_Vector delta = nullptr;
    /** container for relative error calcuation? */
    N_Vector rel_x_newton = nullptr;
    /** container for absolute error calcuation? */
    N_Vector x_newton = nullptr;
};

} // namespace amici
#endif // STEADYSTATEPROBLEM_H
