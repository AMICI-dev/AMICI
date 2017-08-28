#ifndef STEADYSTATEPROBLEM_H
#define STEADYSTATEPROBLEM_H

#include "include/amici_defines.h"
#include <sundials/sundials_nvector.h>

class UserData;
class TempData;
class ReturnData;
class ExpData;
class Solver;
class Model;
class NewtonSolver;

/**
 * @brief The SteadystateProblem class solves a steady-state problem using
 * Newton's method and falls back to integration on failure.
 */

class SteadystateProblem {
  public:

    static int workSteadyStateProblem(UserData *udata, TempData *tdata,
                                      ReturnData *rdata, Solver *solver,
                                      Model *model, int it);

    /**
     * applyNewtonsMethod applies Newtons method to the current state x to
     * find the steady state
     */
    static int applyNewtonsMethod(UserData *udata, ReturnData *rdata,
                                  TempData *tdata, Model *model,
                                  NewtonSolver *newtonSolver, int newton_try);

    static void getNewtonOutput(TempData *tdata, ReturnData *rdata, Model *model,
                               int newton_status, double run_time);

    static int getNewtonSimulation(UserData *udata, TempData *tdata,
                                   ReturnData *rdata, Solver *solver,
                                   Model *model);
    
    static int linsolveSPBCG(UserData *udata, ReturnData *rdata,
                             TempData *tdata, Model *model, int ntry,
                             int nnewt, N_Vector ns_delta);

  private:
    SteadystateProblem();
};

#endif // STEADYSTATEPROBLEM_H
