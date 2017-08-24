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
 * Newton's method
 * and falls back to integration on failure.
 */

class SteadystateProblem {
  public:

    static int workSteadyStateProblem(UserData *udata, TempData *tdata,
                                      ReturnData *rdata, int it, Solver *solver,
                                      Model *model);

    /**
     * applyNewtonsMethod applies Newtons method to the current state x to
     * find the steady state
     */
    static int applyNewtonsMethod(UserData *udata, ReturnData *rdata,
                                  TempData *tdata, int newton_try, Model *model,
                                  NewtonSolver *newtonSolver);

    static int getNewtonOutput(TempData *tdata, ReturnData *rdata,
                               int newton_status, double run_time, int nx);

    static int getNewtonSimulation(UserData *udata, TempData *tdata,
                                   ReturnData *rdata, Solver *solver,
                                   Model *model);

    static int linsolveSPBCG(UserData *udata, ReturnData *rdata,
                             TempData *tdata, int ntry, int nnewt,
                             N_Vector ns_delta, Model *model);

  private:
    SteadystateProblem();
};

#endif // STEADYSTATEPROBLEM_H
