#ifndef FORWARDPROBLEM_H
#define FORWARDPROBLEM_H

#include "include/amici_defines.h"
namespace amici {

class UserData;
class TempData;
class ReturnData;
class ExpData;
class Solver;
class Model;

/**
 * @brief The ForwardProblem class groups all functions for solving the
 * backwards problem.
 * Has only static members.
 */
class ForwardProblem {
  public:
    static void workForwardProblem(const UserData *udata, TempData *tdata,
                                  ReturnData *rdata, const ExpData *edata,
                                  Model *model);

    static void handleEvent(realtype *tlastroot, const UserData *udata,
                           ReturnData *rdata, const ExpData *edata,
                           TempData *tdata, int seflag, Solver *solver,
                           Model *model);

    static void storeJacobianAndDerivativeInReturnData(TempData *tdata,
                                                      ReturnData *rdata,
                                                      Model *model);

    static void getEventOutput(const UserData *udata, ReturnData *rdata,
                              const ExpData *edata, TempData *tdata,
                              Model *model);

    static void prepEventSensis(int ie, ReturnData *rdata, const ExpData *edata,
                               TempData *tdata, Model *model);

    static void getEventSensisFSA(int ie, ReturnData *rdata,
                                 const ExpData *edata, TempData *tdata,
                                 Model *model);

    static void handleDataPoint(int it, const UserData *udata, ReturnData *rdata,
                               const ExpData *edata, TempData *tdata,
                               Solver *solver, Model *model);

    static void getDataOutput(int it, const UserData *udata, ReturnData *rdata,
                             const ExpData *edata, TempData *tdata,
                             Solver *solver, Model *model);

    static void prepDataSensis(int it, ReturnData *rdata, const ExpData *edata,
                              TempData *tdata, Model *model);

    static void getDataSensisFSA(int it, const UserData *udata,
                                ReturnData *rdata, const ExpData *edata,
                                TempData *tdata, Solver *solver, Model *model);

    static void applyEventBolus(TempData *tdata, Model *model);

    static void applyEventSensiBolusFSA(TempData *tdata, Model *model);

    static void updateHeaviside(TempData *tdata, const int ne);

  private:
    ForwardProblem();
};


} // namespace amici

#endif // FORWARDPROBLEM_H
