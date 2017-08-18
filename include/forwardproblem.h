#ifndef FORWARDPROBLEM_H
#define FORWARDPROBLEM_H

#include "include/amici_defines.h"

class UserData;
class TempData;
class ReturnData;
class ExpData;
class Solver;
class Model;

/**
 * @brief The ForwardProblem class groups all functions for solving the backwards problem.
 * Has only static members.
 */
class ForwardProblem
{
public:

    static int workForwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, Solver *solver, Model *model);

    static int handleEvent(realtype *tlastroot, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, int seflag, Solver *solver, Model *model);

    static int storeJacobianAndDerivativeInReturnData(UserData *udata, TempData *tdata,  ReturnData *rdata, Model *model);

    static int getEventOutput(realtype *tlastroot, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Model *model);

    static int prepEventSensis(int ie, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Model *model);

    static int getEventSensisFSA(int ie, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Model *model);

    static int handleDataPoint(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Solver *solver, Model *model);

    static int applyEventBolus( UserData *udata, TempData *tdata, Model *model);

    static int applyEventSensiBolusFSA(UserData *udata, TempData *tdata, Model *model);

    static int updateHeaviside(TempData *tdata, int ne);

private:
    ForwardProblem();
};

#endif // FORWARDPROBLEM_H
