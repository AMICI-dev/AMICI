#ifndef BACKWARDPROBLEM_H
#define BACKWARDPROBLEM_H

#include "include/amici_defines.h"

class UserData;
class TempData;
class ReturnData;
class ExpData;
class Solver;
class Model;

class BackwardProblem
{
public:

    static int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, Model *model);

    static int handleEventB(int iroot, TempData *tdata, Model *model);

    static int handleDataPointB(int it, ReturnData *rdata, TempData *tdata, Solver *solver, Model *model);

    static int updateHeavisideB(int iroot, TempData *tdata, int ne);

    static realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, Model *model);

private:
    BackwardProblem();
};

#endif // BACKWARDPROBLEM_H
