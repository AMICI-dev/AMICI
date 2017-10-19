#ifndef _am_model_steadystate_deltasx_h
#define _am_model_steadystate_deltasx_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

using namespace amici;

namespace amici {
class UserData;
class ReturnData;
class TempData;
class ExpData;
}

int deltasx_model_steadystate(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, amici::TempData *tdata);


#endif /* _am_model_steadystate_deltasx_h */
