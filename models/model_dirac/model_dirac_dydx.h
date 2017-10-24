#ifndef _am_model_dirac_dydx_h
#define _am_model_dirac_dydx_h

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

int dydx_model_dirac(realtype t, int it, N_Vector x, amici::TempData *tdata);


#endif /* _am_model_dirac_dydx_h */
