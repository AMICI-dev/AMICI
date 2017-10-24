#ifndef _am_model_jakstat_adjoint_sigma_y_h
#define _am_model_jakstat_adjoint_sigma_y_h

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

int sigma_y_model_jakstat_adjoint(realtype t, amici::TempData *tdata);


#endif /* _am_model_jakstat_adjoint_sigma_y_h */
