#ifndef _am_model_steadystate_sigma_y_h
#define _am_model_steadystate_sigma_y_h

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

void sigma_y_model_steadystate(realtype t, amici::TempData *tdata);


#endif /* _am_model_steadystate_sigma_y_h */
