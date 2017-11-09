#ifndef _am_model_steadystate_dJydsigma_h
#define _am_model_steadystate_dJydsigma_h

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

void dJydsigma_model_steadystate(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);


#endif /* _am_model_steadystate_dJydsigma_h */
