#ifndef _am_model_robertson_Jy_h
#define _am_model_robertson_Jy_h

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

void Jy_model_robertson(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);


#endif /* _am_model_robertson_Jy_h */
