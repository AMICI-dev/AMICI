#ifndef _am_model_jakstat_adjoint_o2_JSparse_h
#define _am_model_jakstat_adjoint_o2_JSparse_h

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

void JSparse_model_jakstat_adjoint_o2(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


#endif /* _am_model_jakstat_adjoint_o2_JSparse_h */
