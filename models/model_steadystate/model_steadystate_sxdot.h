#ifndef _am_model_steadystate_sxdot_h
#define _am_model_steadystate_sxdot_h

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

void sxdot_model_steadystate(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


#endif /* _am_model_steadystate_sxdot_h */
