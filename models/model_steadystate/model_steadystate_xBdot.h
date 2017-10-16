#ifndef _am_model_steadystate_xBdot_h
#define _am_model_steadystate_xBdot_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int xBdot_model_steadystate(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, void *user_data);


#endif /* _am_model_steadystate_xBdot_h */
