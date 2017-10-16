#ifndef _am_model_jakstat_adjoint_xdot_h
#define _am_model_jakstat_adjoint_xdot_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int xdot_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data);


#endif /* _am_model_jakstat_adjoint_xdot_h */
