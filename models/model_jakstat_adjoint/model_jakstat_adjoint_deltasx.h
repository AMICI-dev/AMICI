#ifndef _am_model_jakstat_adjoint_deltasx_h
#define _am_model_jakstat_adjoint_deltasx_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int deltasx_model_jakstat_adjoint(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, TempData *tdata);


#endif /* _am_model_jakstat_adjoint_deltasx_h */
