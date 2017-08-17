#ifndef _am_model_steadystate_deltax_h
#define _am_model_steadystate_deltax_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int deltax_model_steadystate(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, TempData *tdata);


#endif /* _am_model_steadystate_deltax_h */
