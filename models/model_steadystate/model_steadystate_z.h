#ifndef _am_model_steadystate_z_h
#define _am_model_steadystate_z_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int z_model_steadystate(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata);


#endif /* _am_model_steadystate_z_h */
