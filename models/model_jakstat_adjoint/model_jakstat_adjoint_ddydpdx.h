#ifndef _am_model_jakstat_adjoint_ddydpdx_h
#define _am_model_jakstat_adjoint_ddydpdx_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int ddydpdx_model_jakstat_adjoint(realtype t, int it, N_Vector x, TempData *tdata);


#endif /* _am_model_jakstat_adjoint_ddydpdx_h */
