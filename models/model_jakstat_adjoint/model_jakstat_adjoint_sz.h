#ifndef _am_model_jakstat_adjoint_sz_h
#define _am_model_jakstat_adjoint_sz_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int sz_model_jakstat_adjoint(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata);


#endif /* _am_model_jakstat_adjoint_sz_h */
