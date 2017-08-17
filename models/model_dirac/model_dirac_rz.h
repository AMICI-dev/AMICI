#ifndef _am_model_dirac_rz_h
#define _am_model_dirac_rz_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int rz_model_dirac(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata);


#endif /* _am_model_dirac_rz_h */
