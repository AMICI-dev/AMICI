#ifndef _am_model_dirac_dsigma_ydp_h
#define _am_model_dirac_dsigma_ydp_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int dsigma_ydp_model_dirac(realtype t, TempData *tdata);


#endif /* _am_model_dirac_dsigma_ydp_h */
