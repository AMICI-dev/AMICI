#ifndef _am_model_events_dJydsigma_h
#define _am_model_events_dJydsigma_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int dJydsigma_model_events(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata);


#endif /* _am_model_events_dJydsigma_h */
