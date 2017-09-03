#ifndef _am_model_nested_events_dydx_h
#define _am_model_nested_events_dydx_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int dydx_model_nested_events(realtype t, int it, N_Vector x, TempData *tdata);


#endif /* _am_model_nested_events_dydx_h */
