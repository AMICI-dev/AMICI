#ifndef _am_model_nested_events_dzdx_h
#define _am_model_nested_events_dzdx_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int dzdx_model_nested_events(realtype t, int ie, N_Vector x, TempData *tdata);


#endif /* _am_model_nested_events_dzdx_h */
