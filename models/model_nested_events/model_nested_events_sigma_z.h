#ifndef _am_model_nested_events_sigma_z_h
#define _am_model_nested_events_sigma_z_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int sigma_z_model_nested_events(realtype t, int ie, TempData *tdata);


#endif /* _am_model_nested_events_sigma_z_h */
