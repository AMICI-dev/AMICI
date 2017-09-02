#ifndef _am_model_events_stau_h
#define _am_model_events_stau_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int stau_model_events(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata);


#endif /* _am_model_events_stau_h */
