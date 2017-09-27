#ifndef _am_model_nested_events_qBdot_h
#define _am_model_nested_events_qBdot_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int qBdot_model_nested_events(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data);


#endif /* _am_model_nested_events_qBdot_h */
