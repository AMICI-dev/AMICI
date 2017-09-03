#ifndef _am_model_nested_events_Jv_h
#define _am_model_nested_events_Jv_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int Jv_model_nested_events(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp);


#endif /* _am_model_nested_events_Jv_h */
