#ifndef _am_model_jakstat_adjoint_o2_JvB_h
#define _am_model_jakstat_adjoint_o2_JvB_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int JvB_model_jakstat_adjoint_o2(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB);


#endif /* _am_model_jakstat_adjoint_o2_JvB_h */
