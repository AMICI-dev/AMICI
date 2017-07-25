#ifndef _am_model_dirac_xBdot_h
#define _am_model_dirac_xBdot_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int xBdot_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data);


#endif /* _am_model_dirac_xBdot_h */
