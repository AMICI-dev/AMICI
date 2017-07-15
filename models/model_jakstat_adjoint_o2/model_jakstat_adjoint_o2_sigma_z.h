#ifndef _am_model_jakstat_adjoint_o2_sigma_z_h
#define _am_model_jakstat_adjoint_o2_sigma_z_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int sigma_z_model_jakstat_adjoint_o2(realtype t, int ie, void *user_data, TempData *tdata);


#endif /* _am_model_jakstat_adjoint_o2_sigma_z_h */
