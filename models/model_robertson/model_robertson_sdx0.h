#ifndef _am_model_robertson_sdx0_h
#define _am_model_robertson_sdx0_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

using namespace amici;

namespace amici {
class UserData;
class ReturnData;
class TempData;
class ExpData;
}

int sdx0_model_robertson(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data);


#endif /* _am_model_robertson_sdx0_h */
