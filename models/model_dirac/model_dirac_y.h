#ifndef _am_model_dirac_y_h
#define _am_model_dirac_y_h

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

void y_model_dirac(realtype t, int it, N_Vector x, void *user_data, amici::ReturnData *rdata);


#endif /* _am_model_dirac_y_h */
