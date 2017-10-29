#ifndef _am_model_dirac_stau_h
#define _am_model_dirac_stau_h

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

void stau_model_dirac(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata);


#endif /* _am_model_dirac_stau_h */
