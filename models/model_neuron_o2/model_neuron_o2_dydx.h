#ifndef _am_model_neuron_o2_dydx_h
#define _am_model_neuron_o2_dydx_h

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

void dydx_model_neuron_o2(realtype t, int it, N_Vector x, amici::TempData *tdata);


#endif /* _am_model_neuron_o2_dydx_h */
