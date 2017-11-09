#ifndef _am_model_neuron_o2_dJzdz_h
#define _am_model_neuron_o2_dJzdz_h

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

void dJzdz_model_neuron_o2(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);


#endif /* _am_model_neuron_o2_dJzdz_h */
