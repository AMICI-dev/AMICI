#ifndef _am_model_neuron_srz_h
#define _am_model_neuron_srz_h

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

void srz_model_neuron(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata);


#endif /* _am_model_neuron_srz_h */
