#ifndef _am_model_neuron_sigma_z_h
#define _am_model_neuron_sigma_z_h

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

int sigma_z_model_neuron(realtype t, int ie, amici::TempData *tdata);


#endif /* _am_model_neuron_sigma_z_h */
