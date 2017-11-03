#ifndef _am_model_neuron_deltax_h
#define _am_model_neuron_deltax_h

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

void deltax_model_neuron(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata);


#endif /* _am_model_neuron_deltax_h */
