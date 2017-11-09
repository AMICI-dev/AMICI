#ifndef _am_model_neuron_qBdot_h
#define _am_model_neuron_qBdot_h

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

void qBdot_model_neuron(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data);


#endif /* _am_model_neuron_qBdot_h */
