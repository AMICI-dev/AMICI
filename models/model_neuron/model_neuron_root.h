#ifndef _am_model_neuron_root_h
#define _am_model_neuron_root_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int root_model_neuron(realtype t, N_Vector x, realtype *root, void *user_data);


#endif /* _am_model_neuron_root_h */
