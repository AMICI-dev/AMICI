#ifndef _am_model_neuron_JDiag_h
#define _am_model_neuron_JDiag_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int JDiag_model_neuron(realtype t, N_Vector JDiag, N_Vector x, void *user_data);


#endif /* _am_model_neuron_JDiag_h */
