#ifndef _am_model_neuron_o2_sigma_y_h
#define _am_model_neuron_o2_sigma_y_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int sigma_y_model_neuron_o2(realtype t, TempData *tdata);


#endif /* _am_model_neuron_o2_sigma_y_h */
