#ifndef _am_model_neuron_o2_dsigma_zdp_h
#define _am_model_neuron_o2_dsigma_zdp_h

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int dsigma_zdp_model_neuron_o2(realtype t, int ie, TempData *tdata);


#endif /* _am_model_neuron_o2_dsigma_zdp_h */
