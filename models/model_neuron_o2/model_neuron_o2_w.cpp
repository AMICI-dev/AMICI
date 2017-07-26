
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

int w_model_neuron_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(udata->w,0,sizeof(realtype)*2);
  udata->w[0] = x_tmp[0]*(2.0/2.5E1);
  udata->w[1] = udata->w[0]+5.0;
return(status);

}


