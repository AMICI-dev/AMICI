
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_neuron_o2_w.h"

int y_model_neuron_o2(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_neuron_o2(t,x,NULL,user_data);
  rdata->y[it + udata->nt*0] = x_tmp[0];
  rdata->y[it + udata->nt*1] = x_tmp[2];
  rdata->y[it + udata->nt*2] = x_tmp[4];
  rdata->y[it + udata->nt*3] = x_tmp[6];
  rdata->y[it + udata->nt*4] = x_tmp[8];
return(status);

}


