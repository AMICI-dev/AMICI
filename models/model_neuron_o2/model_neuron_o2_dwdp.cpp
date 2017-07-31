
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

int dwdp_model_neuron_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(udata->dwdp,0,sizeof(realtype)*0);
status = w_model_neuron_o2(t,x,NULL,user_data);
return(status);

}


