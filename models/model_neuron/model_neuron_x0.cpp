
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_w.h"

int x0_model_neuron(N_Vector x0, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x0_tmp = N_VGetArrayPointer(x0);
memset(x0_tmp,0,sizeof(realtype)*2);
realtype t = udata->tstart;
  x0_tmp[0] = udata->k[0];
  x0_tmp[1] = udata->k[0]*udata->p[1];
return(status);

}


