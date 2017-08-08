
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_neuron_w.h"

int rz_model_neuron(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_neuron(t,x,NULL,user_data);
  rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0] = x_tmp[0]-3.0E1;
return(status);

}


