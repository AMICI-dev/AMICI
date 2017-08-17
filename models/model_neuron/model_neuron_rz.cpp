
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_neuron_w.h"

int rz_model_neuron(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_neuron(t,x,NULL,tdata);
  rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0] = x_tmp[0]-3.0E1;
return(status);

}


