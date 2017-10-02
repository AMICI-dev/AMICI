
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

int dydx_model_neuron_o2(realtype t, int it, N_Vector x, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
status = w_model_neuron_o2(t,x,NULL,tdata);
  tdata->dydx[0+0*5] = 1.0;
  tdata->dydx[1+2*5] = 1.0;
  tdata->dydx[2+4*5] = 1.0;
  tdata->dydx[3+6*5] = 1.0;
  tdata->dydx[4+8*5] = 1.0;
return(status);

}


