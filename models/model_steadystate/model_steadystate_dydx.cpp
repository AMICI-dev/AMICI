
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include "model_steadystate_w.h"

int dydx_model_steadystate(realtype t, int it, N_Vector x, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,tdata);
  tdata->dydx[0+0*3] = 1.0;
  tdata->dydx[1+1*3] = 1.0;
  tdata->dydx[2+2*3] = 1.0;
return(status);

}


