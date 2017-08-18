
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int dydx_model_jakstat_adjoint(realtype t, int it, N_Vector x, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
  tdata->dydx[0+1*3] = tdata->p[13]/tdata->p[4];
  tdata->dydx[0+2*3] = (tdata->p[13]*2.0)/tdata->p[4];
  tdata->dydx[1+0*3] = tdata->p[12]/tdata->p[4];
  tdata->dydx[1+1*3] = tdata->p[12]/tdata->p[4];
  tdata->dydx[1+2*3] = (tdata->p[12]*2.0)/tdata->p[4];
return(status);

}


