
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

int dydx_model_events(realtype t, int it, N_Vector x, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,tdata);
  tdata->dydx[0+0*1] = tdata->p[3];
  tdata->dydx[0+1*1] = tdata->p[3];
  tdata->dydx[0+2*1] = tdata->p[3];
return(status);

}


