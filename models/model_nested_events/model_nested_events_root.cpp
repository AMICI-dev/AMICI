
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_nested_events_w.h"

using namespace amici;

void root_model_nested_events(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
w_model_nested_events(t,x,NULL,tdata);
  root[0] = -t+tdata->p[2];
  root[1] = x_tmp[0]-1.0;
  root[2] = t-tdata->p[2];
return;

}


