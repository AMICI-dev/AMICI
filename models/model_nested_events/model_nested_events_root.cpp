
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_nested_events_w.h"

int root_model_nested_events(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_nested_events(t,x,NULL,tdata);
  root[0] = -t+tdata->p[2];
  root[1] = x_tmp[0]-1.0;
  root[2] = t-tdata->p[2];
return(status);

}


