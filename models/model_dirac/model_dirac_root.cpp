
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_dirac_w.h"

int root_model_dirac(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
status = w_model_dirac(t,x,NULL,tdata);
  root[0] = t-tdata->p[1];
  root[1] = -t+tdata->p[1];
return(status);

}


