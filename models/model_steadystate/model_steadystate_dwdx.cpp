
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

using namespace amici;

void dwdx_model_steadystate(realtype t, N_Vector x, N_Vector dx, void *user_data) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
memset(tdata->dwdx,0,sizeof(realtype)*2);
w_model_steadystate(t,x,NULL,tdata);
  tdata->dwdx[0] = x_tmp[0]*2.0;
  tdata->dwdx[1] = tdata->p[3];
return;

}


