
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_robertson_dwdx.h"
#include "model_robertson_w.h"

using namespace amici;

void dfdx_model_robertson(realtype t, N_Vector x, N_Vector dx, void *user_data) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
memset(tdata->dfdx,0,sizeof(realtype)*9);
dwdx_model_robertson(t,x,dx,user_data);
  tdata->dfdx[0+0*3] = -tdata->p[0];
  tdata->dfdx[0+1*3] = tdata->dwdx[0];
  tdata->dfdx[0+2*3] = tdata->dwdx[1];
  tdata->dfdx[1+0*3] = tdata->p[0];
  tdata->dfdx[1+1*3] = -tdata->dwdx[0]-tdata->p[2]*x_tmp[1]*2.0;
  tdata->dfdx[1+2*3] = -tdata->dwdx[1];
  tdata->dfdx[2+0*3] = 1.0;
  tdata->dfdx[2+1*3] = 1.0;
  tdata->dfdx[2+2*3] = 1.0;
return;

}


