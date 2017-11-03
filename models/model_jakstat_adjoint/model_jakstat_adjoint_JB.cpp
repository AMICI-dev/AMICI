
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

using namespace amici;

void JB_model_jakstat_adjoint(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *xB_tmp = nullptr;
if(xB)
    xB_tmp = N_VGetArrayPointer(xB);
realtype *dxB_tmp = nullptr;
if(dxB)
    dxB_tmp = N_VGetArrayPointer(dxB);
realtype *xBdot_tmp = nullptr;
if(xBdot)
    xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(JB->data,0,sizeof(realtype)*81);
w_model_jakstat_adjoint(t,x,NULL,tdata);
dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  JB->data[0+0*9] = tdata->p[0]*tdata->w[0];
  JB->data[0+1*9] = -tdata->p[0]*tdata->w[0];
  JB->data[1+1*9] = tdata->p[1]*tdata->dwdx[0]*2.0;
  JB->data[1+2*9] = -tdata->p[1]*tdata->dwdx[0];
  JB->data[2+2*9] = tdata->p[2];
  JB->data[2+3*9] = -(udata->k[0]*tdata->p[2])/udata->k[1];
  JB->data[3+3*9] = tdata->p[3];
  JB->data[3+4*9] = tdata->p[3]*-2.0;
  JB->data[4+4*9] = tdata->p[3];
  JB->data[4+5*9] = -tdata->p[3];
  JB->data[5+5*9] = tdata->p[3];
  JB->data[5+6*9] = -tdata->p[3];
  JB->data[6+6*9] = tdata->p[3];
  JB->data[6+7*9] = -tdata->p[3];
  JB->data[7+7*9] = tdata->p[3];
  JB->data[7+8*9] = -tdata->p[3];
  JB->data[8+0*9] = -(udata->k[1]*tdata->p[3])/udata->k[0];
  JB->data[8+8*9] = tdata->p[3];
return;

}


