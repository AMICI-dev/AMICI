
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

using namespace amici;

void JB_model_neuron_o2(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
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
  memset(JB->data,0,sizeof(realtype)*100);
w_model_neuron_o2(t,x,NULL,tdata);
dwdx_model_neuron_o2(t,x,NULL,user_data);
  JB->data[0+0*10] = x_tmp[0]*(-2.0/2.5E1)-5.0;
  JB->data[0+1*10] = -tdata->p[0]*tdata->p[1];
  JB->data[1+0*10] = 1.0;
  JB->data[1+1*10] = tdata->p[0];
  JB->data[2+0*10] = -x_tmp[2]*tdata->dwdx[1];
  JB->data[2+1*10] = -tdata->p[1];
  JB->data[2+2*10] = -tdata->w[1];
  JB->data[2+3*10] = -tdata->p[0]*tdata->p[1];
  JB->data[3+1*10] = 1.0;
  JB->data[3+2*10] = 1.0;
  JB->data[3+3*10] = tdata->p[0];
  JB->data[4+0*10] = -x_tmp[4]*tdata->dwdx[1];
  JB->data[4+1*10] = -tdata->p[0];
  JB->data[4+4*10] = -tdata->w[1];
  JB->data[4+5*10] = -tdata->p[0]*tdata->p[1];
  JB->data[5+4*10] = 1.0;
  JB->data[5+5*10] = tdata->p[0];
  JB->data[6+0*10] = -x_tmp[6]*tdata->dwdx[1];
  JB->data[6+6*10] = -tdata->w[1];
  JB->data[6+7*10] = -tdata->p[0]*tdata->p[1];
  JB->data[7+6*10] = 1.0;
  JB->data[7+7*10] = tdata->p[0];
  JB->data[8+0*10] = -x_tmp[8]*tdata->dwdx[1];
  JB->data[8+8*10] = -tdata->w[1];
  JB->data[8+9*10] = -tdata->p[0]*tdata->p[1];
  JB->data[9+8*10] = 1.0;
  JB->data[9+9*10] = tdata->p[0];
return;

}


