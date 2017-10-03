
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_dwdx.h"
#include "model_events_w.h"

int JSparseB_model_events(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
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
realtype *xB_tmp = nullptr;
if(xB)
    xB_tmp = N_VGetArrayPointer(xB);
realtype *dxB_tmp = nullptr;
if(dxB)
    dxB_tmp = N_VGetArrayPointer(dxB);
realtype *xBdot_tmp = nullptr;
if(xBdot)
    xBdot_tmp = N_VGetArrayPointer(xBdot);
  SparseSetMatToZero(JB);
  JB->indexvals[0] = 0;
  JB->indexvals[1] = 0;
  JB->indexvals[2] = 1;
  JB->indexvals[3] = 2;
  JB->indexptrs[0] = 0;
  JB->indexptrs[1] = 1;
  JB->indexptrs[2] = 3;
  JB->indexptrs[3] = 4;
status = w_model_events(t,x,NULL,tdata);
status = dwdx_model_events(t,x,NULL,user_data);
  JB->data[0] = tdata->h[3]*tdata->p[0];
  JB->data[1] = -tdata->p[1]*exp(t*(-1.0/1.0E1));
  JB->data[2] = tdata->p[2];
  JB->data[3] = 1.0;
return(status);

}


