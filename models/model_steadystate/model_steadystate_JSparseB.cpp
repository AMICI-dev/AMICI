
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_dwdx.h"
#include "model_steadystate_w.h"

int JSparseB_model_steadystate(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  SparseSetMatToZero(JB);
  JB->indexvals[0] = 0;
  JB->indexvals[1] = 1;
  JB->indexvals[2] = 2;
  JB->indexvals[3] = 0;
  JB->indexvals[4] = 1;
  JB->indexvals[5] = 2;
  JB->indexvals[6] = 0;
  JB->indexvals[7] = 1;
  JB->indexvals[8] = 2;
  JB->indexptrs[0] = 0;
  JB->indexptrs[1] = 3;
  JB->indexptrs[2] = 6;
  JB->indexptrs[3] = 9;
status = w_model_steadystate(t,x,NULL,tdata);
status = dwdx_model_steadystate(t,x,NULL,user_data);
  JB->data[0] = udata->p[1]*x_tmp[1]+udata->p[0]*tdata->dwdx[0]*2.0;
  JB->data[1] = udata->p[2]*-2.0+udata->p[1]*x_tmp[0];
  JB->data[2] = -tdata->dwdx[1];
  JB->data[3] = udata->p[1]*x_tmp[1]-udata->p[0]*tdata->dwdx[0];
  JB->data[4] = udata->p[2]+udata->p[1]*x_tmp[0];
  JB->data[5] = -tdata->dwdx[1];
  JB->data[6] = -udata->p[1]*x_tmp[1];
  JB->data[7] = -udata->p[1]*x_tmp[0];
  JB->data[8] = udata->k[3]+tdata->dwdx[1];
return(status);

}


