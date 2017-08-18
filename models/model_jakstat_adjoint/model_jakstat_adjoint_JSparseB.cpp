
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int JSparseB_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  SparseSetMatToZero(JB);
  JB->indexvals[0] = 0;
  JB->indexvals[1] = 8;
  JB->indexvals[2] = 0;
  JB->indexvals[3] = 1;
  JB->indexvals[4] = 1;
  JB->indexvals[5] = 2;
  JB->indexvals[6] = 2;
  JB->indexvals[7] = 3;
  JB->indexvals[8] = 3;
  JB->indexvals[9] = 4;
  JB->indexvals[10] = 4;
  JB->indexvals[11] = 5;
  JB->indexvals[12] = 5;
  JB->indexvals[13] = 6;
  JB->indexvals[14] = 6;
  JB->indexvals[15] = 7;
  JB->indexvals[16] = 7;
  JB->indexvals[17] = 8;
  JB->indexptrs[0] = 0;
  JB->indexptrs[1] = 2;
  JB->indexptrs[2] = 4;
  JB->indexptrs[3] = 6;
  JB->indexptrs[4] = 8;
  JB->indexptrs[5] = 10;
  JB->indexptrs[6] = 12;
  JB->indexptrs[7] = 14;
  JB->indexptrs[8] = 16;
  JB->indexptrs[9] = 18;
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  JB->data[0] = udata->p[0]*tdata->w[0];
  JB->data[1] = -(udata->k[1]*udata->p[3])/udata->k[0];
  JB->data[2] = -udata->p[0]*tdata->w[0];
  JB->data[3] = udata->p[1]*tdata->dwdx[0]*2.0;
  JB->data[4] = -udata->p[1]*tdata->dwdx[0];
  JB->data[5] = udata->p[2];
  JB->data[6] = -(udata->k[0]*udata->p[2])/udata->k[1];
  JB->data[7] = udata->p[3];
  JB->data[8] = udata->p[3]*-2.0;
  JB->data[9] = udata->p[3];
  JB->data[10] = -udata->p[3];
  JB->data[11] = udata->p[3];
  JB->data[12] = -udata->p[3];
  JB->data[13] = udata->p[3];
  JB->data[14] = -udata->p[3];
  JB->data[15] = udata->p[3];
  JB->data[16] = -udata->p[3];
  JB->data[17] = udata->p[3];
return(status);

}


