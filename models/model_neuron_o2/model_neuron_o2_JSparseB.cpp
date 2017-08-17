
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

int JSparseB_model_neuron_o2(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
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
  JB->indexvals[2] = 0;
  JB->indexvals[3] = 1;
  JB->indexvals[4] = 0;
  JB->indexvals[5] = 2;
  JB->indexvals[6] = 3;
  JB->indexvals[7] = 0;
  JB->indexvals[8] = 1;
  JB->indexvals[9] = 2;
  JB->indexvals[10] = 3;
  JB->indexvals[11] = 0;
  JB->indexvals[12] = 4;
  JB->indexvals[13] = 5;
  JB->indexvals[14] = 0;
  JB->indexvals[15] = 4;
  JB->indexvals[16] = 5;
  JB->indexvals[17] = 0;
  JB->indexvals[18] = 6;
  JB->indexvals[19] = 7;
  JB->indexvals[20] = 6;
  JB->indexvals[21] = 7;
  JB->indexvals[22] = 0;
  JB->indexvals[23] = 8;
  JB->indexvals[24] = 9;
  JB->indexvals[25] = 8;
  JB->indexvals[26] = 9;
  JB->indexptrs[0] = 0;
  JB->indexptrs[1] = 2;
  JB->indexptrs[2] = 4;
  JB->indexptrs[3] = 7;
  JB->indexptrs[4] = 11;
  JB->indexptrs[5] = 14;
  JB->indexptrs[6] = 17;
  JB->indexptrs[7] = 20;
  JB->indexptrs[8] = 22;
  JB->indexptrs[9] = 25;
  JB->indexptrs[10] = 27;
status = w_model_neuron_o2(t,x,NULL,tdata);
status = dwdx_model_neuron_o2(t,x,NULL,user_data);
  JB->data[0] = x_tmp[0]*(-2.0/2.5E1)-5.0;
  JB->data[1] = 1.0;
  JB->data[2] = -udata->p[0]*udata->p[1];
  JB->data[3] = udata->p[0];
  JB->data[5] = -udata->w[1];
  JB->data[6] = 1.0;
  JB->data[9] = -udata->p[0]*udata->p[1];
  JB->data[10] = udata->p[0];
  JB->data[12] = -udata->w[1];
  JB->data[13] = 1.0;
  JB->data[15] = -udata->p[0]*udata->p[1];
  JB->data[16] = udata->p[0];
  JB->data[18] = -udata->w[1];
  JB->data[19] = 1.0;
  JB->data[20] = -udata->p[0]*udata->p[1];
  JB->data[21] = udata->p[0];
  JB->data[23] = -udata->w[1];
  JB->data[24] = 1.0;
  JB->data[25] = -udata->p[0]*udata->p[1];
  JB->data[26] = udata->p[0];
return(status);

}


