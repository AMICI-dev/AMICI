
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int JSparseB_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
UserData *udata = (UserData*) user_data;
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
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  JB->data[0] = p[0]*w_tmp[0];
  JB->data[1] = -(k[1]*p[3])/k[0];
  JB->data[2] = -p[0]*w_tmp[0];
  JB->data[3] = dwdx_tmp[0]*p[1]*2.0;
  JB->data[4] = -dwdx_tmp[0]*p[1];
  JB->data[5] = p[2];
  JB->data[6] = -(k[0]*p[2])/k[1];
  JB->data[7] = p[3];
  JB->data[8] = p[3]*-2.0;
  JB->data[9] = p[3];
  JB->data[10] = -p[3];
  JB->data[11] = p[3];
  JB->data[12] = -p[3];
  JB->data[13] = p[3];
  JB->data[14] = -p[3];
  JB->data[15] = p[3];
  JB->data[16] = -p[3];
  JB->data[17] = p[3];
return(status);

}


