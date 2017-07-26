
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_dwdx.h"
#include "model_steadystate_w.h"

int JB_model_steadystate(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(JB->data,0,sizeof(realtype)*9);
status = w_model_steadystate(t,x,NULL,user_data);
status = dwdx_model_steadystate(t,x,NULL,user_data);
  JB->data[0+0*3] = udata->p[1]*x_tmp[1]+udata->p[0]*udata->dwdx[0]*2.0;
  JB->data[0+1*3] = udata->p[1]*x_tmp[1]-udata->p[0]*udata->dwdx[0];
  JB->data[0+2*3] = -udata->p[1]*x_tmp[1];
  JB->data[1+0*3] = udata->p[2]*-2.0+udata->p[1]*x_tmp[0];
  JB->data[1+1*3] = udata->p[2]+udata->p[1]*x_tmp[0];
  JB->data[1+2*3] = -udata->p[1]*x_tmp[0];
  JB->data[2+0*3] = -udata->dwdx[1];
  JB->data[2+1*3] = -udata->dwdx[1];
  JB->data[2+2*3] = udata->k[3]+udata->dwdx[1];
return(status);

}


