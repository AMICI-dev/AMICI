
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

int JB_model_neuron_o2(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(JB->data,0,sizeof(realtype)*100);
status = w_model_neuron_o2(t,x,NULL,user_data);
status = dwdx_model_neuron_o2(t,x,NULL,user_data);
  JB->data[0+0*10] = x_tmp[0]*(-2.0/2.5E1)-5.0;
  JB->data[0+1*10] = -udata->p[0]*udata->p[1];
  JB->data[1+0*10] = 1.0;
  JB->data[1+1*10] = udata->p[0];
  JB->data[2+0*10] = -x_tmp[2]*udata->dwdx[1];
  JB->data[2+1*10] = -udata->p[1];
  JB->data[2+2*10] = -udata->w[1];
  JB->data[2+3*10] = -udata->p[0]*udata->p[1];
  JB->data[3+1*10] = 1.0;
  JB->data[3+2*10] = 1.0;
  JB->data[3+3*10] = udata->p[0];
  JB->data[4+0*10] = -x_tmp[4]*udata->dwdx[1];
  JB->data[4+1*10] = -udata->p[0];
  JB->data[4+4*10] = -udata->w[1];
  JB->data[4+5*10] = -udata->p[0]*udata->p[1];
  JB->data[5+4*10] = 1.0;
  JB->data[5+5*10] = udata->p[0];
  JB->data[6+0*10] = -x_tmp[6]*udata->dwdx[1];
  JB->data[6+6*10] = -udata->w[1];
  JB->data[6+7*10] = -udata->p[0]*udata->p[1];
  JB->data[7+6*10] = 1.0;
  JB->data[7+7*10] = udata->p[0];
  JB->data[8+0*10] = -x_tmp[8]*udata->dwdx[1];
  JB->data[8+8*10] = -udata->w[1];
  JB->data[8+9*10] = -udata->p[0]*udata->p[1];
  JB->data[9+8*10] = 1.0;
  JB->data[9+9*10] = udata->p[0];
return(status);

}


