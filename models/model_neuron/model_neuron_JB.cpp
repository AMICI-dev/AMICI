
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_dwdx.h"
#include "model_neuron_w.h"

int JB_model_neuron(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(JB->data,0,sizeof(realtype)*4);
status = w_model_neuron(t,x,NULL,tdata);
status = dwdx_model_neuron(t,x,NULL,user_data);
  JB->data[0+0*2] = x_tmp[0]*(-2.0/2.5E1)-5.0;
  JB->data[0+1*2] = -udata->p[0]*udata->p[1];
  JB->data[1+0*2] = 1.0;
  JB->data[1+1*2] = udata->p[0];
return(status);

}


