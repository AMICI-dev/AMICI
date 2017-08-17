
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

int JvB_model_neuron_o2(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
realtype *vB_tmp = N_VGetArrayPointer(vB);
realtype *JvB_tmp = N_VGetArrayPointer(JvB);
memset(JvB_tmp,0,sizeof(realtype)*10);
status = w_model_neuron_o2(t,x,NULL,tdata);
  JvB_tmp[0] = -vB_tmp[0]*(x_tmp[0]*(2.0/2.5E1)+5.0)-udata->p[0]*udata->p[1]*vB_tmp[1];
  JvB_tmp[1] = vB_tmp[0]+udata->p[0]*vB_tmp[1];
  JvB_tmp[2] = -udata->p[1]*vB_tmp[1]-udata->w[1]*vB_tmp[2]-udata->p[0]*udata->p[1]*vB_tmp[3]-x_tmp[2]*vB_tmp[0]*udata->dwdx[1];
  JvB_tmp[3] = vB_tmp[1]+vB_tmp[2]+udata->p[0]*vB_tmp[3];
  JvB_tmp[4] = -udata->p[0]*vB_tmp[1]-udata->w[1]*vB_tmp[4]-udata->p[0]*udata->p[1]*vB_tmp[5]-x_tmp[4]*vB_tmp[0]*udata->dwdx[1];
  JvB_tmp[5] = vB_tmp[4]+udata->p[0]*vB_tmp[5];
  JvB_tmp[6] = -udata->w[1]*vB_tmp[6]-udata->p[0]*udata->p[1]*vB_tmp[7]-x_tmp[6]*vB_tmp[0]*udata->dwdx[1];
  JvB_tmp[7] = vB_tmp[6]+udata->p[0]*vB_tmp[7];
  JvB_tmp[8] = -udata->w[1]*vB_tmp[8]-udata->p[0]*udata->p[1]*vB_tmp[9]-x_tmp[8]*vB_tmp[0]*udata->dwdx[1];
  JvB_tmp[9] = vB_tmp[8]+udata->p[0]*vB_tmp[9];
return(status);

}


