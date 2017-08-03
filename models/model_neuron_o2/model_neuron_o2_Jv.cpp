
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

int Jv_model_neuron_o2(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *v_tmp = N_VGetArrayPointer(v);
realtype *Jv_tmp = N_VGetArrayPointer(Jv);
memset(Jv_tmp,0,sizeof(realtype)*10);
status = w_model_neuron_o2(t,x,NULL,user_data);
  Jv_tmp[0] = -v_tmp[1]+v_tmp[0]*(x_tmp[0]*(2.0/2.5E1)+5.0);
  Jv_tmp[1] = -udata->p[0]*v_tmp[1]+udata->p[0]*udata->p[1]*v_tmp[0];
  Jv_tmp[2] = -v_tmp[3]+v_tmp[2]*udata->w[1]+v_tmp[0]*x_tmp[2]*udata->dwdx[1];
  Jv_tmp[3] = -v_tmp[1]+udata->p[1]*v_tmp[0]-udata->p[0]*v_tmp[3]+udata->p[0]*udata->p[1]*v_tmp[2];
  Jv_tmp[4] = -v_tmp[5]+udata->w[1]*v_tmp[4]+v_tmp[0]*x_tmp[4]*udata->dwdx[1];
  Jv_tmp[5] = udata->p[0]*v_tmp[0]-udata->p[0]*v_tmp[5]+udata->p[0]*udata->p[1]*v_tmp[4];
  Jv_tmp[6] = -v_tmp[7]+udata->w[1]*v_tmp[6]+v_tmp[0]*x_tmp[6]*udata->dwdx[1];
  Jv_tmp[7] = -udata->p[0]*v_tmp[7]+udata->p[0]*udata->p[1]*v_tmp[6];
  Jv_tmp[8] = -v_tmp[9]+udata->w[1]*v_tmp[8]+v_tmp[0]*x_tmp[8]*udata->dwdx[1];
  Jv_tmp[9] = -udata->p[0]*v_tmp[9]+udata->p[0]*udata->p[1]*v_tmp[8];
return(status);

}


