
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_JSparse.h"
#include "model_neuron_o2_dxdotdp.h"
#include "model_neuron_o2_w.h"

int sxdot_model_neuron_o2(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp = N_VGetArrayPointer(sx);
realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*10);
if(ip == 0) {
    status = JSparse_model_neuron_o2(t,x,xdot,udata->J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_neuron_o2(t,x,NULL,user_data);
}
  sxdot_tmp[0] = udata->dxdotdp[0 + ip*udata->nx]+udata->J->data[0]*sx_tmp[0]+udata->J->data[8]*sx_tmp[1];
  sxdot_tmp[1] = udata->dxdotdp[1 + ip*udata->nx]+udata->J->data[1]*sx_tmp[0]+udata->J->data[9]*sx_tmp[1];
  sxdot_tmp[2] = udata->dxdotdp[2 + ip*udata->nx]+udata->J->data[2]*sx_tmp[0]+udata->J->data[11]*sx_tmp[2]+udata->J->data[13]*sx_tmp[3];
  sxdot_tmp[3] = udata->dxdotdp[3 + ip*udata->nx]+udata->J->data[3]*sx_tmp[0]+udata->J->data[10]*sx_tmp[1]+udata->J->data[12]*sx_tmp[2]+udata->J->data[14]*sx_tmp[3];
  sxdot_tmp[4] = udata->dxdotdp[4 + ip*udata->nx]+udata->J->data[4]*sx_tmp[0]+udata->J->data[15]*sx_tmp[4]+udata->J->data[17]*sx_tmp[5];
  sxdot_tmp[5] = udata->dxdotdp[5 + ip*udata->nx]+udata->J->data[5]*sx_tmp[0]+udata->J->data[16]*sx_tmp[4]+udata->J->data[18]*sx_tmp[5];
  sxdot_tmp[6] = udata->dxdotdp[6 + ip*udata->nx]+udata->J->data[6]*sx_tmp[0]+udata->J->data[19]*sx_tmp[6]+udata->J->data[21]*sx_tmp[7];
  sxdot_tmp[7] = udata->dxdotdp[7 + ip*udata->nx]+udata->J->data[20]*sx_tmp[6]+udata->J->data[22]*sx_tmp[7];
  sxdot_tmp[8] = udata->dxdotdp[8 + ip*udata->nx]+udata->J->data[7]*sx_tmp[0]+udata->J->data[23]*sx_tmp[8]+udata->J->data[25]*sx_tmp[9];
  sxdot_tmp[9] = udata->dxdotdp[9 + ip*udata->nx]+udata->J->data[24]*sx_tmp[8]+udata->J->data[26]*sx_tmp[9];
return(status);

}


