
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_JSparse.h"
#include "model_neuron_o2_dxdotdp.h"
#include "model_neuron_o2_w.h"

using namespace amici;

int sxdot_model_neuron_o2(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *sx_tmp = nullptr;
if(sx)
    sx_tmp = N_VGetArrayPointer(sx);
realtype *sdx_tmp = nullptr;
if(sdx)
    sdx_tmp = N_VGetArrayPointer(sdx);
realtype *sxdot_tmp = nullptr;
if(sxdot)
    sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*10);
if(ip == 0) {
    status = JSparse_model_neuron_o2(t,0.0,x,NULL,xdot,tdata->J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_neuron_o2(t,x,NULL,user_data);
}
  sxdot_tmp[0] = tdata->dxdotdp[0 + ip*model->nx]+tdata->J->data[0]*sx_tmp[0]+tdata->J->data[8]*sx_tmp[1];
  sxdot_tmp[1] = tdata->dxdotdp[1 + ip*model->nx]+tdata->J->data[1]*sx_tmp[0]+tdata->J->data[9]*sx_tmp[1];
  sxdot_tmp[2] = tdata->dxdotdp[2 + ip*model->nx]+tdata->J->data[2]*sx_tmp[0]+tdata->J->data[11]*sx_tmp[2]+tdata->J->data[13]*sx_tmp[3];
  sxdot_tmp[3] = tdata->dxdotdp[3 + ip*model->nx]+tdata->J->data[3]*sx_tmp[0]+tdata->J->data[10]*sx_tmp[1]+tdata->J->data[12]*sx_tmp[2]+tdata->J->data[14]*sx_tmp[3];
  sxdot_tmp[4] = tdata->dxdotdp[4 + ip*model->nx]+tdata->J->data[4]*sx_tmp[0]+tdata->J->data[15]*sx_tmp[4]+tdata->J->data[17]*sx_tmp[5];
  sxdot_tmp[5] = tdata->dxdotdp[5 + ip*model->nx]+tdata->J->data[5]*sx_tmp[0]+tdata->J->data[16]*sx_tmp[4]+tdata->J->data[18]*sx_tmp[5];
  sxdot_tmp[6] = tdata->dxdotdp[6 + ip*model->nx]+tdata->J->data[6]*sx_tmp[0]+tdata->J->data[19]*sx_tmp[6]+tdata->J->data[21]*sx_tmp[7];
  sxdot_tmp[7] = tdata->dxdotdp[7 + ip*model->nx]+tdata->J->data[20]*sx_tmp[6]+tdata->J->data[22]*sx_tmp[7];
  sxdot_tmp[8] = tdata->dxdotdp[8 + ip*model->nx]+tdata->J->data[7]*sx_tmp[0]+tdata->J->data[23]*sx_tmp[8]+tdata->J->data[25]*sx_tmp[9];
  sxdot_tmp[9] = tdata->dxdotdp[9 + ip*model->nx]+tdata->J->data[24]*sx_tmp[8]+tdata->J->data[26]*sx_tmp[9];
return(status);

}


