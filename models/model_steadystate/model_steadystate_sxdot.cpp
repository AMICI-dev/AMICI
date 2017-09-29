
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_steadystate_JSparse.h"
#include "model_steadystate_dxdotdp.h"
#include "model_steadystate_w.h"

int sxdot_model_steadystate(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp = N_VGetArrayPointer(sx);
realtype *sdx_tmp = N_VGetArrayPointer(sdx);
realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*3);
if(ip == 0) {
    status = JSparse_model_steadystate(t,x,xdot,tdata->J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_steadystate(t,x,NULL,user_data);
}
  sxdot_tmp[0] = tdata->dxdotdp[0 + ip*model->nx]+tdata->J->data[0]*sx_tmp[0]+tdata->J->data[3]*sx_tmp[1]+tdata->J->data[6]*sx_tmp[2];
  sxdot_tmp[1] = tdata->dxdotdp[1 + ip*model->nx]+tdata->J->data[1]*sx_tmp[0]+tdata->J->data[4]*sx_tmp[1]+tdata->J->data[7]*sx_tmp[2];
  sxdot_tmp[2] = tdata->dxdotdp[2 + ip*model->nx]+tdata->J->data[2]*sx_tmp[0]+tdata->J->data[5]*sx_tmp[1]+tdata->J->data[8]*sx_tmp[2];
return(status);

}


