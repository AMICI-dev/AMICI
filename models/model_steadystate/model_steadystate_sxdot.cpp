
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_JSparse.h"
#include "model_steadystate_dxdotdp.h"
#include "model_steadystate_w.h"

int sxdot_model_steadystate(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp = N_VGetArrayPointer(sx);
realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*3);
if(ip == 0) {
    status = JSparse_model_steadystate(t,x,xdot,udata->J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_steadystate(t,x,NULL,user_data);
}
  sxdot_tmp[0] = udata->dxdotdp[0 + ip*udata->nx]+udata->J->data[0]*sx_tmp[0]+udata->J->data[3]*sx_tmp[1]+udata->J->data[6]*sx_tmp[2];
  sxdot_tmp[1] = udata->dxdotdp[1 + ip*udata->nx]+udata->J->data[1]*sx_tmp[0]+udata->J->data[4]*sx_tmp[1]+udata->J->data[7]*sx_tmp[2];
  sxdot_tmp[2] = udata->dxdotdp[2 + ip*udata->nx]+udata->J->data[2]*sx_tmp[0]+udata->J->data[5]*sx_tmp[1]+udata->J->data[8]*sx_tmp[2];
return(status);

}


