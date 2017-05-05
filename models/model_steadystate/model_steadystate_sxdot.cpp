
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
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
    status = JSparse_model_steadystate(t,x,xdot,tmp_J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_steadystate(t,tmp_dxdotdp,x,NULL,user_data);
}
  sxdot_tmp[0] = tmp_dxdotdp[0 + ip*3]+sx_tmp[0]*tmp_J->data[0]+sx_tmp[1]*tmp_J->data[3]+sx_tmp[2]*tmp_J->data[6];
  sxdot_tmp[1] = tmp_dxdotdp[1 + ip*3]+sx_tmp[0]*tmp_J->data[1]+sx_tmp[1]*tmp_J->data[4]+sx_tmp[2]*tmp_J->data[7];
  sxdot_tmp[2] = tmp_dxdotdp[2 + ip*3]+sx_tmp[0]*tmp_J->data[2]+sx_tmp[1]*tmp_J->data[5]+sx_tmp[2]*tmp_J->data[8];
return(status);

}


