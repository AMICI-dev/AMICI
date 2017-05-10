
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_JSparse.h"
#include "model_jakstat_adjoint_dxdotdp.h"
#include "model_jakstat_adjoint_w.h"

int sxdot_model_jakstat_adjoint(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp = N_VGetArrayPointer(sx);
realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*9);
if(ip == 0) {
    status = JSparse_model_jakstat_adjoint(t,x,xdot,tmp_J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_jakstat_adjoint(t,tmp_dxdotdp,x,NULL,user_data);
}
  sxdot_tmp[0] = tmp_dxdotdp[0 + ip*9]+sx_tmp[0]*tmp_J->data[0]+sx_tmp[8]*tmp_J->data[16];
  sxdot_tmp[1] = tmp_dxdotdp[1 + ip*9]+sx_tmp[0]*tmp_J->data[1]+sx_tmp[1]*tmp_J->data[2];
  sxdot_tmp[2] = tmp_dxdotdp[2 + ip*9]+sx_tmp[1]*tmp_J->data[3]+sx_tmp[2]*tmp_J->data[4];
  sxdot_tmp[3] = tmp_dxdotdp[3 + ip*9]+sx_tmp[2]*tmp_J->data[5]+sx_tmp[3]*tmp_J->data[6];
  sxdot_tmp[4] = tmp_dxdotdp[4 + ip*9]+sx_tmp[3]*tmp_J->data[7]+sx_tmp[4]*tmp_J->data[8];
  sxdot_tmp[5] = tmp_dxdotdp[5 + ip*9]+sx_tmp[4]*tmp_J->data[9]+sx_tmp[5]*tmp_J->data[10];
  sxdot_tmp[6] = tmp_dxdotdp[6 + ip*9]+sx_tmp[5]*tmp_J->data[11]+sx_tmp[6]*tmp_J->data[12];
  sxdot_tmp[7] = tmp_dxdotdp[7 + ip*9]+sx_tmp[6]*tmp_J->data[13]+sx_tmp[7]*tmp_J->data[14];
  sxdot_tmp[8] = tmp_dxdotdp[8 + ip*9]+sx_tmp[7]*tmp_J->data[15]+sx_tmp[8]*tmp_J->data[17];
return(status);

}


