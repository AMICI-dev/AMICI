
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int dwdp_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(udata->dwdp,0,sizeof(realtype)*30);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  udata->dwdp[0] = am_Dspline_pos(4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[1] = am_DDspline_pos(4,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[2] = am_DDspline_pos(6,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[3] = am_DDspline_pos(8,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[4] = am_DDspline_pos(10,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[5] = am_DDspline_pos(12,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[6] = am_Dspline_pos(6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[7] = am_DDspline_pos(6,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[8] = am_DDspline_pos(6,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[9] = am_DDspline_pos(8,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[10] = am_DDspline_pos(10,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[11] = am_DDspline_pos(12,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[12] = am_Dspline_pos(8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[13] = am_DDspline_pos(8,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[14] = am_DDspline_pos(8,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[15] = am_DDspline_pos(8,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[16] = am_DDspline_pos(10,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[17] = am_DDspline_pos(12,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[18] = am_Dspline_pos(10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[19] = am_DDspline_pos(10,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[20] = am_DDspline_pos(10,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[21] = am_DDspline_pos(10,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[22] = am_DDspline_pos(10,10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[23] = am_DDspline_pos(12,10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[24] = am_Dspline_pos(12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[25] = am_DDspline_pos(12,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[26] = am_DDspline_pos(12,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[27] = am_DDspline_pos(12,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[28] = am_DDspline_pos(12,10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->dwdp[29] = am_DDspline_pos(12,12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
return(status);

}


