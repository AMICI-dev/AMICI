
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int dwdp_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dwdp_tmp,0,sizeof(realtype)*30);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  dwdp_tmp[0] = am_Dspline_pos(4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[1] = am_DDspline_pos(4,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[2] = am_DDspline_pos(6,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[3] = am_DDspline_pos(8,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[4] = am_DDspline_pos(10,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[5] = am_DDspline_pos(12,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[6] = am_Dspline_pos(6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[7] = am_DDspline_pos(6,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[8] = am_DDspline_pos(6,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[9] = am_DDspline_pos(8,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[10] = am_DDspline_pos(10,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[11] = am_DDspline_pos(12,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[12] = am_Dspline_pos(8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[13] = am_DDspline_pos(8,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[14] = am_DDspline_pos(8,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[15] = am_DDspline_pos(8,8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[16] = am_DDspline_pos(10,8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[17] = am_DDspline_pos(12,8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[18] = am_Dspline_pos(10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[19] = am_DDspline_pos(10,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[20] = am_DDspline_pos(10,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[21] = am_DDspline_pos(10,8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[22] = am_DDspline_pos(10,10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[23] = am_DDspline_pos(12,10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[24] = am_Dspline_pos(12,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[25] = am_DDspline_pos(12,4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[26] = am_DDspline_pos(12,6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[27] = am_DDspline_pos(12,8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[28] = am_DDspline_pos(12,10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[29] = am_DDspline_pos(12,12,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
return(status);

}


