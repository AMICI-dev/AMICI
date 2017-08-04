
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_jakstat_adjoint_o2_w.h"

int dydp_model_jakstat_adjoint_o2(realtype t, int it, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 4: {
  tdata->dydp[ip*udata->ny + 0] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  tdata->dydp[ip*udata->ny + 1] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);
  tdata->dydp[ip*udata->ny + 3] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[10]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[11]*2.0;
  tdata->dydp[ip*udata->ny + 4] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[9]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[10]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[11]*2.0;
  tdata->dydp[ip*udata->ny + 6] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[19]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[20]*2.0;
  tdata->dydp[ip*udata->ny + 7] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[18]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[19]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[20]*2.0;
  tdata->dydp[ip*udata->ny + 9] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[28]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[29]*2.0;
  tdata->dydp[ip*udata->ny + 10] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[27]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[28]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[29]*2.0;
  tdata->dydp[ip*udata->ny + 12] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[37]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[38]*2.0;
  tdata->dydp[ip*udata->ny + 13] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[36]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[37]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[38]*2.0;
  tdata->dydp[ip*udata->ny + 15] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[46]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[47]*2.0+1.0/(udata->p[4]*udata->p[4]*udata->p[4])*udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0)*2.0;
  tdata->dydp[ip*udata->ny + 16] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[45]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[46]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[47]*2.0+1.0/(udata->p[4]*udata->p[4]*udata->p[4])*udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)*2.0;
  tdata->dydp[ip*udata->ny + 18] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[55]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[56]*2.0;
  tdata->dydp[ip*udata->ny + 19] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[54]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[55]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[56]*2.0;
  tdata->dydp[ip*udata->ny + 21] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[64]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[65]*2.0;
  tdata->dydp[ip*udata->ny + 22] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[63]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[64]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[65]*2.0;
  tdata->dydp[ip*udata->ny + 24] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[73]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[74]*2.0;
  tdata->dydp[ip*udata->ny + 25] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[72]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[73]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[74]*2.0;
  tdata->dydp[ip*udata->ny + 27] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[82]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[83]*2.0;
  tdata->dydp[ip*udata->ny + 28] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[81]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[82]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[83]*2.0;
  tdata->dydp[ip*udata->ny + 30] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[91]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[92]*2.0;
  tdata->dydp[ip*udata->ny + 31] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[90]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[91]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[92]*2.0;
  tdata->dydp[ip*udata->ny + 33] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[100]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[101]*2.0;
  tdata->dydp[ip*udata->ny + 34] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[99]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[100]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[101]*2.0;
  tdata->dydp[ip*udata->ny + 36] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[109]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[110]*2.0;
  tdata->dydp[ip*udata->ny + 37] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[108]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[109]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[110]*2.0;
  tdata->dydp[ip*udata->ny + 39] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[118]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[119]*2.0;
  tdata->dydp[ip*udata->ny + 40] = -1.0/(udata->p[4]*udata->p[4])*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[117]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[118]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[119]*2.0;
  tdata->dydp[ip*udata->ny + 42] = -1.0/(udata->p[4]*udata->p[4])*(x_tmp[1]+x_tmp[2]*2.0)-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[127]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[128]*2.0;
  tdata->dydp[ip*udata->ny + 43] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[126]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[127]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[128]*2.0;
  tdata->dydp[ip*udata->ny + 45] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[136]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[137]*2.0;
  tdata->dydp[ip*udata->ny + 46] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[135]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[136]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[137]*2.0;
  tdata->dydp[ip*udata->ny + 48] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[145]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[146]*2.0;
  tdata->dydp[ip*udata->ny + 49] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[144]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[145]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[146]*2.0;
  tdata->dydp[ip*udata->ny + 51] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[154]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*x_tmp[155]*2.0;
  tdata->dydp[ip*udata->ny + 52] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[153]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[154]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*x_tmp[155]*2.0;

  } break;

  case 5: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 20] = am_DDspline_pos(4,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 23] = am_DDspline_pos(6,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 26] = am_DDspline_pos(8,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 29] = am_DDspline_pos(10,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 32] = am_DDspline_pos(12,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 6: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 20] = am_DDspline_pos(6,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 23] = am_DDspline_pos(6,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 26] = am_DDspline_pos(8,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 29] = am_DDspline_pos(10,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 32] = am_DDspline_pos(12,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 7: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 20] = am_DDspline_pos(8,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 23] = am_DDspline_pos(8,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 26] = am_DDspline_pos(8,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 29] = am_DDspline_pos(10,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 32] = am_DDspline_pos(12,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 8: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 20] = am_DDspline_pos(10,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 23] = am_DDspline_pos(10,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 26] = am_DDspline_pos(10,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 29] = am_DDspline_pos(10,10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 32] = am_DDspline_pos(12,10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 9: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 20] = am_DDspline_pos(12,4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 23] = am_DDspline_pos(12,6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 26] = am_DDspline_pos(12,8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 29] = am_DDspline_pos(12,10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  tdata->dydp[ip*udata->ny + 32] = am_DDspline_pos(12,12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 10: {
  tdata->dydp[ip*udata->ny + 1] = 1.0;

  } break;

  case 11: {
  tdata->dydp[ip*udata->ny + 0] = 1.0;

  } break;

  case 12: {
  tdata->dydp[ip*udata->ny + 1] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 4] = x_tmp[9]/udata->p[4]+x_tmp[10]/udata->p[4]+(x_tmp[11]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 7] = x_tmp[18]/udata->p[4]+x_tmp[19]/udata->p[4]+(x_tmp[20]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 10] = x_tmp[27]/udata->p[4]+x_tmp[28]/udata->p[4]+(x_tmp[29]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 13] = x_tmp[36]/udata->p[4]+x_tmp[37]/udata->p[4]+(x_tmp[38]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 16] = -1.0/(udata->p[4]*udata->p[4])*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)+x_tmp[45]/udata->p[4]+x_tmp[46]/udata->p[4]+(x_tmp[47]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 19] = x_tmp[54]/udata->p[4]+x_tmp[55]/udata->p[4]+(x_tmp[56]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 22] = x_tmp[63]/udata->p[4]+x_tmp[64]/udata->p[4]+(x_tmp[65]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 25] = x_tmp[72]/udata->p[4]+x_tmp[73]/udata->p[4]+(x_tmp[74]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 28] = x_tmp[81]/udata->p[4]+x_tmp[82]/udata->p[4]+(x_tmp[83]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 31] = x_tmp[90]/udata->p[4]+x_tmp[91]/udata->p[4]+(x_tmp[92]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 34] = x_tmp[99]/udata->p[4]+x_tmp[100]/udata->p[4]+(x_tmp[101]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 37] = x_tmp[108]/udata->p[4]+x_tmp[109]/udata->p[4]+(x_tmp[110]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 40] = x_tmp[117]/udata->p[4]+x_tmp[118]/udata->p[4]+(x_tmp[119]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 43] = x_tmp[126]/udata->p[4]+x_tmp[127]/udata->p[4]+(x_tmp[128]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 46] = x_tmp[135]/udata->p[4]+x_tmp[136]/udata->p[4]+(x_tmp[137]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 49] = x_tmp[144]/udata->p[4]+x_tmp[145]/udata->p[4]+(x_tmp[146]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 52] = x_tmp[153]/udata->p[4]+x_tmp[154]/udata->p[4]+(x_tmp[155]*2.0)/udata->p[4];

  } break;

  case 13: {
  tdata->dydp[ip*udata->ny + 0] = (x_tmp[1]+x_tmp[2]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 3] = x_tmp[10]/udata->p[4]+(x_tmp[11]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 6] = x_tmp[19]/udata->p[4]+(x_tmp[20]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 9] = x_tmp[28]/udata->p[4]+(x_tmp[29]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 12] = x_tmp[37]/udata->p[4]+(x_tmp[38]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 15] = x_tmp[46]/udata->p[4]+(x_tmp[47]*2.0)/udata->p[4]-1.0/(udata->p[4]*udata->p[4])*(x_tmp[1]+x_tmp[2]*2.0);
  tdata->dydp[ip*udata->ny + 18] = x_tmp[55]/udata->p[4]+(x_tmp[56]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 21] = x_tmp[64]/udata->p[4]+(x_tmp[65]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 24] = x_tmp[73]/udata->p[4]+(x_tmp[74]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 27] = x_tmp[82]/udata->p[4]+(x_tmp[83]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 30] = x_tmp[91]/udata->p[4]+(x_tmp[92]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 33] = x_tmp[100]/udata->p[4]+(x_tmp[101]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 36] = x_tmp[109]/udata->p[4]+(x_tmp[110]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 39] = x_tmp[118]/udata->p[4]+(x_tmp[119]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 42] = x_tmp[127]/udata->p[4]+(x_tmp[128]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 45] = x_tmp[136]/udata->p[4]+(x_tmp[137]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 48] = x_tmp[145]/udata->p[4]+(x_tmp[146]*2.0)/udata->p[4];
  tdata->dydp[ip*udata->ny + 51] = x_tmp[154]/udata->p[4]+(x_tmp[155]*2.0)/udata->p[4];

  } break;

}
}
return(status);

}


