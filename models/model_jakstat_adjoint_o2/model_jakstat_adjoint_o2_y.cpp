
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int y_model_jakstat_adjoint_o2(realtype t, int it, realtype *y, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  y[it + udata->nt*(0)] = udata->p[11]+(udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0))/udata->p[4];
  y[it + udata->nt*(1)] = udata->p[10]+(udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0))/udata->p[4];
  y[it + udata->nt*(2)] = am_spline_pos(t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  y[it + udata->nt*(3)] = (udata->p[13]*x_tmp[10])/udata->p[4]+(udata->p[13]*x_tmp[11]*2.0)/udata->p[4];
  y[it + udata->nt*(4)] = (udata->p[12]*x_tmp[9])/udata->p[4]+(udata->p[12]*x_tmp[10])/udata->p[4]+(udata->p[12]*x_tmp[11]*2.0)/udata->p[4];
  y[it + udata->nt*(6)] = (udata->p[13]*x_tmp[19])/udata->p[4]+(udata->p[13]*x_tmp[20]*2.0)/udata->p[4];
  y[it + udata->nt*(7)] = (udata->p[12]*x_tmp[18])/udata->p[4]+(udata->p[12]*x_tmp[19])/udata->p[4]+(udata->p[12]*x_tmp[20]*2.0)/udata->p[4];
  y[it + udata->nt*(9)] = (udata->p[13]*x_tmp[28])/udata->p[4]+(udata->p[13]*x_tmp[29]*2.0)/udata->p[4];
  y[it + udata->nt*(10)] = (udata->p[12]*x_tmp[27])/udata->p[4]+(udata->p[12]*x_tmp[28])/udata->p[4]+(udata->p[12]*x_tmp[29]*2.0)/udata->p[4];
  y[it + udata->nt*(12)] = (udata->p[13]*x_tmp[37])/udata->p[4]+(udata->p[13]*x_tmp[38]*2.0)/udata->p[4];
  y[it + udata->nt*(13)] = (udata->p[12]*x_tmp[36])/udata->p[4]+(udata->p[12]*x_tmp[37])/udata->p[4]+(udata->p[12]*x_tmp[38]*2.0)/udata->p[4];
  y[it + udata->nt*(15)] = (udata->p[13]*x_tmp[46])/udata->p[4]+(udata->p[13]*x_tmp[47]*2.0)/udata->p[4]-1.0/(udata->p[4]*udata->p[4])*udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  y[it + udata->nt*(16)] = (udata->p[12]*x_tmp[45])/udata->p[4]+(udata->p[12]*x_tmp[46])/udata->p[4]+(udata->p[12]*x_tmp[47]*2.0)/udata->p[4]-1.0/(udata->p[4]*udata->p[4])*udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);
  y[it + udata->nt*(18)] = (udata->p[13]*x_tmp[55])/udata->p[4]+(udata->p[13]*x_tmp[56]*2.0)/udata->p[4];
  y[it + udata->nt*(19)] = (udata->p[12]*x_tmp[54])/udata->p[4]+(udata->p[12]*x_tmp[55])/udata->p[4]+(udata->p[12]*x_tmp[56]*2.0)/udata->p[4];
  y[it + udata->nt*(20)] = am_Dspline_pos(4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  y[it + udata->nt*(21)] = (udata->p[13]*x_tmp[64])/udata->p[4]+(udata->p[13]*x_tmp[65]*2.0)/udata->p[4];
  y[it + udata->nt*(22)] = (udata->p[12]*x_tmp[63])/udata->p[4]+(udata->p[12]*x_tmp[64])/udata->p[4]+(udata->p[12]*x_tmp[65]*2.0)/udata->p[4];
  y[it + udata->nt*(23)] = am_Dspline_pos(6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  y[it + udata->nt*(24)] = (udata->p[13]*x_tmp[73])/udata->p[4]+(udata->p[13]*x_tmp[74]*2.0)/udata->p[4];
  y[it + udata->nt*(25)] = (udata->p[12]*x_tmp[72])/udata->p[4]+(udata->p[12]*x_tmp[73])/udata->p[4]+(udata->p[12]*x_tmp[74]*2.0)/udata->p[4];
  y[it + udata->nt*(26)] = am_Dspline_pos(8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  y[it + udata->nt*(27)] = (udata->p[13]*x_tmp[82])/udata->p[4]+(udata->p[13]*x_tmp[83]*2.0)/udata->p[4];
  y[it + udata->nt*(28)] = (udata->p[12]*x_tmp[81])/udata->p[4]+(udata->p[12]*x_tmp[82])/udata->p[4]+(udata->p[12]*x_tmp[83]*2.0)/udata->p[4];
  y[it + udata->nt*(29)] = am_Dspline_pos(10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  y[it + udata->nt*(30)] = (udata->p[13]*x_tmp[91])/udata->p[4]+(udata->p[13]*x_tmp[92]*2.0)/udata->p[4];
  y[it + udata->nt*(31)] = (udata->p[12]*x_tmp[90])/udata->p[4]+(udata->p[12]*x_tmp[91])/udata->p[4]+(udata->p[12]*x_tmp[92]*2.0)/udata->p[4];
  y[it + udata->nt*(32)] = am_Dspline_pos(12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  y[it + udata->nt*(33)] = (udata->p[13]*x_tmp[100])/udata->p[4]+(udata->p[13]*x_tmp[101]*2.0)/udata->p[4];
  y[it + udata->nt*(34)] = (udata->p[12]*x_tmp[99])/udata->p[4]+(udata->p[12]*x_tmp[100])/udata->p[4]+(udata->p[12]*x_tmp[101]*2.0)/udata->p[4]+1.0;
  y[it + udata->nt*(36)] = (udata->p[13]*x_tmp[109])/udata->p[4]+(udata->p[13]*x_tmp[110]*2.0)/udata->p[4]+1.0;
  y[it + udata->nt*(37)] = (udata->p[12]*x_tmp[108])/udata->p[4]+(udata->p[12]*x_tmp[109])/udata->p[4]+(udata->p[12]*x_tmp[110]*2.0)/udata->p[4];
  y[it + udata->nt*(39)] = (udata->p[13]*x_tmp[118])/udata->p[4]+(udata->p[13]*x_tmp[119]*2.0)/udata->p[4];
  y[it + udata->nt*(40)] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/udata->p[4]+(udata->p[12]*x_tmp[117])/udata->p[4]+(udata->p[12]*x_tmp[118])/udata->p[4]+(udata->p[12]*x_tmp[119]*2.0)/udata->p[4];
  y[it + udata->nt*(42)] = (x_tmp[1]+x_tmp[2]*2.0)/udata->p[4]+(udata->p[13]*x_tmp[127])/udata->p[4]+(udata->p[13]*x_tmp[128]*2.0)/udata->p[4];
  y[it + udata->nt*(43)] = (udata->p[12]*x_tmp[126])/udata->p[4]+(udata->p[12]*x_tmp[127])/udata->p[4]+(udata->p[12]*x_tmp[128]*2.0)/udata->p[4];
  y[it + udata->nt*(45)] = (udata->p[13]*x_tmp[136])/udata->p[4]+(udata->p[13]*x_tmp[137]*2.0)/udata->p[4];
  y[it + udata->nt*(46)] = (udata->p[12]*x_tmp[135])/udata->p[4]+(udata->p[12]*x_tmp[136])/udata->p[4]+(udata->p[12]*x_tmp[137]*2.0)/udata->p[4];
  y[it + udata->nt*(48)] = (udata->p[13]*x_tmp[145])/udata->p[4]+(udata->p[13]*x_tmp[146]*2.0)/udata->p[4];
  y[it + udata->nt*(49)] = (udata->p[12]*x_tmp[144])/udata->p[4]+(udata->p[12]*x_tmp[145])/udata->p[4]+(udata->p[12]*x_tmp[146]*2.0)/udata->p[4];
  y[it + udata->nt*(51)] = (udata->p[13]*x_tmp[154])/udata->p[4]+(udata->p[13]*x_tmp[155]*2.0)/udata->p[4];
  y[it + udata->nt*(52)] = (udata->p[12]*x_tmp[153])/udata->p[4]+(udata->p[12]*x_tmp[154])/udata->p[4]+(udata->p[12]*x_tmp[155]*2.0)/udata->p[4];
return(status);

}


