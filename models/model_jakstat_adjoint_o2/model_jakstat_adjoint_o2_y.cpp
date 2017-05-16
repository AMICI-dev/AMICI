
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int y_model_jakstat_adjoint_o2(realtype t, int it, realtype *y, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  y[it+nt*(0)] = p[11]+(p[13]*(x_tmp[1]+x_tmp[2]*2.0))/p[4];
  y[it+nt*(1)] = p[10]+(p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0))/p[4];
  y[it+nt*(2)] = am_spline_pos(t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  y[it+nt*(3)] = (p[13]*x_tmp[10])/p[4]+(p[13]*x_tmp[11]*2.0)/p[4];
  y[it+nt*(4)] = (p[12]*x_tmp[9])/p[4]+(p[12]*x_tmp[10])/p[4]+(p[12]*x_tmp[11]*2.0)/p[4];
  y[it+nt*(6)] = (p[13]*x_tmp[19])/p[4]+(p[13]*x_tmp[20]*2.0)/p[4];
  y[it+nt*(7)] = (p[12]*x_tmp[18])/p[4]+(p[12]*x_tmp[19])/p[4]+(p[12]*x_tmp[20]*2.0)/p[4];
  y[it+nt*(9)] = (p[13]*x_tmp[28])/p[4]+(p[13]*x_tmp[29]*2.0)/p[4];
  y[it+nt*(10)] = (p[12]*x_tmp[27])/p[4]+(p[12]*x_tmp[28])/p[4]+(p[12]*x_tmp[29]*2.0)/p[4];
  y[it+nt*(12)] = (p[13]*x_tmp[37])/p[4]+(p[13]*x_tmp[38]*2.0)/p[4];
  y[it+nt*(13)] = (p[12]*x_tmp[36])/p[4]+(p[12]*x_tmp[37])/p[4]+(p[12]*x_tmp[38]*2.0)/p[4];
  y[it+nt*(15)] = (p[13]*x_tmp[46])/p[4]+(p[13]*x_tmp[47]*2.0)/p[4]-1.0/(p[4]*p[4])*p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  y[it+nt*(16)] = (p[12]*x_tmp[45])/p[4]+(p[12]*x_tmp[46])/p[4]+(p[12]*x_tmp[47]*2.0)/p[4]-1.0/(p[4]*p[4])*p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);
  y[it+nt*(18)] = (p[13]*x_tmp[55])/p[4]+(p[13]*x_tmp[56]*2.0)/p[4];
  y[it+nt*(19)] = (p[12]*x_tmp[54])/p[4]+(p[12]*x_tmp[55])/p[4]+(p[12]*x_tmp[56]*2.0)/p[4];
  y[it+nt*(20)] = am_Dspline_pos(4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  y[it+nt*(21)] = (p[13]*x_tmp[64])/p[4]+(p[13]*x_tmp[65]*2.0)/p[4];
  y[it+nt*(22)] = (p[12]*x_tmp[63])/p[4]+(p[12]*x_tmp[64])/p[4]+(p[12]*x_tmp[65]*2.0)/p[4];
  y[it+nt*(23)] = am_Dspline_pos(6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  y[it+nt*(24)] = (p[13]*x_tmp[73])/p[4]+(p[13]*x_tmp[74]*2.0)/p[4];
  y[it+nt*(25)] = (p[12]*x_tmp[72])/p[4]+(p[12]*x_tmp[73])/p[4]+(p[12]*x_tmp[74]*2.0)/p[4];
  y[it+nt*(26)] = am_Dspline_pos(8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  y[it+nt*(27)] = (p[13]*x_tmp[82])/p[4]+(p[13]*x_tmp[83]*2.0)/p[4];
  y[it+nt*(28)] = (p[12]*x_tmp[81])/p[4]+(p[12]*x_tmp[82])/p[4]+(p[12]*x_tmp[83]*2.0)/p[4];
  y[it+nt*(29)] = am_Dspline_pos(10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  y[it+nt*(30)] = (p[13]*x_tmp[91])/p[4]+(p[13]*x_tmp[92]*2.0)/p[4];
  y[it+nt*(31)] = (p[12]*x_tmp[90])/p[4]+(p[12]*x_tmp[91])/p[4]+(p[12]*x_tmp[92]*2.0)/p[4];
  y[it+nt*(32)] = am_Dspline_pos(12,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  y[it+nt*(33)] = (p[13]*x_tmp[100])/p[4]+(p[13]*x_tmp[101]*2.0)/p[4];
  y[it+nt*(34)] = (p[12]*x_tmp[99])/p[4]+(p[12]*x_tmp[100])/p[4]+(p[12]*x_tmp[101]*2.0)/p[4]+1.0;
  y[it+nt*(36)] = (p[13]*x_tmp[109])/p[4]+(p[13]*x_tmp[110]*2.0)/p[4]+1.0;
  y[it+nt*(37)] = (p[12]*x_tmp[108])/p[4]+(p[12]*x_tmp[109])/p[4]+(p[12]*x_tmp[110]*2.0)/p[4];
  y[it+nt*(39)] = (p[13]*x_tmp[118])/p[4]+(p[13]*x_tmp[119]*2.0)/p[4];
  y[it+nt*(40)] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/p[4]+(p[12]*x_tmp[117])/p[4]+(p[12]*x_tmp[118])/p[4]+(p[12]*x_tmp[119]*2.0)/p[4];
  y[it+nt*(42)] = (x_tmp[1]+x_tmp[2]*2.0)/p[4]+(p[13]*x_tmp[127])/p[4]+(p[13]*x_tmp[128]*2.0)/p[4];
  y[it+nt*(43)] = (p[12]*x_tmp[126])/p[4]+(p[12]*x_tmp[127])/p[4]+(p[12]*x_tmp[128]*2.0)/p[4];
  y[it+nt*(45)] = (p[13]*x_tmp[136])/p[4]+(p[13]*x_tmp[137]*2.0)/p[4];
  y[it+nt*(46)] = (p[12]*x_tmp[135])/p[4]+(p[12]*x_tmp[136])/p[4]+(p[12]*x_tmp[137]*2.0)/p[4];
  y[it+nt*(48)] = (p[13]*x_tmp[145])/p[4]+(p[13]*x_tmp[146]*2.0)/p[4];
  y[it+nt*(49)] = (p[12]*x_tmp[144])/p[4]+(p[12]*x_tmp[145])/p[4]+(p[12]*x_tmp[146]*2.0)/p[4];
  y[it+nt*(51)] = (p[13]*x_tmp[154])/p[4]+(p[13]*x_tmp[155]*2.0)/p[4];
  y[it+nt*(52)] = (p[12]*x_tmp[153])/p[4]+(p[12]*x_tmp[154])/p[4]+(p[12]*x_tmp[155]*2.0)/p[4];
return(status);

}


