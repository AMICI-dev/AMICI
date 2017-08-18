
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_jakstat_adjoint_o2_w.h"

int y_model_jakstat_adjoint_o2(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
  rdata->y[it + udata->nt*0] = tdata->p[11]+(tdata->p[13]*(x_tmp[1]+x_tmp[2]*2.0))/tdata->p[4];
  rdata->y[it + udata->nt*1] = tdata->p[10]+(tdata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0))/tdata->p[4];
  rdata->y[it + udata->nt*2] = am_spline_pos(t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  rdata->y[it + udata->nt*3] = (tdata->p[13]*x_tmp[10])/tdata->p[4]+(tdata->p[13]*x_tmp[11]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*4] = (tdata->p[12]*x_tmp[9])/tdata->p[4]+(tdata->p[12]*x_tmp[10])/tdata->p[4]+(tdata->p[12]*x_tmp[11]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*6] = (tdata->p[13]*x_tmp[19])/tdata->p[4]+(tdata->p[13]*x_tmp[20]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*7] = (tdata->p[12]*x_tmp[18])/tdata->p[4]+(tdata->p[12]*x_tmp[19])/tdata->p[4]+(tdata->p[12]*x_tmp[20]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*9] = (tdata->p[13]*x_tmp[28])/tdata->p[4]+(tdata->p[13]*x_tmp[29]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*10] = (tdata->p[12]*x_tmp[27])/tdata->p[4]+(tdata->p[12]*x_tmp[28])/tdata->p[4]+(tdata->p[12]*x_tmp[29]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*12] = (tdata->p[13]*x_tmp[37])/tdata->p[4]+(tdata->p[13]*x_tmp[38]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*13] = (tdata->p[12]*x_tmp[36])/tdata->p[4]+(tdata->p[12]*x_tmp[37])/tdata->p[4]+(tdata->p[12]*x_tmp[38]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*15] = (tdata->p[13]*x_tmp[46])/tdata->p[4]+(tdata->p[13]*x_tmp[47]*2.0)/tdata->p[4]-1.0/(tdata->p[4]*tdata->p[4])*tdata->p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  rdata->y[it + udata->nt*16] = (tdata->p[12]*x_tmp[45])/tdata->p[4]+(tdata->p[12]*x_tmp[46])/tdata->p[4]+(tdata->p[12]*x_tmp[47]*2.0)/tdata->p[4]-1.0/(tdata->p[4]*tdata->p[4])*tdata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);
  rdata->y[it + udata->nt*18] = (tdata->p[13]*x_tmp[55])/tdata->p[4]+(tdata->p[13]*x_tmp[56]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*19] = (tdata->p[12]*x_tmp[54])/tdata->p[4]+(tdata->p[12]*x_tmp[55])/tdata->p[4]+(tdata->p[12]*x_tmp[56]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*20] = am_Dspline_pos(4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  rdata->y[it + udata->nt*21] = (tdata->p[13]*x_tmp[64])/tdata->p[4]+(tdata->p[13]*x_tmp[65]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*22] = (tdata->p[12]*x_tmp[63])/tdata->p[4]+(tdata->p[12]*x_tmp[64])/tdata->p[4]+(tdata->p[12]*x_tmp[65]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*23] = am_Dspline_pos(6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  rdata->y[it + udata->nt*24] = (tdata->p[13]*x_tmp[73])/tdata->p[4]+(tdata->p[13]*x_tmp[74]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*25] = (tdata->p[12]*x_tmp[72])/tdata->p[4]+(tdata->p[12]*x_tmp[73])/tdata->p[4]+(tdata->p[12]*x_tmp[74]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*26] = am_Dspline_pos(8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  rdata->y[it + udata->nt*27] = (tdata->p[13]*x_tmp[82])/tdata->p[4]+(tdata->p[13]*x_tmp[83]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*28] = (tdata->p[12]*x_tmp[81])/tdata->p[4]+(tdata->p[12]*x_tmp[82])/tdata->p[4]+(tdata->p[12]*x_tmp[83]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*29] = am_Dspline_pos(10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  rdata->y[it + udata->nt*30] = (tdata->p[13]*x_tmp[91])/tdata->p[4]+(tdata->p[13]*x_tmp[92]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*31] = (tdata->p[12]*x_tmp[90])/tdata->p[4]+(tdata->p[12]*x_tmp[91])/tdata->p[4]+(tdata->p[12]*x_tmp[92]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*32] = am_Dspline_pos(12,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  rdata->y[it + udata->nt*33] = (tdata->p[13]*x_tmp[100])/tdata->p[4]+(tdata->p[13]*x_tmp[101]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*34] = (tdata->p[12]*x_tmp[99])/tdata->p[4]+(tdata->p[12]*x_tmp[100])/tdata->p[4]+(tdata->p[12]*x_tmp[101]*2.0)/tdata->p[4]+1.0;
  rdata->y[it + udata->nt*36] = (tdata->p[13]*x_tmp[109])/tdata->p[4]+(tdata->p[13]*x_tmp[110]*2.0)/tdata->p[4]+1.0;
  rdata->y[it + udata->nt*37] = (tdata->p[12]*x_tmp[108])/tdata->p[4]+(tdata->p[12]*x_tmp[109])/tdata->p[4]+(tdata->p[12]*x_tmp[110]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*39] = (tdata->p[13]*x_tmp[118])/tdata->p[4]+(tdata->p[13]*x_tmp[119]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*40] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/tdata->p[4]+(tdata->p[12]*x_tmp[117])/tdata->p[4]+(tdata->p[12]*x_tmp[118])/tdata->p[4]+(tdata->p[12]*x_tmp[119]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*42] = (x_tmp[1]+x_tmp[2]*2.0)/tdata->p[4]+(tdata->p[13]*x_tmp[127])/tdata->p[4]+(tdata->p[13]*x_tmp[128]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*43] = (tdata->p[12]*x_tmp[126])/tdata->p[4]+(tdata->p[12]*x_tmp[127])/tdata->p[4]+(tdata->p[12]*x_tmp[128]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*45] = (tdata->p[13]*x_tmp[136])/tdata->p[4]+(tdata->p[13]*x_tmp[137]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*46] = (tdata->p[12]*x_tmp[135])/tdata->p[4]+(tdata->p[12]*x_tmp[136])/tdata->p[4]+(tdata->p[12]*x_tmp[137]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*48] = (tdata->p[13]*x_tmp[145])/tdata->p[4]+(tdata->p[13]*x_tmp[146]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*49] = (tdata->p[12]*x_tmp[144])/tdata->p[4]+(tdata->p[12]*x_tmp[145])/tdata->p[4]+(tdata->p[12]*x_tmp[146]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*51] = (tdata->p[13]*x_tmp[154])/tdata->p[4]+(tdata->p[13]*x_tmp[155]*2.0)/tdata->p[4];
  rdata->y[it + udata->nt*52] = (tdata->p[12]*x_tmp[153])/tdata->p[4]+(tdata->p[12]*x_tmp[154])/tdata->p[4]+(tdata->p[12]*x_tmp[155]*2.0)/tdata->p[4];
return(status);

}


