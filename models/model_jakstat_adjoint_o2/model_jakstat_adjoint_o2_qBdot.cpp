
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_dwdp.h"
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

int qBdot_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *xB_tmp = nullptr;
if(xB)
    xB_tmp = N_VGetArrayPointer(xB);
realtype *dxB_tmp = nullptr;
if(dxB)
    dxB_tmp = N_VGetArrayPointer(dxB);
realtype *qBdot_tmp = nullptr;
if(qBdot)
    qBdot_tmp = N_VGetArrayPointer(qBdot);
int ip;
memset(qBdot_tmp,0,sizeof(realtype)*udata->nplist*model->nJ);
status = dwdp_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  qBdot_tmp[ip + udata->nplist*0] = -tdata->w[0]*x_tmp[0]*xB_tmp[1]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*1] = tdata->w[0]*x_tmp[9]*xB_tmp[0]-tdata->w[0]*x_tmp[0]*xB_tmp[10]-tdata->w[0]*x_tmp[9]*xB_tmp[1]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[9];
  qBdot_tmp[ip + udata->nplist*2] = tdata->w[0]*x_tmp[18]*xB_tmp[0]-tdata->w[0]*x_tmp[0]*xB_tmp[19]-tdata->w[0]*x_tmp[18]*xB_tmp[1]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[18];
  qBdot_tmp[ip + udata->nplist*3] = tdata->w[0]*x_tmp[27]*xB_tmp[0]-tdata->w[0]*x_tmp[0]*xB_tmp[28]-tdata->w[0]*x_tmp[27]*xB_tmp[1]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[27];
  qBdot_tmp[ip + udata->nplist*4] = tdata->w[0]*x_tmp[36]*xB_tmp[0]-tdata->w[0]*x_tmp[0]*xB_tmp[37]-tdata->w[0]*x_tmp[36]*xB_tmp[1]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[36];
  qBdot_tmp[ip + udata->nplist*5] = tdata->w[0]*x_tmp[45]*xB_tmp[0]-tdata->w[0]*x_tmp[0]*xB_tmp[46]-tdata->w[0]*x_tmp[45]*xB_tmp[1]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[45];
  qBdot_tmp[ip + udata->nplist*6] = xB_tmp[0]*(x_tmp[0]*tdata->w[5]+tdata->w[0]*x_tmp[54])-xB_tmp[1]*(x_tmp[0]*tdata->w[5]+tdata->w[0]*x_tmp[54])-tdata->w[0]*x_tmp[0]*xB_tmp[55]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[54];
  qBdot_tmp[ip + udata->nplist*7] = xB_tmp[0]*(x_tmp[0]*tdata->w[6]+tdata->w[0]*x_tmp[63])-xB_tmp[1]*(x_tmp[0]*tdata->w[6]+tdata->w[0]*x_tmp[63])-tdata->w[0]*x_tmp[0]*xB_tmp[64]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[63];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[0]*(x_tmp[0]*tdata->w[7]+tdata->w[0]*x_tmp[72])-xB_tmp[1]*(x_tmp[0]*tdata->w[7]+tdata->w[0]*x_tmp[72])-tdata->w[0]*x_tmp[0]*xB_tmp[73]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[72];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[0]*(x_tmp[0]*tdata->w[8]+tdata->w[0]*x_tmp[81])-xB_tmp[1]*(x_tmp[0]*tdata->w[8]+tdata->w[0]*x_tmp[81])-tdata->w[0]*x_tmp[0]*xB_tmp[82]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[81];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[0]*(x_tmp[0]*tdata->w[9]+tdata->w[0]*x_tmp[90])-xB_tmp[1]*(x_tmp[0]*tdata->w[9]+tdata->w[0]*x_tmp[90])-tdata->w[0]*x_tmp[0]*xB_tmp[91]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[90];
  qBdot_tmp[ip + udata->nplist*11] = tdata->w[0]*xB_tmp[0]*x_tmp[99]-tdata->w[0]*x_tmp[0]*xB_tmp[100]-tdata->w[0]*xB_tmp[1]*x_tmp[99]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[99];
  qBdot_tmp[ip + udata->nplist*12] = tdata->w[0]*xB_tmp[0]*x_tmp[108]-tdata->w[0]*x_tmp[0]*xB_tmp[109]-tdata->w[0]*xB_tmp[1]*x_tmp[108]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[108];
  qBdot_tmp[ip + udata->nplist*13] = tdata->w[0]*xB_tmp[0]*x_tmp[117]-tdata->w[0]*x_tmp[0]*xB_tmp[118]-tdata->w[0]*xB_tmp[1]*x_tmp[117]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[117];
  qBdot_tmp[ip + udata->nplist*14] = tdata->w[0]*xB_tmp[0]*x_tmp[126]-tdata->w[0]*x_tmp[0]*xB_tmp[127]-tdata->w[0]*xB_tmp[1]*x_tmp[126]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[126];
  qBdot_tmp[ip + udata->nplist*15] = tdata->w[0]*xB_tmp[0]*x_tmp[135]-tdata->w[0]*x_tmp[0]*xB_tmp[136]-tdata->w[0]*xB_tmp[1]*x_tmp[135]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[135];
  qBdot_tmp[ip + udata->nplist*16] = tdata->w[0]*xB_tmp[0]*x_tmp[144]-tdata->w[0]*x_tmp[0]*xB_tmp[145]-tdata->w[0]*xB_tmp[1]*x_tmp[144]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[144];
  qBdot_tmp[ip + udata->nplist*17] = tdata->w[0]*xB_tmp[0]*x_tmp[153]-tdata->w[0]*x_tmp[0]*xB_tmp[154]-tdata->w[0]*xB_tmp[1]*x_tmp[153]+udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2]*xB_tmp[153];

  } break;

  case 1: {
  qBdot_tmp[ip + udata->nplist*0] = tdata->w[1]*xB_tmp[1]*2.0-tdata->w[1]*xB_tmp[2];
  qBdot_tmp[ip + udata->nplist*1] = tdata->w[1]*xB_tmp[10]*2.0-tdata->w[1]*xB_tmp[11]+x_tmp[1]*x_tmp[10]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[10]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*2] = tdata->w[1]*xB_tmp[19]*2.0-tdata->w[1]*xB_tmp[20]+x_tmp[1]*x_tmp[19]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[19]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*3] = tdata->w[1]*xB_tmp[28]*2.0-tdata->w[1]*xB_tmp[29]+x_tmp[1]*x_tmp[28]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[28]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*4] = tdata->w[1]*xB_tmp[37]*2.0-tdata->w[1]*xB_tmp[38]+x_tmp[1]*x_tmp[37]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[37]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*5] = tdata->w[1]*xB_tmp[46]*2.0-tdata->w[1]*xB_tmp[47]+x_tmp[1]*x_tmp[46]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[46]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*6] = tdata->w[1]*xB_tmp[55]*2.0-tdata->w[1]*xB_tmp[56]+x_tmp[1]*x_tmp[55]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[55]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*7] = tdata->w[1]*xB_tmp[64]*2.0-tdata->w[1]*xB_tmp[65]+x_tmp[1]*x_tmp[64]*xB_tmp[1]*4.0-x_tmp[1]*x_tmp[64]*xB_tmp[2]*2.0;
  qBdot_tmp[ip + udata->nplist*8] = tdata->w[1]*xB_tmp[73]*2.0-tdata->w[1]*xB_tmp[74]+x_tmp[1]*xB_tmp[1]*x_tmp[73]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[73]*2.0;
  qBdot_tmp[ip + udata->nplist*9] = tdata->w[1]*xB_tmp[82]*2.0-tdata->w[1]*xB_tmp[83]+x_tmp[1]*xB_tmp[1]*x_tmp[82]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[82]*2.0;
  qBdot_tmp[ip + udata->nplist*10] = tdata->w[1]*xB_tmp[91]*2.0-tdata->w[1]*xB_tmp[92]+x_tmp[1]*xB_tmp[1]*x_tmp[91]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[91]*2.0;
  qBdot_tmp[ip + udata->nplist*11] = tdata->w[1]*xB_tmp[100]*2.0-tdata->w[1]*xB_tmp[101]+x_tmp[1]*xB_tmp[1]*x_tmp[100]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[100]*2.0;
  qBdot_tmp[ip + udata->nplist*12] = tdata->w[1]*xB_tmp[109]*2.0-tdata->w[1]*xB_tmp[110]+x_tmp[1]*xB_tmp[1]*x_tmp[109]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[109]*2.0;
  qBdot_tmp[ip + udata->nplist*13] = tdata->w[1]*xB_tmp[118]*2.0-tdata->w[1]*xB_tmp[119]+x_tmp[1]*xB_tmp[1]*x_tmp[118]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[118]*2.0;
  qBdot_tmp[ip + udata->nplist*14] = tdata->w[1]*xB_tmp[127]*2.0-tdata->w[1]*xB_tmp[128]+x_tmp[1]*xB_tmp[1]*x_tmp[127]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[127]*2.0;
  qBdot_tmp[ip + udata->nplist*15] = tdata->w[1]*xB_tmp[136]*2.0-tdata->w[1]*xB_tmp[137]+x_tmp[1]*xB_tmp[1]*x_tmp[136]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[136]*2.0;
  qBdot_tmp[ip + udata->nplist*16] = tdata->w[1]*xB_tmp[145]*2.0-tdata->w[1]*xB_tmp[146]+x_tmp[1]*xB_tmp[1]*x_tmp[145]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[145]*2.0;
  qBdot_tmp[ip + udata->nplist*17] = tdata->w[1]*xB_tmp[154]*2.0-tdata->w[1]*xB_tmp[155]+x_tmp[1]*xB_tmp[1]*x_tmp[154]*4.0-x_tmp[1]*xB_tmp[2]*x_tmp[154]*2.0;

  } break;

  case 2: {
  qBdot_tmp[ip + udata->nplist*0] = x_tmp[2]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*1] = x_tmp[2]*xB_tmp[11]+x_tmp[11]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[12]-udata->k[0]*tdata->w[3]*x_tmp[11]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*2] = x_tmp[2]*xB_tmp[20]+x_tmp[20]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[21]-udata->k[0]*tdata->w[3]*x_tmp[20]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*3] = x_tmp[2]*xB_tmp[29]+x_tmp[29]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[30]-udata->k[0]*tdata->w[3]*x_tmp[29]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*4] = x_tmp[2]*xB_tmp[38]+x_tmp[38]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[39]-udata->k[0]*tdata->w[3]*x_tmp[38]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*5] = x_tmp[2]*xB_tmp[47]+x_tmp[47]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[48]-udata->k[0]*tdata->w[3]*x_tmp[47]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*6] = x_tmp[2]*xB_tmp[56]+x_tmp[56]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[57]-udata->k[0]*tdata->w[3]*x_tmp[56]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*7] = x_tmp[2]*xB_tmp[65]+x_tmp[65]*xB_tmp[2]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[66]-udata->k[0]*tdata->w[3]*x_tmp[65]*xB_tmp[3];
  qBdot_tmp[ip + udata->nplist*8] = x_tmp[2]*xB_tmp[74]+xB_tmp[2]*x_tmp[74]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[75]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[74];
  qBdot_tmp[ip + udata->nplist*9] = x_tmp[2]*xB_tmp[83]+xB_tmp[2]*x_tmp[83]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[84]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[83];
  qBdot_tmp[ip + udata->nplist*10] = x_tmp[2]*xB_tmp[92]+xB_tmp[2]*x_tmp[92]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[93]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[92];
  qBdot_tmp[ip + udata->nplist*11] = x_tmp[2]*xB_tmp[101]+xB_tmp[2]*x_tmp[101]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[102]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[101];
  qBdot_tmp[ip + udata->nplist*12] = x_tmp[2]*xB_tmp[110]+xB_tmp[2]*x_tmp[110]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[111]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[110];
  qBdot_tmp[ip + udata->nplist*13] = x_tmp[2]*xB_tmp[119]+xB_tmp[2]*x_tmp[119]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[120]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[119];
  qBdot_tmp[ip + udata->nplist*14] = x_tmp[2]*xB_tmp[128]+xB_tmp[2]*x_tmp[128]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[129]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[128];
  qBdot_tmp[ip + udata->nplist*15] = x_tmp[2]*xB_tmp[137]+xB_tmp[2]*x_tmp[137]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[138]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[137];
  qBdot_tmp[ip + udata->nplist*16] = x_tmp[2]*xB_tmp[146]+xB_tmp[2]*x_tmp[146]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[147]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[146];
  qBdot_tmp[ip + udata->nplist*17] = x_tmp[2]*xB_tmp[155]+xB_tmp[2]*x_tmp[155]-udata->k[0]*tdata->w[3]*x_tmp[2]*xB_tmp[156]-udata->k[0]*tdata->w[3]*xB_tmp[3]*x_tmp[155];

  } break;

  case 3: {
  qBdot_tmp[ip + udata->nplist*0] = -xB_tmp[4]*(tdata->w[4]-x_tmp[4])-xB_tmp[5]*(x_tmp[4]-x_tmp[5])-xB_tmp[6]*(x_tmp[5]-x_tmp[6])-xB_tmp[7]*(x_tmp[6]-x_tmp[7])-xB_tmp[8]*(x_tmp[7]-x_tmp[8])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[3]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*1] = x_tmp[12]*xB_tmp[3]-xB_tmp[13]*(tdata->w[4]-x_tmp[4])-xB_tmp[14]*(x_tmp[4]-x_tmp[5])-xB_tmp[15]*(x_tmp[5]-x_tmp[6])-xB_tmp[16]*(x_tmp[6]-x_tmp[7])-xB_tmp[5]*(x_tmp[13]-x_tmp[14])-xB_tmp[17]*(x_tmp[7]-x_tmp[8])-xB_tmp[6]*(x_tmp[14]-x_tmp[15])-xB_tmp[7]*(x_tmp[15]-x_tmp[16])-xB_tmp[8]*(x_tmp[16]-x_tmp[17])-xB_tmp[4]*(x_tmp[12]*2.0-x_tmp[13])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[12]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[9]-udata->k[1]*tdata->w[2]*x_tmp[17]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*2] = x_tmp[21]*xB_tmp[3]-xB_tmp[22]*(tdata->w[4]-x_tmp[4])-xB_tmp[23]*(x_tmp[4]-x_tmp[5])-xB_tmp[24]*(x_tmp[5]-x_tmp[6])-xB_tmp[25]*(x_tmp[6]-x_tmp[7])-xB_tmp[26]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[22]-x_tmp[23])-xB_tmp[6]*(x_tmp[23]-x_tmp[24])-xB_tmp[7]*(x_tmp[24]-x_tmp[25])-xB_tmp[8]*(x_tmp[25]-x_tmp[26])-xB_tmp[4]*(x_tmp[21]*2.0-x_tmp[22])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[21]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[18]-udata->k[1]*tdata->w[2]*x_tmp[26]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*3] = x_tmp[30]*xB_tmp[3]-xB_tmp[31]*(tdata->w[4]-x_tmp[4])-xB_tmp[32]*(x_tmp[4]-x_tmp[5])-xB_tmp[33]*(x_tmp[5]-x_tmp[6])-xB_tmp[34]*(x_tmp[6]-x_tmp[7])-xB_tmp[35]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[31]-x_tmp[32])-xB_tmp[6]*(x_tmp[32]-x_tmp[33])-xB_tmp[7]*(x_tmp[33]-x_tmp[34])-xB_tmp[8]*(x_tmp[34]-x_tmp[35])-xB_tmp[4]*(x_tmp[30]*2.0-x_tmp[31])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[30]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[27]-udata->k[1]*tdata->w[2]*x_tmp[35]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*4] = x_tmp[39]*xB_tmp[3]-xB_tmp[40]*(tdata->w[4]-x_tmp[4])-xB_tmp[41]*(x_tmp[4]-x_tmp[5])-xB_tmp[42]*(x_tmp[5]-x_tmp[6])-xB_tmp[43]*(x_tmp[6]-x_tmp[7])-xB_tmp[44]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[40]-x_tmp[41])-xB_tmp[6]*(x_tmp[41]-x_tmp[42])-xB_tmp[7]*(x_tmp[42]-x_tmp[43])-xB_tmp[8]*(x_tmp[43]-x_tmp[44])-xB_tmp[4]*(x_tmp[39]*2.0-x_tmp[40])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[39]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[36]-udata->k[1]*tdata->w[2]*x_tmp[44]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*5] = x_tmp[48]*xB_tmp[3]-xB_tmp[49]*(tdata->w[4]-x_tmp[4])-xB_tmp[50]*(x_tmp[4]-x_tmp[5])-xB_tmp[51]*(x_tmp[5]-x_tmp[6])-xB_tmp[52]*(x_tmp[6]-x_tmp[7])-xB_tmp[53]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[49]-x_tmp[50])-xB_tmp[6]*(x_tmp[50]-x_tmp[51])-xB_tmp[7]*(x_tmp[51]-x_tmp[52])-xB_tmp[8]*(x_tmp[52]-x_tmp[53])-xB_tmp[4]*(x_tmp[48]*2.0-x_tmp[49])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[48]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[45]-udata->k[1]*tdata->w[2]*x_tmp[53]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*6] = x_tmp[57]*xB_tmp[3]-xB_tmp[58]*(tdata->w[4]-x_tmp[4])-xB_tmp[59]*(x_tmp[4]-x_tmp[5])-xB_tmp[60]*(x_tmp[5]-x_tmp[6])-xB_tmp[61]*(x_tmp[6]-x_tmp[7])-xB_tmp[62]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[58]-x_tmp[59])-xB_tmp[6]*(x_tmp[59]-x_tmp[60])-xB_tmp[7]*(x_tmp[60]-x_tmp[61])-xB_tmp[8]*(x_tmp[61]-x_tmp[62])-xB_tmp[4]*(x_tmp[57]*2.0-x_tmp[58])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[57]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[54]-udata->k[1]*tdata->w[2]*x_tmp[62]*xB_tmp[0];
  qBdot_tmp[ip + udata->nplist*7] = x_tmp[66]*xB_tmp[3]-xB_tmp[67]*(tdata->w[4]-x_tmp[4])-xB_tmp[68]*(x_tmp[4]-x_tmp[5])-xB_tmp[69]*(x_tmp[5]-x_tmp[6])-xB_tmp[70]*(x_tmp[6]-x_tmp[7])-xB_tmp[71]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[67]-x_tmp[68])-xB_tmp[6]*(x_tmp[68]-x_tmp[69])-xB_tmp[7]*(x_tmp[69]-x_tmp[70])-xB_tmp[8]*(x_tmp[70]-x_tmp[71])-xB_tmp[4]*(x_tmp[66]*2.0-x_tmp[67])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[66]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[63]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[71];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[3]*x_tmp[75]-xB_tmp[76]*(tdata->w[4]-x_tmp[4])-xB_tmp[77]*(x_tmp[4]-x_tmp[5])-xB_tmp[78]*(x_tmp[5]-x_tmp[6])-xB_tmp[79]*(x_tmp[6]-x_tmp[7])-xB_tmp[80]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[76]-x_tmp[77])-xB_tmp[6]*(x_tmp[77]-x_tmp[78])-xB_tmp[7]*(x_tmp[78]-x_tmp[79])-xB_tmp[8]*(x_tmp[79]-x_tmp[80])-xB_tmp[4]*(x_tmp[75]*2.0-x_tmp[76])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[75]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[72]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[80];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[3]*x_tmp[84]-xB_tmp[85]*(tdata->w[4]-x_tmp[4])-xB_tmp[86]*(x_tmp[4]-x_tmp[5])-xB_tmp[87]*(x_tmp[5]-x_tmp[6])-xB_tmp[88]*(x_tmp[6]-x_tmp[7])-xB_tmp[89]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[85]-x_tmp[86])-xB_tmp[6]*(x_tmp[86]-x_tmp[87])-xB_tmp[7]*(x_tmp[87]-x_tmp[88])-xB_tmp[8]*(x_tmp[88]-x_tmp[89])-xB_tmp[4]*(x_tmp[84]*2.0-x_tmp[85])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[84]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[81]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[89];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[3]*x_tmp[93]-xB_tmp[94]*(tdata->w[4]-x_tmp[4])-xB_tmp[95]*(x_tmp[4]-x_tmp[5])-xB_tmp[96]*(x_tmp[5]-x_tmp[6])-xB_tmp[97]*(x_tmp[6]-x_tmp[7])-xB_tmp[98]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[94]-x_tmp[95])-xB_tmp[6]*(x_tmp[95]-x_tmp[96])-xB_tmp[7]*(x_tmp[96]-x_tmp[97])-xB_tmp[8]*(x_tmp[97]-x_tmp[98])-xB_tmp[4]*(x_tmp[93]*2.0-x_tmp[94])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[93]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[90]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[98];
  qBdot_tmp[ip + udata->nplist*11] = xB_tmp[3]*x_tmp[102]-xB_tmp[103]*(tdata->w[4]-x_tmp[4])-xB_tmp[104]*(x_tmp[4]-x_tmp[5])-xB_tmp[105]*(x_tmp[5]-x_tmp[6])-xB_tmp[106]*(x_tmp[6]-x_tmp[7])-xB_tmp[107]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[103]-x_tmp[104])-xB_tmp[6]*(x_tmp[104]-x_tmp[105])-xB_tmp[7]*(x_tmp[105]-x_tmp[106])-xB_tmp[8]*(x_tmp[106]-x_tmp[107])-xB_tmp[4]*(x_tmp[102]*2.0-x_tmp[103])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[102]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[99]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[107];
  qBdot_tmp[ip + udata->nplist*12] = xB_tmp[3]*x_tmp[111]-xB_tmp[112]*(tdata->w[4]-x_tmp[4])-xB_tmp[113]*(x_tmp[4]-x_tmp[5])-xB_tmp[114]*(x_tmp[5]-x_tmp[6])-xB_tmp[115]*(x_tmp[6]-x_tmp[7])-xB_tmp[116]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[112]-x_tmp[113])-xB_tmp[6]*(x_tmp[113]-x_tmp[114])-xB_tmp[7]*(x_tmp[114]-x_tmp[115])-xB_tmp[8]*(x_tmp[115]-x_tmp[116])-xB_tmp[4]*(x_tmp[111]*2.0-x_tmp[112])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[111]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[108]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[116];
  qBdot_tmp[ip + udata->nplist*13] = xB_tmp[3]*x_tmp[120]-xB_tmp[121]*(tdata->w[4]-x_tmp[4])-xB_tmp[122]*(x_tmp[4]-x_tmp[5])-xB_tmp[123]*(x_tmp[5]-x_tmp[6])-xB_tmp[124]*(x_tmp[6]-x_tmp[7])-xB_tmp[125]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[121]-x_tmp[122])-xB_tmp[6]*(x_tmp[122]-x_tmp[123])-xB_tmp[7]*(x_tmp[123]-x_tmp[124])-xB_tmp[8]*(x_tmp[124]-x_tmp[125])-xB_tmp[4]*(x_tmp[120]*2.0-x_tmp[121])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[120]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[117]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[125];
  qBdot_tmp[ip + udata->nplist*14] = xB_tmp[3]*x_tmp[129]-xB_tmp[130]*(tdata->w[4]-x_tmp[4])-xB_tmp[131]*(x_tmp[4]-x_tmp[5])-xB_tmp[132]*(x_tmp[5]-x_tmp[6])-xB_tmp[133]*(x_tmp[6]-x_tmp[7])-xB_tmp[134]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[130]-x_tmp[131])-xB_tmp[6]*(x_tmp[131]-x_tmp[132])-xB_tmp[7]*(x_tmp[132]-x_tmp[133])-xB_tmp[8]*(x_tmp[133]-x_tmp[134])-xB_tmp[4]*(x_tmp[129]*2.0-x_tmp[130])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[129]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[126]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[134];
  qBdot_tmp[ip + udata->nplist*15] = xB_tmp[3]*x_tmp[138]-xB_tmp[139]*(tdata->w[4]-x_tmp[4])-xB_tmp[140]*(x_tmp[4]-x_tmp[5])-xB_tmp[141]*(x_tmp[5]-x_tmp[6])-xB_tmp[142]*(x_tmp[6]-x_tmp[7])-xB_tmp[143]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[139]-x_tmp[140])-xB_tmp[6]*(x_tmp[140]-x_tmp[141])-xB_tmp[7]*(x_tmp[141]-x_tmp[142])-xB_tmp[8]*(x_tmp[142]-x_tmp[143])-xB_tmp[4]*(x_tmp[138]*2.0-x_tmp[139])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[138]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[135]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[143];
  qBdot_tmp[ip + udata->nplist*16] = xB_tmp[3]*x_tmp[147]-xB_tmp[148]*(tdata->w[4]-x_tmp[4])-xB_tmp[149]*(x_tmp[4]-x_tmp[5])-xB_tmp[150]*(x_tmp[5]-x_tmp[6])-xB_tmp[151]*(x_tmp[6]-x_tmp[7])-xB_tmp[152]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[148]-x_tmp[149])-xB_tmp[6]*(x_tmp[149]-x_tmp[150])-xB_tmp[7]*(x_tmp[150]-x_tmp[151])-xB_tmp[8]*(x_tmp[151]-x_tmp[152])-xB_tmp[4]*(x_tmp[147]*2.0-x_tmp[148])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[147]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[144]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[152];
  qBdot_tmp[ip + udata->nplist*17] = xB_tmp[3]*x_tmp[156]-xB_tmp[157]*(tdata->w[4]-x_tmp[4])-xB_tmp[158]*(x_tmp[4]-x_tmp[5])-xB_tmp[159]*(x_tmp[5]-x_tmp[6])-xB_tmp[160]*(x_tmp[6]-x_tmp[7])-xB_tmp[161]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[157]-x_tmp[158])-xB_tmp[6]*(x_tmp[158]-x_tmp[159])-xB_tmp[7]*(x_tmp[159]-x_tmp[160])-xB_tmp[8]*(x_tmp[160]-x_tmp[161])-xB_tmp[4]*(x_tmp[156]*2.0-x_tmp[157])+udata->k[1]*tdata->w[3]*x_tmp[3]*xB_tmp[156]-udata->k[1]*tdata->w[2]*x_tmp[8]*xB_tmp[153]-udata->k[1]*tdata->w[2]*xB_tmp[0]*x_tmp[161];

  } break;

  case 5: {
  qBdot_tmp[ip + udata->nplist*0] = -tdata->p[0]*x_tmp[0]*xB_tmp[1]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[0]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*1] = xB_tmp[0]*tdata->dwdp[0]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-xB_tmp[1]*tdata->dwdp[0]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-tdata->p[0]*x_tmp[0]*xB_tmp[10]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[9]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*2] = tdata->p[0]*x_tmp[18]*xB_tmp[0]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[19]*tdata->dwdp[0]-tdata->p[0]*x_tmp[18]*xB_tmp[1]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[18]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*3] = tdata->p[0]*x_tmp[27]*xB_tmp[0]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[28]*tdata->dwdp[0]-tdata->p[0]*x_tmp[27]*xB_tmp[1]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[27]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*4] = tdata->p[0]*x_tmp[36]*xB_tmp[0]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[37]*tdata->dwdp[0]-tdata->p[0]*x_tmp[36]*xB_tmp[1]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[36]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*5] = tdata->p[0]*x_tmp[45]*xB_tmp[0]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[46]*tdata->dwdp[0]-tdata->p[0]*x_tmp[45]*xB_tmp[1]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[45]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*6] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[1]+tdata->p[0]*x_tmp[54]*tdata->dwdp[0])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[1]+tdata->p[0]*x_tmp[54]*tdata->dwdp[0])-tdata->p[0]*x_tmp[0]*xB_tmp[55]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[54]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*7] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[2]+tdata->p[0]*x_tmp[63]*tdata->dwdp[0])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[2]+tdata->p[0]*x_tmp[63]*tdata->dwdp[0])-tdata->p[0]*x_tmp[0]*xB_tmp[64]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[63]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[3]+tdata->p[0]*x_tmp[72]*tdata->dwdp[0])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[3]+tdata->p[0]*x_tmp[72]*tdata->dwdp[0])-tdata->p[0]*x_tmp[0]*xB_tmp[73]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[72]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[4]+tdata->p[0]*x_tmp[81]*tdata->dwdp[0])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[4]+tdata->p[0]*x_tmp[81]*tdata->dwdp[0])-tdata->p[0]*x_tmp[0]*xB_tmp[82]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[81]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[5]+tdata->p[0]*x_tmp[90]*tdata->dwdp[0])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[5]+tdata->p[0]*x_tmp[90]*tdata->dwdp[0])-tdata->p[0]*x_tmp[0]*xB_tmp[91]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[90]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*11] = tdata->p[0]*xB_tmp[0]*x_tmp[99]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[100]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[99]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[99]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*12] = tdata->p[0]*xB_tmp[0]*x_tmp[108]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[109]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[108]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[108]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*13] = tdata->p[0]*xB_tmp[0]*x_tmp[117]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[118]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[117]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[117]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*14] = tdata->p[0]*xB_tmp[0]*x_tmp[126]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[127]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[126]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[126]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*15] = tdata->p[0]*xB_tmp[0]*x_tmp[135]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[136]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[135]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[135]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*16] = tdata->p[0]*xB_tmp[0]*x_tmp[144]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[145]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[144]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[144]*tdata->dwdp[0];
  qBdot_tmp[ip + udata->nplist*17] = tdata->p[0]*xB_tmp[0]*x_tmp[153]*tdata->dwdp[0]-tdata->p[0]*x_tmp[0]*xB_tmp[154]*tdata->dwdp[0]-tdata->p[0]*xB_tmp[1]*x_tmp[153]*tdata->dwdp[0]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[153]*tdata->dwdp[0];

  } break;

  case 6: {
  qBdot_tmp[ip + udata->nplist*0] = -tdata->p[0]*x_tmp[0]*xB_tmp[1]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[0]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*1] = xB_tmp[0]*tdata->dwdp[6]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-xB_tmp[1]*tdata->dwdp[6]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-tdata->p[0]*x_tmp[0]*xB_tmp[10]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[9]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*2] = tdata->p[0]*x_tmp[18]*xB_tmp[0]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[19]*tdata->dwdp[6]-tdata->p[0]*x_tmp[18]*xB_tmp[1]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[18]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*3] = tdata->p[0]*x_tmp[27]*xB_tmp[0]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[28]*tdata->dwdp[6]-tdata->p[0]*x_tmp[27]*xB_tmp[1]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[27]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*4] = tdata->p[0]*x_tmp[36]*xB_tmp[0]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[37]*tdata->dwdp[6]-tdata->p[0]*x_tmp[36]*xB_tmp[1]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[36]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*5] = tdata->p[0]*x_tmp[45]*xB_tmp[0]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[46]*tdata->dwdp[6]-tdata->p[0]*x_tmp[45]*xB_tmp[1]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[45]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*6] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[7]+tdata->p[0]*x_tmp[54]*tdata->dwdp[6])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[7]+tdata->p[0]*x_tmp[54]*tdata->dwdp[6])-tdata->p[0]*x_tmp[0]*xB_tmp[55]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[54]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*7] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[8]+tdata->p[0]*x_tmp[63]*tdata->dwdp[6])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[8]+tdata->p[0]*x_tmp[63]*tdata->dwdp[6])-tdata->p[0]*x_tmp[0]*xB_tmp[64]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[63]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[9]+tdata->p[0]*x_tmp[72]*tdata->dwdp[6])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[9]+tdata->p[0]*x_tmp[72]*tdata->dwdp[6])-tdata->p[0]*x_tmp[0]*xB_tmp[73]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[72]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[10]+tdata->p[0]*x_tmp[81]*tdata->dwdp[6])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[10]+tdata->p[0]*x_tmp[81]*tdata->dwdp[6])-tdata->p[0]*x_tmp[0]*xB_tmp[82]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[81]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[11]+tdata->p[0]*x_tmp[90]*tdata->dwdp[6])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[11]+tdata->p[0]*x_tmp[90]*tdata->dwdp[6])-tdata->p[0]*x_tmp[0]*xB_tmp[91]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[90]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*11] = tdata->p[0]*xB_tmp[0]*x_tmp[99]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[100]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[99]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[99]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*12] = tdata->p[0]*xB_tmp[0]*x_tmp[108]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[109]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[108]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[108]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*13] = tdata->p[0]*xB_tmp[0]*x_tmp[117]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[118]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[117]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[117]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*14] = tdata->p[0]*xB_tmp[0]*x_tmp[126]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[127]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[126]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[126]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*15] = tdata->p[0]*xB_tmp[0]*x_tmp[135]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[136]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[135]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[135]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*16] = tdata->p[0]*xB_tmp[0]*x_tmp[144]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[145]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[144]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[144]*tdata->dwdp[6];
  qBdot_tmp[ip + udata->nplist*17] = tdata->p[0]*xB_tmp[0]*x_tmp[153]*tdata->dwdp[6]-tdata->p[0]*x_tmp[0]*xB_tmp[154]*tdata->dwdp[6]-tdata->p[0]*xB_tmp[1]*x_tmp[153]*tdata->dwdp[6]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[153]*tdata->dwdp[6];

  } break;

  case 7: {
  qBdot_tmp[ip + udata->nplist*0] = -tdata->p[0]*x_tmp[0]*xB_tmp[1]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[0]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*1] = xB_tmp[0]*tdata->dwdp[12]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-xB_tmp[1]*tdata->dwdp[12]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-tdata->p[0]*x_tmp[0]*xB_tmp[10]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[9]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*2] = tdata->p[0]*x_tmp[18]*xB_tmp[0]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[19]*tdata->dwdp[12]-tdata->p[0]*x_tmp[18]*xB_tmp[1]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[18]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*3] = tdata->p[0]*x_tmp[27]*xB_tmp[0]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[28]*tdata->dwdp[12]-tdata->p[0]*x_tmp[27]*xB_tmp[1]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[27]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*4] = tdata->p[0]*x_tmp[36]*xB_tmp[0]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[37]*tdata->dwdp[12]-tdata->p[0]*x_tmp[36]*xB_tmp[1]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[36]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*5] = tdata->p[0]*x_tmp[45]*xB_tmp[0]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[46]*tdata->dwdp[12]-tdata->p[0]*x_tmp[45]*xB_tmp[1]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[45]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*6] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[13]+tdata->p[0]*x_tmp[54]*tdata->dwdp[12])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[13]+tdata->p[0]*x_tmp[54]*tdata->dwdp[12])-tdata->p[0]*x_tmp[0]*xB_tmp[55]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[54]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*7] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[14]+tdata->p[0]*x_tmp[63]*tdata->dwdp[12])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[14]+tdata->p[0]*x_tmp[63]*tdata->dwdp[12])-tdata->p[0]*x_tmp[0]*xB_tmp[64]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[63]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[15]+tdata->p[0]*x_tmp[72]*tdata->dwdp[12])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[15]+tdata->p[0]*x_tmp[72]*tdata->dwdp[12])-tdata->p[0]*x_tmp[0]*xB_tmp[73]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[72]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[16]+tdata->p[0]*x_tmp[81]*tdata->dwdp[12])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[16]+tdata->p[0]*x_tmp[81]*tdata->dwdp[12])-tdata->p[0]*x_tmp[0]*xB_tmp[82]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[81]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[17]+tdata->p[0]*x_tmp[90]*tdata->dwdp[12])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[17]+tdata->p[0]*x_tmp[90]*tdata->dwdp[12])-tdata->p[0]*x_tmp[0]*xB_tmp[91]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[90]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*11] = tdata->p[0]*xB_tmp[0]*x_tmp[99]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[100]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[99]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[99]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*12] = tdata->p[0]*xB_tmp[0]*x_tmp[108]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[109]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[108]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[108]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*13] = tdata->p[0]*xB_tmp[0]*x_tmp[117]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[118]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[117]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[117]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*14] = tdata->p[0]*xB_tmp[0]*x_tmp[126]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[127]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[126]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[126]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*15] = tdata->p[0]*xB_tmp[0]*x_tmp[135]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[136]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[135]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[135]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*16] = tdata->p[0]*xB_tmp[0]*x_tmp[144]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[145]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[144]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[144]*tdata->dwdp[12];
  qBdot_tmp[ip + udata->nplist*17] = tdata->p[0]*xB_tmp[0]*x_tmp[153]*tdata->dwdp[12]-tdata->p[0]*x_tmp[0]*xB_tmp[154]*tdata->dwdp[12]-tdata->p[0]*xB_tmp[1]*x_tmp[153]*tdata->dwdp[12]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[153]*tdata->dwdp[12];

  } break;

  case 8: {
  qBdot_tmp[ip + udata->nplist*0] = -tdata->p[0]*x_tmp[0]*xB_tmp[1]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[0]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*1] = xB_tmp[0]*tdata->dwdp[18]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-xB_tmp[1]*tdata->dwdp[18]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-tdata->p[0]*x_tmp[0]*xB_tmp[10]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[9]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*2] = tdata->p[0]*x_tmp[18]*xB_tmp[0]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[19]*tdata->dwdp[18]-tdata->p[0]*x_tmp[18]*xB_tmp[1]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[18]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*3] = tdata->p[0]*x_tmp[27]*xB_tmp[0]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[28]*tdata->dwdp[18]-tdata->p[0]*x_tmp[27]*xB_tmp[1]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[27]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*4] = tdata->p[0]*x_tmp[36]*xB_tmp[0]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[37]*tdata->dwdp[18]-tdata->p[0]*x_tmp[36]*xB_tmp[1]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[36]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*5] = tdata->p[0]*x_tmp[45]*xB_tmp[0]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[46]*tdata->dwdp[18]-tdata->p[0]*x_tmp[45]*xB_tmp[1]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[45]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*6] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[19]+tdata->p[0]*x_tmp[54]*tdata->dwdp[18])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[19]+tdata->p[0]*x_tmp[54]*tdata->dwdp[18])-tdata->p[0]*x_tmp[0]*xB_tmp[55]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[54]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*7] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[20]+tdata->p[0]*x_tmp[63]*tdata->dwdp[18])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[20]+tdata->p[0]*x_tmp[63]*tdata->dwdp[18])-tdata->p[0]*x_tmp[0]*xB_tmp[64]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[63]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[21]+tdata->p[0]*x_tmp[72]*tdata->dwdp[18])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[21]+tdata->p[0]*x_tmp[72]*tdata->dwdp[18])-tdata->p[0]*x_tmp[0]*xB_tmp[73]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[72]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[22]+tdata->p[0]*x_tmp[81]*tdata->dwdp[18])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[22]+tdata->p[0]*x_tmp[81]*tdata->dwdp[18])-tdata->p[0]*x_tmp[0]*xB_tmp[82]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[81]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[23]+tdata->p[0]*x_tmp[90]*tdata->dwdp[18])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[23]+tdata->p[0]*x_tmp[90]*tdata->dwdp[18])-tdata->p[0]*x_tmp[0]*xB_tmp[91]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[90]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*11] = tdata->p[0]*xB_tmp[0]*x_tmp[99]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[100]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[99]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[99]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*12] = tdata->p[0]*xB_tmp[0]*x_tmp[108]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[109]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[108]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[108]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*13] = tdata->p[0]*xB_tmp[0]*x_tmp[117]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[118]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[117]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[117]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*14] = tdata->p[0]*xB_tmp[0]*x_tmp[126]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[127]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[126]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[126]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*15] = tdata->p[0]*xB_tmp[0]*x_tmp[135]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[136]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[135]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[135]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*16] = tdata->p[0]*xB_tmp[0]*x_tmp[144]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[145]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[144]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[144]*tdata->dwdp[18];
  qBdot_tmp[ip + udata->nplist*17] = tdata->p[0]*xB_tmp[0]*x_tmp[153]*tdata->dwdp[18]-tdata->p[0]*x_tmp[0]*xB_tmp[154]*tdata->dwdp[18]-tdata->p[0]*xB_tmp[1]*x_tmp[153]*tdata->dwdp[18]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[153]*tdata->dwdp[18];

  } break;

  case 9: {
  qBdot_tmp[ip + udata->nplist*0] = -tdata->p[0]*x_tmp[0]*xB_tmp[1]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[0]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*1] = xB_tmp[0]*tdata->dwdp[24]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-xB_tmp[1]*tdata->dwdp[24]*(x_tmp[0]+tdata->p[0]*x_tmp[9])-tdata->p[0]*x_tmp[0]*xB_tmp[10]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[9]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*2] = tdata->p[0]*x_tmp[18]*xB_tmp[0]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[19]*tdata->dwdp[24]-tdata->p[0]*x_tmp[18]*xB_tmp[1]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[18]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*3] = tdata->p[0]*x_tmp[27]*xB_tmp[0]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[28]*tdata->dwdp[24]-tdata->p[0]*x_tmp[27]*xB_tmp[1]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[27]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*4] = tdata->p[0]*x_tmp[36]*xB_tmp[0]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[37]*tdata->dwdp[24]-tdata->p[0]*x_tmp[36]*xB_tmp[1]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[36]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*5] = tdata->p[0]*x_tmp[45]*xB_tmp[0]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[46]*tdata->dwdp[24]-tdata->p[0]*x_tmp[45]*xB_tmp[1]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[45]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*6] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[25]+tdata->p[0]*x_tmp[54]*tdata->dwdp[24])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[25]+tdata->p[0]*x_tmp[54]*tdata->dwdp[24])-tdata->p[0]*x_tmp[0]*xB_tmp[55]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[54]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*7] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[26]+tdata->p[0]*x_tmp[63]*tdata->dwdp[24])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[26]+tdata->p[0]*x_tmp[63]*tdata->dwdp[24])-tdata->p[0]*x_tmp[0]*xB_tmp[64]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[63]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*8] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[27]+tdata->p[0]*x_tmp[72]*tdata->dwdp[24])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[27]+tdata->p[0]*x_tmp[72]*tdata->dwdp[24])-tdata->p[0]*x_tmp[0]*xB_tmp[73]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[72]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*9] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[28]+tdata->p[0]*x_tmp[81]*tdata->dwdp[24])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[28]+tdata->p[0]*x_tmp[81]*tdata->dwdp[24])-tdata->p[0]*x_tmp[0]*xB_tmp[82]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[81]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*10] = xB_tmp[0]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[29]+tdata->p[0]*x_tmp[90]*tdata->dwdp[24])-xB_tmp[1]*(tdata->p[0]*x_tmp[0]*tdata->dwdp[29]+tdata->p[0]*x_tmp[90]*tdata->dwdp[24])-tdata->p[0]*x_tmp[0]*xB_tmp[91]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[90]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*11] = tdata->p[0]*xB_tmp[0]*x_tmp[99]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[100]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[99]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[99]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*12] = tdata->p[0]*xB_tmp[0]*x_tmp[108]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[109]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[108]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[108]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*13] = tdata->p[0]*xB_tmp[0]*x_tmp[117]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[118]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[117]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[117]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*14] = tdata->p[0]*xB_tmp[0]*x_tmp[126]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[127]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[126]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[126]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*15] = tdata->p[0]*xB_tmp[0]*x_tmp[135]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[136]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[135]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[135]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*16] = tdata->p[0]*xB_tmp[0]*x_tmp[144]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[145]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[144]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[144]*tdata->dwdp[24];
  qBdot_tmp[ip + udata->nplist*17] = tdata->p[0]*xB_tmp[0]*x_tmp[153]*tdata->dwdp[24]-tdata->p[0]*x_tmp[0]*xB_tmp[154]*tdata->dwdp[24]-tdata->p[0]*xB_tmp[1]*x_tmp[153]*tdata->dwdp[24]+udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*xB_tmp[153]*tdata->dwdp[24];

  } break;

}
}
for(ip = 0; ip<udata->nplist*model->nJ; ip++) {
   if(amiIsNaN(qBdot_tmp[ip])) {
       qBdot_tmp[ip] = 0;       if(!tdata->nan_qBdot) {
           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_qBdot = TRUE;
       }
   }   if(amiIsInf(qBdot_tmp[ip])) {
       warnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


