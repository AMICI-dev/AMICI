
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_dwdp.h"
#include "model_jakstat_adjoint_o2_w.h"

int qBdot_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
int ip;
memset(qBdot_tmp,0,sizeof(realtype)*nplist*ng);
status = dwdp_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
  case 0: {
  qBdot_tmp[ip+nplist*0] = -w_tmp[0]*xB_tmp[1]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[0]*x_tmp[0];
  qBdot_tmp[ip+nplist*1] = w_tmp[0]*xB_tmp[0]*x_tmp[9]-w_tmp[0]*xB_tmp[1]*x_tmp[9]-w_tmp[0]*xB_tmp[10]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[9]*x_tmp[0];
  qBdot_tmp[ip+nplist*2] = w_tmp[0]*xB_tmp[0]*x_tmp[18]-w_tmp[0]*xB_tmp[1]*x_tmp[18]-w_tmp[0]*xB_tmp[19]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[18]*x_tmp[0];
  qBdot_tmp[ip+nplist*3] = w_tmp[0]*xB_tmp[0]*x_tmp[27]-w_tmp[0]*xB_tmp[1]*x_tmp[27]-w_tmp[0]*xB_tmp[28]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[27]*x_tmp[0];
  qBdot_tmp[ip+nplist*4] = w_tmp[0]*xB_tmp[0]*x_tmp[36]-w_tmp[0]*xB_tmp[1]*x_tmp[36]-w_tmp[0]*xB_tmp[37]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[36]*x_tmp[0];
  qBdot_tmp[ip+nplist*5] = w_tmp[0]*xB_tmp[0]*x_tmp[45]-w_tmp[0]*xB_tmp[1]*x_tmp[45]-w_tmp[0]*xB_tmp[46]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[45]*x_tmp[0];
  qBdot_tmp[ip+nplist*6] = xB_tmp[0]*(w_tmp[5]*x_tmp[0]+w_tmp[0]*x_tmp[54])-xB_tmp[1]*(w_tmp[5]*x_tmp[0]+w_tmp[0]*x_tmp[54])-w_tmp[0]*xB_tmp[55]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[54]*x_tmp[0];
  qBdot_tmp[ip+nplist*7] = xB_tmp[0]*(w_tmp[6]*x_tmp[0]+w_tmp[0]*x_tmp[63])-xB_tmp[1]*(w_tmp[6]*x_tmp[0]+w_tmp[0]*x_tmp[63])-w_tmp[0]*xB_tmp[64]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[63]*x_tmp[0];
  qBdot_tmp[ip+nplist*8] = xB_tmp[0]*(w_tmp[7]*x_tmp[0]+w_tmp[0]*x_tmp[72])-xB_tmp[1]*(w_tmp[7]*x_tmp[0]+w_tmp[0]*x_tmp[72])-w_tmp[0]*xB_tmp[73]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[72]*x_tmp[0];
  qBdot_tmp[ip+nplist*9] = xB_tmp[0]*(w_tmp[8]*x_tmp[0]+w_tmp[0]*x_tmp[81])-xB_tmp[1]*(w_tmp[8]*x_tmp[0]+w_tmp[0]*x_tmp[81])-w_tmp[0]*xB_tmp[82]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[81]*x_tmp[0];
  qBdot_tmp[ip+nplist*10] = xB_tmp[0]*(w_tmp[9]*x_tmp[0]+w_tmp[0]*x_tmp[90])-xB_tmp[1]*(w_tmp[9]*x_tmp[0]+w_tmp[0]*x_tmp[90])-w_tmp[0]*xB_tmp[91]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[90]*x_tmp[0];
  qBdot_tmp[ip+nplist*11] = w_tmp[0]*xB_tmp[0]*x_tmp[99]-w_tmp[0]*xB_tmp[1]*x_tmp[99]-w_tmp[0]*xB_tmp[100]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[99]*x_tmp[0];
  qBdot_tmp[ip+nplist*12] = w_tmp[0]*xB_tmp[0]*x_tmp[108]-w_tmp[0]*xB_tmp[1]*x_tmp[108]-w_tmp[0]*xB_tmp[109]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[108]*x_tmp[0];
  qBdot_tmp[ip+nplist*13] = w_tmp[0]*xB_tmp[0]*x_tmp[117]-w_tmp[0]*xB_tmp[1]*x_tmp[117]-w_tmp[0]*xB_tmp[118]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[117]*x_tmp[0];
  qBdot_tmp[ip+nplist*14] = w_tmp[0]*xB_tmp[0]*x_tmp[126]-w_tmp[0]*xB_tmp[1]*x_tmp[126]-w_tmp[0]*xB_tmp[127]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[126]*x_tmp[0];
  qBdot_tmp[ip+nplist*15] = w_tmp[0]*xB_tmp[0]*x_tmp[135]-w_tmp[0]*xB_tmp[1]*x_tmp[135]-w_tmp[0]*xB_tmp[136]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[135]*x_tmp[0];
  qBdot_tmp[ip+nplist*16] = w_tmp[0]*xB_tmp[0]*x_tmp[144]-w_tmp[0]*xB_tmp[1]*x_tmp[144]-w_tmp[0]*xB_tmp[145]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[144]*x_tmp[0];
  qBdot_tmp[ip+nplist*17] = w_tmp[0]*xB_tmp[0]*x_tmp[153]-w_tmp[0]*xB_tmp[1]*x_tmp[153]-w_tmp[0]*xB_tmp[154]*x_tmp[0]+k[0]*w_tmp[0]*w_tmp[2]*xB_tmp[153]*x_tmp[0];

  } break;

  case 1: {
  qBdot_tmp[ip+nplist*0] = w_tmp[1]*xB_tmp[1]*2.0-w_tmp[1]*xB_tmp[2];
  qBdot_tmp[ip+nplist*1] = w_tmp[1]*xB_tmp[10]*2.0-w_tmp[1]*xB_tmp[11]+xB_tmp[1]*x_tmp[1]*x_tmp[10]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[10]*2.0;
  qBdot_tmp[ip+nplist*2] = w_tmp[1]*xB_tmp[19]*2.0-w_tmp[1]*xB_tmp[20]+xB_tmp[1]*x_tmp[1]*x_tmp[19]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[19]*2.0;
  qBdot_tmp[ip+nplist*3] = w_tmp[1]*xB_tmp[28]*2.0-w_tmp[1]*xB_tmp[29]+xB_tmp[1]*x_tmp[1]*x_tmp[28]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[28]*2.0;
  qBdot_tmp[ip+nplist*4] = w_tmp[1]*xB_tmp[37]*2.0-w_tmp[1]*xB_tmp[38]+xB_tmp[1]*x_tmp[1]*x_tmp[37]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[37]*2.0;
  qBdot_tmp[ip+nplist*5] = w_tmp[1]*xB_tmp[46]*2.0-w_tmp[1]*xB_tmp[47]+xB_tmp[1]*x_tmp[1]*x_tmp[46]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[46]*2.0;
  qBdot_tmp[ip+nplist*6] = w_tmp[1]*xB_tmp[55]*2.0-w_tmp[1]*xB_tmp[56]+xB_tmp[1]*x_tmp[1]*x_tmp[55]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[55]*2.0;
  qBdot_tmp[ip+nplist*7] = w_tmp[1]*xB_tmp[64]*2.0-w_tmp[1]*xB_tmp[65]+xB_tmp[1]*x_tmp[1]*x_tmp[64]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[64]*2.0;
  qBdot_tmp[ip+nplist*8] = w_tmp[1]*xB_tmp[73]*2.0-w_tmp[1]*xB_tmp[74]+xB_tmp[1]*x_tmp[1]*x_tmp[73]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[73]*2.0;
  qBdot_tmp[ip+nplist*9] = w_tmp[1]*xB_tmp[82]*2.0-w_tmp[1]*xB_tmp[83]+xB_tmp[1]*x_tmp[1]*x_tmp[82]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[82]*2.0;
  qBdot_tmp[ip+nplist*10] = w_tmp[1]*xB_tmp[91]*2.0-w_tmp[1]*xB_tmp[92]+xB_tmp[1]*x_tmp[1]*x_tmp[91]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[91]*2.0;
  qBdot_tmp[ip+nplist*11] = w_tmp[1]*xB_tmp[100]*2.0-w_tmp[1]*xB_tmp[101]+xB_tmp[1]*x_tmp[1]*x_tmp[100]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[100]*2.0;
  qBdot_tmp[ip+nplist*12] = w_tmp[1]*xB_tmp[109]*2.0-w_tmp[1]*xB_tmp[110]+xB_tmp[1]*x_tmp[1]*x_tmp[109]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[109]*2.0;
  qBdot_tmp[ip+nplist*13] = w_tmp[1]*xB_tmp[118]*2.0-w_tmp[1]*xB_tmp[119]+xB_tmp[1]*x_tmp[1]*x_tmp[118]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[118]*2.0;
  qBdot_tmp[ip+nplist*14] = w_tmp[1]*xB_tmp[127]*2.0-w_tmp[1]*xB_tmp[128]+xB_tmp[1]*x_tmp[1]*x_tmp[127]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[127]*2.0;
  qBdot_tmp[ip+nplist*15] = w_tmp[1]*xB_tmp[136]*2.0-w_tmp[1]*xB_tmp[137]+xB_tmp[1]*x_tmp[1]*x_tmp[136]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[136]*2.0;
  qBdot_tmp[ip+nplist*16] = w_tmp[1]*xB_tmp[145]*2.0-w_tmp[1]*xB_tmp[146]+xB_tmp[1]*x_tmp[1]*x_tmp[145]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[145]*2.0;
  qBdot_tmp[ip+nplist*17] = w_tmp[1]*xB_tmp[154]*2.0-w_tmp[1]*xB_tmp[155]+xB_tmp[1]*x_tmp[1]*x_tmp[154]*4.0-xB_tmp[2]*x_tmp[1]*x_tmp[154]*2.0;

  } break;

  case 2: {
  qBdot_tmp[ip+nplist*0] = xB_tmp[2]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[2];
  qBdot_tmp[ip+nplist*1] = xB_tmp[2]*x_tmp[11]+xB_tmp[11]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[11]-k[0]*w_tmp[3]*xB_tmp[12]*x_tmp[2];
  qBdot_tmp[ip+nplist*2] = xB_tmp[2]*x_tmp[20]+xB_tmp[20]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[20]-k[0]*w_tmp[3]*xB_tmp[21]*x_tmp[2];
  qBdot_tmp[ip+nplist*3] = xB_tmp[2]*x_tmp[29]+xB_tmp[29]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[29]-k[0]*w_tmp[3]*xB_tmp[30]*x_tmp[2];
  qBdot_tmp[ip+nplist*4] = xB_tmp[2]*x_tmp[38]+xB_tmp[38]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[38]-k[0]*w_tmp[3]*xB_tmp[39]*x_tmp[2];
  qBdot_tmp[ip+nplist*5] = xB_tmp[2]*x_tmp[47]+xB_tmp[47]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[47]-k[0]*w_tmp[3]*xB_tmp[48]*x_tmp[2];
  qBdot_tmp[ip+nplist*6] = xB_tmp[2]*x_tmp[56]+xB_tmp[56]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[56]-k[0]*w_tmp[3]*xB_tmp[57]*x_tmp[2];
  qBdot_tmp[ip+nplist*7] = xB_tmp[2]*x_tmp[65]+xB_tmp[65]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[65]-k[0]*w_tmp[3]*xB_tmp[66]*x_tmp[2];
  qBdot_tmp[ip+nplist*8] = xB_tmp[2]*x_tmp[74]+xB_tmp[74]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[74]-k[0]*w_tmp[3]*xB_tmp[75]*x_tmp[2];
  qBdot_tmp[ip+nplist*9] = xB_tmp[2]*x_tmp[83]+xB_tmp[83]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[83]-k[0]*w_tmp[3]*xB_tmp[84]*x_tmp[2];
  qBdot_tmp[ip+nplist*10] = xB_tmp[2]*x_tmp[92]+xB_tmp[92]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[92]-k[0]*w_tmp[3]*xB_tmp[93]*x_tmp[2];
  qBdot_tmp[ip+nplist*11] = xB_tmp[2]*x_tmp[101]+xB_tmp[101]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[101]-k[0]*w_tmp[3]*xB_tmp[102]*x_tmp[2];
  qBdot_tmp[ip+nplist*12] = xB_tmp[2]*x_tmp[110]+xB_tmp[110]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[110]-k[0]*w_tmp[3]*xB_tmp[111]*x_tmp[2];
  qBdot_tmp[ip+nplist*13] = xB_tmp[2]*x_tmp[119]+xB_tmp[119]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[119]-k[0]*w_tmp[3]*xB_tmp[120]*x_tmp[2];
  qBdot_tmp[ip+nplist*14] = xB_tmp[2]*x_tmp[128]+xB_tmp[128]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[128]-k[0]*w_tmp[3]*xB_tmp[129]*x_tmp[2];
  qBdot_tmp[ip+nplist*15] = xB_tmp[2]*x_tmp[137]+xB_tmp[137]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[137]-k[0]*w_tmp[3]*xB_tmp[138]*x_tmp[2];
  qBdot_tmp[ip+nplist*16] = xB_tmp[2]*x_tmp[146]+xB_tmp[146]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[146]-k[0]*w_tmp[3]*xB_tmp[147]*x_tmp[2];
  qBdot_tmp[ip+nplist*17] = xB_tmp[2]*x_tmp[155]+xB_tmp[155]*x_tmp[2]-k[0]*w_tmp[3]*xB_tmp[3]*x_tmp[155]-k[0]*w_tmp[3]*xB_tmp[156]*x_tmp[2];

  } break;

  case 3: {
  qBdot_tmp[ip+nplist*0] = -xB_tmp[4]*(w_tmp[4]-x_tmp[4])-xB_tmp[5]*(x_tmp[4]-x_tmp[5])-xB_tmp[6]*(x_tmp[5]-x_tmp[6])-xB_tmp[7]*(x_tmp[6]-x_tmp[7])-xB_tmp[8]*(x_tmp[7]-x_tmp[8])+k[1]*w_tmp[3]*xB_tmp[3]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[8];
  qBdot_tmp[ip+nplist*1] = xB_tmp[3]*x_tmp[12]-xB_tmp[13]*(w_tmp[4]-x_tmp[4])-xB_tmp[14]*(x_tmp[4]-x_tmp[5])-xB_tmp[15]*(x_tmp[5]-x_tmp[6])-xB_tmp[16]*(x_tmp[6]-x_tmp[7])-xB_tmp[5]*(x_tmp[13]-x_tmp[14])-xB_tmp[17]*(x_tmp[7]-x_tmp[8])-xB_tmp[6]*(x_tmp[14]-x_tmp[15])-xB_tmp[7]*(x_tmp[15]-x_tmp[16])-xB_tmp[8]*(x_tmp[16]-x_tmp[17])-xB_tmp[4]*(x_tmp[12]*2.0-x_tmp[13])+k[1]*w_tmp[3]*xB_tmp[12]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[17]-k[1]*w_tmp[2]*xB_tmp[9]*x_tmp[8];
  qBdot_tmp[ip+nplist*2] = xB_tmp[3]*x_tmp[21]-xB_tmp[22]*(w_tmp[4]-x_tmp[4])-xB_tmp[23]*(x_tmp[4]-x_tmp[5])-xB_tmp[24]*(x_tmp[5]-x_tmp[6])-xB_tmp[25]*(x_tmp[6]-x_tmp[7])-xB_tmp[26]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[22]-x_tmp[23])-xB_tmp[6]*(x_tmp[23]-x_tmp[24])-xB_tmp[7]*(x_tmp[24]-x_tmp[25])-xB_tmp[8]*(x_tmp[25]-x_tmp[26])-xB_tmp[4]*(x_tmp[21]*2.0-x_tmp[22])+k[1]*w_tmp[3]*xB_tmp[21]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[26]-k[1]*w_tmp[2]*xB_tmp[18]*x_tmp[8];
  qBdot_tmp[ip+nplist*3] = xB_tmp[3]*x_tmp[30]-xB_tmp[31]*(w_tmp[4]-x_tmp[4])-xB_tmp[32]*(x_tmp[4]-x_tmp[5])-xB_tmp[33]*(x_tmp[5]-x_tmp[6])-xB_tmp[34]*(x_tmp[6]-x_tmp[7])-xB_tmp[35]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[31]-x_tmp[32])-xB_tmp[6]*(x_tmp[32]-x_tmp[33])-xB_tmp[7]*(x_tmp[33]-x_tmp[34])-xB_tmp[8]*(x_tmp[34]-x_tmp[35])-xB_tmp[4]*(x_tmp[30]*2.0-x_tmp[31])+k[1]*w_tmp[3]*xB_tmp[30]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[35]-k[1]*w_tmp[2]*xB_tmp[27]*x_tmp[8];
  qBdot_tmp[ip+nplist*4] = xB_tmp[3]*x_tmp[39]-xB_tmp[40]*(w_tmp[4]-x_tmp[4])-xB_tmp[41]*(x_tmp[4]-x_tmp[5])-xB_tmp[42]*(x_tmp[5]-x_tmp[6])-xB_tmp[43]*(x_tmp[6]-x_tmp[7])-xB_tmp[44]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[40]-x_tmp[41])-xB_tmp[6]*(x_tmp[41]-x_tmp[42])-xB_tmp[7]*(x_tmp[42]-x_tmp[43])-xB_tmp[8]*(x_tmp[43]-x_tmp[44])-xB_tmp[4]*(x_tmp[39]*2.0-x_tmp[40])+k[1]*w_tmp[3]*xB_tmp[39]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[44]-k[1]*w_tmp[2]*xB_tmp[36]*x_tmp[8];
  qBdot_tmp[ip+nplist*5] = xB_tmp[3]*x_tmp[48]-xB_tmp[49]*(w_tmp[4]-x_tmp[4])-xB_tmp[50]*(x_tmp[4]-x_tmp[5])-xB_tmp[51]*(x_tmp[5]-x_tmp[6])-xB_tmp[52]*(x_tmp[6]-x_tmp[7])-xB_tmp[53]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[49]-x_tmp[50])-xB_tmp[6]*(x_tmp[50]-x_tmp[51])-xB_tmp[7]*(x_tmp[51]-x_tmp[52])-xB_tmp[8]*(x_tmp[52]-x_tmp[53])-xB_tmp[4]*(x_tmp[48]*2.0-x_tmp[49])+k[1]*w_tmp[3]*xB_tmp[48]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[53]-k[1]*w_tmp[2]*xB_tmp[45]*x_tmp[8];
  qBdot_tmp[ip+nplist*6] = xB_tmp[3]*x_tmp[57]-xB_tmp[58]*(w_tmp[4]-x_tmp[4])-xB_tmp[59]*(x_tmp[4]-x_tmp[5])-xB_tmp[60]*(x_tmp[5]-x_tmp[6])-xB_tmp[61]*(x_tmp[6]-x_tmp[7])-xB_tmp[62]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[58]-x_tmp[59])-xB_tmp[6]*(x_tmp[59]-x_tmp[60])-xB_tmp[7]*(x_tmp[60]-x_tmp[61])-xB_tmp[8]*(x_tmp[61]-x_tmp[62])-xB_tmp[4]*(x_tmp[57]*2.0-x_tmp[58])+k[1]*w_tmp[3]*xB_tmp[57]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[62]-k[1]*w_tmp[2]*xB_tmp[54]*x_tmp[8];
  qBdot_tmp[ip+nplist*7] = xB_tmp[3]*x_tmp[66]-xB_tmp[67]*(w_tmp[4]-x_tmp[4])-xB_tmp[68]*(x_tmp[4]-x_tmp[5])-xB_tmp[69]*(x_tmp[5]-x_tmp[6])-xB_tmp[70]*(x_tmp[6]-x_tmp[7])-xB_tmp[71]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[67]-x_tmp[68])-xB_tmp[6]*(x_tmp[68]-x_tmp[69])-xB_tmp[7]*(x_tmp[69]-x_tmp[70])-xB_tmp[8]*(x_tmp[70]-x_tmp[71])-xB_tmp[4]*(x_tmp[66]*2.0-x_tmp[67])+k[1]*w_tmp[3]*xB_tmp[66]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[71]-k[1]*w_tmp[2]*xB_tmp[63]*x_tmp[8];
  qBdot_tmp[ip+nplist*8] = xB_tmp[3]*x_tmp[75]-xB_tmp[76]*(w_tmp[4]-x_tmp[4])-xB_tmp[77]*(x_tmp[4]-x_tmp[5])-xB_tmp[78]*(x_tmp[5]-x_tmp[6])-xB_tmp[79]*(x_tmp[6]-x_tmp[7])-xB_tmp[80]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[76]-x_tmp[77])-xB_tmp[6]*(x_tmp[77]-x_tmp[78])-xB_tmp[7]*(x_tmp[78]-x_tmp[79])-xB_tmp[8]*(x_tmp[79]-x_tmp[80])-xB_tmp[4]*(x_tmp[75]*2.0-x_tmp[76])+k[1]*w_tmp[3]*xB_tmp[75]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[80]-k[1]*w_tmp[2]*xB_tmp[72]*x_tmp[8];
  qBdot_tmp[ip+nplist*9] = xB_tmp[3]*x_tmp[84]-xB_tmp[85]*(w_tmp[4]-x_tmp[4])-xB_tmp[86]*(x_tmp[4]-x_tmp[5])-xB_tmp[87]*(x_tmp[5]-x_tmp[6])-xB_tmp[88]*(x_tmp[6]-x_tmp[7])-xB_tmp[89]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[85]-x_tmp[86])-xB_tmp[6]*(x_tmp[86]-x_tmp[87])-xB_tmp[7]*(x_tmp[87]-x_tmp[88])-xB_tmp[8]*(x_tmp[88]-x_tmp[89])-xB_tmp[4]*(x_tmp[84]*2.0-x_tmp[85])+k[1]*w_tmp[3]*xB_tmp[84]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[89]-k[1]*w_tmp[2]*xB_tmp[81]*x_tmp[8];
  qBdot_tmp[ip+nplist*10] = xB_tmp[3]*x_tmp[93]-xB_tmp[94]*(w_tmp[4]-x_tmp[4])-xB_tmp[95]*(x_tmp[4]-x_tmp[5])-xB_tmp[96]*(x_tmp[5]-x_tmp[6])-xB_tmp[97]*(x_tmp[6]-x_tmp[7])-xB_tmp[98]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[94]-x_tmp[95])-xB_tmp[6]*(x_tmp[95]-x_tmp[96])-xB_tmp[7]*(x_tmp[96]-x_tmp[97])-xB_tmp[8]*(x_tmp[97]-x_tmp[98])-xB_tmp[4]*(x_tmp[93]*2.0-x_tmp[94])+k[1]*w_tmp[3]*xB_tmp[93]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[98]-k[1]*w_tmp[2]*xB_tmp[90]*x_tmp[8];
  qBdot_tmp[ip+nplist*11] = xB_tmp[3]*x_tmp[102]-xB_tmp[103]*(w_tmp[4]-x_tmp[4])-xB_tmp[104]*(x_tmp[4]-x_tmp[5])-xB_tmp[105]*(x_tmp[5]-x_tmp[6])-xB_tmp[106]*(x_tmp[6]-x_tmp[7])-xB_tmp[107]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[103]-x_tmp[104])-xB_tmp[6]*(x_tmp[104]-x_tmp[105])-xB_tmp[7]*(x_tmp[105]-x_tmp[106])-xB_tmp[8]*(x_tmp[106]-x_tmp[107])-xB_tmp[4]*(x_tmp[102]*2.0-x_tmp[103])+k[1]*w_tmp[3]*xB_tmp[102]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[107]-k[1]*w_tmp[2]*xB_tmp[99]*x_tmp[8];
  qBdot_tmp[ip+nplist*12] = xB_tmp[3]*x_tmp[111]-xB_tmp[112]*(w_tmp[4]-x_tmp[4])-xB_tmp[113]*(x_tmp[4]-x_tmp[5])-xB_tmp[114]*(x_tmp[5]-x_tmp[6])-xB_tmp[115]*(x_tmp[6]-x_tmp[7])-xB_tmp[116]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[112]-x_tmp[113])-xB_tmp[6]*(x_tmp[113]-x_tmp[114])-xB_tmp[7]*(x_tmp[114]-x_tmp[115])-xB_tmp[8]*(x_tmp[115]-x_tmp[116])-xB_tmp[4]*(x_tmp[111]*2.0-x_tmp[112])+k[1]*w_tmp[3]*xB_tmp[111]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[116]-k[1]*w_tmp[2]*xB_tmp[108]*x_tmp[8];
  qBdot_tmp[ip+nplist*13] = xB_tmp[3]*x_tmp[120]-xB_tmp[121]*(w_tmp[4]-x_tmp[4])-xB_tmp[122]*(x_tmp[4]-x_tmp[5])-xB_tmp[123]*(x_tmp[5]-x_tmp[6])-xB_tmp[124]*(x_tmp[6]-x_tmp[7])-xB_tmp[125]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[121]-x_tmp[122])-xB_tmp[6]*(x_tmp[122]-x_tmp[123])-xB_tmp[7]*(x_tmp[123]-x_tmp[124])-xB_tmp[8]*(x_tmp[124]-x_tmp[125])-xB_tmp[4]*(x_tmp[120]*2.0-x_tmp[121])+k[1]*w_tmp[3]*xB_tmp[120]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[125]-k[1]*w_tmp[2]*xB_tmp[117]*x_tmp[8];
  qBdot_tmp[ip+nplist*14] = xB_tmp[3]*x_tmp[129]-xB_tmp[130]*(w_tmp[4]-x_tmp[4])-xB_tmp[131]*(x_tmp[4]-x_tmp[5])-xB_tmp[132]*(x_tmp[5]-x_tmp[6])-xB_tmp[133]*(x_tmp[6]-x_tmp[7])-xB_tmp[134]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[130]-x_tmp[131])-xB_tmp[6]*(x_tmp[131]-x_tmp[132])-xB_tmp[7]*(x_tmp[132]-x_tmp[133])-xB_tmp[8]*(x_tmp[133]-x_tmp[134])-xB_tmp[4]*(x_tmp[129]*2.0-x_tmp[130])+k[1]*w_tmp[3]*xB_tmp[129]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[134]-k[1]*w_tmp[2]*xB_tmp[126]*x_tmp[8];
  qBdot_tmp[ip+nplist*15] = xB_tmp[3]*x_tmp[138]-xB_tmp[139]*(w_tmp[4]-x_tmp[4])-xB_tmp[140]*(x_tmp[4]-x_tmp[5])-xB_tmp[141]*(x_tmp[5]-x_tmp[6])-xB_tmp[142]*(x_tmp[6]-x_tmp[7])-xB_tmp[143]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[139]-x_tmp[140])-xB_tmp[6]*(x_tmp[140]-x_tmp[141])-xB_tmp[7]*(x_tmp[141]-x_tmp[142])-xB_tmp[8]*(x_tmp[142]-x_tmp[143])-xB_tmp[4]*(x_tmp[138]*2.0-x_tmp[139])+k[1]*w_tmp[3]*xB_tmp[138]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[143]-k[1]*w_tmp[2]*xB_tmp[135]*x_tmp[8];
  qBdot_tmp[ip+nplist*16] = xB_tmp[3]*x_tmp[147]-xB_tmp[148]*(w_tmp[4]-x_tmp[4])-xB_tmp[149]*(x_tmp[4]-x_tmp[5])-xB_tmp[150]*(x_tmp[5]-x_tmp[6])-xB_tmp[151]*(x_tmp[6]-x_tmp[7])-xB_tmp[152]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[148]-x_tmp[149])-xB_tmp[6]*(x_tmp[149]-x_tmp[150])-xB_tmp[7]*(x_tmp[150]-x_tmp[151])-xB_tmp[8]*(x_tmp[151]-x_tmp[152])-xB_tmp[4]*(x_tmp[147]*2.0-x_tmp[148])+k[1]*w_tmp[3]*xB_tmp[147]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[152]-k[1]*w_tmp[2]*xB_tmp[144]*x_tmp[8];
  qBdot_tmp[ip+nplist*17] = xB_tmp[3]*x_tmp[156]-xB_tmp[157]*(w_tmp[4]-x_tmp[4])-xB_tmp[158]*(x_tmp[4]-x_tmp[5])-xB_tmp[159]*(x_tmp[5]-x_tmp[6])-xB_tmp[160]*(x_tmp[6]-x_tmp[7])-xB_tmp[161]*(x_tmp[7]-x_tmp[8])-xB_tmp[5]*(x_tmp[157]-x_tmp[158])-xB_tmp[6]*(x_tmp[158]-x_tmp[159])-xB_tmp[7]*(x_tmp[159]-x_tmp[160])-xB_tmp[8]*(x_tmp[160]-x_tmp[161])-xB_tmp[4]*(x_tmp[156]*2.0-x_tmp[157])+k[1]*w_tmp[3]*xB_tmp[156]*x_tmp[3]-k[1]*w_tmp[2]*xB_tmp[0]*x_tmp[161]-k[1]*w_tmp[2]*xB_tmp[153]*x_tmp[8];

  } break;

  case 5: {
  qBdot_tmp[ip+nplist*0] = -dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[0]*x_tmp[0];
  qBdot_tmp[ip+nplist*1] = dwdp_tmp[0]*xB_tmp[0]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[0]*xB_tmp[1]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[0]*p[0]*xB_tmp[10]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[9]*x_tmp[0];
  qBdot_tmp[ip+nplist*2] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[18]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[18]-dwdp_tmp[0]*p[0]*xB_tmp[19]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[18]*x_tmp[0];
  qBdot_tmp[ip+nplist*3] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[27]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[27]-dwdp_tmp[0]*p[0]*xB_tmp[28]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[27]*x_tmp[0];
  qBdot_tmp[ip+nplist*4] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[36]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[36]-dwdp_tmp[0]*p[0]*xB_tmp[37]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[36]*x_tmp[0];
  qBdot_tmp[ip+nplist*5] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[45]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[45]-dwdp_tmp[0]*p[0]*xB_tmp[46]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[45]*x_tmp[0];
  qBdot_tmp[ip+nplist*6] = xB_tmp[0]*(dwdp_tmp[1]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[54])-xB_tmp[1]*(dwdp_tmp[1]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[54])-dwdp_tmp[0]*p[0]*xB_tmp[55]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[54]*x_tmp[0];
  qBdot_tmp[ip+nplist*7] = xB_tmp[0]*(dwdp_tmp[2]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[63])-xB_tmp[1]*(dwdp_tmp[2]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[63])-dwdp_tmp[0]*p[0]*xB_tmp[64]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[63]*x_tmp[0];
  qBdot_tmp[ip+nplist*8] = xB_tmp[0]*(dwdp_tmp[3]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[72])-xB_tmp[1]*(dwdp_tmp[3]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[72])-dwdp_tmp[0]*p[0]*xB_tmp[73]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[72]*x_tmp[0];
  qBdot_tmp[ip+nplist*9] = xB_tmp[0]*(dwdp_tmp[4]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[81])-xB_tmp[1]*(dwdp_tmp[4]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[81])-dwdp_tmp[0]*p[0]*xB_tmp[82]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[81]*x_tmp[0];
  qBdot_tmp[ip+nplist*10] = xB_tmp[0]*(dwdp_tmp[5]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[90])-xB_tmp[1]*(dwdp_tmp[5]*p[0]*x_tmp[0]+dwdp_tmp[0]*p[0]*x_tmp[90])-dwdp_tmp[0]*p[0]*xB_tmp[91]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[90]*x_tmp[0];
  qBdot_tmp[ip+nplist*11] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[99]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[99]-dwdp_tmp[0]*p[0]*xB_tmp[100]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[99]*x_tmp[0];
  qBdot_tmp[ip+nplist*12] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[108]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[108]-dwdp_tmp[0]*p[0]*xB_tmp[109]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[108]*x_tmp[0];
  qBdot_tmp[ip+nplist*13] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[117]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[117]-dwdp_tmp[0]*p[0]*xB_tmp[118]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[117]*x_tmp[0];
  qBdot_tmp[ip+nplist*14] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[126]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[126]-dwdp_tmp[0]*p[0]*xB_tmp[127]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[126]*x_tmp[0];
  qBdot_tmp[ip+nplist*15] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[135]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[135]-dwdp_tmp[0]*p[0]*xB_tmp[136]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[135]*x_tmp[0];
  qBdot_tmp[ip+nplist*16] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[144]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[144]-dwdp_tmp[0]*p[0]*xB_tmp[145]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[144]*x_tmp[0];
  qBdot_tmp[ip+nplist*17] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[153]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[153]-dwdp_tmp[0]*p[0]*xB_tmp[154]*x_tmp[0]+dwdp_tmp[0]*k[0]*p[0]*w_tmp[2]*xB_tmp[153]*x_tmp[0];

  } break;

  case 6: {
  qBdot_tmp[ip+nplist*0] = -dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[0]*x_tmp[0];
  qBdot_tmp[ip+nplist*1] = dwdp_tmp[6]*xB_tmp[0]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[6]*xB_tmp[1]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[6]*p[0]*xB_tmp[10]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[9]*x_tmp[0];
  qBdot_tmp[ip+nplist*2] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[18]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[18]-dwdp_tmp[6]*p[0]*xB_tmp[19]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[18]*x_tmp[0];
  qBdot_tmp[ip+nplist*3] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[27]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[27]-dwdp_tmp[6]*p[0]*xB_tmp[28]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[27]*x_tmp[0];
  qBdot_tmp[ip+nplist*4] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[36]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[36]-dwdp_tmp[6]*p[0]*xB_tmp[37]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[36]*x_tmp[0];
  qBdot_tmp[ip+nplist*5] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[45]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[45]-dwdp_tmp[6]*p[0]*xB_tmp[46]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[45]*x_tmp[0];
  qBdot_tmp[ip+nplist*6] = xB_tmp[0]*(dwdp_tmp[7]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[54])-xB_tmp[1]*(dwdp_tmp[7]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[54])-dwdp_tmp[6]*p[0]*xB_tmp[55]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[54]*x_tmp[0];
  qBdot_tmp[ip+nplist*7] = xB_tmp[0]*(dwdp_tmp[8]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[63])-xB_tmp[1]*(dwdp_tmp[8]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[63])-dwdp_tmp[6]*p[0]*xB_tmp[64]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[63]*x_tmp[0];
  qBdot_tmp[ip+nplist*8] = xB_tmp[0]*(dwdp_tmp[9]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[72])-xB_tmp[1]*(dwdp_tmp[9]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[72])-dwdp_tmp[6]*p[0]*xB_tmp[73]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[72]*x_tmp[0];
  qBdot_tmp[ip+nplist*9] = xB_tmp[0]*(dwdp_tmp[10]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[81])-xB_tmp[1]*(dwdp_tmp[10]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[81])-dwdp_tmp[6]*p[0]*xB_tmp[82]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[81]*x_tmp[0];
  qBdot_tmp[ip+nplist*10] = xB_tmp[0]*(dwdp_tmp[11]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[90])-xB_tmp[1]*(dwdp_tmp[11]*p[0]*x_tmp[0]+dwdp_tmp[6]*p[0]*x_tmp[90])-dwdp_tmp[6]*p[0]*xB_tmp[91]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[90]*x_tmp[0];
  qBdot_tmp[ip+nplist*11] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[99]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[99]-dwdp_tmp[6]*p[0]*xB_tmp[100]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[99]*x_tmp[0];
  qBdot_tmp[ip+nplist*12] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[108]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[108]-dwdp_tmp[6]*p[0]*xB_tmp[109]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[108]*x_tmp[0];
  qBdot_tmp[ip+nplist*13] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[117]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[117]-dwdp_tmp[6]*p[0]*xB_tmp[118]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[117]*x_tmp[0];
  qBdot_tmp[ip+nplist*14] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[126]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[126]-dwdp_tmp[6]*p[0]*xB_tmp[127]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[126]*x_tmp[0];
  qBdot_tmp[ip+nplist*15] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[135]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[135]-dwdp_tmp[6]*p[0]*xB_tmp[136]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[135]*x_tmp[0];
  qBdot_tmp[ip+nplist*16] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[144]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[144]-dwdp_tmp[6]*p[0]*xB_tmp[145]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[144]*x_tmp[0];
  qBdot_tmp[ip+nplist*17] = dwdp_tmp[6]*p[0]*xB_tmp[0]*x_tmp[153]-dwdp_tmp[6]*p[0]*xB_tmp[1]*x_tmp[153]-dwdp_tmp[6]*p[0]*xB_tmp[154]*x_tmp[0]+dwdp_tmp[6]*k[0]*p[0]*w_tmp[2]*xB_tmp[153]*x_tmp[0];

  } break;

  case 7: {
  qBdot_tmp[ip+nplist*0] = -dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[0]*x_tmp[0];
  qBdot_tmp[ip+nplist*1] = dwdp_tmp[12]*xB_tmp[0]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[12]*xB_tmp[1]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[12]*p[0]*xB_tmp[10]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[9]*x_tmp[0];
  qBdot_tmp[ip+nplist*2] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[18]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[18]-dwdp_tmp[12]*p[0]*xB_tmp[19]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[18]*x_tmp[0];
  qBdot_tmp[ip+nplist*3] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[27]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[27]-dwdp_tmp[12]*p[0]*xB_tmp[28]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[27]*x_tmp[0];
  qBdot_tmp[ip+nplist*4] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[36]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[36]-dwdp_tmp[12]*p[0]*xB_tmp[37]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[36]*x_tmp[0];
  qBdot_tmp[ip+nplist*5] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[45]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[45]-dwdp_tmp[12]*p[0]*xB_tmp[46]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[45]*x_tmp[0];
  qBdot_tmp[ip+nplist*6] = xB_tmp[0]*(dwdp_tmp[13]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[54])-xB_tmp[1]*(dwdp_tmp[13]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[54])-dwdp_tmp[12]*p[0]*xB_tmp[55]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[54]*x_tmp[0];
  qBdot_tmp[ip+nplist*7] = xB_tmp[0]*(dwdp_tmp[14]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[63])-xB_tmp[1]*(dwdp_tmp[14]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[63])-dwdp_tmp[12]*p[0]*xB_tmp[64]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[63]*x_tmp[0];
  qBdot_tmp[ip+nplist*8] = xB_tmp[0]*(dwdp_tmp[15]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[72])-xB_tmp[1]*(dwdp_tmp[15]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[72])-dwdp_tmp[12]*p[0]*xB_tmp[73]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[72]*x_tmp[0];
  qBdot_tmp[ip+nplist*9] = xB_tmp[0]*(dwdp_tmp[16]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[81])-xB_tmp[1]*(dwdp_tmp[16]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[81])-dwdp_tmp[12]*p[0]*xB_tmp[82]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[81]*x_tmp[0];
  qBdot_tmp[ip+nplist*10] = xB_tmp[0]*(dwdp_tmp[17]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[90])-xB_tmp[1]*(dwdp_tmp[17]*p[0]*x_tmp[0]+dwdp_tmp[12]*p[0]*x_tmp[90])-dwdp_tmp[12]*p[0]*xB_tmp[91]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[90]*x_tmp[0];
  qBdot_tmp[ip+nplist*11] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[99]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[99]-dwdp_tmp[12]*p[0]*xB_tmp[100]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[99]*x_tmp[0];
  qBdot_tmp[ip+nplist*12] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[108]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[108]-dwdp_tmp[12]*p[0]*xB_tmp[109]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[108]*x_tmp[0];
  qBdot_tmp[ip+nplist*13] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[117]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[117]-dwdp_tmp[12]*p[0]*xB_tmp[118]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[117]*x_tmp[0];
  qBdot_tmp[ip+nplist*14] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[126]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[126]-dwdp_tmp[12]*p[0]*xB_tmp[127]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[126]*x_tmp[0];
  qBdot_tmp[ip+nplist*15] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[135]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[135]-dwdp_tmp[12]*p[0]*xB_tmp[136]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[135]*x_tmp[0];
  qBdot_tmp[ip+nplist*16] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[144]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[144]-dwdp_tmp[12]*p[0]*xB_tmp[145]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[144]*x_tmp[0];
  qBdot_tmp[ip+nplist*17] = dwdp_tmp[12]*p[0]*xB_tmp[0]*x_tmp[153]-dwdp_tmp[12]*p[0]*xB_tmp[1]*x_tmp[153]-dwdp_tmp[12]*p[0]*xB_tmp[154]*x_tmp[0]+dwdp_tmp[12]*k[0]*p[0]*w_tmp[2]*xB_tmp[153]*x_tmp[0];

  } break;

  case 8: {
  qBdot_tmp[ip+nplist*0] = -dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[0]*x_tmp[0];
  qBdot_tmp[ip+nplist*1] = dwdp_tmp[18]*xB_tmp[0]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[18]*xB_tmp[1]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[18]*p[0]*xB_tmp[10]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[9]*x_tmp[0];
  qBdot_tmp[ip+nplist*2] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[18]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[18]-dwdp_tmp[18]*p[0]*xB_tmp[19]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[18]*x_tmp[0];
  qBdot_tmp[ip+nplist*3] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[27]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[27]-dwdp_tmp[18]*p[0]*xB_tmp[28]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[27]*x_tmp[0];
  qBdot_tmp[ip+nplist*4] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[36]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[36]-dwdp_tmp[18]*p[0]*xB_tmp[37]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[36]*x_tmp[0];
  qBdot_tmp[ip+nplist*5] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[45]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[45]-dwdp_tmp[18]*p[0]*xB_tmp[46]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[45]*x_tmp[0];
  qBdot_tmp[ip+nplist*6] = xB_tmp[0]*(dwdp_tmp[19]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[54])-xB_tmp[1]*(dwdp_tmp[19]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[54])-dwdp_tmp[18]*p[0]*xB_tmp[55]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[54]*x_tmp[0];
  qBdot_tmp[ip+nplist*7] = xB_tmp[0]*(dwdp_tmp[20]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[63])-xB_tmp[1]*(dwdp_tmp[20]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[63])-dwdp_tmp[18]*p[0]*xB_tmp[64]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[63]*x_tmp[0];
  qBdot_tmp[ip+nplist*8] = xB_tmp[0]*(dwdp_tmp[21]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[72])-xB_tmp[1]*(dwdp_tmp[21]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[72])-dwdp_tmp[18]*p[0]*xB_tmp[73]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[72]*x_tmp[0];
  qBdot_tmp[ip+nplist*9] = xB_tmp[0]*(dwdp_tmp[22]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[81])-xB_tmp[1]*(dwdp_tmp[22]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[81])-dwdp_tmp[18]*p[0]*xB_tmp[82]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[81]*x_tmp[0];
  qBdot_tmp[ip+nplist*10] = xB_tmp[0]*(dwdp_tmp[23]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[90])-xB_tmp[1]*(dwdp_tmp[23]*p[0]*x_tmp[0]+dwdp_tmp[18]*p[0]*x_tmp[90])-dwdp_tmp[18]*p[0]*xB_tmp[91]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[90]*x_tmp[0];
  qBdot_tmp[ip+nplist*11] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[99]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[99]-dwdp_tmp[18]*p[0]*xB_tmp[100]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[99]*x_tmp[0];
  qBdot_tmp[ip+nplist*12] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[108]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[108]-dwdp_tmp[18]*p[0]*xB_tmp[109]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[108]*x_tmp[0];
  qBdot_tmp[ip+nplist*13] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[117]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[117]-dwdp_tmp[18]*p[0]*xB_tmp[118]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[117]*x_tmp[0];
  qBdot_tmp[ip+nplist*14] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[126]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[126]-dwdp_tmp[18]*p[0]*xB_tmp[127]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[126]*x_tmp[0];
  qBdot_tmp[ip+nplist*15] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[135]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[135]-dwdp_tmp[18]*p[0]*xB_tmp[136]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[135]*x_tmp[0];
  qBdot_tmp[ip+nplist*16] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[144]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[144]-dwdp_tmp[18]*p[0]*xB_tmp[145]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[144]*x_tmp[0];
  qBdot_tmp[ip+nplist*17] = dwdp_tmp[18]*p[0]*xB_tmp[0]*x_tmp[153]-dwdp_tmp[18]*p[0]*xB_tmp[1]*x_tmp[153]-dwdp_tmp[18]*p[0]*xB_tmp[154]*x_tmp[0]+dwdp_tmp[18]*k[0]*p[0]*w_tmp[2]*xB_tmp[153]*x_tmp[0];

  } break;

  case 9: {
  qBdot_tmp[ip+nplist*0] = -dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[0]*x_tmp[0];
  qBdot_tmp[ip+nplist*1] = dwdp_tmp[24]*xB_tmp[0]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[24]*xB_tmp[1]*(x_tmp[0]+p[0]*x_tmp[9])-dwdp_tmp[24]*p[0]*xB_tmp[10]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[9]*x_tmp[0];
  qBdot_tmp[ip+nplist*2] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[18]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[18]-dwdp_tmp[24]*p[0]*xB_tmp[19]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[18]*x_tmp[0];
  qBdot_tmp[ip+nplist*3] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[27]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[27]-dwdp_tmp[24]*p[0]*xB_tmp[28]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[27]*x_tmp[0];
  qBdot_tmp[ip+nplist*4] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[36]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[36]-dwdp_tmp[24]*p[0]*xB_tmp[37]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[36]*x_tmp[0];
  qBdot_tmp[ip+nplist*5] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[45]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[45]-dwdp_tmp[24]*p[0]*xB_tmp[46]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[45]*x_tmp[0];
  qBdot_tmp[ip+nplist*6] = xB_tmp[0]*(dwdp_tmp[25]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[54])-xB_tmp[1]*(dwdp_tmp[25]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[54])-dwdp_tmp[24]*p[0]*xB_tmp[55]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[54]*x_tmp[0];
  qBdot_tmp[ip+nplist*7] = xB_tmp[0]*(dwdp_tmp[26]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[63])-xB_tmp[1]*(dwdp_tmp[26]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[63])-dwdp_tmp[24]*p[0]*xB_tmp[64]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[63]*x_tmp[0];
  qBdot_tmp[ip+nplist*8] = xB_tmp[0]*(dwdp_tmp[27]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[72])-xB_tmp[1]*(dwdp_tmp[27]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[72])-dwdp_tmp[24]*p[0]*xB_tmp[73]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[72]*x_tmp[0];
  qBdot_tmp[ip+nplist*9] = xB_tmp[0]*(dwdp_tmp[28]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[81])-xB_tmp[1]*(dwdp_tmp[28]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[81])-dwdp_tmp[24]*p[0]*xB_tmp[82]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[81]*x_tmp[0];
  qBdot_tmp[ip+nplist*10] = xB_tmp[0]*(dwdp_tmp[29]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[90])-xB_tmp[1]*(dwdp_tmp[29]*p[0]*x_tmp[0]+dwdp_tmp[24]*p[0]*x_tmp[90])-dwdp_tmp[24]*p[0]*xB_tmp[91]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[90]*x_tmp[0];
  qBdot_tmp[ip+nplist*11] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[99]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[99]-dwdp_tmp[24]*p[0]*xB_tmp[100]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[99]*x_tmp[0];
  qBdot_tmp[ip+nplist*12] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[108]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[108]-dwdp_tmp[24]*p[0]*xB_tmp[109]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[108]*x_tmp[0];
  qBdot_tmp[ip+nplist*13] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[117]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[117]-dwdp_tmp[24]*p[0]*xB_tmp[118]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[117]*x_tmp[0];
  qBdot_tmp[ip+nplist*14] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[126]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[126]-dwdp_tmp[24]*p[0]*xB_tmp[127]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[126]*x_tmp[0];
  qBdot_tmp[ip+nplist*15] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[135]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[135]-dwdp_tmp[24]*p[0]*xB_tmp[136]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[135]*x_tmp[0];
  qBdot_tmp[ip+nplist*16] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[144]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[144]-dwdp_tmp[24]*p[0]*xB_tmp[145]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[144]*x_tmp[0];
  qBdot_tmp[ip+nplist*17] = dwdp_tmp[24]*p[0]*xB_tmp[0]*x_tmp[153]-dwdp_tmp[24]*p[0]*xB_tmp[1]*x_tmp[153]-dwdp_tmp[24]*p[0]*xB_tmp[154]*x_tmp[0]+dwdp_tmp[24]*k[0]*p[0]*w_tmp[2]*xB_tmp[153]*x_tmp[0];

  } break;

}
}
for(ip = 0; ip<nplist*ng; ip++) {
   if(amiIsNaN(qBdot_tmp[ip])) {
       qBdot_tmp[ip] = 0;       if(!udata->am_nan_qBdot) {
           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_qBdot = TRUE;
       }
   }   if(amiIsInf(qBdot_tmp[ip])) {
       warnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


