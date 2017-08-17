
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_dwdx.h"
#include "model_jakstat_adjoint_o2_w.h"

int xBdot_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
int ix;
memset(xBdot_tmp,0,sizeof(realtype)*162);
status = w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
status = dwdx_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  xBdot_tmp[0] = -udata->p[0]*tdata->w[0]*xB_tmp[1]+udata->k[0]*udata->p[0]*tdata->w[0]*tdata->w[2]*xB_tmp[0];
  xBdot_tmp[1] = udata->p[1]*xB_tmp[1]*tdata->dwdx[0]*2.0-udata->p[1]*xB_tmp[2]*tdata->dwdx[0];
  xBdot_tmp[2] = udata->p[2]*xB_tmp[2]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[3];
  xBdot_tmp[3] = -udata->p[3]*xB_tmp[4]*tdata->dwdx[1]+udata->k[1]*udata->p[3]*tdata->w[3]*xB_tmp[3];
  xBdot_tmp[4] = udata->p[3]*xB_tmp[4]-udata->p[3]*xB_tmp[5];
  xBdot_tmp[5] = udata->p[3]*xB_tmp[5]-udata->p[3]*xB_tmp[6];
  xBdot_tmp[6] = udata->p[3]*xB_tmp[6]-udata->p[3]*xB_tmp[7];
  xBdot_tmp[7] = udata->p[3]*xB_tmp[7]-udata->p[3]*xB_tmp[8];
  xBdot_tmp[8] = udata->p[3]*xB_tmp[8]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[0];
  xBdot_tmp[9] = tdata->w[0]*xB_tmp[0]-tdata->w[0]*xB_tmp[1]+udata->p[0]*tdata->w[0]*xB_tmp[9]-udata->p[0]*tdata->w[0]*xB_tmp[10];
  xBdot_tmp[10] = udata->p[1]*x_tmp[1]*xB_tmp[10]*4.0+udata->p[1]*x_tmp[10]*xB_tmp[1]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[11]*2.0-udata->p[1]*x_tmp[10]*xB_tmp[2]*2.0;
  xBdot_tmp[11] = udata->p[2]*xB_tmp[11]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[12];
  xBdot_tmp[12] = udata->p[3]*xB_tmp[12]-udata->p[3]*xB_tmp[13]*2.0;
  xBdot_tmp[13] = udata->p[3]*xB_tmp[13]-udata->p[3]*xB_tmp[14];
  xBdot_tmp[14] = udata->p[3]*xB_tmp[14]-udata->p[3]*xB_tmp[15];
  xBdot_tmp[15] = udata->p[3]*xB_tmp[15]-udata->p[3]*xB_tmp[16];
  xBdot_tmp[16] = udata->p[3]*xB_tmp[16]-udata->p[3]*xB_tmp[17];
  xBdot_tmp[17] = udata->p[3]*xB_tmp[17]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[9];
  xBdot_tmp[18] = udata->p[0]*tdata->w[0]*xB_tmp[18]-udata->p[0]*tdata->w[0]*xB_tmp[19];
  xBdot_tmp[19] = xB_tmp[1]*(tdata->dwdx[0]*2.0+udata->p[1]*x_tmp[19]*4.0)-xB_tmp[2]*(tdata->dwdx[0]+udata->p[1]*x_tmp[19]*2.0)+udata->p[1]*x_tmp[1]*xB_tmp[19]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[20]*2.0;
  xBdot_tmp[20] = udata->p[2]*xB_tmp[20]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[21];
  xBdot_tmp[21] = udata->p[3]*xB_tmp[21]-udata->p[3]*xB_tmp[22]*2.0;
  xBdot_tmp[22] = udata->p[3]*xB_tmp[22]-udata->p[3]*xB_tmp[23];
  xBdot_tmp[23] = udata->p[3]*xB_tmp[23]-udata->p[3]*xB_tmp[24];
  xBdot_tmp[24] = udata->p[3]*xB_tmp[24]-udata->p[3]*xB_tmp[25];
  xBdot_tmp[25] = udata->p[3]*xB_tmp[25]-udata->p[3]*xB_tmp[26];
  xBdot_tmp[26] = udata->p[3]*xB_tmp[26]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[18];
  xBdot_tmp[27] = udata->p[0]*tdata->w[0]*xB_tmp[27]-udata->p[0]*tdata->w[0]*xB_tmp[28];
  xBdot_tmp[28] = udata->p[1]*x_tmp[1]*xB_tmp[28]*4.0+udata->p[1]*x_tmp[28]*xB_tmp[1]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[29]*2.0-udata->p[1]*x_tmp[28]*xB_tmp[2]*2.0;
  xBdot_tmp[29] = xB_tmp[2]+udata->p[2]*xB_tmp[29]-udata->k[0]*tdata->w[3]*xB_tmp[3]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[30];
  xBdot_tmp[30] = udata->p[3]*xB_tmp[30]-udata->p[3]*xB_tmp[31]*2.0;
  xBdot_tmp[31] = udata->p[3]*xB_tmp[31]-udata->p[3]*xB_tmp[32];
  xBdot_tmp[32] = udata->p[3]*xB_tmp[32]-udata->p[3]*xB_tmp[33];
  xBdot_tmp[33] = udata->p[3]*xB_tmp[33]-udata->p[3]*xB_tmp[34];
  xBdot_tmp[34] = udata->p[3]*xB_tmp[34]-udata->p[3]*xB_tmp[35];
  xBdot_tmp[35] = udata->p[3]*xB_tmp[35]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[27];
  xBdot_tmp[36] = udata->p[0]*tdata->w[0]*xB_tmp[36]-udata->p[0]*tdata->w[0]*xB_tmp[37];
  xBdot_tmp[37] = udata->p[1]*x_tmp[1]*xB_tmp[37]*4.0+udata->p[1]*x_tmp[37]*xB_tmp[1]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[38]*2.0-udata->p[1]*x_tmp[37]*xB_tmp[2]*2.0;
  xBdot_tmp[38] = udata->p[2]*xB_tmp[38]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[39];
  xBdot_tmp[39] = xB_tmp[3]+udata->p[3]*xB_tmp[39]-udata->p[3]*xB_tmp[40]*2.0-xB_tmp[4]*tdata->dwdx[1];
  xBdot_tmp[40] = xB_tmp[4]-xB_tmp[5]+udata->p[3]*xB_tmp[40]-udata->p[3]*xB_tmp[41];
  xBdot_tmp[41] = xB_tmp[5]-xB_tmp[6]+udata->p[3]*xB_tmp[41]-udata->p[3]*xB_tmp[42];
  xBdot_tmp[42] = xB_tmp[6]-xB_tmp[7]+udata->p[3]*xB_tmp[42]-udata->p[3]*xB_tmp[43];
  xBdot_tmp[43] = xB_tmp[7]-xB_tmp[8]+udata->p[3]*xB_tmp[43]-udata->p[3]*xB_tmp[44];
  xBdot_tmp[44] = xB_tmp[8]+udata->p[3]*xB_tmp[44]-udata->k[1]*tdata->w[2]*xB_tmp[0]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[36];
  xBdot_tmp[45] = udata->p[0]*tdata->w[0]*xB_tmp[45]-udata->p[0]*tdata->w[0]*xB_tmp[46];
  xBdot_tmp[46] = udata->p[1]*x_tmp[1]*xB_tmp[46]*4.0+udata->p[1]*x_tmp[46]*xB_tmp[1]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[47]*2.0-udata->p[1]*x_tmp[46]*xB_tmp[2]*2.0;
  xBdot_tmp[47] = udata->p[2]*xB_tmp[47]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[48];
  xBdot_tmp[48] = udata->p[3]*xB_tmp[48]-udata->p[3]*xB_tmp[49]*2.0;
  xBdot_tmp[49] = udata->p[3]*xB_tmp[49]-udata->p[3]*xB_tmp[50];
  xBdot_tmp[50] = udata->p[3]*xB_tmp[50]-udata->p[3]*xB_tmp[51];
  xBdot_tmp[51] = udata->p[3]*xB_tmp[51]-udata->p[3]*xB_tmp[52];
  xBdot_tmp[52] = udata->p[3]*xB_tmp[52]-udata->p[3]*xB_tmp[53];
  xBdot_tmp[53] = udata->p[3]*xB_tmp[53]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[45];
  xBdot_tmp[54] = udata->p[0]*tdata->w[5]*xB_tmp[0]-udata->p[0]*tdata->w[5]*xB_tmp[1]+udata->p[0]*tdata->w[0]*xB_tmp[54]-udata->p[0]*tdata->w[0]*xB_tmp[55];
  xBdot_tmp[55] = udata->p[1]*x_tmp[1]*xB_tmp[55]*4.0+udata->p[1]*x_tmp[55]*xB_tmp[1]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[56]*2.0-udata->p[1]*x_tmp[55]*xB_tmp[2]*2.0;
  xBdot_tmp[56] = udata->p[2]*xB_tmp[56]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[57];
  xBdot_tmp[57] = udata->p[3]*xB_tmp[57]-udata->p[3]*xB_tmp[58]*2.0;
  xBdot_tmp[58] = udata->p[3]*xB_tmp[58]-udata->p[3]*xB_tmp[59];
  xBdot_tmp[59] = udata->p[3]*xB_tmp[59]-udata->p[3]*xB_tmp[60];
  xBdot_tmp[60] = udata->p[3]*xB_tmp[60]-udata->p[3]*xB_tmp[61];
  xBdot_tmp[61] = udata->p[3]*xB_tmp[61]-udata->p[3]*xB_tmp[62];
  xBdot_tmp[62] = udata->p[3]*xB_tmp[62]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[54];
  xBdot_tmp[63] = udata->p[0]*tdata->w[6]*xB_tmp[0]-udata->p[0]*tdata->w[6]*xB_tmp[1]+udata->p[0]*tdata->w[0]*xB_tmp[63]-udata->p[0]*tdata->w[0]*xB_tmp[64];
  xBdot_tmp[64] = udata->p[1]*x_tmp[1]*xB_tmp[64]*4.0+udata->p[1]*x_tmp[64]*xB_tmp[1]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[65]*2.0-udata->p[1]*x_tmp[64]*xB_tmp[2]*2.0;
  xBdot_tmp[65] = udata->p[2]*xB_tmp[65]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[66];
  xBdot_tmp[66] = udata->p[3]*xB_tmp[66]-udata->p[3]*xB_tmp[67]*2.0;
  xBdot_tmp[67] = udata->p[3]*xB_tmp[67]-udata->p[3]*xB_tmp[68];
  xBdot_tmp[68] = udata->p[3]*xB_tmp[68]-udata->p[3]*xB_tmp[69];
  xBdot_tmp[69] = udata->p[3]*xB_tmp[69]-udata->p[3]*xB_tmp[70];
  xBdot_tmp[70] = udata->p[3]*xB_tmp[70]-udata->p[3]*xB_tmp[71];
  xBdot_tmp[71] = udata->p[3]*xB_tmp[71]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[63];
  xBdot_tmp[72] = udata->p[0]*tdata->w[7]*xB_tmp[0]-udata->p[0]*tdata->w[7]*xB_tmp[1]+udata->p[0]*tdata->w[0]*xB_tmp[72]-udata->p[0]*tdata->w[0]*xB_tmp[73];
  xBdot_tmp[73] = udata->p[1]*x_tmp[1]*xB_tmp[73]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[73]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[74]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[73]*2.0;
  xBdot_tmp[74] = udata->p[2]*xB_tmp[74]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[75];
  xBdot_tmp[75] = udata->p[3]*xB_tmp[75]-udata->p[3]*xB_tmp[76]*2.0;
  xBdot_tmp[76] = udata->p[3]*xB_tmp[76]-udata->p[3]*xB_tmp[77];
  xBdot_tmp[77] = udata->p[3]*xB_tmp[77]-udata->p[3]*xB_tmp[78];
  xBdot_tmp[78] = udata->p[3]*xB_tmp[78]-udata->p[3]*xB_tmp[79];
  xBdot_tmp[79] = udata->p[3]*xB_tmp[79]-udata->p[3]*xB_tmp[80];
  xBdot_tmp[80] = udata->p[3]*xB_tmp[80]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[72];
  xBdot_tmp[81] = udata->p[0]*tdata->w[8]*xB_tmp[0]-udata->p[0]*tdata->w[8]*xB_tmp[1]+udata->p[0]*tdata->w[0]*xB_tmp[81]-udata->p[0]*tdata->w[0]*xB_tmp[82];
  xBdot_tmp[82] = udata->p[1]*x_tmp[1]*xB_tmp[82]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[82]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[83]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[82]*2.0;
  xBdot_tmp[83] = udata->p[2]*xB_tmp[83]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[84];
  xBdot_tmp[84] = udata->p[3]*xB_tmp[84]-udata->p[3]*xB_tmp[85]*2.0;
  xBdot_tmp[85] = udata->p[3]*xB_tmp[85]-udata->p[3]*xB_tmp[86];
  xBdot_tmp[86] = udata->p[3]*xB_tmp[86]-udata->p[3]*xB_tmp[87];
  xBdot_tmp[87] = udata->p[3]*xB_tmp[87]-udata->p[3]*xB_tmp[88];
  xBdot_tmp[88] = udata->p[3]*xB_tmp[88]-udata->p[3]*xB_tmp[89];
  xBdot_tmp[89] = udata->p[3]*xB_tmp[89]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[81];
  xBdot_tmp[90] = udata->p[0]*tdata->w[9]*xB_tmp[0]-udata->p[0]*tdata->w[9]*xB_tmp[1]+udata->p[0]*tdata->w[0]*xB_tmp[90]-udata->p[0]*tdata->w[0]*xB_tmp[91];
  xBdot_tmp[91] = udata->p[1]*x_tmp[1]*xB_tmp[91]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[91]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[92]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[91]*2.0;
  xBdot_tmp[92] = udata->p[2]*xB_tmp[92]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[93];
  xBdot_tmp[93] = udata->p[3]*xB_tmp[93]-udata->p[3]*xB_tmp[94]*2.0;
  xBdot_tmp[94] = udata->p[3]*xB_tmp[94]-udata->p[3]*xB_tmp[95];
  xBdot_tmp[95] = udata->p[3]*xB_tmp[95]-udata->p[3]*xB_tmp[96];
  xBdot_tmp[96] = udata->p[3]*xB_tmp[96]-udata->p[3]*xB_tmp[97];
  xBdot_tmp[97] = udata->p[3]*xB_tmp[97]-udata->p[3]*xB_tmp[98];
  xBdot_tmp[98] = udata->p[3]*xB_tmp[98]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[90];
  xBdot_tmp[99] = udata->p[0]*tdata->w[0]*xB_tmp[99]-udata->p[0]*tdata->w[0]*xB_tmp[100];
  xBdot_tmp[100] = udata->p[1]*x_tmp[1]*xB_tmp[100]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[100]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[101]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[100]*2.0;
  xBdot_tmp[101] = udata->p[2]*xB_tmp[101]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[102];
  xBdot_tmp[102] = udata->p[3]*xB_tmp[102]-udata->p[3]*xB_tmp[103]*2.0;
  xBdot_tmp[103] = udata->p[3]*xB_tmp[103]-udata->p[3]*xB_tmp[104];
  xBdot_tmp[104] = udata->p[3]*xB_tmp[104]-udata->p[3]*xB_tmp[105];
  xBdot_tmp[105] = udata->p[3]*xB_tmp[105]-udata->p[3]*xB_tmp[106];
  xBdot_tmp[106] = udata->p[3]*xB_tmp[106]-udata->p[3]*xB_tmp[107];
  xBdot_tmp[107] = udata->p[3]*xB_tmp[107]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[99];
  xBdot_tmp[108] = udata->p[0]*tdata->w[0]*xB_tmp[108]-udata->p[0]*tdata->w[0]*xB_tmp[109];
  xBdot_tmp[109] = udata->p[1]*x_tmp[1]*xB_tmp[109]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[109]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[110]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[109]*2.0;
  xBdot_tmp[110] = udata->p[2]*xB_tmp[110]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[111];
  xBdot_tmp[111] = udata->p[3]*xB_tmp[111]-udata->p[3]*xB_tmp[112]*2.0;
  xBdot_tmp[112] = udata->p[3]*xB_tmp[112]-udata->p[3]*xB_tmp[113];
  xBdot_tmp[113] = udata->p[3]*xB_tmp[113]-udata->p[3]*xB_tmp[114];
  xBdot_tmp[114] = udata->p[3]*xB_tmp[114]-udata->p[3]*xB_tmp[115];
  xBdot_tmp[115] = udata->p[3]*xB_tmp[115]-udata->p[3]*xB_tmp[116];
  xBdot_tmp[116] = udata->p[3]*xB_tmp[116]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[108];
  xBdot_tmp[117] = udata->p[0]*tdata->w[0]*xB_tmp[117]-udata->p[0]*tdata->w[0]*xB_tmp[118];
  xBdot_tmp[118] = udata->p[1]*x_tmp[1]*xB_tmp[118]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[118]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[119]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[118]*2.0;
  xBdot_tmp[119] = udata->p[2]*xB_tmp[119]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[120];
  xBdot_tmp[120] = udata->p[3]*xB_tmp[120]-udata->p[3]*xB_tmp[121]*2.0;
  xBdot_tmp[121] = udata->p[3]*xB_tmp[121]-udata->p[3]*xB_tmp[122];
  xBdot_tmp[122] = udata->p[3]*xB_tmp[122]-udata->p[3]*xB_tmp[123];
  xBdot_tmp[123] = udata->p[3]*xB_tmp[123]-udata->p[3]*xB_tmp[124];
  xBdot_tmp[124] = udata->p[3]*xB_tmp[124]-udata->p[3]*xB_tmp[125];
  xBdot_tmp[125] = udata->p[3]*xB_tmp[125]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[117];
  xBdot_tmp[126] = udata->p[0]*tdata->w[0]*xB_tmp[126]-udata->p[0]*tdata->w[0]*xB_tmp[127];
  xBdot_tmp[127] = udata->p[1]*x_tmp[1]*xB_tmp[127]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[127]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[128]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[127]*2.0;
  xBdot_tmp[128] = udata->p[2]*xB_tmp[128]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[129];
  xBdot_tmp[129] = udata->p[3]*xB_tmp[129]-udata->p[3]*xB_tmp[130]*2.0;
  xBdot_tmp[130] = udata->p[3]*xB_tmp[130]-udata->p[3]*xB_tmp[131];
  xBdot_tmp[131] = udata->p[3]*xB_tmp[131]-udata->p[3]*xB_tmp[132];
  xBdot_tmp[132] = udata->p[3]*xB_tmp[132]-udata->p[3]*xB_tmp[133];
  xBdot_tmp[133] = udata->p[3]*xB_tmp[133]-udata->p[3]*xB_tmp[134];
  xBdot_tmp[134] = udata->p[3]*xB_tmp[134]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[126];
  xBdot_tmp[135] = udata->p[0]*tdata->w[0]*xB_tmp[135]-udata->p[0]*tdata->w[0]*xB_tmp[136];
  xBdot_tmp[136] = udata->p[1]*x_tmp[1]*xB_tmp[136]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[136]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[137]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[136]*2.0;
  xBdot_tmp[137] = udata->p[2]*xB_tmp[137]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[138];
  xBdot_tmp[138] = udata->p[3]*xB_tmp[138]-udata->p[3]*xB_tmp[139]*2.0;
  xBdot_tmp[139] = udata->p[3]*xB_tmp[139]-udata->p[3]*xB_tmp[140];
  xBdot_tmp[140] = udata->p[3]*xB_tmp[140]-udata->p[3]*xB_tmp[141];
  xBdot_tmp[141] = udata->p[3]*xB_tmp[141]-udata->p[3]*xB_tmp[142];
  xBdot_tmp[142] = udata->p[3]*xB_tmp[142]-udata->p[3]*xB_tmp[143];
  xBdot_tmp[143] = udata->p[3]*xB_tmp[143]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[135];
  xBdot_tmp[144] = udata->p[0]*tdata->w[0]*xB_tmp[144]-udata->p[0]*tdata->w[0]*xB_tmp[145];
  xBdot_tmp[145] = udata->p[1]*x_tmp[1]*xB_tmp[145]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[145]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[146]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[145]*2.0;
  xBdot_tmp[146] = udata->p[2]*xB_tmp[146]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[147];
  xBdot_tmp[147] = udata->p[3]*xB_tmp[147]-udata->p[3]*xB_tmp[148]*2.0;
  xBdot_tmp[148] = udata->p[3]*xB_tmp[148]-udata->p[3]*xB_tmp[149];
  xBdot_tmp[149] = udata->p[3]*xB_tmp[149]-udata->p[3]*xB_tmp[150];
  xBdot_tmp[150] = udata->p[3]*xB_tmp[150]-udata->p[3]*xB_tmp[151];
  xBdot_tmp[151] = udata->p[3]*xB_tmp[151]-udata->p[3]*xB_tmp[152];
  xBdot_tmp[152] = udata->p[3]*xB_tmp[152]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[144];
  xBdot_tmp[153] = udata->p[0]*tdata->w[0]*xB_tmp[153]-udata->p[0]*tdata->w[0]*xB_tmp[154];
  xBdot_tmp[154] = udata->p[1]*x_tmp[1]*xB_tmp[154]*4.0+udata->p[1]*xB_tmp[1]*x_tmp[154]*4.0-udata->p[1]*x_tmp[1]*xB_tmp[155]*2.0-udata->p[1]*xB_tmp[2]*x_tmp[154]*2.0;
  xBdot_tmp[155] = udata->p[2]*xB_tmp[155]-udata->k[0]*udata->p[2]*tdata->w[3]*xB_tmp[156];
  xBdot_tmp[156] = udata->p[3]*xB_tmp[156]-udata->p[3]*xB_tmp[157]*2.0;
  xBdot_tmp[157] = udata->p[3]*xB_tmp[157]-udata->p[3]*xB_tmp[158];
  xBdot_tmp[158] = udata->p[3]*xB_tmp[158]-udata->p[3]*xB_tmp[159];
  xBdot_tmp[159] = udata->p[3]*xB_tmp[159]-udata->p[3]*xB_tmp[160];
  xBdot_tmp[160] = udata->p[3]*xB_tmp[160]-udata->p[3]*xB_tmp[161];
  xBdot_tmp[161] = udata->p[3]*xB_tmp[161]-udata->k[1]*udata->p[3]*tdata->w[2]*xB_tmp[153];
for(ix = 0; ix<162; ix++) {
   if(amiIsNaN(xBdot_tmp[ix])) {
       xBdot_tmp[ix] = 0;       if(!tdata->nan_xBdot) {
           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_xBdot = TRUE;
       }
   }   if(amiIsInf(xBdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


