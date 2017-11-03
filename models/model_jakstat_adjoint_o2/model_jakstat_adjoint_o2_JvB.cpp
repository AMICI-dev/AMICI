
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

void JvB_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB, realtype cj, void *user_data, N_Vector tmpB1, N_Vector tmpB2) {
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
realtype *xBdot_tmp = nullptr;
if(xBdot)
    xBdot_tmp = N_VGetArrayPointer(xBdot);
realtype *vB_tmp = nullptr;
if(vB)
    vB_tmp = N_VGetArrayPointer(vB);
realtype *JvB_tmp = nullptr;
if(JvB)
    JvB_tmp = N_VGetArrayPointer(JvB);
memset(JvB_tmp,0,sizeof(realtype)*162);
w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
  JvB_tmp[0] = -tdata->p[0]*tdata->w[0]*vB_tmp[1]+udata->k[0]*tdata->p[0]*tdata->w[0]*tdata->w[2]*vB_tmp[0];
  JvB_tmp[1] = tdata->p[1]*vB_tmp[1]*tdata->dwdx[0]*2.0-tdata->p[1]*vB_tmp[2]*tdata->dwdx[0];
  JvB_tmp[2] = tdata->p[2]*vB_tmp[2]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[3];
  JvB_tmp[3] = -tdata->p[3]*vB_tmp[4]*tdata->dwdx[1]+udata->k[1]*tdata->p[3]*tdata->w[3]*vB_tmp[3];
  JvB_tmp[4] = tdata->p[3]*vB_tmp[4]-tdata->p[3]*vB_tmp[5];
  JvB_tmp[5] = tdata->p[3]*vB_tmp[5]-tdata->p[3]*vB_tmp[6];
  JvB_tmp[6] = tdata->p[3]*vB_tmp[6]-tdata->p[3]*vB_tmp[7];
  JvB_tmp[7] = tdata->p[3]*vB_tmp[7]-tdata->p[3]*vB_tmp[8];
  JvB_tmp[8] = tdata->p[3]*vB_tmp[8]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[0];
  JvB_tmp[9] = tdata->w[0]*vB_tmp[0]-tdata->w[0]*vB_tmp[1]+tdata->p[0]*tdata->w[0]*vB_tmp[9]-tdata->p[0]*tdata->w[0]*vB_tmp[10];
  JvB_tmp[10] = tdata->p[1]*x_tmp[1]*vB_tmp[10]*4.0+tdata->p[1]*x_tmp[10]*vB_tmp[1]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[11]*2.0-tdata->p[1]*x_tmp[10]*vB_tmp[2]*2.0;
  JvB_tmp[11] = tdata->p[2]*vB_tmp[11]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[12];
  JvB_tmp[12] = tdata->p[3]*vB_tmp[12]-tdata->p[3]*vB_tmp[13]*2.0;
  JvB_tmp[13] = tdata->p[3]*vB_tmp[13]-tdata->p[3]*vB_tmp[14];
  JvB_tmp[14] = tdata->p[3]*vB_tmp[14]-tdata->p[3]*vB_tmp[15];
  JvB_tmp[15] = tdata->p[3]*vB_tmp[15]-tdata->p[3]*vB_tmp[16];
  JvB_tmp[16] = tdata->p[3]*vB_tmp[16]-tdata->p[3]*vB_tmp[17];
  JvB_tmp[17] = tdata->p[3]*vB_tmp[17]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[9];
  JvB_tmp[18] = tdata->p[0]*tdata->w[0]*vB_tmp[18]-tdata->p[0]*tdata->w[0]*vB_tmp[19];
  JvB_tmp[19] = vB_tmp[1]*(tdata->dwdx[0]*2.0+tdata->p[1]*x_tmp[19]*4.0)-vB_tmp[2]*(tdata->dwdx[0]+tdata->p[1]*x_tmp[19]*2.0)+tdata->p[1]*x_tmp[1]*vB_tmp[19]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[20]*2.0;
  JvB_tmp[20] = tdata->p[2]*vB_tmp[20]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[21];
  JvB_tmp[21] = tdata->p[3]*vB_tmp[21]-tdata->p[3]*vB_tmp[22]*2.0;
  JvB_tmp[22] = tdata->p[3]*vB_tmp[22]-tdata->p[3]*vB_tmp[23];
  JvB_tmp[23] = tdata->p[3]*vB_tmp[23]-tdata->p[3]*vB_tmp[24];
  JvB_tmp[24] = tdata->p[3]*vB_tmp[24]-tdata->p[3]*vB_tmp[25];
  JvB_tmp[25] = tdata->p[3]*vB_tmp[25]-tdata->p[3]*vB_tmp[26];
  JvB_tmp[26] = tdata->p[3]*vB_tmp[26]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[18];
  JvB_tmp[27] = tdata->p[0]*tdata->w[0]*vB_tmp[27]-tdata->p[0]*tdata->w[0]*vB_tmp[28];
  JvB_tmp[28] = tdata->p[1]*x_tmp[1]*vB_tmp[28]*4.0+tdata->p[1]*x_tmp[28]*vB_tmp[1]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[29]*2.0-tdata->p[1]*x_tmp[28]*vB_tmp[2]*2.0;
  JvB_tmp[29] = vB_tmp[2]+tdata->p[2]*vB_tmp[29]-udata->k[0]*tdata->w[3]*vB_tmp[3]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[30];
  JvB_tmp[30] = tdata->p[3]*vB_tmp[30]-tdata->p[3]*vB_tmp[31]*2.0;
  JvB_tmp[31] = tdata->p[3]*vB_tmp[31]-tdata->p[3]*vB_tmp[32];
  JvB_tmp[32] = tdata->p[3]*vB_tmp[32]-tdata->p[3]*vB_tmp[33];
  JvB_tmp[33] = tdata->p[3]*vB_tmp[33]-tdata->p[3]*vB_tmp[34];
  JvB_tmp[34] = tdata->p[3]*vB_tmp[34]-tdata->p[3]*vB_tmp[35];
  JvB_tmp[35] = tdata->p[3]*vB_tmp[35]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[27];
  JvB_tmp[36] = tdata->p[0]*tdata->w[0]*vB_tmp[36]-tdata->p[0]*tdata->w[0]*vB_tmp[37];
  JvB_tmp[37] = tdata->p[1]*x_tmp[1]*vB_tmp[37]*4.0+tdata->p[1]*x_tmp[37]*vB_tmp[1]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[38]*2.0-tdata->p[1]*x_tmp[37]*vB_tmp[2]*2.0;
  JvB_tmp[38] = tdata->p[2]*vB_tmp[38]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[39];
  JvB_tmp[39] = vB_tmp[3]+tdata->p[3]*vB_tmp[39]-tdata->p[3]*vB_tmp[40]*2.0-vB_tmp[4]*tdata->dwdx[1];
  JvB_tmp[40] = vB_tmp[4]-vB_tmp[5]+tdata->p[3]*vB_tmp[40]-tdata->p[3]*vB_tmp[41];
  JvB_tmp[41] = vB_tmp[5]-vB_tmp[6]+tdata->p[3]*vB_tmp[41]-tdata->p[3]*vB_tmp[42];
  JvB_tmp[42] = vB_tmp[6]-vB_tmp[7]+tdata->p[3]*vB_tmp[42]-tdata->p[3]*vB_tmp[43];
  JvB_tmp[43] = vB_tmp[7]-vB_tmp[8]+tdata->p[3]*vB_tmp[43]-tdata->p[3]*vB_tmp[44];
  JvB_tmp[44] = vB_tmp[8]+tdata->p[3]*vB_tmp[44]-udata->k[1]*tdata->w[2]*vB_tmp[0]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[36];
  JvB_tmp[45] = tdata->p[0]*tdata->w[0]*vB_tmp[45]-tdata->p[0]*tdata->w[0]*vB_tmp[46];
  JvB_tmp[46] = tdata->p[1]*x_tmp[1]*vB_tmp[46]*4.0+tdata->p[1]*x_tmp[46]*vB_tmp[1]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[47]*2.0-tdata->p[1]*x_tmp[46]*vB_tmp[2]*2.0;
  JvB_tmp[47] = tdata->p[2]*vB_tmp[47]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[48];
  JvB_tmp[48] = tdata->p[3]*vB_tmp[48]-tdata->p[3]*vB_tmp[49]*2.0;
  JvB_tmp[49] = tdata->p[3]*vB_tmp[49]-tdata->p[3]*vB_tmp[50];
  JvB_tmp[50] = tdata->p[3]*vB_tmp[50]-tdata->p[3]*vB_tmp[51];
  JvB_tmp[51] = tdata->p[3]*vB_tmp[51]-tdata->p[3]*vB_tmp[52];
  JvB_tmp[52] = tdata->p[3]*vB_tmp[52]-tdata->p[3]*vB_tmp[53];
  JvB_tmp[53] = tdata->p[3]*vB_tmp[53]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[45];
  JvB_tmp[54] = tdata->p[0]*tdata->w[5]*vB_tmp[0]-tdata->p[0]*tdata->w[5]*vB_tmp[1]+tdata->p[0]*tdata->w[0]*vB_tmp[54]-tdata->p[0]*tdata->w[0]*vB_tmp[55];
  JvB_tmp[55] = tdata->p[1]*x_tmp[1]*vB_tmp[55]*4.0+tdata->p[1]*x_tmp[55]*vB_tmp[1]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[56]*2.0-tdata->p[1]*x_tmp[55]*vB_tmp[2]*2.0;
  JvB_tmp[56] = tdata->p[2]*vB_tmp[56]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[57];
  JvB_tmp[57] = tdata->p[3]*vB_tmp[57]-tdata->p[3]*vB_tmp[58]*2.0;
  JvB_tmp[58] = tdata->p[3]*vB_tmp[58]-tdata->p[3]*vB_tmp[59];
  JvB_tmp[59] = tdata->p[3]*vB_tmp[59]-tdata->p[3]*vB_tmp[60];
  JvB_tmp[60] = tdata->p[3]*vB_tmp[60]-tdata->p[3]*vB_tmp[61];
  JvB_tmp[61] = tdata->p[3]*vB_tmp[61]-tdata->p[3]*vB_tmp[62];
  JvB_tmp[62] = tdata->p[3]*vB_tmp[62]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[54];
  JvB_tmp[63] = tdata->p[0]*tdata->w[6]*vB_tmp[0]-tdata->p[0]*tdata->w[6]*vB_tmp[1]+tdata->p[0]*tdata->w[0]*vB_tmp[63]-tdata->p[0]*tdata->w[0]*vB_tmp[64];
  JvB_tmp[64] = tdata->p[1]*x_tmp[1]*vB_tmp[64]*4.0+tdata->p[1]*x_tmp[64]*vB_tmp[1]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[65]*2.0-tdata->p[1]*x_tmp[64]*vB_tmp[2]*2.0;
  JvB_tmp[65] = tdata->p[2]*vB_tmp[65]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[66];
  JvB_tmp[66] = tdata->p[3]*vB_tmp[66]-tdata->p[3]*vB_tmp[67]*2.0;
  JvB_tmp[67] = tdata->p[3]*vB_tmp[67]-tdata->p[3]*vB_tmp[68];
  JvB_tmp[68] = tdata->p[3]*vB_tmp[68]-tdata->p[3]*vB_tmp[69];
  JvB_tmp[69] = tdata->p[3]*vB_tmp[69]-tdata->p[3]*vB_tmp[70];
  JvB_tmp[70] = tdata->p[3]*vB_tmp[70]-tdata->p[3]*vB_tmp[71];
  JvB_tmp[71] = tdata->p[3]*vB_tmp[71]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[63];
  JvB_tmp[72] = tdata->p[0]*tdata->w[7]*vB_tmp[0]-tdata->p[0]*tdata->w[7]*vB_tmp[1]+tdata->p[0]*tdata->w[0]*vB_tmp[72]-tdata->p[0]*tdata->w[0]*vB_tmp[73];
  JvB_tmp[73] = tdata->p[1]*x_tmp[1]*vB_tmp[73]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[73]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[74]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[73]*2.0;
  JvB_tmp[74] = tdata->p[2]*vB_tmp[74]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[75];
  JvB_tmp[75] = tdata->p[3]*vB_tmp[75]-tdata->p[3]*vB_tmp[76]*2.0;
  JvB_tmp[76] = tdata->p[3]*vB_tmp[76]-tdata->p[3]*vB_tmp[77];
  JvB_tmp[77] = tdata->p[3]*vB_tmp[77]-tdata->p[3]*vB_tmp[78];
  JvB_tmp[78] = tdata->p[3]*vB_tmp[78]-tdata->p[3]*vB_tmp[79];
  JvB_tmp[79] = tdata->p[3]*vB_tmp[79]-tdata->p[3]*vB_tmp[80];
  JvB_tmp[80] = tdata->p[3]*vB_tmp[80]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[72];
  JvB_tmp[81] = tdata->p[0]*tdata->w[8]*vB_tmp[0]-tdata->p[0]*tdata->w[8]*vB_tmp[1]+tdata->p[0]*tdata->w[0]*vB_tmp[81]-tdata->p[0]*tdata->w[0]*vB_tmp[82];
  JvB_tmp[82] = tdata->p[1]*x_tmp[1]*vB_tmp[82]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[82]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[83]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[82]*2.0;
  JvB_tmp[83] = tdata->p[2]*vB_tmp[83]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[84];
  JvB_tmp[84] = tdata->p[3]*vB_tmp[84]-tdata->p[3]*vB_tmp[85]*2.0;
  JvB_tmp[85] = tdata->p[3]*vB_tmp[85]-tdata->p[3]*vB_tmp[86];
  JvB_tmp[86] = tdata->p[3]*vB_tmp[86]-tdata->p[3]*vB_tmp[87];
  JvB_tmp[87] = tdata->p[3]*vB_tmp[87]-tdata->p[3]*vB_tmp[88];
  JvB_tmp[88] = tdata->p[3]*vB_tmp[88]-tdata->p[3]*vB_tmp[89];
  JvB_tmp[89] = tdata->p[3]*vB_tmp[89]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[81];
  JvB_tmp[90] = tdata->p[0]*tdata->w[9]*vB_tmp[0]-tdata->p[0]*tdata->w[9]*vB_tmp[1]+tdata->p[0]*tdata->w[0]*vB_tmp[90]-tdata->p[0]*tdata->w[0]*vB_tmp[91];
  JvB_tmp[91] = tdata->p[1]*x_tmp[1]*vB_tmp[91]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[91]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[92]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[91]*2.0;
  JvB_tmp[92] = tdata->p[2]*vB_tmp[92]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[93];
  JvB_tmp[93] = tdata->p[3]*vB_tmp[93]-tdata->p[3]*vB_tmp[94]*2.0;
  JvB_tmp[94] = tdata->p[3]*vB_tmp[94]-tdata->p[3]*vB_tmp[95];
  JvB_tmp[95] = tdata->p[3]*vB_tmp[95]-tdata->p[3]*vB_tmp[96];
  JvB_tmp[96] = tdata->p[3]*vB_tmp[96]-tdata->p[3]*vB_tmp[97];
  JvB_tmp[97] = tdata->p[3]*vB_tmp[97]-tdata->p[3]*vB_tmp[98];
  JvB_tmp[98] = tdata->p[3]*vB_tmp[98]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[90];
  JvB_tmp[99] = tdata->p[0]*tdata->w[0]*vB_tmp[99]-tdata->p[0]*tdata->w[0]*vB_tmp[100];
  JvB_tmp[100] = tdata->p[1]*x_tmp[1]*vB_tmp[100]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[100]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[101]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[100]*2.0;
  JvB_tmp[101] = tdata->p[2]*vB_tmp[101]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[102];
  JvB_tmp[102] = tdata->p[3]*vB_tmp[102]-tdata->p[3]*vB_tmp[103]*2.0;
  JvB_tmp[103] = tdata->p[3]*vB_tmp[103]-tdata->p[3]*vB_tmp[104];
  JvB_tmp[104] = tdata->p[3]*vB_tmp[104]-tdata->p[3]*vB_tmp[105];
  JvB_tmp[105] = tdata->p[3]*vB_tmp[105]-tdata->p[3]*vB_tmp[106];
  JvB_tmp[106] = tdata->p[3]*vB_tmp[106]-tdata->p[3]*vB_tmp[107];
  JvB_tmp[107] = tdata->p[3]*vB_tmp[107]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[99];
  JvB_tmp[108] = tdata->p[0]*tdata->w[0]*vB_tmp[108]-tdata->p[0]*tdata->w[0]*vB_tmp[109];
  JvB_tmp[109] = tdata->p[1]*x_tmp[1]*vB_tmp[109]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[109]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[110]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[109]*2.0;
  JvB_tmp[110] = tdata->p[2]*vB_tmp[110]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[111];
  JvB_tmp[111] = tdata->p[3]*vB_tmp[111]-tdata->p[3]*vB_tmp[112]*2.0;
  JvB_tmp[112] = tdata->p[3]*vB_tmp[112]-tdata->p[3]*vB_tmp[113];
  JvB_tmp[113] = tdata->p[3]*vB_tmp[113]-tdata->p[3]*vB_tmp[114];
  JvB_tmp[114] = tdata->p[3]*vB_tmp[114]-tdata->p[3]*vB_tmp[115];
  JvB_tmp[115] = tdata->p[3]*vB_tmp[115]-tdata->p[3]*vB_tmp[116];
  JvB_tmp[116] = tdata->p[3]*vB_tmp[116]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[108];
  JvB_tmp[117] = tdata->p[0]*tdata->w[0]*vB_tmp[117]-tdata->p[0]*tdata->w[0]*vB_tmp[118];
  JvB_tmp[118] = tdata->p[1]*x_tmp[1]*vB_tmp[118]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[118]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[119]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[118]*2.0;
  JvB_tmp[119] = tdata->p[2]*vB_tmp[119]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[120];
  JvB_tmp[120] = tdata->p[3]*vB_tmp[120]-tdata->p[3]*vB_tmp[121]*2.0;
  JvB_tmp[121] = tdata->p[3]*vB_tmp[121]-tdata->p[3]*vB_tmp[122];
  JvB_tmp[122] = tdata->p[3]*vB_tmp[122]-tdata->p[3]*vB_tmp[123];
  JvB_tmp[123] = tdata->p[3]*vB_tmp[123]-tdata->p[3]*vB_tmp[124];
  JvB_tmp[124] = tdata->p[3]*vB_tmp[124]-tdata->p[3]*vB_tmp[125];
  JvB_tmp[125] = tdata->p[3]*vB_tmp[125]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[117];
  JvB_tmp[126] = tdata->p[0]*tdata->w[0]*vB_tmp[126]-tdata->p[0]*tdata->w[0]*vB_tmp[127];
  JvB_tmp[127] = tdata->p[1]*x_tmp[1]*vB_tmp[127]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[127]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[128]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[127]*2.0;
  JvB_tmp[128] = tdata->p[2]*vB_tmp[128]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[129];
  JvB_tmp[129] = tdata->p[3]*vB_tmp[129]-tdata->p[3]*vB_tmp[130]*2.0;
  JvB_tmp[130] = tdata->p[3]*vB_tmp[130]-tdata->p[3]*vB_tmp[131];
  JvB_tmp[131] = tdata->p[3]*vB_tmp[131]-tdata->p[3]*vB_tmp[132];
  JvB_tmp[132] = tdata->p[3]*vB_tmp[132]-tdata->p[3]*vB_tmp[133];
  JvB_tmp[133] = tdata->p[3]*vB_tmp[133]-tdata->p[3]*vB_tmp[134];
  JvB_tmp[134] = tdata->p[3]*vB_tmp[134]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[126];
  JvB_tmp[135] = tdata->p[0]*tdata->w[0]*vB_tmp[135]-tdata->p[0]*tdata->w[0]*vB_tmp[136];
  JvB_tmp[136] = tdata->p[1]*x_tmp[1]*vB_tmp[136]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[136]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[137]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[136]*2.0;
  JvB_tmp[137] = tdata->p[2]*vB_tmp[137]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[138];
  JvB_tmp[138] = tdata->p[3]*vB_tmp[138]-tdata->p[3]*vB_tmp[139]*2.0;
  JvB_tmp[139] = tdata->p[3]*vB_tmp[139]-tdata->p[3]*vB_tmp[140];
  JvB_tmp[140] = tdata->p[3]*vB_tmp[140]-tdata->p[3]*vB_tmp[141];
  JvB_tmp[141] = tdata->p[3]*vB_tmp[141]-tdata->p[3]*vB_tmp[142];
  JvB_tmp[142] = tdata->p[3]*vB_tmp[142]-tdata->p[3]*vB_tmp[143];
  JvB_tmp[143] = tdata->p[3]*vB_tmp[143]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[135];
  JvB_tmp[144] = tdata->p[0]*tdata->w[0]*vB_tmp[144]-tdata->p[0]*tdata->w[0]*vB_tmp[145];
  JvB_tmp[145] = tdata->p[1]*x_tmp[1]*vB_tmp[145]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[145]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[146]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[145]*2.0;
  JvB_tmp[146] = tdata->p[2]*vB_tmp[146]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[147];
  JvB_tmp[147] = tdata->p[3]*vB_tmp[147]-tdata->p[3]*vB_tmp[148]*2.0;
  JvB_tmp[148] = tdata->p[3]*vB_tmp[148]-tdata->p[3]*vB_tmp[149];
  JvB_tmp[149] = tdata->p[3]*vB_tmp[149]-tdata->p[3]*vB_tmp[150];
  JvB_tmp[150] = tdata->p[3]*vB_tmp[150]-tdata->p[3]*vB_tmp[151];
  JvB_tmp[151] = tdata->p[3]*vB_tmp[151]-tdata->p[3]*vB_tmp[152];
  JvB_tmp[152] = tdata->p[3]*vB_tmp[152]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[144];
  JvB_tmp[153] = tdata->p[0]*tdata->w[0]*vB_tmp[153]-tdata->p[0]*tdata->w[0]*vB_tmp[154];
  JvB_tmp[154] = tdata->p[1]*x_tmp[1]*vB_tmp[154]*4.0+tdata->p[1]*vB_tmp[1]*x_tmp[154]*4.0-tdata->p[1]*x_tmp[1]*vB_tmp[155]*2.0-tdata->p[1]*vB_tmp[2]*x_tmp[154]*2.0;
  JvB_tmp[155] = tdata->p[2]*vB_tmp[155]-udata->k[0]*tdata->p[2]*tdata->w[3]*vB_tmp[156];
  JvB_tmp[156] = tdata->p[3]*vB_tmp[156]-tdata->p[3]*vB_tmp[157]*2.0;
  JvB_tmp[157] = tdata->p[3]*vB_tmp[157]-tdata->p[3]*vB_tmp[158];
  JvB_tmp[158] = tdata->p[3]*vB_tmp[158]-tdata->p[3]*vB_tmp[159];
  JvB_tmp[159] = tdata->p[3]*vB_tmp[159]-tdata->p[3]*vB_tmp[160];
  JvB_tmp[160] = tdata->p[3]*vB_tmp[160]-tdata->p[3]*vB_tmp[161];
  JvB_tmp[161] = tdata->p[3]*vB_tmp[161]-udata->k[1]*tdata->p[3]*tdata->w[2]*vB_tmp[153];
return;

}


