
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_dwdx.h"
#include "model_jakstat_adjoint_o2_w.h"

int xBdot_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
int ix;
memset(xBdot_tmp,0,sizeof(realtype)*162);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
status = dwdx_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  xBdot_tmp[0] = -p[0]*w_tmp[0]*xB_tmp[1]+k[0]*p[0]*w_tmp[0]*w_tmp[2]*xB_tmp[0];
  xBdot_tmp[1] = dwdx_tmp[0]*p[1]*xB_tmp[1]*2.0-dwdx_tmp[0]*p[1]*xB_tmp[2];
  xBdot_tmp[2] = p[2]*xB_tmp[2]-k[0]*p[2]*w_tmp[3]*xB_tmp[3];
  xBdot_tmp[3] = -dwdx_tmp[1]*p[3]*xB_tmp[4]+k[1]*p[3]*w_tmp[3]*xB_tmp[3];
  xBdot_tmp[4] = p[3]*xB_tmp[4]-p[3]*xB_tmp[5];
  xBdot_tmp[5] = p[3]*xB_tmp[5]-p[3]*xB_tmp[6];
  xBdot_tmp[6] = p[3]*xB_tmp[6]-p[3]*xB_tmp[7];
  xBdot_tmp[7] = p[3]*xB_tmp[7]-p[3]*xB_tmp[8];
  xBdot_tmp[8] = p[3]*xB_tmp[8]-k[1]*p[3]*w_tmp[2]*xB_tmp[0];
  xBdot_tmp[9] = w_tmp[0]*xB_tmp[0]-w_tmp[0]*xB_tmp[1]+p[0]*w_tmp[0]*xB_tmp[9]-p[0]*w_tmp[0]*xB_tmp[10];
  xBdot_tmp[10] = p[1]*xB_tmp[1]*x_tmp[10]*4.0+p[1]*xB_tmp[10]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[10]*2.0-p[1]*xB_tmp[11]*x_tmp[1]*2.0;
  xBdot_tmp[11] = p[2]*xB_tmp[11]-k[0]*p[2]*w_tmp[3]*xB_tmp[12];
  xBdot_tmp[12] = p[3]*xB_tmp[12]-p[3]*xB_tmp[13]*2.0;
  xBdot_tmp[13] = p[3]*xB_tmp[13]-p[3]*xB_tmp[14];
  xBdot_tmp[14] = p[3]*xB_tmp[14]-p[3]*xB_tmp[15];
  xBdot_tmp[15] = p[3]*xB_tmp[15]-p[3]*xB_tmp[16];
  xBdot_tmp[16] = p[3]*xB_tmp[16]-p[3]*xB_tmp[17];
  xBdot_tmp[17] = p[3]*xB_tmp[17]-k[1]*p[3]*w_tmp[2]*xB_tmp[9];
  xBdot_tmp[18] = p[0]*w_tmp[0]*xB_tmp[18]-p[0]*w_tmp[0]*xB_tmp[19];
  xBdot_tmp[19] = xB_tmp[1]*(dwdx_tmp[0]*2.0+p[1]*x_tmp[19]*4.0)-xB_tmp[2]*(dwdx_tmp[0]+p[1]*x_tmp[19]*2.0)+p[1]*xB_tmp[19]*x_tmp[1]*4.0-p[1]*xB_tmp[20]*x_tmp[1]*2.0;
  xBdot_tmp[20] = p[2]*xB_tmp[20]-k[0]*p[2]*w_tmp[3]*xB_tmp[21];
  xBdot_tmp[21] = p[3]*xB_tmp[21]-p[3]*xB_tmp[22]*2.0;
  xBdot_tmp[22] = p[3]*xB_tmp[22]-p[3]*xB_tmp[23];
  xBdot_tmp[23] = p[3]*xB_tmp[23]-p[3]*xB_tmp[24];
  xBdot_tmp[24] = p[3]*xB_tmp[24]-p[3]*xB_tmp[25];
  xBdot_tmp[25] = p[3]*xB_tmp[25]-p[3]*xB_tmp[26];
  xBdot_tmp[26] = p[3]*xB_tmp[26]-k[1]*p[3]*w_tmp[2]*xB_tmp[18];
  xBdot_tmp[27] = p[0]*w_tmp[0]*xB_tmp[27]-p[0]*w_tmp[0]*xB_tmp[28];
  xBdot_tmp[28] = p[1]*xB_tmp[1]*x_tmp[28]*4.0+p[1]*xB_tmp[28]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[28]*2.0-p[1]*xB_tmp[29]*x_tmp[1]*2.0;
  xBdot_tmp[29] = xB_tmp[2]+p[2]*xB_tmp[29]-k[0]*w_tmp[3]*xB_tmp[3]-k[0]*p[2]*w_tmp[3]*xB_tmp[30];
  xBdot_tmp[30] = p[3]*xB_tmp[30]-p[3]*xB_tmp[31]*2.0;
  xBdot_tmp[31] = p[3]*xB_tmp[31]-p[3]*xB_tmp[32];
  xBdot_tmp[32] = p[3]*xB_tmp[32]-p[3]*xB_tmp[33];
  xBdot_tmp[33] = p[3]*xB_tmp[33]-p[3]*xB_tmp[34];
  xBdot_tmp[34] = p[3]*xB_tmp[34]-p[3]*xB_tmp[35];
  xBdot_tmp[35] = p[3]*xB_tmp[35]-k[1]*p[3]*w_tmp[2]*xB_tmp[27];
  xBdot_tmp[36] = p[0]*w_tmp[0]*xB_tmp[36]-p[0]*w_tmp[0]*xB_tmp[37];
  xBdot_tmp[37] = p[1]*xB_tmp[1]*x_tmp[37]*4.0+p[1]*xB_tmp[37]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[37]*2.0-p[1]*xB_tmp[38]*x_tmp[1]*2.0;
  xBdot_tmp[38] = p[2]*xB_tmp[38]-k[0]*p[2]*w_tmp[3]*xB_tmp[39];
  xBdot_tmp[39] = xB_tmp[3]-dwdx_tmp[1]*xB_tmp[4]+p[3]*xB_tmp[39]-p[3]*xB_tmp[40]*2.0;
  xBdot_tmp[40] = xB_tmp[4]-xB_tmp[5]+p[3]*xB_tmp[40]-p[3]*xB_tmp[41];
  xBdot_tmp[41] = xB_tmp[5]-xB_tmp[6]+p[3]*xB_tmp[41]-p[3]*xB_tmp[42];
  xBdot_tmp[42] = xB_tmp[6]-xB_tmp[7]+p[3]*xB_tmp[42]-p[3]*xB_tmp[43];
  xBdot_tmp[43] = xB_tmp[7]-xB_tmp[8]+p[3]*xB_tmp[43]-p[3]*xB_tmp[44];
  xBdot_tmp[44] = xB_tmp[8]+p[3]*xB_tmp[44]-k[1]*w_tmp[2]*xB_tmp[0]-k[1]*p[3]*w_tmp[2]*xB_tmp[36];
  xBdot_tmp[45] = p[0]*w_tmp[0]*xB_tmp[45]-p[0]*w_tmp[0]*xB_tmp[46];
  xBdot_tmp[46] = p[1]*xB_tmp[1]*x_tmp[46]*4.0+p[1]*xB_tmp[46]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[46]*2.0-p[1]*xB_tmp[47]*x_tmp[1]*2.0;
  xBdot_tmp[47] = p[2]*xB_tmp[47]-k[0]*p[2]*w_tmp[3]*xB_tmp[48];
  xBdot_tmp[48] = p[3]*xB_tmp[48]-p[3]*xB_tmp[49]*2.0;
  xBdot_tmp[49] = p[3]*xB_tmp[49]-p[3]*xB_tmp[50];
  xBdot_tmp[50] = p[3]*xB_tmp[50]-p[3]*xB_tmp[51];
  xBdot_tmp[51] = p[3]*xB_tmp[51]-p[3]*xB_tmp[52];
  xBdot_tmp[52] = p[3]*xB_tmp[52]-p[3]*xB_tmp[53];
  xBdot_tmp[53] = p[3]*xB_tmp[53]-k[1]*p[3]*w_tmp[2]*xB_tmp[45];
  xBdot_tmp[54] = p[0]*w_tmp[5]*xB_tmp[0]-p[0]*w_tmp[5]*xB_tmp[1]+p[0]*w_tmp[0]*xB_tmp[54]-p[0]*w_tmp[0]*xB_tmp[55];
  xBdot_tmp[55] = p[1]*xB_tmp[1]*x_tmp[55]*4.0+p[1]*xB_tmp[55]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[55]*2.0-p[1]*xB_tmp[56]*x_tmp[1]*2.0;
  xBdot_tmp[56] = p[2]*xB_tmp[56]-k[0]*p[2]*w_tmp[3]*xB_tmp[57];
  xBdot_tmp[57] = p[3]*xB_tmp[57]-p[3]*xB_tmp[58]*2.0;
  xBdot_tmp[58] = p[3]*xB_tmp[58]-p[3]*xB_tmp[59];
  xBdot_tmp[59] = p[3]*xB_tmp[59]-p[3]*xB_tmp[60];
  xBdot_tmp[60] = p[3]*xB_tmp[60]-p[3]*xB_tmp[61];
  xBdot_tmp[61] = p[3]*xB_tmp[61]-p[3]*xB_tmp[62];
  xBdot_tmp[62] = p[3]*xB_tmp[62]-k[1]*p[3]*w_tmp[2]*xB_tmp[54];
  xBdot_tmp[63] = p[0]*w_tmp[6]*xB_tmp[0]-p[0]*w_tmp[6]*xB_tmp[1]+p[0]*w_tmp[0]*xB_tmp[63]-p[0]*w_tmp[0]*xB_tmp[64];
  xBdot_tmp[64] = p[1]*xB_tmp[1]*x_tmp[64]*4.0+p[1]*xB_tmp[64]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[64]*2.0-p[1]*xB_tmp[65]*x_tmp[1]*2.0;
  xBdot_tmp[65] = p[2]*xB_tmp[65]-k[0]*p[2]*w_tmp[3]*xB_tmp[66];
  xBdot_tmp[66] = p[3]*xB_tmp[66]-p[3]*xB_tmp[67]*2.0;
  xBdot_tmp[67] = p[3]*xB_tmp[67]-p[3]*xB_tmp[68];
  xBdot_tmp[68] = p[3]*xB_tmp[68]-p[3]*xB_tmp[69];
  xBdot_tmp[69] = p[3]*xB_tmp[69]-p[3]*xB_tmp[70];
  xBdot_tmp[70] = p[3]*xB_tmp[70]-p[3]*xB_tmp[71];
  xBdot_tmp[71] = p[3]*xB_tmp[71]-k[1]*p[3]*w_tmp[2]*xB_tmp[63];
  xBdot_tmp[72] = p[0]*w_tmp[7]*xB_tmp[0]-p[0]*w_tmp[7]*xB_tmp[1]+p[0]*w_tmp[0]*xB_tmp[72]-p[0]*w_tmp[0]*xB_tmp[73];
  xBdot_tmp[73] = p[1]*xB_tmp[1]*x_tmp[73]*4.0+p[1]*xB_tmp[73]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[73]*2.0-p[1]*xB_tmp[74]*x_tmp[1]*2.0;
  xBdot_tmp[74] = p[2]*xB_tmp[74]-k[0]*p[2]*w_tmp[3]*xB_tmp[75];
  xBdot_tmp[75] = p[3]*xB_tmp[75]-p[3]*xB_tmp[76]*2.0;
  xBdot_tmp[76] = p[3]*xB_tmp[76]-p[3]*xB_tmp[77];
  xBdot_tmp[77] = p[3]*xB_tmp[77]-p[3]*xB_tmp[78];
  xBdot_tmp[78] = p[3]*xB_tmp[78]-p[3]*xB_tmp[79];
  xBdot_tmp[79] = p[3]*xB_tmp[79]-p[3]*xB_tmp[80];
  xBdot_tmp[80] = p[3]*xB_tmp[80]-k[1]*p[3]*w_tmp[2]*xB_tmp[72];
  xBdot_tmp[81] = p[0]*w_tmp[8]*xB_tmp[0]-p[0]*w_tmp[8]*xB_tmp[1]+p[0]*w_tmp[0]*xB_tmp[81]-p[0]*w_tmp[0]*xB_tmp[82];
  xBdot_tmp[82] = p[1]*xB_tmp[1]*x_tmp[82]*4.0+p[1]*xB_tmp[82]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[82]*2.0-p[1]*xB_tmp[83]*x_tmp[1]*2.0;
  xBdot_tmp[83] = p[2]*xB_tmp[83]-k[0]*p[2]*w_tmp[3]*xB_tmp[84];
  xBdot_tmp[84] = p[3]*xB_tmp[84]-p[3]*xB_tmp[85]*2.0;
  xBdot_tmp[85] = p[3]*xB_tmp[85]-p[3]*xB_tmp[86];
  xBdot_tmp[86] = p[3]*xB_tmp[86]-p[3]*xB_tmp[87];
  xBdot_tmp[87] = p[3]*xB_tmp[87]-p[3]*xB_tmp[88];
  xBdot_tmp[88] = p[3]*xB_tmp[88]-p[3]*xB_tmp[89];
  xBdot_tmp[89] = p[3]*xB_tmp[89]-k[1]*p[3]*w_tmp[2]*xB_tmp[81];
  xBdot_tmp[90] = p[0]*w_tmp[9]*xB_tmp[0]-p[0]*w_tmp[9]*xB_tmp[1]+p[0]*w_tmp[0]*xB_tmp[90]-p[0]*w_tmp[0]*xB_tmp[91];
  xBdot_tmp[91] = p[1]*xB_tmp[1]*x_tmp[91]*4.0+p[1]*xB_tmp[91]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[91]*2.0-p[1]*xB_tmp[92]*x_tmp[1]*2.0;
  xBdot_tmp[92] = p[2]*xB_tmp[92]-k[0]*p[2]*w_tmp[3]*xB_tmp[93];
  xBdot_tmp[93] = p[3]*xB_tmp[93]-p[3]*xB_tmp[94]*2.0;
  xBdot_tmp[94] = p[3]*xB_tmp[94]-p[3]*xB_tmp[95];
  xBdot_tmp[95] = p[3]*xB_tmp[95]-p[3]*xB_tmp[96];
  xBdot_tmp[96] = p[3]*xB_tmp[96]-p[3]*xB_tmp[97];
  xBdot_tmp[97] = p[3]*xB_tmp[97]-p[3]*xB_tmp[98];
  xBdot_tmp[98] = p[3]*xB_tmp[98]-k[1]*p[3]*w_tmp[2]*xB_tmp[90];
  xBdot_tmp[99] = p[0]*w_tmp[0]*xB_tmp[99]-p[0]*w_tmp[0]*xB_tmp[100];
  xBdot_tmp[100] = p[1]*xB_tmp[1]*x_tmp[100]*4.0+p[1]*xB_tmp[100]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[100]*2.0-p[1]*xB_tmp[101]*x_tmp[1]*2.0;
  xBdot_tmp[101] = p[2]*xB_tmp[101]-k[0]*p[2]*w_tmp[3]*xB_tmp[102];
  xBdot_tmp[102] = p[3]*xB_tmp[102]-p[3]*xB_tmp[103]*2.0;
  xBdot_tmp[103] = p[3]*xB_tmp[103]-p[3]*xB_tmp[104];
  xBdot_tmp[104] = p[3]*xB_tmp[104]-p[3]*xB_tmp[105];
  xBdot_tmp[105] = p[3]*xB_tmp[105]-p[3]*xB_tmp[106];
  xBdot_tmp[106] = p[3]*xB_tmp[106]-p[3]*xB_tmp[107];
  xBdot_tmp[107] = p[3]*xB_tmp[107]-k[1]*p[3]*w_tmp[2]*xB_tmp[99];
  xBdot_tmp[108] = p[0]*w_tmp[0]*xB_tmp[108]-p[0]*w_tmp[0]*xB_tmp[109];
  xBdot_tmp[109] = p[1]*xB_tmp[1]*x_tmp[109]*4.0+p[1]*xB_tmp[109]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[109]*2.0-p[1]*xB_tmp[110]*x_tmp[1]*2.0;
  xBdot_tmp[110] = p[2]*xB_tmp[110]-k[0]*p[2]*w_tmp[3]*xB_tmp[111];
  xBdot_tmp[111] = p[3]*xB_tmp[111]-p[3]*xB_tmp[112]*2.0;
  xBdot_tmp[112] = p[3]*xB_tmp[112]-p[3]*xB_tmp[113];
  xBdot_tmp[113] = p[3]*xB_tmp[113]-p[3]*xB_tmp[114];
  xBdot_tmp[114] = p[3]*xB_tmp[114]-p[3]*xB_tmp[115];
  xBdot_tmp[115] = p[3]*xB_tmp[115]-p[3]*xB_tmp[116];
  xBdot_tmp[116] = p[3]*xB_tmp[116]-k[1]*p[3]*w_tmp[2]*xB_tmp[108];
  xBdot_tmp[117] = p[0]*w_tmp[0]*xB_tmp[117]-p[0]*w_tmp[0]*xB_tmp[118];
  xBdot_tmp[118] = p[1]*xB_tmp[1]*x_tmp[118]*4.0+p[1]*xB_tmp[118]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[118]*2.0-p[1]*xB_tmp[119]*x_tmp[1]*2.0;
  xBdot_tmp[119] = p[2]*xB_tmp[119]-k[0]*p[2]*w_tmp[3]*xB_tmp[120];
  xBdot_tmp[120] = p[3]*xB_tmp[120]-p[3]*xB_tmp[121]*2.0;
  xBdot_tmp[121] = p[3]*xB_tmp[121]-p[3]*xB_tmp[122];
  xBdot_tmp[122] = p[3]*xB_tmp[122]-p[3]*xB_tmp[123];
  xBdot_tmp[123] = p[3]*xB_tmp[123]-p[3]*xB_tmp[124];
  xBdot_tmp[124] = p[3]*xB_tmp[124]-p[3]*xB_tmp[125];
  xBdot_tmp[125] = p[3]*xB_tmp[125]-k[1]*p[3]*w_tmp[2]*xB_tmp[117];
  xBdot_tmp[126] = p[0]*w_tmp[0]*xB_tmp[126]-p[0]*w_tmp[0]*xB_tmp[127];
  xBdot_tmp[127] = p[1]*xB_tmp[1]*x_tmp[127]*4.0+p[1]*xB_tmp[127]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[127]*2.0-p[1]*xB_tmp[128]*x_tmp[1]*2.0;
  xBdot_tmp[128] = p[2]*xB_tmp[128]-k[0]*p[2]*w_tmp[3]*xB_tmp[129];
  xBdot_tmp[129] = p[3]*xB_tmp[129]-p[3]*xB_tmp[130]*2.0;
  xBdot_tmp[130] = p[3]*xB_tmp[130]-p[3]*xB_tmp[131];
  xBdot_tmp[131] = p[3]*xB_tmp[131]-p[3]*xB_tmp[132];
  xBdot_tmp[132] = p[3]*xB_tmp[132]-p[3]*xB_tmp[133];
  xBdot_tmp[133] = p[3]*xB_tmp[133]-p[3]*xB_tmp[134];
  xBdot_tmp[134] = p[3]*xB_tmp[134]-k[1]*p[3]*w_tmp[2]*xB_tmp[126];
  xBdot_tmp[135] = p[0]*w_tmp[0]*xB_tmp[135]-p[0]*w_tmp[0]*xB_tmp[136];
  xBdot_tmp[136] = p[1]*xB_tmp[1]*x_tmp[136]*4.0+p[1]*xB_tmp[136]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[136]*2.0-p[1]*xB_tmp[137]*x_tmp[1]*2.0;
  xBdot_tmp[137] = p[2]*xB_tmp[137]-k[0]*p[2]*w_tmp[3]*xB_tmp[138];
  xBdot_tmp[138] = p[3]*xB_tmp[138]-p[3]*xB_tmp[139]*2.0;
  xBdot_tmp[139] = p[3]*xB_tmp[139]-p[3]*xB_tmp[140];
  xBdot_tmp[140] = p[3]*xB_tmp[140]-p[3]*xB_tmp[141];
  xBdot_tmp[141] = p[3]*xB_tmp[141]-p[3]*xB_tmp[142];
  xBdot_tmp[142] = p[3]*xB_tmp[142]-p[3]*xB_tmp[143];
  xBdot_tmp[143] = p[3]*xB_tmp[143]-k[1]*p[3]*w_tmp[2]*xB_tmp[135];
  xBdot_tmp[144] = p[0]*w_tmp[0]*xB_tmp[144]-p[0]*w_tmp[0]*xB_tmp[145];
  xBdot_tmp[145] = p[1]*xB_tmp[1]*x_tmp[145]*4.0+p[1]*xB_tmp[145]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[145]*2.0-p[1]*xB_tmp[146]*x_tmp[1]*2.0;
  xBdot_tmp[146] = p[2]*xB_tmp[146]-k[0]*p[2]*w_tmp[3]*xB_tmp[147];
  xBdot_tmp[147] = p[3]*xB_tmp[147]-p[3]*xB_tmp[148]*2.0;
  xBdot_tmp[148] = p[3]*xB_tmp[148]-p[3]*xB_tmp[149];
  xBdot_tmp[149] = p[3]*xB_tmp[149]-p[3]*xB_tmp[150];
  xBdot_tmp[150] = p[3]*xB_tmp[150]-p[3]*xB_tmp[151];
  xBdot_tmp[151] = p[3]*xB_tmp[151]-p[3]*xB_tmp[152];
  xBdot_tmp[152] = p[3]*xB_tmp[152]-k[1]*p[3]*w_tmp[2]*xB_tmp[144];
  xBdot_tmp[153] = p[0]*w_tmp[0]*xB_tmp[153]-p[0]*w_tmp[0]*xB_tmp[154];
  xBdot_tmp[154] = p[1]*xB_tmp[1]*x_tmp[154]*4.0+p[1]*xB_tmp[154]*x_tmp[1]*4.0-p[1]*xB_tmp[2]*x_tmp[154]*2.0-p[1]*xB_tmp[155]*x_tmp[1]*2.0;
  xBdot_tmp[155] = p[2]*xB_tmp[155]-k[0]*p[2]*w_tmp[3]*xB_tmp[156];
  xBdot_tmp[156] = p[3]*xB_tmp[156]-p[3]*xB_tmp[157]*2.0;
  xBdot_tmp[157] = p[3]*xB_tmp[157]-p[3]*xB_tmp[158];
  xBdot_tmp[158] = p[3]*xB_tmp[158]-p[3]*xB_tmp[159];
  xBdot_tmp[159] = p[3]*xB_tmp[159]-p[3]*xB_tmp[160];
  xBdot_tmp[160] = p[3]*xB_tmp[160]-p[3]*xB_tmp[161];
  xBdot_tmp[161] = p[3]*xB_tmp[161]-k[1]*p[3]*w_tmp[2]*xB_tmp[153];
for(ix = 0; ix<162; ix++) {
   if(amiIsNaN(xBdot_tmp[ix])) {
       xBdot_tmp[ix] = 0;       if(!udata->am_nan_xBdot) {
           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_xBdot = TRUE;
       }
   }   if(amiIsInf(xBdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


