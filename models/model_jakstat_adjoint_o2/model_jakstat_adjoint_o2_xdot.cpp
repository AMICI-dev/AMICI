
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int xdot_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(xdot_tmp,0,sizeof(realtype)*162);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  xdot_tmp[0] = w_tmp[2]*(k[1]*p[3]*x_tmp[8]-k[0]*p[0]*w_tmp[0]*x_tmp[0]);
  xdot_tmp[1] = p[1]*w_tmp[1]*-2.0+p[0]*w_tmp[0]*x_tmp[0];
  xdot_tmp[2] = p[1]*w_tmp[1]-p[2]*x_tmp[2];
  xdot_tmp[3] = w_tmp[3]*(k[0]*p[2]*x_tmp[2]-k[1]*p[3]*x_tmp[3]);
  xdot_tmp[4] = p[3]*(w_tmp[4]-x_tmp[4]);
  xdot_tmp[5] = p[3]*(x_tmp[4]-x_tmp[5]);
  xdot_tmp[6] = p[3]*(x_tmp[5]-x_tmp[6]);
  xdot_tmp[7] = p[3]*(x_tmp[6]-x_tmp[7]);
  xdot_tmp[8] = p[3]*(x_tmp[7]-x_tmp[8]);
  xdot_tmp[9] = -w_tmp[0]*x_tmp[0]-p[0]*w_tmp[0]*x_tmp[9]+k[1]*p[3]*w_tmp[2]*x_tmp[17];
  xdot_tmp[10] = w_tmp[0]*x_tmp[0]+p[0]*w_tmp[0]*x_tmp[9]-p[1]*x_tmp[1]*x_tmp[10]*4.0;
  xdot_tmp[11] = -p[2]*x_tmp[11]+p[1]*x_tmp[1]*x_tmp[10]*2.0;
  xdot_tmp[12] = -p[3]*x_tmp[12]+k[0]*p[2]*w_tmp[3]*x_tmp[11];
  xdot_tmp[13] = p[3]*x_tmp[12]*2.0-p[3]*x_tmp[13];
  xdot_tmp[14] = p[3]*x_tmp[13]-p[3]*x_tmp[14];
  xdot_tmp[15] = p[3]*x_tmp[14]-p[3]*x_tmp[15];
  xdot_tmp[16] = p[3]*x_tmp[15]-p[3]*x_tmp[16];
  xdot_tmp[17] = p[3]*x_tmp[16]-p[3]*x_tmp[17];
  xdot_tmp[18] = -p[0]*w_tmp[0]*x_tmp[18]+k[1]*p[3]*w_tmp[2]*x_tmp[26];
  xdot_tmp[19] = w_tmp[1]*-2.0+p[0]*w_tmp[0]*x_tmp[18]-p[1]*x_tmp[1]*x_tmp[19]*4.0;
  xdot_tmp[20] = w_tmp[1]-p[2]*x_tmp[20]+p[1]*x_tmp[1]*x_tmp[19]*2.0;
  xdot_tmp[21] = -p[3]*x_tmp[21]+k[0]*p[2]*w_tmp[3]*x_tmp[20];
  xdot_tmp[22] = p[3]*x_tmp[21]*2.0-p[3]*x_tmp[22];
  xdot_tmp[23] = p[3]*x_tmp[22]-p[3]*x_tmp[23];
  xdot_tmp[24] = p[3]*x_tmp[23]-p[3]*x_tmp[24];
  xdot_tmp[25] = p[3]*x_tmp[24]-p[3]*x_tmp[25];
  xdot_tmp[26] = p[3]*x_tmp[25]-p[3]*x_tmp[26];
  xdot_tmp[27] = -p[0]*w_tmp[0]*x_tmp[27]+k[1]*p[3]*w_tmp[2]*x_tmp[35];
  xdot_tmp[28] = p[0]*w_tmp[0]*x_tmp[27]-p[1]*x_tmp[1]*x_tmp[28]*4.0;
  xdot_tmp[29] = -x_tmp[2]-p[2]*x_tmp[29]+p[1]*x_tmp[1]*x_tmp[28]*2.0;
  xdot_tmp[30] = -p[3]*x_tmp[30]+k[0]*w_tmp[3]*x_tmp[2]+k[0]*p[2]*w_tmp[3]*x_tmp[29];
  xdot_tmp[31] = p[3]*x_tmp[30]*2.0-p[3]*x_tmp[31];
  xdot_tmp[32] = p[3]*x_tmp[31]-p[3]*x_tmp[32];
  xdot_tmp[33] = p[3]*x_tmp[32]-p[3]*x_tmp[33];
  xdot_tmp[34] = p[3]*x_tmp[33]-p[3]*x_tmp[34];
  xdot_tmp[35] = p[3]*x_tmp[34]-p[3]*x_tmp[35];
  xdot_tmp[36] = k[1]*w_tmp[2]*x_tmp[8]-p[0]*w_tmp[0]*x_tmp[36]+k[1]*p[3]*w_tmp[2]*x_tmp[44];
  xdot_tmp[37] = p[0]*w_tmp[0]*x_tmp[36]-p[1]*x_tmp[1]*x_tmp[37]*4.0;
  xdot_tmp[38] = -p[2]*x_tmp[38]+p[1]*x_tmp[1]*x_tmp[37]*2.0;
  xdot_tmp[39] = -x_tmp[3]-p[3]*x_tmp[39]+k[0]*p[2]*w_tmp[3]*x_tmp[38];
  xdot_tmp[40] = w_tmp[4]-x_tmp[4]+p[3]*x_tmp[39]*2.0-p[3]*x_tmp[40];
  xdot_tmp[41] = x_tmp[4]-x_tmp[5]+p[3]*x_tmp[40]-p[3]*x_tmp[41];
  xdot_tmp[42] = x_tmp[5]-x_tmp[6]+p[3]*x_tmp[41]-p[3]*x_tmp[42];
  xdot_tmp[43] = x_tmp[6]-x_tmp[7]+p[3]*x_tmp[42]-p[3]*x_tmp[43];
  xdot_tmp[44] = x_tmp[7]-x_tmp[8]+p[3]*x_tmp[43]-p[3]*x_tmp[44];
  xdot_tmp[45] = -p[0]*w_tmp[0]*x_tmp[45]+k[1]*p[3]*w_tmp[2]*x_tmp[53];
  xdot_tmp[46] = p[0]*w_tmp[0]*x_tmp[45]-p[1]*x_tmp[1]*x_tmp[46]*4.0;
  xdot_tmp[47] = -p[2]*x_tmp[47]+p[1]*x_tmp[1]*x_tmp[46]*2.0;
  xdot_tmp[48] = -p[3]*x_tmp[48]+k[0]*p[2]*w_tmp[3]*x_tmp[47];
  xdot_tmp[49] = p[3]*x_tmp[48]*2.0-p[3]*x_tmp[49];
  xdot_tmp[50] = p[3]*x_tmp[49]-p[3]*x_tmp[50];
  xdot_tmp[51] = p[3]*x_tmp[50]-p[3]*x_tmp[51];
  xdot_tmp[52] = p[3]*x_tmp[51]-p[3]*x_tmp[52];
  xdot_tmp[53] = p[3]*x_tmp[52]-p[3]*x_tmp[53];
  xdot_tmp[54] = -p[0]*w_tmp[5]*x_tmp[0]-p[0]*w_tmp[0]*x_tmp[54]+k[1]*p[3]*w_tmp[2]*x_tmp[62];
  xdot_tmp[55] = p[0]*w_tmp[5]*x_tmp[0]+p[0]*w_tmp[0]*x_tmp[54]-p[1]*x_tmp[1]*x_tmp[55]*4.0;
  xdot_tmp[56] = -p[2]*x_tmp[56]+p[1]*x_tmp[1]*x_tmp[55]*2.0;
  xdot_tmp[57] = -p[3]*x_tmp[57]+k[0]*p[2]*w_tmp[3]*x_tmp[56];
  xdot_tmp[58] = p[3]*x_tmp[57]*2.0-p[3]*x_tmp[58];
  xdot_tmp[59] = p[3]*x_tmp[58]-p[3]*x_tmp[59];
  xdot_tmp[60] = p[3]*x_tmp[59]-p[3]*x_tmp[60];
  xdot_tmp[61] = p[3]*x_tmp[60]-p[3]*x_tmp[61];
  xdot_tmp[62] = p[3]*x_tmp[61]-p[3]*x_tmp[62];
  xdot_tmp[63] = -p[0]*w_tmp[6]*x_tmp[0]-p[0]*w_tmp[0]*x_tmp[63]+k[1]*p[3]*w_tmp[2]*x_tmp[71];
  xdot_tmp[64] = p[0]*w_tmp[6]*x_tmp[0]+p[0]*w_tmp[0]*x_tmp[63]-p[1]*x_tmp[1]*x_tmp[64]*4.0;
  xdot_tmp[65] = -p[2]*x_tmp[65]+p[1]*x_tmp[1]*x_tmp[64]*2.0;
  xdot_tmp[66] = -p[3]*x_tmp[66]+k[0]*p[2]*w_tmp[3]*x_tmp[65];
  xdot_tmp[67] = p[3]*x_tmp[66]*2.0-p[3]*x_tmp[67];
  xdot_tmp[68] = p[3]*x_tmp[67]-p[3]*x_tmp[68];
  xdot_tmp[69] = p[3]*x_tmp[68]-p[3]*x_tmp[69];
  xdot_tmp[70] = p[3]*x_tmp[69]-p[3]*x_tmp[70];
  xdot_tmp[71] = p[3]*x_tmp[70]-p[3]*x_tmp[71];
  xdot_tmp[72] = -p[0]*w_tmp[7]*x_tmp[0]-p[0]*w_tmp[0]*x_tmp[72]+k[1]*p[3]*w_tmp[2]*x_tmp[80];
  xdot_tmp[73] = p[0]*w_tmp[7]*x_tmp[0]+p[0]*w_tmp[0]*x_tmp[72]-p[1]*x_tmp[1]*x_tmp[73]*4.0;
  xdot_tmp[74] = -p[2]*x_tmp[74]+p[1]*x_tmp[1]*x_tmp[73]*2.0;
  xdot_tmp[75] = -p[3]*x_tmp[75]+k[0]*p[2]*w_tmp[3]*x_tmp[74];
  xdot_tmp[76] = p[3]*x_tmp[75]*2.0-p[3]*x_tmp[76];
  xdot_tmp[77] = p[3]*x_tmp[76]-p[3]*x_tmp[77];
  xdot_tmp[78] = p[3]*x_tmp[77]-p[3]*x_tmp[78];
  xdot_tmp[79] = p[3]*x_tmp[78]-p[3]*x_tmp[79];
  xdot_tmp[80] = p[3]*x_tmp[79]-p[3]*x_tmp[80];
  xdot_tmp[81] = -p[0]*w_tmp[8]*x_tmp[0]-p[0]*w_tmp[0]*x_tmp[81]+k[1]*p[3]*w_tmp[2]*x_tmp[89];
  xdot_tmp[82] = p[0]*w_tmp[8]*x_tmp[0]+p[0]*w_tmp[0]*x_tmp[81]-p[1]*x_tmp[1]*x_tmp[82]*4.0;
  xdot_tmp[83] = -p[2]*x_tmp[83]+p[1]*x_tmp[1]*x_tmp[82]*2.0;
  xdot_tmp[84] = -p[3]*x_tmp[84]+k[0]*p[2]*w_tmp[3]*x_tmp[83];
  xdot_tmp[85] = p[3]*x_tmp[84]*2.0-p[3]*x_tmp[85];
  xdot_tmp[86] = p[3]*x_tmp[85]-p[3]*x_tmp[86];
  xdot_tmp[87] = p[3]*x_tmp[86]-p[3]*x_tmp[87];
  xdot_tmp[88] = p[3]*x_tmp[87]-p[3]*x_tmp[88];
  xdot_tmp[89] = p[3]*x_tmp[88]-p[3]*x_tmp[89];
  xdot_tmp[90] = -p[0]*w_tmp[9]*x_tmp[0]-p[0]*w_tmp[0]*x_tmp[90]+k[1]*p[3]*w_tmp[2]*x_tmp[98];
  xdot_tmp[91] = p[0]*w_tmp[9]*x_tmp[0]+p[0]*w_tmp[0]*x_tmp[90]-p[1]*x_tmp[1]*x_tmp[91]*4.0;
  xdot_tmp[92] = -p[2]*x_tmp[92]+p[1]*x_tmp[1]*x_tmp[91]*2.0;
  xdot_tmp[93] = -p[3]*x_tmp[93]+k[0]*p[2]*w_tmp[3]*x_tmp[92];
  xdot_tmp[94] = p[3]*x_tmp[93]*2.0-p[3]*x_tmp[94];
  xdot_tmp[95] = p[3]*x_tmp[94]-p[3]*x_tmp[95];
  xdot_tmp[96] = p[3]*x_tmp[95]-p[3]*x_tmp[96];
  xdot_tmp[97] = p[3]*x_tmp[96]-p[3]*x_tmp[97];
  xdot_tmp[98] = p[3]*x_tmp[97]-p[3]*x_tmp[98];
  xdot_tmp[99] = -p[0]*w_tmp[0]*x_tmp[99]+k[1]*p[3]*w_tmp[2]*x_tmp[107];
  xdot_tmp[100] = p[0]*w_tmp[0]*x_tmp[99]-p[1]*x_tmp[1]*x_tmp[100]*4.0;
  xdot_tmp[101] = -p[2]*x_tmp[101]+p[1]*x_tmp[1]*x_tmp[100]*2.0;
  xdot_tmp[102] = -p[3]*x_tmp[102]+k[0]*p[2]*w_tmp[3]*x_tmp[101];
  xdot_tmp[103] = p[3]*x_tmp[102]*2.0-p[3]*x_tmp[103];
  xdot_tmp[104] = p[3]*x_tmp[103]-p[3]*x_tmp[104];
  xdot_tmp[105] = p[3]*x_tmp[104]-p[3]*x_tmp[105];
  xdot_tmp[106] = p[3]*x_tmp[105]-p[3]*x_tmp[106];
  xdot_tmp[107] = p[3]*x_tmp[106]-p[3]*x_tmp[107];
  xdot_tmp[108] = -p[0]*w_tmp[0]*x_tmp[108]+k[1]*p[3]*w_tmp[2]*x_tmp[116];
  xdot_tmp[109] = p[0]*w_tmp[0]*x_tmp[108]-p[1]*x_tmp[1]*x_tmp[109]*4.0;
  xdot_tmp[110] = -p[2]*x_tmp[110]+p[1]*x_tmp[1]*x_tmp[109]*2.0;
  xdot_tmp[111] = -p[3]*x_tmp[111]+k[0]*p[2]*w_tmp[3]*x_tmp[110];
  xdot_tmp[112] = p[3]*x_tmp[111]*2.0-p[3]*x_tmp[112];
  xdot_tmp[113] = p[3]*x_tmp[112]-p[3]*x_tmp[113];
  xdot_tmp[114] = p[3]*x_tmp[113]-p[3]*x_tmp[114];
  xdot_tmp[115] = p[3]*x_tmp[114]-p[3]*x_tmp[115];
  xdot_tmp[116] = p[3]*x_tmp[115]-p[3]*x_tmp[116];
  xdot_tmp[117] = -p[0]*w_tmp[0]*x_tmp[117]+k[1]*p[3]*w_tmp[2]*x_tmp[125];
  xdot_tmp[118] = p[0]*w_tmp[0]*x_tmp[117]-p[1]*x_tmp[1]*x_tmp[118]*4.0;
  xdot_tmp[119] = -p[2]*x_tmp[119]+p[1]*x_tmp[1]*x_tmp[118]*2.0;
  xdot_tmp[120] = -p[3]*x_tmp[120]+k[0]*p[2]*w_tmp[3]*x_tmp[119];
  xdot_tmp[121] = p[3]*x_tmp[120]*2.0-p[3]*x_tmp[121];
  xdot_tmp[122] = p[3]*x_tmp[121]-p[3]*x_tmp[122];
  xdot_tmp[123] = p[3]*x_tmp[122]-p[3]*x_tmp[123];
  xdot_tmp[124] = p[3]*x_tmp[123]-p[3]*x_tmp[124];
  xdot_tmp[125] = p[3]*x_tmp[124]-p[3]*x_tmp[125];
  xdot_tmp[126] = -p[0]*w_tmp[0]*x_tmp[126]+k[1]*p[3]*w_tmp[2]*x_tmp[134];
  xdot_tmp[127] = p[0]*w_tmp[0]*x_tmp[126]-p[1]*x_tmp[1]*x_tmp[127]*4.0;
  xdot_tmp[128] = -p[2]*x_tmp[128]+p[1]*x_tmp[1]*x_tmp[127]*2.0;
  xdot_tmp[129] = -p[3]*x_tmp[129]+k[0]*p[2]*w_tmp[3]*x_tmp[128];
  xdot_tmp[130] = p[3]*x_tmp[129]*2.0-p[3]*x_tmp[130];
  xdot_tmp[131] = p[3]*x_tmp[130]-p[3]*x_tmp[131];
  xdot_tmp[132] = p[3]*x_tmp[131]-p[3]*x_tmp[132];
  xdot_tmp[133] = p[3]*x_tmp[132]-p[3]*x_tmp[133];
  xdot_tmp[134] = p[3]*x_tmp[133]-p[3]*x_tmp[134];
  xdot_tmp[135] = -p[0]*w_tmp[0]*x_tmp[135]+k[1]*p[3]*w_tmp[2]*x_tmp[143];
  xdot_tmp[136] = p[0]*w_tmp[0]*x_tmp[135]-p[1]*x_tmp[1]*x_tmp[136]*4.0;
  xdot_tmp[137] = -p[2]*x_tmp[137]+p[1]*x_tmp[1]*x_tmp[136]*2.0;
  xdot_tmp[138] = -p[3]*x_tmp[138]+k[0]*p[2]*w_tmp[3]*x_tmp[137];
  xdot_tmp[139] = p[3]*x_tmp[138]*2.0-p[3]*x_tmp[139];
  xdot_tmp[140] = p[3]*x_tmp[139]-p[3]*x_tmp[140];
  xdot_tmp[141] = p[3]*x_tmp[140]-p[3]*x_tmp[141];
  xdot_tmp[142] = p[3]*x_tmp[141]-p[3]*x_tmp[142];
  xdot_tmp[143] = p[3]*x_tmp[142]-p[3]*x_tmp[143];
  xdot_tmp[144] = -p[0]*w_tmp[0]*x_tmp[144]+k[1]*p[3]*w_tmp[2]*x_tmp[152];
  xdot_tmp[145] = p[0]*w_tmp[0]*x_tmp[144]-p[1]*x_tmp[1]*x_tmp[145]*4.0;
  xdot_tmp[146] = -p[2]*x_tmp[146]+p[1]*x_tmp[1]*x_tmp[145]*2.0;
  xdot_tmp[147] = -p[3]*x_tmp[147]+k[0]*p[2]*w_tmp[3]*x_tmp[146];
  xdot_tmp[148] = p[3]*x_tmp[147]*2.0-p[3]*x_tmp[148];
  xdot_tmp[149] = p[3]*x_tmp[148]-p[3]*x_tmp[149];
  xdot_tmp[150] = p[3]*x_tmp[149]-p[3]*x_tmp[150];
  xdot_tmp[151] = p[3]*x_tmp[150]-p[3]*x_tmp[151];
  xdot_tmp[152] = p[3]*x_tmp[151]-p[3]*x_tmp[152];
  xdot_tmp[153] = -p[0]*w_tmp[0]*x_tmp[153]+k[1]*p[3]*w_tmp[2]*x_tmp[161];
  xdot_tmp[154] = p[0]*w_tmp[0]*x_tmp[153]-p[1]*x_tmp[1]*x_tmp[154]*4.0;
  xdot_tmp[155] = -p[2]*x_tmp[155]+p[1]*x_tmp[1]*x_tmp[154]*2.0;
  xdot_tmp[156] = -p[3]*x_tmp[156]+k[0]*p[2]*w_tmp[3]*x_tmp[155];
  xdot_tmp[157] = p[3]*x_tmp[156]*2.0-p[3]*x_tmp[157];
  xdot_tmp[158] = p[3]*x_tmp[157]-p[3]*x_tmp[158];
  xdot_tmp[159] = p[3]*x_tmp[158]-p[3]*x_tmp[159];
  xdot_tmp[160] = p[3]*x_tmp[159]-p[3]*x_tmp[160];
  xdot_tmp[161] = p[3]*x_tmp[160]-p[3]*x_tmp[161];
for(ix = 0; ix<162; ix++) {
   if(amiIsNaN(xdot_tmp[ix])) {
       xdot_tmp[ix] = 0;
       if(!udata->am_nan_xdot) {
           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_xdot = TRUE;
       }
   }
   if(amiIsInf(xdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");
       return(-1);
   }   if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {
       xdot_tmp[ix] = -xdot_tmp[ix];
   }
}
return(status);

}


