
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

void xdot_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(xdot_tmp,0,sizeof(realtype)*162);
w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
  xdot_tmp[0] = tdata->w[2]*(udata->k[1]*tdata->p[3]*x_tmp[8]-udata->k[0]*tdata->p[0]*tdata->w[0]*x_tmp[0]);
  xdot_tmp[1] = tdata->p[1]*tdata->w[1]*-2.0+tdata->p[0]*tdata->w[0]*x_tmp[0];
  xdot_tmp[2] = tdata->p[1]*tdata->w[1]-tdata->p[2]*x_tmp[2];
  xdot_tmp[3] = tdata->w[3]*(udata->k[0]*tdata->p[2]*x_tmp[2]-udata->k[1]*tdata->p[3]*x_tmp[3]);
  xdot_tmp[4] = tdata->p[3]*(tdata->w[4]-x_tmp[4]);
  xdot_tmp[5] = tdata->p[3]*(x_tmp[4]-x_tmp[5]);
  xdot_tmp[6] = tdata->p[3]*(x_tmp[5]-x_tmp[6]);
  xdot_tmp[7] = tdata->p[3]*(x_tmp[6]-x_tmp[7]);
  xdot_tmp[8] = tdata->p[3]*(x_tmp[7]-x_tmp[8]);
  xdot_tmp[9] = -tdata->w[0]*x_tmp[0]-tdata->p[0]*tdata->w[0]*x_tmp[9]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[17];
  xdot_tmp[10] = tdata->w[0]*x_tmp[0]+tdata->p[0]*tdata->w[0]*x_tmp[9]-tdata->p[1]*x_tmp[1]*x_tmp[10]*4.0;
  xdot_tmp[11] = -tdata->p[2]*x_tmp[11]+tdata->p[1]*x_tmp[1]*x_tmp[10]*2.0;
  xdot_tmp[12] = -tdata->p[3]*x_tmp[12]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[11];
  xdot_tmp[13] = tdata->p[3]*x_tmp[12]*2.0-tdata->p[3]*x_tmp[13];
  xdot_tmp[14] = tdata->p[3]*x_tmp[13]-tdata->p[3]*x_tmp[14];
  xdot_tmp[15] = tdata->p[3]*x_tmp[14]-tdata->p[3]*x_tmp[15];
  xdot_tmp[16] = tdata->p[3]*x_tmp[15]-tdata->p[3]*x_tmp[16];
  xdot_tmp[17] = tdata->p[3]*x_tmp[16]-tdata->p[3]*x_tmp[17];
  xdot_tmp[18] = -tdata->p[0]*tdata->w[0]*x_tmp[18]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[26];
  xdot_tmp[19] = tdata->w[1]*-2.0+tdata->p[0]*tdata->w[0]*x_tmp[18]-tdata->p[1]*x_tmp[1]*x_tmp[19]*4.0;
  xdot_tmp[20] = tdata->w[1]-tdata->p[2]*x_tmp[20]+tdata->p[1]*x_tmp[1]*x_tmp[19]*2.0;
  xdot_tmp[21] = -tdata->p[3]*x_tmp[21]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[20];
  xdot_tmp[22] = tdata->p[3]*x_tmp[21]*2.0-tdata->p[3]*x_tmp[22];
  xdot_tmp[23] = tdata->p[3]*x_tmp[22]-tdata->p[3]*x_tmp[23];
  xdot_tmp[24] = tdata->p[3]*x_tmp[23]-tdata->p[3]*x_tmp[24];
  xdot_tmp[25] = tdata->p[3]*x_tmp[24]-tdata->p[3]*x_tmp[25];
  xdot_tmp[26] = tdata->p[3]*x_tmp[25]-tdata->p[3]*x_tmp[26];
  xdot_tmp[27] = -tdata->p[0]*tdata->w[0]*x_tmp[27]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[35];
  xdot_tmp[28] = tdata->p[0]*tdata->w[0]*x_tmp[27]-tdata->p[1]*x_tmp[1]*x_tmp[28]*4.0;
  xdot_tmp[29] = -x_tmp[2]-tdata->p[2]*x_tmp[29]+tdata->p[1]*x_tmp[1]*x_tmp[28]*2.0;
  xdot_tmp[30] = -tdata->p[3]*x_tmp[30]+udata->k[0]*tdata->w[3]*x_tmp[2]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[29];
  xdot_tmp[31] = tdata->p[3]*x_tmp[30]*2.0-tdata->p[3]*x_tmp[31];
  xdot_tmp[32] = tdata->p[3]*x_tmp[31]-tdata->p[3]*x_tmp[32];
  xdot_tmp[33] = tdata->p[3]*x_tmp[32]-tdata->p[3]*x_tmp[33];
  xdot_tmp[34] = tdata->p[3]*x_tmp[33]-tdata->p[3]*x_tmp[34];
  xdot_tmp[35] = tdata->p[3]*x_tmp[34]-tdata->p[3]*x_tmp[35];
  xdot_tmp[36] = udata->k[1]*tdata->w[2]*x_tmp[8]-tdata->p[0]*tdata->w[0]*x_tmp[36]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[44];
  xdot_tmp[37] = tdata->p[0]*tdata->w[0]*x_tmp[36]-tdata->p[1]*x_tmp[1]*x_tmp[37]*4.0;
  xdot_tmp[38] = -tdata->p[2]*x_tmp[38]+tdata->p[1]*x_tmp[1]*x_tmp[37]*2.0;
  xdot_tmp[39] = -x_tmp[3]-tdata->p[3]*x_tmp[39]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[38];
  xdot_tmp[40] = tdata->w[4]-x_tmp[4]+tdata->p[3]*x_tmp[39]*2.0-tdata->p[3]*x_tmp[40];
  xdot_tmp[41] = x_tmp[4]-x_tmp[5]+tdata->p[3]*x_tmp[40]-tdata->p[3]*x_tmp[41];
  xdot_tmp[42] = x_tmp[5]-x_tmp[6]+tdata->p[3]*x_tmp[41]-tdata->p[3]*x_tmp[42];
  xdot_tmp[43] = x_tmp[6]-x_tmp[7]+tdata->p[3]*x_tmp[42]-tdata->p[3]*x_tmp[43];
  xdot_tmp[44] = x_tmp[7]-x_tmp[8]+tdata->p[3]*x_tmp[43]-tdata->p[3]*x_tmp[44];
  xdot_tmp[45] = -tdata->p[0]*tdata->w[0]*x_tmp[45]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[53];
  xdot_tmp[46] = tdata->p[0]*tdata->w[0]*x_tmp[45]-tdata->p[1]*x_tmp[1]*x_tmp[46]*4.0;
  xdot_tmp[47] = -tdata->p[2]*x_tmp[47]+tdata->p[1]*x_tmp[1]*x_tmp[46]*2.0;
  xdot_tmp[48] = -tdata->p[3]*x_tmp[48]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[47];
  xdot_tmp[49] = tdata->p[3]*x_tmp[48]*2.0-tdata->p[3]*x_tmp[49];
  xdot_tmp[50] = tdata->p[3]*x_tmp[49]-tdata->p[3]*x_tmp[50];
  xdot_tmp[51] = tdata->p[3]*x_tmp[50]-tdata->p[3]*x_tmp[51];
  xdot_tmp[52] = tdata->p[3]*x_tmp[51]-tdata->p[3]*x_tmp[52];
  xdot_tmp[53] = tdata->p[3]*x_tmp[52]-tdata->p[3]*x_tmp[53];
  xdot_tmp[54] = -tdata->p[0]*x_tmp[0]*tdata->w[5]-tdata->p[0]*tdata->w[0]*x_tmp[54]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[62];
  xdot_tmp[55] = tdata->p[0]*x_tmp[0]*tdata->w[5]+tdata->p[0]*tdata->w[0]*x_tmp[54]-tdata->p[1]*x_tmp[1]*x_tmp[55]*4.0;
  xdot_tmp[56] = -tdata->p[2]*x_tmp[56]+tdata->p[1]*x_tmp[1]*x_tmp[55]*2.0;
  xdot_tmp[57] = -tdata->p[3]*x_tmp[57]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[56];
  xdot_tmp[58] = tdata->p[3]*x_tmp[57]*2.0-tdata->p[3]*x_tmp[58];
  xdot_tmp[59] = tdata->p[3]*x_tmp[58]-tdata->p[3]*x_tmp[59];
  xdot_tmp[60] = tdata->p[3]*x_tmp[59]-tdata->p[3]*x_tmp[60];
  xdot_tmp[61] = tdata->p[3]*x_tmp[60]-tdata->p[3]*x_tmp[61];
  xdot_tmp[62] = tdata->p[3]*x_tmp[61]-tdata->p[3]*x_tmp[62];
  xdot_tmp[63] = -tdata->p[0]*x_tmp[0]*tdata->w[6]-tdata->p[0]*tdata->w[0]*x_tmp[63]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[71];
  xdot_tmp[64] = tdata->p[0]*x_tmp[0]*tdata->w[6]+tdata->p[0]*tdata->w[0]*x_tmp[63]-tdata->p[1]*x_tmp[1]*x_tmp[64]*4.0;
  xdot_tmp[65] = -tdata->p[2]*x_tmp[65]+tdata->p[1]*x_tmp[1]*x_tmp[64]*2.0;
  xdot_tmp[66] = -tdata->p[3]*x_tmp[66]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[65];
  xdot_tmp[67] = tdata->p[3]*x_tmp[66]*2.0-tdata->p[3]*x_tmp[67];
  xdot_tmp[68] = tdata->p[3]*x_tmp[67]-tdata->p[3]*x_tmp[68];
  xdot_tmp[69] = tdata->p[3]*x_tmp[68]-tdata->p[3]*x_tmp[69];
  xdot_tmp[70] = tdata->p[3]*x_tmp[69]-tdata->p[3]*x_tmp[70];
  xdot_tmp[71] = tdata->p[3]*x_tmp[70]-tdata->p[3]*x_tmp[71];
  xdot_tmp[72] = -tdata->p[0]*x_tmp[0]*tdata->w[7]-tdata->p[0]*tdata->w[0]*x_tmp[72]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[80];
  xdot_tmp[73] = tdata->p[0]*x_tmp[0]*tdata->w[7]+tdata->p[0]*tdata->w[0]*x_tmp[72]-tdata->p[1]*x_tmp[1]*x_tmp[73]*4.0;
  xdot_tmp[74] = -tdata->p[2]*x_tmp[74]+tdata->p[1]*x_tmp[1]*x_tmp[73]*2.0;
  xdot_tmp[75] = -tdata->p[3]*x_tmp[75]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[74];
  xdot_tmp[76] = tdata->p[3]*x_tmp[75]*2.0-tdata->p[3]*x_tmp[76];
  xdot_tmp[77] = tdata->p[3]*x_tmp[76]-tdata->p[3]*x_tmp[77];
  xdot_tmp[78] = tdata->p[3]*x_tmp[77]-tdata->p[3]*x_tmp[78];
  xdot_tmp[79] = tdata->p[3]*x_tmp[78]-tdata->p[3]*x_tmp[79];
  xdot_tmp[80] = tdata->p[3]*x_tmp[79]-tdata->p[3]*x_tmp[80];
  xdot_tmp[81] = -tdata->p[0]*x_tmp[0]*tdata->w[8]-tdata->p[0]*tdata->w[0]*x_tmp[81]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[89];
  xdot_tmp[82] = tdata->p[0]*x_tmp[0]*tdata->w[8]+tdata->p[0]*tdata->w[0]*x_tmp[81]-tdata->p[1]*x_tmp[1]*x_tmp[82]*4.0;
  xdot_tmp[83] = -tdata->p[2]*x_tmp[83]+tdata->p[1]*x_tmp[1]*x_tmp[82]*2.0;
  xdot_tmp[84] = -tdata->p[3]*x_tmp[84]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[83];
  xdot_tmp[85] = tdata->p[3]*x_tmp[84]*2.0-tdata->p[3]*x_tmp[85];
  xdot_tmp[86] = tdata->p[3]*x_tmp[85]-tdata->p[3]*x_tmp[86];
  xdot_tmp[87] = tdata->p[3]*x_tmp[86]-tdata->p[3]*x_tmp[87];
  xdot_tmp[88] = tdata->p[3]*x_tmp[87]-tdata->p[3]*x_tmp[88];
  xdot_tmp[89] = tdata->p[3]*x_tmp[88]-tdata->p[3]*x_tmp[89];
  xdot_tmp[90] = -tdata->p[0]*x_tmp[0]*tdata->w[9]-tdata->p[0]*tdata->w[0]*x_tmp[90]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[98];
  xdot_tmp[91] = tdata->p[0]*x_tmp[0]*tdata->w[9]+tdata->p[0]*tdata->w[0]*x_tmp[90]-tdata->p[1]*x_tmp[1]*x_tmp[91]*4.0;
  xdot_tmp[92] = -tdata->p[2]*x_tmp[92]+tdata->p[1]*x_tmp[1]*x_tmp[91]*2.0;
  xdot_tmp[93] = -tdata->p[3]*x_tmp[93]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[92];
  xdot_tmp[94] = tdata->p[3]*x_tmp[93]*2.0-tdata->p[3]*x_tmp[94];
  xdot_tmp[95] = tdata->p[3]*x_tmp[94]-tdata->p[3]*x_tmp[95];
  xdot_tmp[96] = tdata->p[3]*x_tmp[95]-tdata->p[3]*x_tmp[96];
  xdot_tmp[97] = tdata->p[3]*x_tmp[96]-tdata->p[3]*x_tmp[97];
  xdot_tmp[98] = tdata->p[3]*x_tmp[97]-tdata->p[3]*x_tmp[98];
  xdot_tmp[99] = -tdata->p[0]*tdata->w[0]*x_tmp[99]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[107];
  xdot_tmp[100] = tdata->p[0]*tdata->w[0]*x_tmp[99]-tdata->p[1]*x_tmp[1]*x_tmp[100]*4.0;
  xdot_tmp[101] = -tdata->p[2]*x_tmp[101]+tdata->p[1]*x_tmp[1]*x_tmp[100]*2.0;
  xdot_tmp[102] = -tdata->p[3]*x_tmp[102]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[101];
  xdot_tmp[103] = tdata->p[3]*x_tmp[102]*2.0-tdata->p[3]*x_tmp[103];
  xdot_tmp[104] = tdata->p[3]*x_tmp[103]-tdata->p[3]*x_tmp[104];
  xdot_tmp[105] = tdata->p[3]*x_tmp[104]-tdata->p[3]*x_tmp[105];
  xdot_tmp[106] = tdata->p[3]*x_tmp[105]-tdata->p[3]*x_tmp[106];
  xdot_tmp[107] = tdata->p[3]*x_tmp[106]-tdata->p[3]*x_tmp[107];
  xdot_tmp[108] = -tdata->p[0]*tdata->w[0]*x_tmp[108]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[116];
  xdot_tmp[109] = tdata->p[0]*tdata->w[0]*x_tmp[108]-tdata->p[1]*x_tmp[1]*x_tmp[109]*4.0;
  xdot_tmp[110] = -tdata->p[2]*x_tmp[110]+tdata->p[1]*x_tmp[1]*x_tmp[109]*2.0;
  xdot_tmp[111] = -tdata->p[3]*x_tmp[111]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[110];
  xdot_tmp[112] = tdata->p[3]*x_tmp[111]*2.0-tdata->p[3]*x_tmp[112];
  xdot_tmp[113] = tdata->p[3]*x_tmp[112]-tdata->p[3]*x_tmp[113];
  xdot_tmp[114] = tdata->p[3]*x_tmp[113]-tdata->p[3]*x_tmp[114];
  xdot_tmp[115] = tdata->p[3]*x_tmp[114]-tdata->p[3]*x_tmp[115];
  xdot_tmp[116] = tdata->p[3]*x_tmp[115]-tdata->p[3]*x_tmp[116];
  xdot_tmp[117] = -tdata->p[0]*tdata->w[0]*x_tmp[117]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[125];
  xdot_tmp[118] = tdata->p[0]*tdata->w[0]*x_tmp[117]-tdata->p[1]*x_tmp[1]*x_tmp[118]*4.0;
  xdot_tmp[119] = -tdata->p[2]*x_tmp[119]+tdata->p[1]*x_tmp[1]*x_tmp[118]*2.0;
  xdot_tmp[120] = -tdata->p[3]*x_tmp[120]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[119];
  xdot_tmp[121] = tdata->p[3]*x_tmp[120]*2.0-tdata->p[3]*x_tmp[121];
  xdot_tmp[122] = tdata->p[3]*x_tmp[121]-tdata->p[3]*x_tmp[122];
  xdot_tmp[123] = tdata->p[3]*x_tmp[122]-tdata->p[3]*x_tmp[123];
  xdot_tmp[124] = tdata->p[3]*x_tmp[123]-tdata->p[3]*x_tmp[124];
  xdot_tmp[125] = tdata->p[3]*x_tmp[124]-tdata->p[3]*x_tmp[125];
  xdot_tmp[126] = -tdata->p[0]*tdata->w[0]*x_tmp[126]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[134];
  xdot_tmp[127] = tdata->p[0]*tdata->w[0]*x_tmp[126]-tdata->p[1]*x_tmp[1]*x_tmp[127]*4.0;
  xdot_tmp[128] = -tdata->p[2]*x_tmp[128]+tdata->p[1]*x_tmp[1]*x_tmp[127]*2.0;
  xdot_tmp[129] = -tdata->p[3]*x_tmp[129]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[128];
  xdot_tmp[130] = tdata->p[3]*x_tmp[129]*2.0-tdata->p[3]*x_tmp[130];
  xdot_tmp[131] = tdata->p[3]*x_tmp[130]-tdata->p[3]*x_tmp[131];
  xdot_tmp[132] = tdata->p[3]*x_tmp[131]-tdata->p[3]*x_tmp[132];
  xdot_tmp[133] = tdata->p[3]*x_tmp[132]-tdata->p[3]*x_tmp[133];
  xdot_tmp[134] = tdata->p[3]*x_tmp[133]-tdata->p[3]*x_tmp[134];
  xdot_tmp[135] = -tdata->p[0]*tdata->w[0]*x_tmp[135]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[143];
  xdot_tmp[136] = tdata->p[0]*tdata->w[0]*x_tmp[135]-tdata->p[1]*x_tmp[1]*x_tmp[136]*4.0;
  xdot_tmp[137] = -tdata->p[2]*x_tmp[137]+tdata->p[1]*x_tmp[1]*x_tmp[136]*2.0;
  xdot_tmp[138] = -tdata->p[3]*x_tmp[138]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[137];
  xdot_tmp[139] = tdata->p[3]*x_tmp[138]*2.0-tdata->p[3]*x_tmp[139];
  xdot_tmp[140] = tdata->p[3]*x_tmp[139]-tdata->p[3]*x_tmp[140];
  xdot_tmp[141] = tdata->p[3]*x_tmp[140]-tdata->p[3]*x_tmp[141];
  xdot_tmp[142] = tdata->p[3]*x_tmp[141]-tdata->p[3]*x_tmp[142];
  xdot_tmp[143] = tdata->p[3]*x_tmp[142]-tdata->p[3]*x_tmp[143];
  xdot_tmp[144] = -tdata->p[0]*tdata->w[0]*x_tmp[144]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[152];
  xdot_tmp[145] = tdata->p[0]*tdata->w[0]*x_tmp[144]-tdata->p[1]*x_tmp[1]*x_tmp[145]*4.0;
  xdot_tmp[146] = -tdata->p[2]*x_tmp[146]+tdata->p[1]*x_tmp[1]*x_tmp[145]*2.0;
  xdot_tmp[147] = -tdata->p[3]*x_tmp[147]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[146];
  xdot_tmp[148] = tdata->p[3]*x_tmp[147]*2.0-tdata->p[3]*x_tmp[148];
  xdot_tmp[149] = tdata->p[3]*x_tmp[148]-tdata->p[3]*x_tmp[149];
  xdot_tmp[150] = tdata->p[3]*x_tmp[149]-tdata->p[3]*x_tmp[150];
  xdot_tmp[151] = tdata->p[3]*x_tmp[150]-tdata->p[3]*x_tmp[151];
  xdot_tmp[152] = tdata->p[3]*x_tmp[151]-tdata->p[3]*x_tmp[152];
  xdot_tmp[153] = -tdata->p[0]*tdata->w[0]*x_tmp[153]+udata->k[1]*tdata->p[3]*tdata->w[2]*x_tmp[161];
  xdot_tmp[154] = tdata->p[0]*tdata->w[0]*x_tmp[153]-tdata->p[1]*x_tmp[1]*x_tmp[154]*4.0;
  xdot_tmp[155] = -tdata->p[2]*x_tmp[155]+tdata->p[1]*x_tmp[1]*x_tmp[154]*2.0;
  xdot_tmp[156] = -tdata->p[3]*x_tmp[156]+udata->k[0]*tdata->p[2]*tdata->w[3]*x_tmp[155];
  xdot_tmp[157] = tdata->p[3]*x_tmp[156]*2.0-tdata->p[3]*x_tmp[157];
  xdot_tmp[158] = tdata->p[3]*x_tmp[157]-tdata->p[3]*x_tmp[158];
  xdot_tmp[159] = tdata->p[3]*x_tmp[158]-tdata->p[3]*x_tmp[159];
  xdot_tmp[160] = tdata->p[3]*x_tmp[159]-tdata->p[3]*x_tmp[160];
  xdot_tmp[161] = tdata->p[3]*x_tmp[160]-tdata->p[3]*x_tmp[161];
for(ix = 0; ix<162; ix++) {
   if(amiIsNaN(xdot_tmp[ix])) {
       xdot_tmp[ix] = 0;
       if(!tdata->nan_xdot) {
           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_xdot = TRUE;
       }
   }
   if(amiIsInf(xdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");
       return;
   }   if(udata->qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {
       xdot_tmp[ix] = -xdot_tmp[ix];
   }
}
return;

}


