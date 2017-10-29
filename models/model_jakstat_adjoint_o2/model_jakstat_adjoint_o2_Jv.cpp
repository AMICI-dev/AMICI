
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

void Jv_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) {
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
realtype *v_tmp = nullptr;
if(v)
    v_tmp = N_VGetArrayPointer(v);
realtype *Jv_tmp = nullptr;
if(Jv)
    Jv_tmp = N_VGetArrayPointer(Jv);
memset(Jv_tmp,0,sizeof(realtype)*162);
w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
  Jv_tmp[0] = udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[8]-udata->k[0]*tdata->p[0]*v_tmp[0]*tdata->w[0]*tdata->w[2];
  Jv_tmp[1] = tdata->p[0]*v_tmp[0]*tdata->w[0]-tdata->p[1]*v_tmp[1]*tdata->dwdx[0]*2.0;
  Jv_tmp[2] = -tdata->p[2]*v_tmp[2]+tdata->p[1]*v_tmp[1]*tdata->dwdx[0];
  Jv_tmp[3] = udata->k[0]*tdata->p[2]*v_tmp[2]*tdata->w[3]-udata->k[1]*tdata->p[3]*v_tmp[3]*tdata->w[3];
  Jv_tmp[4] = -tdata->p[3]*v_tmp[4]+tdata->p[3]*v_tmp[3]*tdata->dwdx[1];
  Jv_tmp[5] = tdata->p[3]*v_tmp[4]-tdata->p[3]*v_tmp[5];
  Jv_tmp[6] = tdata->p[3]*v_tmp[5]-tdata->p[3]*v_tmp[6];
  Jv_tmp[7] = tdata->p[3]*v_tmp[6]-tdata->p[3]*v_tmp[7];
  Jv_tmp[8] = tdata->p[3]*v_tmp[7]-tdata->p[3]*v_tmp[8];
  Jv_tmp[9] = -v_tmp[0]*tdata->w[0]-tdata->p[0]*tdata->w[0]*v_tmp[9]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[17];
  Jv_tmp[10] = v_tmp[0]*tdata->w[0]+tdata->p[0]*tdata->w[0]*v_tmp[9]-tdata->p[1]*v_tmp[1]*x_tmp[10]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[10]*4.0;
  Jv_tmp[11] = -tdata->p[2]*v_tmp[11]+tdata->p[1]*v_tmp[1]*x_tmp[10]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[10]*2.0;
  Jv_tmp[12] = -tdata->p[3]*v_tmp[12]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[11];
  Jv_tmp[13] = tdata->p[3]*v_tmp[12]*2.0-tdata->p[3]*v_tmp[13];
  Jv_tmp[14] = tdata->p[3]*v_tmp[13]-tdata->p[3]*v_tmp[14];
  Jv_tmp[15] = tdata->p[3]*v_tmp[14]-tdata->p[3]*v_tmp[15];
  Jv_tmp[16] = tdata->p[3]*v_tmp[15]-tdata->p[3]*v_tmp[16];
  Jv_tmp[17] = tdata->p[3]*v_tmp[16]-tdata->p[3]*v_tmp[17];
  Jv_tmp[18] = -tdata->p[0]*tdata->w[0]*v_tmp[18]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[26];
  Jv_tmp[19] = -v_tmp[1]*(tdata->dwdx[0]*2.0+tdata->p[1]*x_tmp[19]*4.0)+tdata->p[0]*tdata->w[0]*v_tmp[18]-tdata->p[1]*x_tmp[1]*v_tmp[19]*4.0;
  Jv_tmp[20] = -tdata->p[2]*v_tmp[20]+v_tmp[1]*(tdata->dwdx[0]+tdata->p[1]*x_tmp[19]*2.0)+tdata->p[1]*x_tmp[1]*v_tmp[19]*2.0;
  Jv_tmp[21] = -tdata->p[3]*v_tmp[21]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[20];
  Jv_tmp[22] = tdata->p[3]*v_tmp[21]*2.0-tdata->p[3]*v_tmp[22];
  Jv_tmp[23] = tdata->p[3]*v_tmp[22]-tdata->p[3]*v_tmp[23];
  Jv_tmp[24] = tdata->p[3]*v_tmp[23]-tdata->p[3]*v_tmp[24];
  Jv_tmp[25] = tdata->p[3]*v_tmp[24]-tdata->p[3]*v_tmp[25];
  Jv_tmp[26] = tdata->p[3]*v_tmp[25]-tdata->p[3]*v_tmp[26];
  Jv_tmp[27] = -tdata->p[0]*tdata->w[0]*v_tmp[27]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[35];
  Jv_tmp[28] = tdata->p[0]*tdata->w[0]*v_tmp[27]-tdata->p[1]*v_tmp[1]*x_tmp[28]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[28]*4.0;
  Jv_tmp[29] = -v_tmp[2]-tdata->p[2]*v_tmp[29]+tdata->p[1]*v_tmp[1]*x_tmp[28]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[28]*2.0;
  Jv_tmp[30] = -tdata->p[3]*v_tmp[30]+udata->k[0]*v_tmp[2]*tdata->w[3]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[29];
  Jv_tmp[31] = tdata->p[3]*v_tmp[30]*2.0-tdata->p[3]*v_tmp[31];
  Jv_tmp[32] = tdata->p[3]*v_tmp[31]-tdata->p[3]*v_tmp[32];
  Jv_tmp[33] = tdata->p[3]*v_tmp[32]-tdata->p[3]*v_tmp[33];
  Jv_tmp[34] = tdata->p[3]*v_tmp[33]-tdata->p[3]*v_tmp[34];
  Jv_tmp[35] = tdata->p[3]*v_tmp[34]-tdata->p[3]*v_tmp[35];
  Jv_tmp[36] = udata->k[1]*tdata->w[2]*v_tmp[8]-tdata->p[0]*tdata->w[0]*v_tmp[36]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[44];
  Jv_tmp[37] = tdata->p[0]*tdata->w[0]*v_tmp[36]-tdata->p[1]*v_tmp[1]*x_tmp[37]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[37]*4.0;
  Jv_tmp[38] = -tdata->p[2]*v_tmp[38]+tdata->p[1]*v_tmp[1]*x_tmp[37]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[37]*2.0;
  Jv_tmp[39] = -v_tmp[3]-tdata->p[3]*v_tmp[39]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[38];
  Jv_tmp[40] = -v_tmp[4]+tdata->p[3]*v_tmp[39]*2.0-tdata->p[3]*v_tmp[40]+v_tmp[3]*tdata->dwdx[1];
  Jv_tmp[41] = v_tmp[4]-v_tmp[5]+tdata->p[3]*v_tmp[40]-tdata->p[3]*v_tmp[41];
  Jv_tmp[42] = v_tmp[5]-v_tmp[6]+tdata->p[3]*v_tmp[41]-tdata->p[3]*v_tmp[42];
  Jv_tmp[43] = v_tmp[6]-v_tmp[7]+tdata->p[3]*v_tmp[42]-tdata->p[3]*v_tmp[43];
  Jv_tmp[44] = v_tmp[7]-v_tmp[8]+tdata->p[3]*v_tmp[43]-tdata->p[3]*v_tmp[44];
  Jv_tmp[45] = -tdata->p[0]*tdata->w[0]*v_tmp[45]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[53];
  Jv_tmp[46] = tdata->p[0]*tdata->w[0]*v_tmp[45]-tdata->p[1]*v_tmp[1]*x_tmp[46]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[46]*4.0;
  Jv_tmp[47] = -tdata->p[2]*v_tmp[47]+tdata->p[1]*v_tmp[1]*x_tmp[46]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[46]*2.0;
  Jv_tmp[48] = -tdata->p[3]*v_tmp[48]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[47];
  Jv_tmp[49] = tdata->p[3]*v_tmp[48]*2.0-tdata->p[3]*v_tmp[49];
  Jv_tmp[50] = tdata->p[3]*v_tmp[49]-tdata->p[3]*v_tmp[50];
  Jv_tmp[51] = tdata->p[3]*v_tmp[50]-tdata->p[3]*v_tmp[51];
  Jv_tmp[52] = tdata->p[3]*v_tmp[51]-tdata->p[3]*v_tmp[52];
  Jv_tmp[53] = tdata->p[3]*v_tmp[52]-tdata->p[3]*v_tmp[53];
  Jv_tmp[54] = -tdata->p[0]*v_tmp[0]*tdata->w[5]-tdata->p[0]*tdata->w[0]*v_tmp[54]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[62];
  Jv_tmp[55] = tdata->p[0]*v_tmp[0]*tdata->w[5]+tdata->p[0]*tdata->w[0]*v_tmp[54]-tdata->p[1]*v_tmp[1]*x_tmp[55]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[55]*4.0;
  Jv_tmp[56] = -tdata->p[2]*v_tmp[56]+tdata->p[1]*v_tmp[1]*x_tmp[55]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[55]*2.0;
  Jv_tmp[57] = -tdata->p[3]*v_tmp[57]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[56];
  Jv_tmp[58] = tdata->p[3]*v_tmp[57]*2.0-tdata->p[3]*v_tmp[58];
  Jv_tmp[59] = tdata->p[3]*v_tmp[58]-tdata->p[3]*v_tmp[59];
  Jv_tmp[60] = tdata->p[3]*v_tmp[59]-tdata->p[3]*v_tmp[60];
  Jv_tmp[61] = tdata->p[3]*v_tmp[60]-tdata->p[3]*v_tmp[61];
  Jv_tmp[62] = tdata->p[3]*v_tmp[61]-tdata->p[3]*v_tmp[62];
  Jv_tmp[63] = -tdata->p[0]*v_tmp[0]*tdata->w[6]-tdata->p[0]*tdata->w[0]*v_tmp[63]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[71];
  Jv_tmp[64] = tdata->p[0]*v_tmp[0]*tdata->w[6]+tdata->p[0]*tdata->w[0]*v_tmp[63]-tdata->p[1]*v_tmp[1]*x_tmp[64]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[64]*4.0;
  Jv_tmp[65] = -tdata->p[2]*v_tmp[65]+tdata->p[1]*v_tmp[1]*x_tmp[64]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[64]*2.0;
  Jv_tmp[66] = -tdata->p[3]*v_tmp[66]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[65];
  Jv_tmp[67] = tdata->p[3]*v_tmp[66]*2.0-tdata->p[3]*v_tmp[67];
  Jv_tmp[68] = tdata->p[3]*v_tmp[67]-tdata->p[3]*v_tmp[68];
  Jv_tmp[69] = tdata->p[3]*v_tmp[68]-tdata->p[3]*v_tmp[69];
  Jv_tmp[70] = tdata->p[3]*v_tmp[69]-tdata->p[3]*v_tmp[70];
  Jv_tmp[71] = tdata->p[3]*v_tmp[70]-tdata->p[3]*v_tmp[71];
  Jv_tmp[72] = -tdata->p[0]*v_tmp[0]*tdata->w[7]-tdata->p[0]*tdata->w[0]*v_tmp[72]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[80];
  Jv_tmp[73] = tdata->p[0]*v_tmp[0]*tdata->w[7]+tdata->p[0]*tdata->w[0]*v_tmp[72]-tdata->p[1]*v_tmp[1]*x_tmp[73]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[73]*4.0;
  Jv_tmp[74] = -tdata->p[2]*v_tmp[74]+tdata->p[1]*v_tmp[1]*x_tmp[73]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[73]*2.0;
  Jv_tmp[75] = -tdata->p[3]*v_tmp[75]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[74];
  Jv_tmp[76] = tdata->p[3]*v_tmp[75]*2.0-tdata->p[3]*v_tmp[76];
  Jv_tmp[77] = tdata->p[3]*v_tmp[76]-tdata->p[3]*v_tmp[77];
  Jv_tmp[78] = tdata->p[3]*v_tmp[77]-tdata->p[3]*v_tmp[78];
  Jv_tmp[79] = tdata->p[3]*v_tmp[78]-tdata->p[3]*v_tmp[79];
  Jv_tmp[80] = tdata->p[3]*v_tmp[79]-tdata->p[3]*v_tmp[80];
  Jv_tmp[81] = -tdata->p[0]*v_tmp[0]*tdata->w[8]-tdata->p[0]*tdata->w[0]*v_tmp[81]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[89];
  Jv_tmp[82] = tdata->p[0]*v_tmp[0]*tdata->w[8]+tdata->p[0]*tdata->w[0]*v_tmp[81]-tdata->p[1]*v_tmp[1]*x_tmp[82]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[82]*4.0;
  Jv_tmp[83] = -tdata->p[2]*v_tmp[83]+tdata->p[1]*v_tmp[1]*x_tmp[82]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[82]*2.0;
  Jv_tmp[84] = -tdata->p[3]*v_tmp[84]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[83];
  Jv_tmp[85] = tdata->p[3]*v_tmp[84]*2.0-tdata->p[3]*v_tmp[85];
  Jv_tmp[86] = tdata->p[3]*v_tmp[85]-tdata->p[3]*v_tmp[86];
  Jv_tmp[87] = tdata->p[3]*v_tmp[86]-tdata->p[3]*v_tmp[87];
  Jv_tmp[88] = tdata->p[3]*v_tmp[87]-tdata->p[3]*v_tmp[88];
  Jv_tmp[89] = tdata->p[3]*v_tmp[88]-tdata->p[3]*v_tmp[89];
  Jv_tmp[90] = -tdata->p[0]*v_tmp[0]*tdata->w[9]-tdata->p[0]*tdata->w[0]*v_tmp[90]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[98];
  Jv_tmp[91] = tdata->p[0]*v_tmp[0]*tdata->w[9]+tdata->p[0]*tdata->w[0]*v_tmp[90]-tdata->p[1]*v_tmp[1]*x_tmp[91]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[91]*4.0;
  Jv_tmp[92] = -tdata->p[2]*v_tmp[92]+tdata->p[1]*v_tmp[1]*x_tmp[91]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[91]*2.0;
  Jv_tmp[93] = -tdata->p[3]*v_tmp[93]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[92];
  Jv_tmp[94] = tdata->p[3]*v_tmp[93]*2.0-tdata->p[3]*v_tmp[94];
  Jv_tmp[95] = tdata->p[3]*v_tmp[94]-tdata->p[3]*v_tmp[95];
  Jv_tmp[96] = tdata->p[3]*v_tmp[95]-tdata->p[3]*v_tmp[96];
  Jv_tmp[97] = tdata->p[3]*v_tmp[96]-tdata->p[3]*v_tmp[97];
  Jv_tmp[98] = tdata->p[3]*v_tmp[97]-tdata->p[3]*v_tmp[98];
  Jv_tmp[99] = -tdata->p[0]*tdata->w[0]*v_tmp[99]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[107];
  Jv_tmp[100] = tdata->p[0]*tdata->w[0]*v_tmp[99]-tdata->p[1]*v_tmp[1]*x_tmp[100]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[100]*4.0;
  Jv_tmp[101] = -tdata->p[2]*v_tmp[101]+tdata->p[1]*v_tmp[1]*x_tmp[100]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[100]*2.0;
  Jv_tmp[102] = -tdata->p[3]*v_tmp[102]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[101];
  Jv_tmp[103] = tdata->p[3]*v_tmp[102]*2.0-tdata->p[3]*v_tmp[103];
  Jv_tmp[104] = tdata->p[3]*v_tmp[103]-tdata->p[3]*v_tmp[104];
  Jv_tmp[105] = tdata->p[3]*v_tmp[104]-tdata->p[3]*v_tmp[105];
  Jv_tmp[106] = tdata->p[3]*v_tmp[105]-tdata->p[3]*v_tmp[106];
  Jv_tmp[107] = tdata->p[3]*v_tmp[106]-tdata->p[3]*v_tmp[107];
  Jv_tmp[108] = -tdata->p[0]*tdata->w[0]*v_tmp[108]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[116];
  Jv_tmp[109] = tdata->p[0]*tdata->w[0]*v_tmp[108]-tdata->p[1]*v_tmp[1]*x_tmp[109]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[109]*4.0;
  Jv_tmp[110] = -tdata->p[2]*v_tmp[110]+tdata->p[1]*v_tmp[1]*x_tmp[109]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[109]*2.0;
  Jv_tmp[111] = -tdata->p[3]*v_tmp[111]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[110];
  Jv_tmp[112] = tdata->p[3]*v_tmp[111]*2.0-tdata->p[3]*v_tmp[112];
  Jv_tmp[113] = tdata->p[3]*v_tmp[112]-tdata->p[3]*v_tmp[113];
  Jv_tmp[114] = tdata->p[3]*v_tmp[113]-tdata->p[3]*v_tmp[114];
  Jv_tmp[115] = tdata->p[3]*v_tmp[114]-tdata->p[3]*v_tmp[115];
  Jv_tmp[116] = tdata->p[3]*v_tmp[115]-tdata->p[3]*v_tmp[116];
  Jv_tmp[117] = -tdata->p[0]*tdata->w[0]*v_tmp[117]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[125];
  Jv_tmp[118] = tdata->p[0]*tdata->w[0]*v_tmp[117]-tdata->p[1]*v_tmp[1]*x_tmp[118]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[118]*4.0;
  Jv_tmp[119] = -tdata->p[2]*v_tmp[119]+tdata->p[1]*v_tmp[1]*x_tmp[118]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[118]*2.0;
  Jv_tmp[120] = -tdata->p[3]*v_tmp[120]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[119];
  Jv_tmp[121] = tdata->p[3]*v_tmp[120]*2.0-tdata->p[3]*v_tmp[121];
  Jv_tmp[122] = tdata->p[3]*v_tmp[121]-tdata->p[3]*v_tmp[122];
  Jv_tmp[123] = tdata->p[3]*v_tmp[122]-tdata->p[3]*v_tmp[123];
  Jv_tmp[124] = tdata->p[3]*v_tmp[123]-tdata->p[3]*v_tmp[124];
  Jv_tmp[125] = tdata->p[3]*v_tmp[124]-tdata->p[3]*v_tmp[125];
  Jv_tmp[126] = -tdata->p[0]*tdata->w[0]*v_tmp[126]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[134];
  Jv_tmp[127] = tdata->p[0]*tdata->w[0]*v_tmp[126]-tdata->p[1]*v_tmp[1]*x_tmp[127]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[127]*4.0;
  Jv_tmp[128] = -tdata->p[2]*v_tmp[128]+tdata->p[1]*v_tmp[1]*x_tmp[127]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[127]*2.0;
  Jv_tmp[129] = -tdata->p[3]*v_tmp[129]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[128];
  Jv_tmp[130] = tdata->p[3]*v_tmp[129]*2.0-tdata->p[3]*v_tmp[130];
  Jv_tmp[131] = tdata->p[3]*v_tmp[130]-tdata->p[3]*v_tmp[131];
  Jv_tmp[132] = tdata->p[3]*v_tmp[131]-tdata->p[3]*v_tmp[132];
  Jv_tmp[133] = tdata->p[3]*v_tmp[132]-tdata->p[3]*v_tmp[133];
  Jv_tmp[134] = tdata->p[3]*v_tmp[133]-tdata->p[3]*v_tmp[134];
  Jv_tmp[135] = -tdata->p[0]*tdata->w[0]*v_tmp[135]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[143];
  Jv_tmp[136] = tdata->p[0]*tdata->w[0]*v_tmp[135]-tdata->p[1]*v_tmp[1]*x_tmp[136]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[136]*4.0;
  Jv_tmp[137] = -tdata->p[2]*v_tmp[137]+tdata->p[1]*v_tmp[1]*x_tmp[136]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[136]*2.0;
  Jv_tmp[138] = -tdata->p[3]*v_tmp[138]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[137];
  Jv_tmp[139] = tdata->p[3]*v_tmp[138]*2.0-tdata->p[3]*v_tmp[139];
  Jv_tmp[140] = tdata->p[3]*v_tmp[139]-tdata->p[3]*v_tmp[140];
  Jv_tmp[141] = tdata->p[3]*v_tmp[140]-tdata->p[3]*v_tmp[141];
  Jv_tmp[142] = tdata->p[3]*v_tmp[141]-tdata->p[3]*v_tmp[142];
  Jv_tmp[143] = tdata->p[3]*v_tmp[142]-tdata->p[3]*v_tmp[143];
  Jv_tmp[144] = -tdata->p[0]*tdata->w[0]*v_tmp[144]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[152];
  Jv_tmp[145] = tdata->p[0]*tdata->w[0]*v_tmp[144]-tdata->p[1]*v_tmp[1]*x_tmp[145]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[145]*4.0;
  Jv_tmp[146] = -tdata->p[2]*v_tmp[146]+tdata->p[1]*v_tmp[1]*x_tmp[145]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[145]*2.0;
  Jv_tmp[147] = -tdata->p[3]*v_tmp[147]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[146];
  Jv_tmp[148] = tdata->p[3]*v_tmp[147]*2.0-tdata->p[3]*v_tmp[148];
  Jv_tmp[149] = tdata->p[3]*v_tmp[148]-tdata->p[3]*v_tmp[149];
  Jv_tmp[150] = tdata->p[3]*v_tmp[149]-tdata->p[3]*v_tmp[150];
  Jv_tmp[151] = tdata->p[3]*v_tmp[150]-tdata->p[3]*v_tmp[151];
  Jv_tmp[152] = tdata->p[3]*v_tmp[151]-tdata->p[3]*v_tmp[152];
  Jv_tmp[153] = -tdata->p[0]*tdata->w[0]*v_tmp[153]+udata->k[1]*tdata->p[3]*tdata->w[2]*v_tmp[161];
  Jv_tmp[154] = tdata->p[0]*tdata->w[0]*v_tmp[153]-tdata->p[1]*v_tmp[1]*x_tmp[154]*4.0-tdata->p[1]*x_tmp[1]*v_tmp[154]*4.0;
  Jv_tmp[155] = -tdata->p[2]*v_tmp[155]+tdata->p[1]*v_tmp[1]*x_tmp[154]*2.0+tdata->p[1]*x_tmp[1]*v_tmp[154]*2.0;
  Jv_tmp[156] = -tdata->p[3]*v_tmp[156]+udata->k[0]*tdata->p[2]*tdata->w[3]*v_tmp[155];
  Jv_tmp[157] = tdata->p[3]*v_tmp[156]*2.0-tdata->p[3]*v_tmp[157];
  Jv_tmp[158] = tdata->p[3]*v_tmp[157]-tdata->p[3]*v_tmp[158];
  Jv_tmp[159] = tdata->p[3]*v_tmp[158]-tdata->p[3]*v_tmp[159];
  Jv_tmp[160] = tdata->p[3]*v_tmp[159]-tdata->p[3]*v_tmp[160];
  Jv_tmp[161] = tdata->p[3]*v_tmp[160]-tdata->p[3]*v_tmp[161];
return;

}


