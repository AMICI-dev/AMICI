
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_dwdp.h"
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

int dxdotdp_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
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
int ip;
int ix;
memset(tdata->dxdotdp,0,sizeof(realtype)*162*udata->nplist);
status = dwdp_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  tdata->dxdotdp[0 + ip*model->nx] = -udata->k[0]*tdata->w[0]*x_tmp[0]*tdata->w[2];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->w[0]*x_tmp[0];
  tdata->dxdotdp[9 + ip*model->nx] = -tdata->w[0]*x_tmp[9];
  tdata->dxdotdp[10 + ip*model->nx] = tdata->w[0]*x_tmp[9];
  tdata->dxdotdp[18 + ip*model->nx] = -tdata->w[0]*x_tmp[18];
  tdata->dxdotdp[19 + ip*model->nx] = tdata->w[0]*x_tmp[18];
  tdata->dxdotdp[27 + ip*model->nx] = -tdata->w[0]*x_tmp[27];
  tdata->dxdotdp[28 + ip*model->nx] = tdata->w[0]*x_tmp[27];
  tdata->dxdotdp[36 + ip*model->nx] = -tdata->w[0]*x_tmp[36];
  tdata->dxdotdp[37 + ip*model->nx] = tdata->w[0]*x_tmp[36];
  tdata->dxdotdp[45 + ip*model->nx] = -tdata->w[0]*x_tmp[45];
  tdata->dxdotdp[46 + ip*model->nx] = tdata->w[0]*x_tmp[45];
  tdata->dxdotdp[54 + ip*model->nx] = -x_tmp[0]*tdata->w[5]-tdata->w[0]*x_tmp[54];
  tdata->dxdotdp[55 + ip*model->nx] = x_tmp[0]*tdata->w[5]+tdata->w[0]*x_tmp[54];
  tdata->dxdotdp[63 + ip*model->nx] = -x_tmp[0]*tdata->w[6]-tdata->w[0]*x_tmp[63];
  tdata->dxdotdp[64 + ip*model->nx] = x_tmp[0]*tdata->w[6]+tdata->w[0]*x_tmp[63];
  tdata->dxdotdp[72 + ip*model->nx] = -x_tmp[0]*tdata->w[7]-tdata->w[0]*x_tmp[72];
  tdata->dxdotdp[73 + ip*model->nx] = x_tmp[0]*tdata->w[7]+tdata->w[0]*x_tmp[72];
  tdata->dxdotdp[81 + ip*model->nx] = -x_tmp[0]*tdata->w[8]-tdata->w[0]*x_tmp[81];
  tdata->dxdotdp[82 + ip*model->nx] = x_tmp[0]*tdata->w[8]+tdata->w[0]*x_tmp[81];
  tdata->dxdotdp[90 + ip*model->nx] = -x_tmp[0]*tdata->w[9]-tdata->w[0]*x_tmp[90];
  tdata->dxdotdp[91 + ip*model->nx] = x_tmp[0]*tdata->w[9]+tdata->w[0]*x_tmp[90];
  tdata->dxdotdp[99 + ip*model->nx] = -tdata->w[0]*x_tmp[99];
  tdata->dxdotdp[100 + ip*model->nx] = tdata->w[0]*x_tmp[99];
  tdata->dxdotdp[108 + ip*model->nx] = -tdata->w[0]*x_tmp[108];
  tdata->dxdotdp[109 + ip*model->nx] = tdata->w[0]*x_tmp[108];
  tdata->dxdotdp[117 + ip*model->nx] = -tdata->w[0]*x_tmp[117];
  tdata->dxdotdp[118 + ip*model->nx] = tdata->w[0]*x_tmp[117];
  tdata->dxdotdp[126 + ip*model->nx] = -tdata->w[0]*x_tmp[126];
  tdata->dxdotdp[127 + ip*model->nx] = tdata->w[0]*x_tmp[126];
  tdata->dxdotdp[135 + ip*model->nx] = -tdata->w[0]*x_tmp[135];
  tdata->dxdotdp[136 + ip*model->nx] = tdata->w[0]*x_tmp[135];
  tdata->dxdotdp[144 + ip*model->nx] = -tdata->w[0]*x_tmp[144];
  tdata->dxdotdp[145 + ip*model->nx] = tdata->w[0]*x_tmp[144];
  tdata->dxdotdp[153 + ip*model->nx] = -tdata->w[0]*x_tmp[153];
  tdata->dxdotdp[154 + ip*model->nx] = tdata->w[0]*x_tmp[153];

  } break;

  case 1: {
  tdata->dxdotdp[1 + ip*model->nx] = tdata->w[1]*-2.0;
  tdata->dxdotdp[2 + ip*model->nx] = tdata->w[1];
  tdata->dxdotdp[10 + ip*model->nx] = x_tmp[1]*x_tmp[10]*-4.0;
  tdata->dxdotdp[11 + ip*model->nx] = x_tmp[1]*x_tmp[10]*2.0;
  tdata->dxdotdp[19 + ip*model->nx] = x_tmp[1]*x_tmp[19]*-4.0;
  tdata->dxdotdp[20 + ip*model->nx] = x_tmp[1]*x_tmp[19]*2.0;
  tdata->dxdotdp[28 + ip*model->nx] = x_tmp[1]*x_tmp[28]*-4.0;
  tdata->dxdotdp[29 + ip*model->nx] = x_tmp[1]*x_tmp[28]*2.0;
  tdata->dxdotdp[37 + ip*model->nx] = x_tmp[1]*x_tmp[37]*-4.0;
  tdata->dxdotdp[38 + ip*model->nx] = x_tmp[1]*x_tmp[37]*2.0;
  tdata->dxdotdp[46 + ip*model->nx] = x_tmp[1]*x_tmp[46]*-4.0;
  tdata->dxdotdp[47 + ip*model->nx] = x_tmp[1]*x_tmp[46]*2.0;
  tdata->dxdotdp[55 + ip*model->nx] = x_tmp[1]*x_tmp[55]*-4.0;
  tdata->dxdotdp[56 + ip*model->nx] = x_tmp[1]*x_tmp[55]*2.0;
  tdata->dxdotdp[64 + ip*model->nx] = x_tmp[1]*x_tmp[64]*-4.0;
  tdata->dxdotdp[65 + ip*model->nx] = x_tmp[1]*x_tmp[64]*2.0;
  tdata->dxdotdp[73 + ip*model->nx] = x_tmp[1]*x_tmp[73]*-4.0;
  tdata->dxdotdp[74 + ip*model->nx] = x_tmp[1]*x_tmp[73]*2.0;
  tdata->dxdotdp[82 + ip*model->nx] = x_tmp[1]*x_tmp[82]*-4.0;
  tdata->dxdotdp[83 + ip*model->nx] = x_tmp[1]*x_tmp[82]*2.0;
  tdata->dxdotdp[91 + ip*model->nx] = x_tmp[1]*x_tmp[91]*-4.0;
  tdata->dxdotdp[92 + ip*model->nx] = x_tmp[1]*x_tmp[91]*2.0;
  tdata->dxdotdp[100 + ip*model->nx] = x_tmp[1]*x_tmp[100]*-4.0;
  tdata->dxdotdp[101 + ip*model->nx] = x_tmp[1]*x_tmp[100]*2.0;
  tdata->dxdotdp[109 + ip*model->nx] = x_tmp[1]*x_tmp[109]*-4.0;
  tdata->dxdotdp[110 + ip*model->nx] = x_tmp[1]*x_tmp[109]*2.0;
  tdata->dxdotdp[118 + ip*model->nx] = x_tmp[1]*x_tmp[118]*-4.0;
  tdata->dxdotdp[119 + ip*model->nx] = x_tmp[1]*x_tmp[118]*2.0;
  tdata->dxdotdp[127 + ip*model->nx] = x_tmp[1]*x_tmp[127]*-4.0;
  tdata->dxdotdp[128 + ip*model->nx] = x_tmp[1]*x_tmp[127]*2.0;
  tdata->dxdotdp[136 + ip*model->nx] = x_tmp[1]*x_tmp[136]*-4.0;
  tdata->dxdotdp[137 + ip*model->nx] = x_tmp[1]*x_tmp[136]*2.0;
  tdata->dxdotdp[145 + ip*model->nx] = x_tmp[1]*x_tmp[145]*-4.0;
  tdata->dxdotdp[146 + ip*model->nx] = x_tmp[1]*x_tmp[145]*2.0;
  tdata->dxdotdp[154 + ip*model->nx] = x_tmp[1]*x_tmp[154]*-4.0;
  tdata->dxdotdp[155 + ip*model->nx] = x_tmp[1]*x_tmp[154]*2.0;

  } break;

  case 2: {
  tdata->dxdotdp[2 + ip*model->nx] = -x_tmp[2];
  tdata->dxdotdp[3 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[2];
  tdata->dxdotdp[11 + ip*model->nx] = -x_tmp[11];
  tdata->dxdotdp[12 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[11];
  tdata->dxdotdp[20 + ip*model->nx] = -x_tmp[20];
  tdata->dxdotdp[21 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[20];
  tdata->dxdotdp[29 + ip*model->nx] = -x_tmp[29];
  tdata->dxdotdp[30 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[29];
  tdata->dxdotdp[38 + ip*model->nx] = -x_tmp[38];
  tdata->dxdotdp[39 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[38];
  tdata->dxdotdp[47 + ip*model->nx] = -x_tmp[47];
  tdata->dxdotdp[48 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[47];
  tdata->dxdotdp[56 + ip*model->nx] = -x_tmp[56];
  tdata->dxdotdp[57 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[56];
  tdata->dxdotdp[65 + ip*model->nx] = -x_tmp[65];
  tdata->dxdotdp[66 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[65];
  tdata->dxdotdp[74 + ip*model->nx] = -x_tmp[74];
  tdata->dxdotdp[75 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[74];
  tdata->dxdotdp[83 + ip*model->nx] = -x_tmp[83];
  tdata->dxdotdp[84 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[83];
  tdata->dxdotdp[92 + ip*model->nx] = -x_tmp[92];
  tdata->dxdotdp[93 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[92];
  tdata->dxdotdp[101 + ip*model->nx] = -x_tmp[101];
  tdata->dxdotdp[102 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[101];
  tdata->dxdotdp[110 + ip*model->nx] = -x_tmp[110];
  tdata->dxdotdp[111 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[110];
  tdata->dxdotdp[119 + ip*model->nx] = -x_tmp[119];
  tdata->dxdotdp[120 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[119];
  tdata->dxdotdp[128 + ip*model->nx] = -x_tmp[128];
  tdata->dxdotdp[129 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[128];
  tdata->dxdotdp[137 + ip*model->nx] = -x_tmp[137];
  tdata->dxdotdp[138 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[137];
  tdata->dxdotdp[146 + ip*model->nx] = -x_tmp[146];
  tdata->dxdotdp[147 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[146];
  tdata->dxdotdp[155 + ip*model->nx] = -x_tmp[155];
  tdata->dxdotdp[156 + ip*model->nx] = udata->k[0]*tdata->w[3]*x_tmp[155];

  } break;

  case 3: {
  tdata->dxdotdp[0 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[8];
  tdata->dxdotdp[3 + ip*model->nx] = -udata->k[1]*tdata->w[3]*x_tmp[3];
  tdata->dxdotdp[4 + ip*model->nx] = tdata->w[4]-x_tmp[4];
  tdata->dxdotdp[5 + ip*model->nx] = x_tmp[4]-x_tmp[5];
  tdata->dxdotdp[6 + ip*model->nx] = x_tmp[5]-x_tmp[6];
  tdata->dxdotdp[7 + ip*model->nx] = x_tmp[6]-x_tmp[7];
  tdata->dxdotdp[8 + ip*model->nx] = x_tmp[7]-x_tmp[8];
  tdata->dxdotdp[9 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[17];
  tdata->dxdotdp[12 + ip*model->nx] = -x_tmp[12];
  tdata->dxdotdp[13 + ip*model->nx] = x_tmp[12]*2.0-x_tmp[13];
  tdata->dxdotdp[14 + ip*model->nx] = x_tmp[13]-x_tmp[14];
  tdata->dxdotdp[15 + ip*model->nx] = x_tmp[14]-x_tmp[15];
  tdata->dxdotdp[16 + ip*model->nx] = x_tmp[15]-x_tmp[16];
  tdata->dxdotdp[17 + ip*model->nx] = x_tmp[16]-x_tmp[17];
  tdata->dxdotdp[18 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[26];
  tdata->dxdotdp[21 + ip*model->nx] = -x_tmp[21];
  tdata->dxdotdp[22 + ip*model->nx] = x_tmp[21]*2.0-x_tmp[22];
  tdata->dxdotdp[23 + ip*model->nx] = x_tmp[22]-x_tmp[23];
  tdata->dxdotdp[24 + ip*model->nx] = x_tmp[23]-x_tmp[24];
  tdata->dxdotdp[25 + ip*model->nx] = x_tmp[24]-x_tmp[25];
  tdata->dxdotdp[26 + ip*model->nx] = x_tmp[25]-x_tmp[26];
  tdata->dxdotdp[27 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[35];
  tdata->dxdotdp[30 + ip*model->nx] = -x_tmp[30];
  tdata->dxdotdp[31 + ip*model->nx] = x_tmp[30]*2.0-x_tmp[31];
  tdata->dxdotdp[32 + ip*model->nx] = x_tmp[31]-x_tmp[32];
  tdata->dxdotdp[33 + ip*model->nx] = x_tmp[32]-x_tmp[33];
  tdata->dxdotdp[34 + ip*model->nx] = x_tmp[33]-x_tmp[34];
  tdata->dxdotdp[35 + ip*model->nx] = x_tmp[34]-x_tmp[35];
  tdata->dxdotdp[36 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[44];
  tdata->dxdotdp[39 + ip*model->nx] = -x_tmp[39];
  tdata->dxdotdp[40 + ip*model->nx] = x_tmp[39]*2.0-x_tmp[40];
  tdata->dxdotdp[41 + ip*model->nx] = x_tmp[40]-x_tmp[41];
  tdata->dxdotdp[42 + ip*model->nx] = x_tmp[41]-x_tmp[42];
  tdata->dxdotdp[43 + ip*model->nx] = x_tmp[42]-x_tmp[43];
  tdata->dxdotdp[44 + ip*model->nx] = x_tmp[43]-x_tmp[44];
  tdata->dxdotdp[45 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[53];
  tdata->dxdotdp[48 + ip*model->nx] = -x_tmp[48];
  tdata->dxdotdp[49 + ip*model->nx] = x_tmp[48]*2.0-x_tmp[49];
  tdata->dxdotdp[50 + ip*model->nx] = x_tmp[49]-x_tmp[50];
  tdata->dxdotdp[51 + ip*model->nx] = x_tmp[50]-x_tmp[51];
  tdata->dxdotdp[52 + ip*model->nx] = x_tmp[51]-x_tmp[52];
  tdata->dxdotdp[53 + ip*model->nx] = x_tmp[52]-x_tmp[53];
  tdata->dxdotdp[54 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[62];
  tdata->dxdotdp[57 + ip*model->nx] = -x_tmp[57];
  tdata->dxdotdp[58 + ip*model->nx] = x_tmp[57]*2.0-x_tmp[58];
  tdata->dxdotdp[59 + ip*model->nx] = x_tmp[58]-x_tmp[59];
  tdata->dxdotdp[60 + ip*model->nx] = x_tmp[59]-x_tmp[60];
  tdata->dxdotdp[61 + ip*model->nx] = x_tmp[60]-x_tmp[61];
  tdata->dxdotdp[62 + ip*model->nx] = x_tmp[61]-x_tmp[62];
  tdata->dxdotdp[63 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[71];
  tdata->dxdotdp[66 + ip*model->nx] = -x_tmp[66];
  tdata->dxdotdp[67 + ip*model->nx] = x_tmp[66]*2.0-x_tmp[67];
  tdata->dxdotdp[68 + ip*model->nx] = x_tmp[67]-x_tmp[68];
  tdata->dxdotdp[69 + ip*model->nx] = x_tmp[68]-x_tmp[69];
  tdata->dxdotdp[70 + ip*model->nx] = x_tmp[69]-x_tmp[70];
  tdata->dxdotdp[71 + ip*model->nx] = x_tmp[70]-x_tmp[71];
  tdata->dxdotdp[72 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[80];
  tdata->dxdotdp[75 + ip*model->nx] = -x_tmp[75];
  tdata->dxdotdp[76 + ip*model->nx] = x_tmp[75]*2.0-x_tmp[76];
  tdata->dxdotdp[77 + ip*model->nx] = x_tmp[76]-x_tmp[77];
  tdata->dxdotdp[78 + ip*model->nx] = x_tmp[77]-x_tmp[78];
  tdata->dxdotdp[79 + ip*model->nx] = x_tmp[78]-x_tmp[79];
  tdata->dxdotdp[80 + ip*model->nx] = x_tmp[79]-x_tmp[80];
  tdata->dxdotdp[81 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[89];
  tdata->dxdotdp[84 + ip*model->nx] = -x_tmp[84];
  tdata->dxdotdp[85 + ip*model->nx] = x_tmp[84]*2.0-x_tmp[85];
  tdata->dxdotdp[86 + ip*model->nx] = x_tmp[85]-x_tmp[86];
  tdata->dxdotdp[87 + ip*model->nx] = x_tmp[86]-x_tmp[87];
  tdata->dxdotdp[88 + ip*model->nx] = x_tmp[87]-x_tmp[88];
  tdata->dxdotdp[89 + ip*model->nx] = x_tmp[88]-x_tmp[89];
  tdata->dxdotdp[90 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[98];
  tdata->dxdotdp[93 + ip*model->nx] = -x_tmp[93];
  tdata->dxdotdp[94 + ip*model->nx] = x_tmp[93]*2.0-x_tmp[94];
  tdata->dxdotdp[95 + ip*model->nx] = x_tmp[94]-x_tmp[95];
  tdata->dxdotdp[96 + ip*model->nx] = x_tmp[95]-x_tmp[96];
  tdata->dxdotdp[97 + ip*model->nx] = x_tmp[96]-x_tmp[97];
  tdata->dxdotdp[98 + ip*model->nx] = x_tmp[97]-x_tmp[98];
  tdata->dxdotdp[99 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[107];
  tdata->dxdotdp[102 + ip*model->nx] = -x_tmp[102];
  tdata->dxdotdp[103 + ip*model->nx] = x_tmp[102]*2.0-x_tmp[103];
  tdata->dxdotdp[104 + ip*model->nx] = x_tmp[103]-x_tmp[104];
  tdata->dxdotdp[105 + ip*model->nx] = x_tmp[104]-x_tmp[105];
  tdata->dxdotdp[106 + ip*model->nx] = x_tmp[105]-x_tmp[106];
  tdata->dxdotdp[107 + ip*model->nx] = x_tmp[106]-x_tmp[107];
  tdata->dxdotdp[108 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[116];
  tdata->dxdotdp[111 + ip*model->nx] = -x_tmp[111];
  tdata->dxdotdp[112 + ip*model->nx] = x_tmp[111]*2.0-x_tmp[112];
  tdata->dxdotdp[113 + ip*model->nx] = x_tmp[112]-x_tmp[113];
  tdata->dxdotdp[114 + ip*model->nx] = x_tmp[113]-x_tmp[114];
  tdata->dxdotdp[115 + ip*model->nx] = x_tmp[114]-x_tmp[115];
  tdata->dxdotdp[116 + ip*model->nx] = x_tmp[115]-x_tmp[116];
  tdata->dxdotdp[117 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[125];
  tdata->dxdotdp[120 + ip*model->nx] = -x_tmp[120];
  tdata->dxdotdp[121 + ip*model->nx] = x_tmp[120]*2.0-x_tmp[121];
  tdata->dxdotdp[122 + ip*model->nx] = x_tmp[121]-x_tmp[122];
  tdata->dxdotdp[123 + ip*model->nx] = x_tmp[122]-x_tmp[123];
  tdata->dxdotdp[124 + ip*model->nx] = x_tmp[123]-x_tmp[124];
  tdata->dxdotdp[125 + ip*model->nx] = x_tmp[124]-x_tmp[125];
  tdata->dxdotdp[126 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[134];
  tdata->dxdotdp[129 + ip*model->nx] = -x_tmp[129];
  tdata->dxdotdp[130 + ip*model->nx] = x_tmp[129]*2.0-x_tmp[130];
  tdata->dxdotdp[131 + ip*model->nx] = x_tmp[130]-x_tmp[131];
  tdata->dxdotdp[132 + ip*model->nx] = x_tmp[131]-x_tmp[132];
  tdata->dxdotdp[133 + ip*model->nx] = x_tmp[132]-x_tmp[133];
  tdata->dxdotdp[134 + ip*model->nx] = x_tmp[133]-x_tmp[134];
  tdata->dxdotdp[135 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[143];
  tdata->dxdotdp[138 + ip*model->nx] = -x_tmp[138];
  tdata->dxdotdp[139 + ip*model->nx] = x_tmp[138]*2.0-x_tmp[139];
  tdata->dxdotdp[140 + ip*model->nx] = x_tmp[139]-x_tmp[140];
  tdata->dxdotdp[141 + ip*model->nx] = x_tmp[140]-x_tmp[141];
  tdata->dxdotdp[142 + ip*model->nx] = x_tmp[141]-x_tmp[142];
  tdata->dxdotdp[143 + ip*model->nx] = x_tmp[142]-x_tmp[143];
  tdata->dxdotdp[144 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[152];
  tdata->dxdotdp[147 + ip*model->nx] = -x_tmp[147];
  tdata->dxdotdp[148 + ip*model->nx] = x_tmp[147]*2.0-x_tmp[148];
  tdata->dxdotdp[149 + ip*model->nx] = x_tmp[148]-x_tmp[149];
  tdata->dxdotdp[150 + ip*model->nx] = x_tmp[149]-x_tmp[150];
  tdata->dxdotdp[151 + ip*model->nx] = x_tmp[150]-x_tmp[151];
  tdata->dxdotdp[152 + ip*model->nx] = x_tmp[151]-x_tmp[152];
  tdata->dxdotdp[153 + ip*model->nx] = udata->k[1]*tdata->w[2]*x_tmp[161];
  tdata->dxdotdp[156 + ip*model->nx] = -x_tmp[156];
  tdata->dxdotdp[157 + ip*model->nx] = x_tmp[156]*2.0-x_tmp[157];
  tdata->dxdotdp[158 + ip*model->nx] = x_tmp[157]-x_tmp[158];
  tdata->dxdotdp[159 + ip*model->nx] = x_tmp[158]-x_tmp[159];
  tdata->dxdotdp[160 + ip*model->nx] = x_tmp[159]-x_tmp[160];
  tdata->dxdotdp[161 + ip*model->nx] = x_tmp[160]-x_tmp[161];

  } break;

  case 5: {
  tdata->dxdotdp[0 + ip*model->nx] = -udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*tdata->dwdp[0];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[0];
  tdata->dxdotdp[9 + ip*model->nx] = -tdata->dwdp[0]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[10 + ip*model->nx] = tdata->dwdp[0]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[18 + ip*model->nx] = -tdata->p[0]*x_tmp[18]*tdata->dwdp[0];
  tdata->dxdotdp[19 + ip*model->nx] = tdata->p[0]*x_tmp[18]*tdata->dwdp[0];
  tdata->dxdotdp[27 + ip*model->nx] = -tdata->p[0]*x_tmp[27]*tdata->dwdp[0];
  tdata->dxdotdp[28 + ip*model->nx] = tdata->p[0]*x_tmp[27]*tdata->dwdp[0];
  tdata->dxdotdp[36 + ip*model->nx] = -tdata->p[0]*x_tmp[36]*tdata->dwdp[0];
  tdata->dxdotdp[37 + ip*model->nx] = tdata->p[0]*x_tmp[36]*tdata->dwdp[0];
  tdata->dxdotdp[45 + ip*model->nx] = -tdata->p[0]*x_tmp[45]*tdata->dwdp[0];
  tdata->dxdotdp[46 + ip*model->nx] = tdata->p[0]*x_tmp[45]*tdata->dwdp[0];
  tdata->dxdotdp[54 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[1]-tdata->p[0]*x_tmp[54]*tdata->dwdp[0];
  tdata->dxdotdp[55 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[1]+tdata->p[0]*x_tmp[54]*tdata->dwdp[0];
  tdata->dxdotdp[63 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[2]-tdata->p[0]*x_tmp[63]*tdata->dwdp[0];
  tdata->dxdotdp[64 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[2]+tdata->p[0]*x_tmp[63]*tdata->dwdp[0];
  tdata->dxdotdp[72 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[3]-tdata->p[0]*x_tmp[72]*tdata->dwdp[0];
  tdata->dxdotdp[73 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[3]+tdata->p[0]*x_tmp[72]*tdata->dwdp[0];
  tdata->dxdotdp[81 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[4]-tdata->p[0]*x_tmp[81]*tdata->dwdp[0];
  tdata->dxdotdp[82 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[4]+tdata->p[0]*x_tmp[81]*tdata->dwdp[0];
  tdata->dxdotdp[90 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[5]-tdata->p[0]*x_tmp[90]*tdata->dwdp[0];
  tdata->dxdotdp[91 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[5]+tdata->p[0]*x_tmp[90]*tdata->dwdp[0];
  tdata->dxdotdp[99 + ip*model->nx] = -tdata->p[0]*x_tmp[99]*tdata->dwdp[0];
  tdata->dxdotdp[100 + ip*model->nx] = tdata->p[0]*x_tmp[99]*tdata->dwdp[0];
  tdata->dxdotdp[108 + ip*model->nx] = -tdata->p[0]*x_tmp[108]*tdata->dwdp[0];
  tdata->dxdotdp[109 + ip*model->nx] = tdata->p[0]*x_tmp[108]*tdata->dwdp[0];
  tdata->dxdotdp[117 + ip*model->nx] = -tdata->p[0]*x_tmp[117]*tdata->dwdp[0];
  tdata->dxdotdp[118 + ip*model->nx] = tdata->p[0]*x_tmp[117]*tdata->dwdp[0];
  tdata->dxdotdp[126 + ip*model->nx] = -tdata->p[0]*x_tmp[126]*tdata->dwdp[0];
  tdata->dxdotdp[127 + ip*model->nx] = tdata->p[0]*x_tmp[126]*tdata->dwdp[0];
  tdata->dxdotdp[135 + ip*model->nx] = -tdata->p[0]*x_tmp[135]*tdata->dwdp[0];
  tdata->dxdotdp[136 + ip*model->nx] = tdata->p[0]*x_tmp[135]*tdata->dwdp[0];
  tdata->dxdotdp[144 + ip*model->nx] = -tdata->p[0]*x_tmp[144]*tdata->dwdp[0];
  tdata->dxdotdp[145 + ip*model->nx] = tdata->p[0]*x_tmp[144]*tdata->dwdp[0];
  tdata->dxdotdp[153 + ip*model->nx] = -tdata->p[0]*x_tmp[153]*tdata->dwdp[0];
  tdata->dxdotdp[154 + ip*model->nx] = tdata->p[0]*x_tmp[153]*tdata->dwdp[0];

  } break;

  case 6: {
  tdata->dxdotdp[0 + ip*model->nx] = -udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*tdata->dwdp[6];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[6];
  tdata->dxdotdp[9 + ip*model->nx] = -tdata->dwdp[6]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[10 + ip*model->nx] = tdata->dwdp[6]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[18 + ip*model->nx] = -tdata->p[0]*x_tmp[18]*tdata->dwdp[6];
  tdata->dxdotdp[19 + ip*model->nx] = tdata->p[0]*x_tmp[18]*tdata->dwdp[6];
  tdata->dxdotdp[27 + ip*model->nx] = -tdata->p[0]*x_tmp[27]*tdata->dwdp[6];
  tdata->dxdotdp[28 + ip*model->nx] = tdata->p[0]*x_tmp[27]*tdata->dwdp[6];
  tdata->dxdotdp[36 + ip*model->nx] = -tdata->p[0]*x_tmp[36]*tdata->dwdp[6];
  tdata->dxdotdp[37 + ip*model->nx] = tdata->p[0]*x_tmp[36]*tdata->dwdp[6];
  tdata->dxdotdp[45 + ip*model->nx] = -tdata->p[0]*x_tmp[45]*tdata->dwdp[6];
  tdata->dxdotdp[46 + ip*model->nx] = tdata->p[0]*x_tmp[45]*tdata->dwdp[6];
  tdata->dxdotdp[54 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[7]-tdata->p[0]*x_tmp[54]*tdata->dwdp[6];
  tdata->dxdotdp[55 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[7]+tdata->p[0]*x_tmp[54]*tdata->dwdp[6];
  tdata->dxdotdp[63 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[8]-tdata->p[0]*x_tmp[63]*tdata->dwdp[6];
  tdata->dxdotdp[64 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[8]+tdata->p[0]*x_tmp[63]*tdata->dwdp[6];
  tdata->dxdotdp[72 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[9]-tdata->p[0]*x_tmp[72]*tdata->dwdp[6];
  tdata->dxdotdp[73 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[9]+tdata->p[0]*x_tmp[72]*tdata->dwdp[6];
  tdata->dxdotdp[81 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[10]-tdata->p[0]*x_tmp[81]*tdata->dwdp[6];
  tdata->dxdotdp[82 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[10]+tdata->p[0]*x_tmp[81]*tdata->dwdp[6];
  tdata->dxdotdp[90 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[11]-tdata->p[0]*x_tmp[90]*tdata->dwdp[6];
  tdata->dxdotdp[91 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[11]+tdata->p[0]*x_tmp[90]*tdata->dwdp[6];
  tdata->dxdotdp[99 + ip*model->nx] = -tdata->p[0]*x_tmp[99]*tdata->dwdp[6];
  tdata->dxdotdp[100 + ip*model->nx] = tdata->p[0]*x_tmp[99]*tdata->dwdp[6];
  tdata->dxdotdp[108 + ip*model->nx] = -tdata->p[0]*x_tmp[108]*tdata->dwdp[6];
  tdata->dxdotdp[109 + ip*model->nx] = tdata->p[0]*x_tmp[108]*tdata->dwdp[6];
  tdata->dxdotdp[117 + ip*model->nx] = -tdata->p[0]*x_tmp[117]*tdata->dwdp[6];
  tdata->dxdotdp[118 + ip*model->nx] = tdata->p[0]*x_tmp[117]*tdata->dwdp[6];
  tdata->dxdotdp[126 + ip*model->nx] = -tdata->p[0]*x_tmp[126]*tdata->dwdp[6];
  tdata->dxdotdp[127 + ip*model->nx] = tdata->p[0]*x_tmp[126]*tdata->dwdp[6];
  tdata->dxdotdp[135 + ip*model->nx] = -tdata->p[0]*x_tmp[135]*tdata->dwdp[6];
  tdata->dxdotdp[136 + ip*model->nx] = tdata->p[0]*x_tmp[135]*tdata->dwdp[6];
  tdata->dxdotdp[144 + ip*model->nx] = -tdata->p[0]*x_tmp[144]*tdata->dwdp[6];
  tdata->dxdotdp[145 + ip*model->nx] = tdata->p[0]*x_tmp[144]*tdata->dwdp[6];
  tdata->dxdotdp[153 + ip*model->nx] = -tdata->p[0]*x_tmp[153]*tdata->dwdp[6];
  tdata->dxdotdp[154 + ip*model->nx] = tdata->p[0]*x_tmp[153]*tdata->dwdp[6];

  } break;

  case 7: {
  tdata->dxdotdp[0 + ip*model->nx] = -udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*tdata->dwdp[12];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[12];
  tdata->dxdotdp[9 + ip*model->nx] = -tdata->dwdp[12]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[10 + ip*model->nx] = tdata->dwdp[12]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[18 + ip*model->nx] = -tdata->p[0]*x_tmp[18]*tdata->dwdp[12];
  tdata->dxdotdp[19 + ip*model->nx] = tdata->p[0]*x_tmp[18]*tdata->dwdp[12];
  tdata->dxdotdp[27 + ip*model->nx] = -tdata->p[0]*x_tmp[27]*tdata->dwdp[12];
  tdata->dxdotdp[28 + ip*model->nx] = tdata->p[0]*x_tmp[27]*tdata->dwdp[12];
  tdata->dxdotdp[36 + ip*model->nx] = -tdata->p[0]*x_tmp[36]*tdata->dwdp[12];
  tdata->dxdotdp[37 + ip*model->nx] = tdata->p[0]*x_tmp[36]*tdata->dwdp[12];
  tdata->dxdotdp[45 + ip*model->nx] = -tdata->p[0]*x_tmp[45]*tdata->dwdp[12];
  tdata->dxdotdp[46 + ip*model->nx] = tdata->p[0]*x_tmp[45]*tdata->dwdp[12];
  tdata->dxdotdp[54 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[13]-tdata->p[0]*x_tmp[54]*tdata->dwdp[12];
  tdata->dxdotdp[55 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[13]+tdata->p[0]*x_tmp[54]*tdata->dwdp[12];
  tdata->dxdotdp[63 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[14]-tdata->p[0]*x_tmp[63]*tdata->dwdp[12];
  tdata->dxdotdp[64 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[14]+tdata->p[0]*x_tmp[63]*tdata->dwdp[12];
  tdata->dxdotdp[72 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[15]-tdata->p[0]*x_tmp[72]*tdata->dwdp[12];
  tdata->dxdotdp[73 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[15]+tdata->p[0]*x_tmp[72]*tdata->dwdp[12];
  tdata->dxdotdp[81 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[16]-tdata->p[0]*x_tmp[81]*tdata->dwdp[12];
  tdata->dxdotdp[82 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[16]+tdata->p[0]*x_tmp[81]*tdata->dwdp[12];
  tdata->dxdotdp[90 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[17]-tdata->p[0]*x_tmp[90]*tdata->dwdp[12];
  tdata->dxdotdp[91 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[17]+tdata->p[0]*x_tmp[90]*tdata->dwdp[12];
  tdata->dxdotdp[99 + ip*model->nx] = -tdata->p[0]*x_tmp[99]*tdata->dwdp[12];
  tdata->dxdotdp[100 + ip*model->nx] = tdata->p[0]*x_tmp[99]*tdata->dwdp[12];
  tdata->dxdotdp[108 + ip*model->nx] = -tdata->p[0]*x_tmp[108]*tdata->dwdp[12];
  tdata->dxdotdp[109 + ip*model->nx] = tdata->p[0]*x_tmp[108]*tdata->dwdp[12];
  tdata->dxdotdp[117 + ip*model->nx] = -tdata->p[0]*x_tmp[117]*tdata->dwdp[12];
  tdata->dxdotdp[118 + ip*model->nx] = tdata->p[0]*x_tmp[117]*tdata->dwdp[12];
  tdata->dxdotdp[126 + ip*model->nx] = -tdata->p[0]*x_tmp[126]*tdata->dwdp[12];
  tdata->dxdotdp[127 + ip*model->nx] = tdata->p[0]*x_tmp[126]*tdata->dwdp[12];
  tdata->dxdotdp[135 + ip*model->nx] = -tdata->p[0]*x_tmp[135]*tdata->dwdp[12];
  tdata->dxdotdp[136 + ip*model->nx] = tdata->p[0]*x_tmp[135]*tdata->dwdp[12];
  tdata->dxdotdp[144 + ip*model->nx] = -tdata->p[0]*x_tmp[144]*tdata->dwdp[12];
  tdata->dxdotdp[145 + ip*model->nx] = tdata->p[0]*x_tmp[144]*tdata->dwdp[12];
  tdata->dxdotdp[153 + ip*model->nx] = -tdata->p[0]*x_tmp[153]*tdata->dwdp[12];
  tdata->dxdotdp[154 + ip*model->nx] = tdata->p[0]*x_tmp[153]*tdata->dwdp[12];

  } break;

  case 8: {
  tdata->dxdotdp[0 + ip*model->nx] = -udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*tdata->dwdp[18];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[18];
  tdata->dxdotdp[9 + ip*model->nx] = -tdata->dwdp[18]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[10 + ip*model->nx] = tdata->dwdp[18]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[18 + ip*model->nx] = -tdata->p[0]*x_tmp[18]*tdata->dwdp[18];
  tdata->dxdotdp[19 + ip*model->nx] = tdata->p[0]*x_tmp[18]*tdata->dwdp[18];
  tdata->dxdotdp[27 + ip*model->nx] = -tdata->p[0]*x_tmp[27]*tdata->dwdp[18];
  tdata->dxdotdp[28 + ip*model->nx] = tdata->p[0]*x_tmp[27]*tdata->dwdp[18];
  tdata->dxdotdp[36 + ip*model->nx] = -tdata->p[0]*x_tmp[36]*tdata->dwdp[18];
  tdata->dxdotdp[37 + ip*model->nx] = tdata->p[0]*x_tmp[36]*tdata->dwdp[18];
  tdata->dxdotdp[45 + ip*model->nx] = -tdata->p[0]*x_tmp[45]*tdata->dwdp[18];
  tdata->dxdotdp[46 + ip*model->nx] = tdata->p[0]*x_tmp[45]*tdata->dwdp[18];
  tdata->dxdotdp[54 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[19]-tdata->p[0]*x_tmp[54]*tdata->dwdp[18];
  tdata->dxdotdp[55 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[19]+tdata->p[0]*x_tmp[54]*tdata->dwdp[18];
  tdata->dxdotdp[63 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[20]-tdata->p[0]*x_tmp[63]*tdata->dwdp[18];
  tdata->dxdotdp[64 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[20]+tdata->p[0]*x_tmp[63]*tdata->dwdp[18];
  tdata->dxdotdp[72 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[21]-tdata->p[0]*x_tmp[72]*tdata->dwdp[18];
  tdata->dxdotdp[73 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[21]+tdata->p[0]*x_tmp[72]*tdata->dwdp[18];
  tdata->dxdotdp[81 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[22]-tdata->p[0]*x_tmp[81]*tdata->dwdp[18];
  tdata->dxdotdp[82 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[22]+tdata->p[0]*x_tmp[81]*tdata->dwdp[18];
  tdata->dxdotdp[90 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[23]-tdata->p[0]*x_tmp[90]*tdata->dwdp[18];
  tdata->dxdotdp[91 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[23]+tdata->p[0]*x_tmp[90]*tdata->dwdp[18];
  tdata->dxdotdp[99 + ip*model->nx] = -tdata->p[0]*x_tmp[99]*tdata->dwdp[18];
  tdata->dxdotdp[100 + ip*model->nx] = tdata->p[0]*x_tmp[99]*tdata->dwdp[18];
  tdata->dxdotdp[108 + ip*model->nx] = -tdata->p[0]*x_tmp[108]*tdata->dwdp[18];
  tdata->dxdotdp[109 + ip*model->nx] = tdata->p[0]*x_tmp[108]*tdata->dwdp[18];
  tdata->dxdotdp[117 + ip*model->nx] = -tdata->p[0]*x_tmp[117]*tdata->dwdp[18];
  tdata->dxdotdp[118 + ip*model->nx] = tdata->p[0]*x_tmp[117]*tdata->dwdp[18];
  tdata->dxdotdp[126 + ip*model->nx] = -tdata->p[0]*x_tmp[126]*tdata->dwdp[18];
  tdata->dxdotdp[127 + ip*model->nx] = tdata->p[0]*x_tmp[126]*tdata->dwdp[18];
  tdata->dxdotdp[135 + ip*model->nx] = -tdata->p[0]*x_tmp[135]*tdata->dwdp[18];
  tdata->dxdotdp[136 + ip*model->nx] = tdata->p[0]*x_tmp[135]*tdata->dwdp[18];
  tdata->dxdotdp[144 + ip*model->nx] = -tdata->p[0]*x_tmp[144]*tdata->dwdp[18];
  tdata->dxdotdp[145 + ip*model->nx] = tdata->p[0]*x_tmp[144]*tdata->dwdp[18];
  tdata->dxdotdp[153 + ip*model->nx] = -tdata->p[0]*x_tmp[153]*tdata->dwdp[18];
  tdata->dxdotdp[154 + ip*model->nx] = tdata->p[0]*x_tmp[153]*tdata->dwdp[18];

  } break;

  case 9: {
  tdata->dxdotdp[0 + ip*model->nx] = -udata->k[0]*tdata->p[0]*x_tmp[0]*tdata->w[2]*tdata->dwdp[24];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[24];
  tdata->dxdotdp[9 + ip*model->nx] = -tdata->dwdp[24]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[10 + ip*model->nx] = tdata->dwdp[24]*(x_tmp[0]+tdata->p[0]*x_tmp[9]);
  tdata->dxdotdp[18 + ip*model->nx] = -tdata->p[0]*x_tmp[18]*tdata->dwdp[24];
  tdata->dxdotdp[19 + ip*model->nx] = tdata->p[0]*x_tmp[18]*tdata->dwdp[24];
  tdata->dxdotdp[27 + ip*model->nx] = -tdata->p[0]*x_tmp[27]*tdata->dwdp[24];
  tdata->dxdotdp[28 + ip*model->nx] = tdata->p[0]*x_tmp[27]*tdata->dwdp[24];
  tdata->dxdotdp[36 + ip*model->nx] = -tdata->p[0]*x_tmp[36]*tdata->dwdp[24];
  tdata->dxdotdp[37 + ip*model->nx] = tdata->p[0]*x_tmp[36]*tdata->dwdp[24];
  tdata->dxdotdp[45 + ip*model->nx] = -tdata->p[0]*x_tmp[45]*tdata->dwdp[24];
  tdata->dxdotdp[46 + ip*model->nx] = tdata->p[0]*x_tmp[45]*tdata->dwdp[24];
  tdata->dxdotdp[54 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[25]-tdata->p[0]*x_tmp[54]*tdata->dwdp[24];
  tdata->dxdotdp[55 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[25]+tdata->p[0]*x_tmp[54]*tdata->dwdp[24];
  tdata->dxdotdp[63 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[26]-tdata->p[0]*x_tmp[63]*tdata->dwdp[24];
  tdata->dxdotdp[64 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[26]+tdata->p[0]*x_tmp[63]*tdata->dwdp[24];
  tdata->dxdotdp[72 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[27]-tdata->p[0]*x_tmp[72]*tdata->dwdp[24];
  tdata->dxdotdp[73 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[27]+tdata->p[0]*x_tmp[72]*tdata->dwdp[24];
  tdata->dxdotdp[81 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[28]-tdata->p[0]*x_tmp[81]*tdata->dwdp[24];
  tdata->dxdotdp[82 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[28]+tdata->p[0]*x_tmp[81]*tdata->dwdp[24];
  tdata->dxdotdp[90 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[29]-tdata->p[0]*x_tmp[90]*tdata->dwdp[24];
  tdata->dxdotdp[91 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[29]+tdata->p[0]*x_tmp[90]*tdata->dwdp[24];
  tdata->dxdotdp[99 + ip*model->nx] = -tdata->p[0]*x_tmp[99]*tdata->dwdp[24];
  tdata->dxdotdp[100 + ip*model->nx] = tdata->p[0]*x_tmp[99]*tdata->dwdp[24];
  tdata->dxdotdp[108 + ip*model->nx] = -tdata->p[0]*x_tmp[108]*tdata->dwdp[24];
  tdata->dxdotdp[109 + ip*model->nx] = tdata->p[0]*x_tmp[108]*tdata->dwdp[24];
  tdata->dxdotdp[117 + ip*model->nx] = -tdata->p[0]*x_tmp[117]*tdata->dwdp[24];
  tdata->dxdotdp[118 + ip*model->nx] = tdata->p[0]*x_tmp[117]*tdata->dwdp[24];
  tdata->dxdotdp[126 + ip*model->nx] = -tdata->p[0]*x_tmp[126]*tdata->dwdp[24];
  tdata->dxdotdp[127 + ip*model->nx] = tdata->p[0]*x_tmp[126]*tdata->dwdp[24];
  tdata->dxdotdp[135 + ip*model->nx] = -tdata->p[0]*x_tmp[135]*tdata->dwdp[24];
  tdata->dxdotdp[136 + ip*model->nx] = tdata->p[0]*x_tmp[135]*tdata->dwdp[24];
  tdata->dxdotdp[144 + ip*model->nx] = -tdata->p[0]*x_tmp[144]*tdata->dwdp[24];
  tdata->dxdotdp[145 + ip*model->nx] = tdata->p[0]*x_tmp[144]*tdata->dwdp[24];
  tdata->dxdotdp[153 + ip*model->nx] = -tdata->p[0]*x_tmp[153]*tdata->dwdp[24];
  tdata->dxdotdp[154 + ip*model->nx] = tdata->p[0]*x_tmp[153]*tdata->dwdp[24];

  } break;

}
}
for(ip = 0; ip<udata->nplist; ip++) {
   for(ix = 0; ix<model->nx; ix++) {
       if(amiIsNaN(tdata->dxdotdp[ix+ip*model->nx])) {
           tdata->dxdotdp[ix+ip*model->nx] = 0;
           if(!tdata->nan_dxdotdp) {
               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");
               tdata->nan_dxdotdp = TRUE;
           }
       }
       if(amiIsInf(tdata->dxdotdp[ix+ip*model->nx])) {
           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");
           return(-1);
       }
   }
}
return(status);

}


