
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_dwdp.h"
#include "model_jakstat_adjoint_o2_w.h"

int dxdotdp_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(udata->dxdotdp,0,sizeof(realtype)*162*udata->nplist);
status = dwdp_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->k[0]*udata->w[0]*x_tmp[0]*udata->w[2];
  udata->dxdotdp[1 + ip*udata->nx] = udata->w[0]*x_tmp[0];
  udata->dxdotdp[9 + ip*udata->nx] = -udata->w[0]*x_tmp[9];
  udata->dxdotdp[10 + ip*udata->nx] = udata->w[0]*x_tmp[9];
  udata->dxdotdp[18 + ip*udata->nx] = -udata->w[0]*x_tmp[18];
  udata->dxdotdp[19 + ip*udata->nx] = udata->w[0]*x_tmp[18];
  udata->dxdotdp[27 + ip*udata->nx] = -udata->w[0]*x_tmp[27];
  udata->dxdotdp[28 + ip*udata->nx] = udata->w[0]*x_tmp[27];
  udata->dxdotdp[36 + ip*udata->nx] = -udata->w[0]*x_tmp[36];
  udata->dxdotdp[37 + ip*udata->nx] = udata->w[0]*x_tmp[36];
  udata->dxdotdp[45 + ip*udata->nx] = -udata->w[0]*x_tmp[45];
  udata->dxdotdp[46 + ip*udata->nx] = udata->w[0]*x_tmp[45];
  udata->dxdotdp[54 + ip*udata->nx] = -x_tmp[0]*udata->w[5]-udata->w[0]*x_tmp[54];
  udata->dxdotdp[55 + ip*udata->nx] = x_tmp[0]*udata->w[5]+udata->w[0]*x_tmp[54];
  udata->dxdotdp[63 + ip*udata->nx] = -x_tmp[0]*udata->w[6]-udata->w[0]*x_tmp[63];
  udata->dxdotdp[64 + ip*udata->nx] = x_tmp[0]*udata->w[6]+udata->w[0]*x_tmp[63];
  udata->dxdotdp[72 + ip*udata->nx] = -x_tmp[0]*udata->w[7]-udata->w[0]*x_tmp[72];
  udata->dxdotdp[73 + ip*udata->nx] = x_tmp[0]*udata->w[7]+udata->w[0]*x_tmp[72];
  udata->dxdotdp[81 + ip*udata->nx] = -x_tmp[0]*udata->w[8]-udata->w[0]*x_tmp[81];
  udata->dxdotdp[82 + ip*udata->nx] = x_tmp[0]*udata->w[8]+udata->w[0]*x_tmp[81];
  udata->dxdotdp[90 + ip*udata->nx] = -x_tmp[0]*udata->w[9]-udata->w[0]*x_tmp[90];
  udata->dxdotdp[91 + ip*udata->nx] = x_tmp[0]*udata->w[9]+udata->w[0]*x_tmp[90];
  udata->dxdotdp[99 + ip*udata->nx] = -udata->w[0]*x_tmp[99];
  udata->dxdotdp[100 + ip*udata->nx] = udata->w[0]*x_tmp[99];
  udata->dxdotdp[108 + ip*udata->nx] = -udata->w[0]*x_tmp[108];
  udata->dxdotdp[109 + ip*udata->nx] = udata->w[0]*x_tmp[108];
  udata->dxdotdp[117 + ip*udata->nx] = -udata->w[0]*x_tmp[117];
  udata->dxdotdp[118 + ip*udata->nx] = udata->w[0]*x_tmp[117];
  udata->dxdotdp[126 + ip*udata->nx] = -udata->w[0]*x_tmp[126];
  udata->dxdotdp[127 + ip*udata->nx] = udata->w[0]*x_tmp[126];
  udata->dxdotdp[135 + ip*udata->nx] = -udata->w[0]*x_tmp[135];
  udata->dxdotdp[136 + ip*udata->nx] = udata->w[0]*x_tmp[135];
  udata->dxdotdp[144 + ip*udata->nx] = -udata->w[0]*x_tmp[144];
  udata->dxdotdp[145 + ip*udata->nx] = udata->w[0]*x_tmp[144];
  udata->dxdotdp[153 + ip*udata->nx] = -udata->w[0]*x_tmp[153];
  udata->dxdotdp[154 + ip*udata->nx] = udata->w[0]*x_tmp[153];

  } break;

  case 1: {
  udata->dxdotdp[1 + ip*udata->nx] = udata->w[1]*-2.0;
  udata->dxdotdp[2 + ip*udata->nx] = udata->w[1];
  udata->dxdotdp[10 + ip*udata->nx] = x_tmp[1]*x_tmp[10]*-4.0;
  udata->dxdotdp[11 + ip*udata->nx] = x_tmp[1]*x_tmp[10]*2.0;
  udata->dxdotdp[19 + ip*udata->nx] = x_tmp[1]*x_tmp[19]*-4.0;
  udata->dxdotdp[20 + ip*udata->nx] = x_tmp[1]*x_tmp[19]*2.0;
  udata->dxdotdp[28 + ip*udata->nx] = x_tmp[1]*x_tmp[28]*-4.0;
  udata->dxdotdp[29 + ip*udata->nx] = x_tmp[1]*x_tmp[28]*2.0;
  udata->dxdotdp[37 + ip*udata->nx] = x_tmp[1]*x_tmp[37]*-4.0;
  udata->dxdotdp[38 + ip*udata->nx] = x_tmp[1]*x_tmp[37]*2.0;
  udata->dxdotdp[46 + ip*udata->nx] = x_tmp[1]*x_tmp[46]*-4.0;
  udata->dxdotdp[47 + ip*udata->nx] = x_tmp[1]*x_tmp[46]*2.0;
  udata->dxdotdp[55 + ip*udata->nx] = x_tmp[1]*x_tmp[55]*-4.0;
  udata->dxdotdp[56 + ip*udata->nx] = x_tmp[1]*x_tmp[55]*2.0;
  udata->dxdotdp[64 + ip*udata->nx] = x_tmp[1]*x_tmp[64]*-4.0;
  udata->dxdotdp[65 + ip*udata->nx] = x_tmp[1]*x_tmp[64]*2.0;
  udata->dxdotdp[73 + ip*udata->nx] = x_tmp[1]*x_tmp[73]*-4.0;
  udata->dxdotdp[74 + ip*udata->nx] = x_tmp[1]*x_tmp[73]*2.0;
  udata->dxdotdp[82 + ip*udata->nx] = x_tmp[1]*x_tmp[82]*-4.0;
  udata->dxdotdp[83 + ip*udata->nx] = x_tmp[1]*x_tmp[82]*2.0;
  udata->dxdotdp[91 + ip*udata->nx] = x_tmp[1]*x_tmp[91]*-4.0;
  udata->dxdotdp[92 + ip*udata->nx] = x_tmp[1]*x_tmp[91]*2.0;
  udata->dxdotdp[100 + ip*udata->nx] = x_tmp[1]*x_tmp[100]*-4.0;
  udata->dxdotdp[101 + ip*udata->nx] = x_tmp[1]*x_tmp[100]*2.0;
  udata->dxdotdp[109 + ip*udata->nx] = x_tmp[1]*x_tmp[109]*-4.0;
  udata->dxdotdp[110 + ip*udata->nx] = x_tmp[1]*x_tmp[109]*2.0;
  udata->dxdotdp[118 + ip*udata->nx] = x_tmp[1]*x_tmp[118]*-4.0;
  udata->dxdotdp[119 + ip*udata->nx] = x_tmp[1]*x_tmp[118]*2.0;
  udata->dxdotdp[127 + ip*udata->nx] = x_tmp[1]*x_tmp[127]*-4.0;
  udata->dxdotdp[128 + ip*udata->nx] = x_tmp[1]*x_tmp[127]*2.0;
  udata->dxdotdp[136 + ip*udata->nx] = x_tmp[1]*x_tmp[136]*-4.0;
  udata->dxdotdp[137 + ip*udata->nx] = x_tmp[1]*x_tmp[136]*2.0;
  udata->dxdotdp[145 + ip*udata->nx] = x_tmp[1]*x_tmp[145]*-4.0;
  udata->dxdotdp[146 + ip*udata->nx] = x_tmp[1]*x_tmp[145]*2.0;
  udata->dxdotdp[154 + ip*udata->nx] = x_tmp[1]*x_tmp[154]*-4.0;
  udata->dxdotdp[155 + ip*udata->nx] = x_tmp[1]*x_tmp[154]*2.0;

  } break;

  case 2: {
  udata->dxdotdp[2 + ip*udata->nx] = -x_tmp[2];
  udata->dxdotdp[3 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[2];
  udata->dxdotdp[11 + ip*udata->nx] = -x_tmp[11];
  udata->dxdotdp[12 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[11];
  udata->dxdotdp[20 + ip*udata->nx] = -x_tmp[20];
  udata->dxdotdp[21 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[20];
  udata->dxdotdp[29 + ip*udata->nx] = -x_tmp[29];
  udata->dxdotdp[30 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[29];
  udata->dxdotdp[38 + ip*udata->nx] = -x_tmp[38];
  udata->dxdotdp[39 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[38];
  udata->dxdotdp[47 + ip*udata->nx] = -x_tmp[47];
  udata->dxdotdp[48 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[47];
  udata->dxdotdp[56 + ip*udata->nx] = -x_tmp[56];
  udata->dxdotdp[57 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[56];
  udata->dxdotdp[65 + ip*udata->nx] = -x_tmp[65];
  udata->dxdotdp[66 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[65];
  udata->dxdotdp[74 + ip*udata->nx] = -x_tmp[74];
  udata->dxdotdp[75 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[74];
  udata->dxdotdp[83 + ip*udata->nx] = -x_tmp[83];
  udata->dxdotdp[84 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[83];
  udata->dxdotdp[92 + ip*udata->nx] = -x_tmp[92];
  udata->dxdotdp[93 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[92];
  udata->dxdotdp[101 + ip*udata->nx] = -x_tmp[101];
  udata->dxdotdp[102 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[101];
  udata->dxdotdp[110 + ip*udata->nx] = -x_tmp[110];
  udata->dxdotdp[111 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[110];
  udata->dxdotdp[119 + ip*udata->nx] = -x_tmp[119];
  udata->dxdotdp[120 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[119];
  udata->dxdotdp[128 + ip*udata->nx] = -x_tmp[128];
  udata->dxdotdp[129 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[128];
  udata->dxdotdp[137 + ip*udata->nx] = -x_tmp[137];
  udata->dxdotdp[138 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[137];
  udata->dxdotdp[146 + ip*udata->nx] = -x_tmp[146];
  udata->dxdotdp[147 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[146];
  udata->dxdotdp[155 + ip*udata->nx] = -x_tmp[155];
  udata->dxdotdp[156 + ip*udata->nx] = udata->k[0]*udata->w[3]*x_tmp[155];

  } break;

  case 3: {
  udata->dxdotdp[0 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[8];
  udata->dxdotdp[3 + ip*udata->nx] = -udata->k[1]*udata->w[3]*x_tmp[3];
  udata->dxdotdp[4 + ip*udata->nx] = udata->w[4]-x_tmp[4];
  udata->dxdotdp[5 + ip*udata->nx] = x_tmp[4]-x_tmp[5];
  udata->dxdotdp[6 + ip*udata->nx] = x_tmp[5]-x_tmp[6];
  udata->dxdotdp[7 + ip*udata->nx] = x_tmp[6]-x_tmp[7];
  udata->dxdotdp[8 + ip*udata->nx] = x_tmp[7]-x_tmp[8];
  udata->dxdotdp[9 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[17];
  udata->dxdotdp[12 + ip*udata->nx] = -x_tmp[12];
  udata->dxdotdp[13 + ip*udata->nx] = x_tmp[12]*2.0-x_tmp[13];
  udata->dxdotdp[14 + ip*udata->nx] = x_tmp[13]-x_tmp[14];
  udata->dxdotdp[15 + ip*udata->nx] = x_tmp[14]-x_tmp[15];
  udata->dxdotdp[16 + ip*udata->nx] = x_tmp[15]-x_tmp[16];
  udata->dxdotdp[17 + ip*udata->nx] = x_tmp[16]-x_tmp[17];
  udata->dxdotdp[18 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[26];
  udata->dxdotdp[21 + ip*udata->nx] = -x_tmp[21];
  udata->dxdotdp[22 + ip*udata->nx] = x_tmp[21]*2.0-x_tmp[22];
  udata->dxdotdp[23 + ip*udata->nx] = x_tmp[22]-x_tmp[23];
  udata->dxdotdp[24 + ip*udata->nx] = x_tmp[23]-x_tmp[24];
  udata->dxdotdp[25 + ip*udata->nx] = x_tmp[24]-x_tmp[25];
  udata->dxdotdp[26 + ip*udata->nx] = x_tmp[25]-x_tmp[26];
  udata->dxdotdp[27 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[35];
  udata->dxdotdp[30 + ip*udata->nx] = -x_tmp[30];
  udata->dxdotdp[31 + ip*udata->nx] = x_tmp[30]*2.0-x_tmp[31];
  udata->dxdotdp[32 + ip*udata->nx] = x_tmp[31]-x_tmp[32];
  udata->dxdotdp[33 + ip*udata->nx] = x_tmp[32]-x_tmp[33];
  udata->dxdotdp[34 + ip*udata->nx] = x_tmp[33]-x_tmp[34];
  udata->dxdotdp[35 + ip*udata->nx] = x_tmp[34]-x_tmp[35];
  udata->dxdotdp[36 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[44];
  udata->dxdotdp[39 + ip*udata->nx] = -x_tmp[39];
  udata->dxdotdp[40 + ip*udata->nx] = x_tmp[39]*2.0-x_tmp[40];
  udata->dxdotdp[41 + ip*udata->nx] = x_tmp[40]-x_tmp[41];
  udata->dxdotdp[42 + ip*udata->nx] = x_tmp[41]-x_tmp[42];
  udata->dxdotdp[43 + ip*udata->nx] = x_tmp[42]-x_tmp[43];
  udata->dxdotdp[44 + ip*udata->nx] = x_tmp[43]-x_tmp[44];
  udata->dxdotdp[45 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[53];
  udata->dxdotdp[48 + ip*udata->nx] = -x_tmp[48];
  udata->dxdotdp[49 + ip*udata->nx] = x_tmp[48]*2.0-x_tmp[49];
  udata->dxdotdp[50 + ip*udata->nx] = x_tmp[49]-x_tmp[50];
  udata->dxdotdp[51 + ip*udata->nx] = x_tmp[50]-x_tmp[51];
  udata->dxdotdp[52 + ip*udata->nx] = x_tmp[51]-x_tmp[52];
  udata->dxdotdp[53 + ip*udata->nx] = x_tmp[52]-x_tmp[53];
  udata->dxdotdp[54 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[62];
  udata->dxdotdp[57 + ip*udata->nx] = -x_tmp[57];
  udata->dxdotdp[58 + ip*udata->nx] = x_tmp[57]*2.0-x_tmp[58];
  udata->dxdotdp[59 + ip*udata->nx] = x_tmp[58]-x_tmp[59];
  udata->dxdotdp[60 + ip*udata->nx] = x_tmp[59]-x_tmp[60];
  udata->dxdotdp[61 + ip*udata->nx] = x_tmp[60]-x_tmp[61];
  udata->dxdotdp[62 + ip*udata->nx] = x_tmp[61]-x_tmp[62];
  udata->dxdotdp[63 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[71];
  udata->dxdotdp[66 + ip*udata->nx] = -x_tmp[66];
  udata->dxdotdp[67 + ip*udata->nx] = x_tmp[66]*2.0-x_tmp[67];
  udata->dxdotdp[68 + ip*udata->nx] = x_tmp[67]-x_tmp[68];
  udata->dxdotdp[69 + ip*udata->nx] = x_tmp[68]-x_tmp[69];
  udata->dxdotdp[70 + ip*udata->nx] = x_tmp[69]-x_tmp[70];
  udata->dxdotdp[71 + ip*udata->nx] = x_tmp[70]-x_tmp[71];
  udata->dxdotdp[72 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[80];
  udata->dxdotdp[75 + ip*udata->nx] = -x_tmp[75];
  udata->dxdotdp[76 + ip*udata->nx] = x_tmp[75]*2.0-x_tmp[76];
  udata->dxdotdp[77 + ip*udata->nx] = x_tmp[76]-x_tmp[77];
  udata->dxdotdp[78 + ip*udata->nx] = x_tmp[77]-x_tmp[78];
  udata->dxdotdp[79 + ip*udata->nx] = x_tmp[78]-x_tmp[79];
  udata->dxdotdp[80 + ip*udata->nx] = x_tmp[79]-x_tmp[80];
  udata->dxdotdp[81 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[89];
  udata->dxdotdp[84 + ip*udata->nx] = -x_tmp[84];
  udata->dxdotdp[85 + ip*udata->nx] = x_tmp[84]*2.0-x_tmp[85];
  udata->dxdotdp[86 + ip*udata->nx] = x_tmp[85]-x_tmp[86];
  udata->dxdotdp[87 + ip*udata->nx] = x_tmp[86]-x_tmp[87];
  udata->dxdotdp[88 + ip*udata->nx] = x_tmp[87]-x_tmp[88];
  udata->dxdotdp[89 + ip*udata->nx] = x_tmp[88]-x_tmp[89];
  udata->dxdotdp[90 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[98];
  udata->dxdotdp[93 + ip*udata->nx] = -x_tmp[93];
  udata->dxdotdp[94 + ip*udata->nx] = x_tmp[93]*2.0-x_tmp[94];
  udata->dxdotdp[95 + ip*udata->nx] = x_tmp[94]-x_tmp[95];
  udata->dxdotdp[96 + ip*udata->nx] = x_tmp[95]-x_tmp[96];
  udata->dxdotdp[97 + ip*udata->nx] = x_tmp[96]-x_tmp[97];
  udata->dxdotdp[98 + ip*udata->nx] = x_tmp[97]-x_tmp[98];
  udata->dxdotdp[99 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[107];
  udata->dxdotdp[102 + ip*udata->nx] = -x_tmp[102];
  udata->dxdotdp[103 + ip*udata->nx] = x_tmp[102]*2.0-x_tmp[103];
  udata->dxdotdp[104 + ip*udata->nx] = x_tmp[103]-x_tmp[104];
  udata->dxdotdp[105 + ip*udata->nx] = x_tmp[104]-x_tmp[105];
  udata->dxdotdp[106 + ip*udata->nx] = x_tmp[105]-x_tmp[106];
  udata->dxdotdp[107 + ip*udata->nx] = x_tmp[106]-x_tmp[107];
  udata->dxdotdp[108 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[116];
  udata->dxdotdp[111 + ip*udata->nx] = -x_tmp[111];
  udata->dxdotdp[112 + ip*udata->nx] = x_tmp[111]*2.0-x_tmp[112];
  udata->dxdotdp[113 + ip*udata->nx] = x_tmp[112]-x_tmp[113];
  udata->dxdotdp[114 + ip*udata->nx] = x_tmp[113]-x_tmp[114];
  udata->dxdotdp[115 + ip*udata->nx] = x_tmp[114]-x_tmp[115];
  udata->dxdotdp[116 + ip*udata->nx] = x_tmp[115]-x_tmp[116];
  udata->dxdotdp[117 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[125];
  udata->dxdotdp[120 + ip*udata->nx] = -x_tmp[120];
  udata->dxdotdp[121 + ip*udata->nx] = x_tmp[120]*2.0-x_tmp[121];
  udata->dxdotdp[122 + ip*udata->nx] = x_tmp[121]-x_tmp[122];
  udata->dxdotdp[123 + ip*udata->nx] = x_tmp[122]-x_tmp[123];
  udata->dxdotdp[124 + ip*udata->nx] = x_tmp[123]-x_tmp[124];
  udata->dxdotdp[125 + ip*udata->nx] = x_tmp[124]-x_tmp[125];
  udata->dxdotdp[126 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[134];
  udata->dxdotdp[129 + ip*udata->nx] = -x_tmp[129];
  udata->dxdotdp[130 + ip*udata->nx] = x_tmp[129]*2.0-x_tmp[130];
  udata->dxdotdp[131 + ip*udata->nx] = x_tmp[130]-x_tmp[131];
  udata->dxdotdp[132 + ip*udata->nx] = x_tmp[131]-x_tmp[132];
  udata->dxdotdp[133 + ip*udata->nx] = x_tmp[132]-x_tmp[133];
  udata->dxdotdp[134 + ip*udata->nx] = x_tmp[133]-x_tmp[134];
  udata->dxdotdp[135 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[143];
  udata->dxdotdp[138 + ip*udata->nx] = -x_tmp[138];
  udata->dxdotdp[139 + ip*udata->nx] = x_tmp[138]*2.0-x_tmp[139];
  udata->dxdotdp[140 + ip*udata->nx] = x_tmp[139]-x_tmp[140];
  udata->dxdotdp[141 + ip*udata->nx] = x_tmp[140]-x_tmp[141];
  udata->dxdotdp[142 + ip*udata->nx] = x_tmp[141]-x_tmp[142];
  udata->dxdotdp[143 + ip*udata->nx] = x_tmp[142]-x_tmp[143];
  udata->dxdotdp[144 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[152];
  udata->dxdotdp[147 + ip*udata->nx] = -x_tmp[147];
  udata->dxdotdp[148 + ip*udata->nx] = x_tmp[147]*2.0-x_tmp[148];
  udata->dxdotdp[149 + ip*udata->nx] = x_tmp[148]-x_tmp[149];
  udata->dxdotdp[150 + ip*udata->nx] = x_tmp[149]-x_tmp[150];
  udata->dxdotdp[151 + ip*udata->nx] = x_tmp[150]-x_tmp[151];
  udata->dxdotdp[152 + ip*udata->nx] = x_tmp[151]-x_tmp[152];
  udata->dxdotdp[153 + ip*udata->nx] = udata->k[1]*udata->w[2]*x_tmp[161];
  udata->dxdotdp[156 + ip*udata->nx] = -x_tmp[156];
  udata->dxdotdp[157 + ip*udata->nx] = x_tmp[156]*2.0-x_tmp[157];
  udata->dxdotdp[158 + ip*udata->nx] = x_tmp[157]-x_tmp[158];
  udata->dxdotdp[159 + ip*udata->nx] = x_tmp[158]-x_tmp[159];
  udata->dxdotdp[160 + ip*udata->nx] = x_tmp[159]-x_tmp[160];
  udata->dxdotdp[161 + ip*udata->nx] = x_tmp[160]-x_tmp[161];

  } break;

  case 5: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->k[0]*udata->p[0]*x_tmp[0]*udata->w[2]*udata->dwdp[0];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[0];
  udata->dxdotdp[9 + ip*udata->nx] = -udata->dwdp[0]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[10 + ip*udata->nx] = udata->dwdp[0]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[18 + ip*udata->nx] = -udata->p[0]*x_tmp[18]*udata->dwdp[0];
  udata->dxdotdp[19 + ip*udata->nx] = udata->p[0]*x_tmp[18]*udata->dwdp[0];
  udata->dxdotdp[27 + ip*udata->nx] = -udata->p[0]*x_tmp[27]*udata->dwdp[0];
  udata->dxdotdp[28 + ip*udata->nx] = udata->p[0]*x_tmp[27]*udata->dwdp[0];
  udata->dxdotdp[36 + ip*udata->nx] = -udata->p[0]*x_tmp[36]*udata->dwdp[0];
  udata->dxdotdp[37 + ip*udata->nx] = udata->p[0]*x_tmp[36]*udata->dwdp[0];
  udata->dxdotdp[45 + ip*udata->nx] = -udata->p[0]*x_tmp[45]*udata->dwdp[0];
  udata->dxdotdp[46 + ip*udata->nx] = udata->p[0]*x_tmp[45]*udata->dwdp[0];
  udata->dxdotdp[54 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[1]-udata->p[0]*x_tmp[54]*udata->dwdp[0];
  udata->dxdotdp[55 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[1]+udata->p[0]*x_tmp[54]*udata->dwdp[0];
  udata->dxdotdp[63 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[2]-udata->p[0]*x_tmp[63]*udata->dwdp[0];
  udata->dxdotdp[64 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[2]+udata->p[0]*x_tmp[63]*udata->dwdp[0];
  udata->dxdotdp[72 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[3]-udata->p[0]*x_tmp[72]*udata->dwdp[0];
  udata->dxdotdp[73 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[3]+udata->p[0]*x_tmp[72]*udata->dwdp[0];
  udata->dxdotdp[81 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[4]-udata->p[0]*x_tmp[81]*udata->dwdp[0];
  udata->dxdotdp[82 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[4]+udata->p[0]*x_tmp[81]*udata->dwdp[0];
  udata->dxdotdp[90 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[5]-udata->p[0]*x_tmp[90]*udata->dwdp[0];
  udata->dxdotdp[91 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[5]+udata->p[0]*x_tmp[90]*udata->dwdp[0];
  udata->dxdotdp[99 + ip*udata->nx] = -udata->p[0]*x_tmp[99]*udata->dwdp[0];
  udata->dxdotdp[100 + ip*udata->nx] = udata->p[0]*x_tmp[99]*udata->dwdp[0];
  udata->dxdotdp[108 + ip*udata->nx] = -udata->p[0]*x_tmp[108]*udata->dwdp[0];
  udata->dxdotdp[109 + ip*udata->nx] = udata->p[0]*x_tmp[108]*udata->dwdp[0];
  udata->dxdotdp[117 + ip*udata->nx] = -udata->p[0]*x_tmp[117]*udata->dwdp[0];
  udata->dxdotdp[118 + ip*udata->nx] = udata->p[0]*x_tmp[117]*udata->dwdp[0];
  udata->dxdotdp[126 + ip*udata->nx] = -udata->p[0]*x_tmp[126]*udata->dwdp[0];
  udata->dxdotdp[127 + ip*udata->nx] = udata->p[0]*x_tmp[126]*udata->dwdp[0];
  udata->dxdotdp[135 + ip*udata->nx] = -udata->p[0]*x_tmp[135]*udata->dwdp[0];
  udata->dxdotdp[136 + ip*udata->nx] = udata->p[0]*x_tmp[135]*udata->dwdp[0];
  udata->dxdotdp[144 + ip*udata->nx] = -udata->p[0]*x_tmp[144]*udata->dwdp[0];
  udata->dxdotdp[145 + ip*udata->nx] = udata->p[0]*x_tmp[144]*udata->dwdp[0];
  udata->dxdotdp[153 + ip*udata->nx] = -udata->p[0]*x_tmp[153]*udata->dwdp[0];
  udata->dxdotdp[154 + ip*udata->nx] = udata->p[0]*x_tmp[153]*udata->dwdp[0];

  } break;

  case 6: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->k[0]*udata->p[0]*x_tmp[0]*udata->w[2]*udata->dwdp[6];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[6];
  udata->dxdotdp[9 + ip*udata->nx] = -udata->dwdp[6]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[10 + ip*udata->nx] = udata->dwdp[6]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[18 + ip*udata->nx] = -udata->p[0]*x_tmp[18]*udata->dwdp[6];
  udata->dxdotdp[19 + ip*udata->nx] = udata->p[0]*x_tmp[18]*udata->dwdp[6];
  udata->dxdotdp[27 + ip*udata->nx] = -udata->p[0]*x_tmp[27]*udata->dwdp[6];
  udata->dxdotdp[28 + ip*udata->nx] = udata->p[0]*x_tmp[27]*udata->dwdp[6];
  udata->dxdotdp[36 + ip*udata->nx] = -udata->p[0]*x_tmp[36]*udata->dwdp[6];
  udata->dxdotdp[37 + ip*udata->nx] = udata->p[0]*x_tmp[36]*udata->dwdp[6];
  udata->dxdotdp[45 + ip*udata->nx] = -udata->p[0]*x_tmp[45]*udata->dwdp[6];
  udata->dxdotdp[46 + ip*udata->nx] = udata->p[0]*x_tmp[45]*udata->dwdp[6];
  udata->dxdotdp[54 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[7]-udata->p[0]*x_tmp[54]*udata->dwdp[6];
  udata->dxdotdp[55 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[7]+udata->p[0]*x_tmp[54]*udata->dwdp[6];
  udata->dxdotdp[63 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[8]-udata->p[0]*x_tmp[63]*udata->dwdp[6];
  udata->dxdotdp[64 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[8]+udata->p[0]*x_tmp[63]*udata->dwdp[6];
  udata->dxdotdp[72 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[9]-udata->p[0]*x_tmp[72]*udata->dwdp[6];
  udata->dxdotdp[73 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[9]+udata->p[0]*x_tmp[72]*udata->dwdp[6];
  udata->dxdotdp[81 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[10]-udata->p[0]*x_tmp[81]*udata->dwdp[6];
  udata->dxdotdp[82 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[10]+udata->p[0]*x_tmp[81]*udata->dwdp[6];
  udata->dxdotdp[90 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[11]-udata->p[0]*x_tmp[90]*udata->dwdp[6];
  udata->dxdotdp[91 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[11]+udata->p[0]*x_tmp[90]*udata->dwdp[6];
  udata->dxdotdp[99 + ip*udata->nx] = -udata->p[0]*x_tmp[99]*udata->dwdp[6];
  udata->dxdotdp[100 + ip*udata->nx] = udata->p[0]*x_tmp[99]*udata->dwdp[6];
  udata->dxdotdp[108 + ip*udata->nx] = -udata->p[0]*x_tmp[108]*udata->dwdp[6];
  udata->dxdotdp[109 + ip*udata->nx] = udata->p[0]*x_tmp[108]*udata->dwdp[6];
  udata->dxdotdp[117 + ip*udata->nx] = -udata->p[0]*x_tmp[117]*udata->dwdp[6];
  udata->dxdotdp[118 + ip*udata->nx] = udata->p[0]*x_tmp[117]*udata->dwdp[6];
  udata->dxdotdp[126 + ip*udata->nx] = -udata->p[0]*x_tmp[126]*udata->dwdp[6];
  udata->dxdotdp[127 + ip*udata->nx] = udata->p[0]*x_tmp[126]*udata->dwdp[6];
  udata->dxdotdp[135 + ip*udata->nx] = -udata->p[0]*x_tmp[135]*udata->dwdp[6];
  udata->dxdotdp[136 + ip*udata->nx] = udata->p[0]*x_tmp[135]*udata->dwdp[6];
  udata->dxdotdp[144 + ip*udata->nx] = -udata->p[0]*x_tmp[144]*udata->dwdp[6];
  udata->dxdotdp[145 + ip*udata->nx] = udata->p[0]*x_tmp[144]*udata->dwdp[6];
  udata->dxdotdp[153 + ip*udata->nx] = -udata->p[0]*x_tmp[153]*udata->dwdp[6];
  udata->dxdotdp[154 + ip*udata->nx] = udata->p[0]*x_tmp[153]*udata->dwdp[6];

  } break;

  case 7: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->k[0]*udata->p[0]*x_tmp[0]*udata->w[2]*udata->dwdp[12];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[12];
  udata->dxdotdp[9 + ip*udata->nx] = -udata->dwdp[12]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[10 + ip*udata->nx] = udata->dwdp[12]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[18 + ip*udata->nx] = -udata->p[0]*x_tmp[18]*udata->dwdp[12];
  udata->dxdotdp[19 + ip*udata->nx] = udata->p[0]*x_tmp[18]*udata->dwdp[12];
  udata->dxdotdp[27 + ip*udata->nx] = -udata->p[0]*x_tmp[27]*udata->dwdp[12];
  udata->dxdotdp[28 + ip*udata->nx] = udata->p[0]*x_tmp[27]*udata->dwdp[12];
  udata->dxdotdp[36 + ip*udata->nx] = -udata->p[0]*x_tmp[36]*udata->dwdp[12];
  udata->dxdotdp[37 + ip*udata->nx] = udata->p[0]*x_tmp[36]*udata->dwdp[12];
  udata->dxdotdp[45 + ip*udata->nx] = -udata->p[0]*x_tmp[45]*udata->dwdp[12];
  udata->dxdotdp[46 + ip*udata->nx] = udata->p[0]*x_tmp[45]*udata->dwdp[12];
  udata->dxdotdp[54 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[13]-udata->p[0]*x_tmp[54]*udata->dwdp[12];
  udata->dxdotdp[55 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[13]+udata->p[0]*x_tmp[54]*udata->dwdp[12];
  udata->dxdotdp[63 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[14]-udata->p[0]*x_tmp[63]*udata->dwdp[12];
  udata->dxdotdp[64 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[14]+udata->p[0]*x_tmp[63]*udata->dwdp[12];
  udata->dxdotdp[72 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[15]-udata->p[0]*x_tmp[72]*udata->dwdp[12];
  udata->dxdotdp[73 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[15]+udata->p[0]*x_tmp[72]*udata->dwdp[12];
  udata->dxdotdp[81 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[16]-udata->p[0]*x_tmp[81]*udata->dwdp[12];
  udata->dxdotdp[82 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[16]+udata->p[0]*x_tmp[81]*udata->dwdp[12];
  udata->dxdotdp[90 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[17]-udata->p[0]*x_tmp[90]*udata->dwdp[12];
  udata->dxdotdp[91 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[17]+udata->p[0]*x_tmp[90]*udata->dwdp[12];
  udata->dxdotdp[99 + ip*udata->nx] = -udata->p[0]*x_tmp[99]*udata->dwdp[12];
  udata->dxdotdp[100 + ip*udata->nx] = udata->p[0]*x_tmp[99]*udata->dwdp[12];
  udata->dxdotdp[108 + ip*udata->nx] = -udata->p[0]*x_tmp[108]*udata->dwdp[12];
  udata->dxdotdp[109 + ip*udata->nx] = udata->p[0]*x_tmp[108]*udata->dwdp[12];
  udata->dxdotdp[117 + ip*udata->nx] = -udata->p[0]*x_tmp[117]*udata->dwdp[12];
  udata->dxdotdp[118 + ip*udata->nx] = udata->p[0]*x_tmp[117]*udata->dwdp[12];
  udata->dxdotdp[126 + ip*udata->nx] = -udata->p[0]*x_tmp[126]*udata->dwdp[12];
  udata->dxdotdp[127 + ip*udata->nx] = udata->p[0]*x_tmp[126]*udata->dwdp[12];
  udata->dxdotdp[135 + ip*udata->nx] = -udata->p[0]*x_tmp[135]*udata->dwdp[12];
  udata->dxdotdp[136 + ip*udata->nx] = udata->p[0]*x_tmp[135]*udata->dwdp[12];
  udata->dxdotdp[144 + ip*udata->nx] = -udata->p[0]*x_tmp[144]*udata->dwdp[12];
  udata->dxdotdp[145 + ip*udata->nx] = udata->p[0]*x_tmp[144]*udata->dwdp[12];
  udata->dxdotdp[153 + ip*udata->nx] = -udata->p[0]*x_tmp[153]*udata->dwdp[12];
  udata->dxdotdp[154 + ip*udata->nx] = udata->p[0]*x_tmp[153]*udata->dwdp[12];

  } break;

  case 8: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->k[0]*udata->p[0]*x_tmp[0]*udata->w[2]*udata->dwdp[18];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[18];
  udata->dxdotdp[9 + ip*udata->nx] = -udata->dwdp[18]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[10 + ip*udata->nx] = udata->dwdp[18]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[18 + ip*udata->nx] = -udata->p[0]*x_tmp[18]*udata->dwdp[18];
  udata->dxdotdp[19 + ip*udata->nx] = udata->p[0]*x_tmp[18]*udata->dwdp[18];
  udata->dxdotdp[27 + ip*udata->nx] = -udata->p[0]*x_tmp[27]*udata->dwdp[18];
  udata->dxdotdp[28 + ip*udata->nx] = udata->p[0]*x_tmp[27]*udata->dwdp[18];
  udata->dxdotdp[36 + ip*udata->nx] = -udata->p[0]*x_tmp[36]*udata->dwdp[18];
  udata->dxdotdp[37 + ip*udata->nx] = udata->p[0]*x_tmp[36]*udata->dwdp[18];
  udata->dxdotdp[45 + ip*udata->nx] = -udata->p[0]*x_tmp[45]*udata->dwdp[18];
  udata->dxdotdp[46 + ip*udata->nx] = udata->p[0]*x_tmp[45]*udata->dwdp[18];
  udata->dxdotdp[54 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[19]-udata->p[0]*x_tmp[54]*udata->dwdp[18];
  udata->dxdotdp[55 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[19]+udata->p[0]*x_tmp[54]*udata->dwdp[18];
  udata->dxdotdp[63 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[20]-udata->p[0]*x_tmp[63]*udata->dwdp[18];
  udata->dxdotdp[64 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[20]+udata->p[0]*x_tmp[63]*udata->dwdp[18];
  udata->dxdotdp[72 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[21]-udata->p[0]*x_tmp[72]*udata->dwdp[18];
  udata->dxdotdp[73 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[21]+udata->p[0]*x_tmp[72]*udata->dwdp[18];
  udata->dxdotdp[81 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[22]-udata->p[0]*x_tmp[81]*udata->dwdp[18];
  udata->dxdotdp[82 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[22]+udata->p[0]*x_tmp[81]*udata->dwdp[18];
  udata->dxdotdp[90 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[23]-udata->p[0]*x_tmp[90]*udata->dwdp[18];
  udata->dxdotdp[91 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[23]+udata->p[0]*x_tmp[90]*udata->dwdp[18];
  udata->dxdotdp[99 + ip*udata->nx] = -udata->p[0]*x_tmp[99]*udata->dwdp[18];
  udata->dxdotdp[100 + ip*udata->nx] = udata->p[0]*x_tmp[99]*udata->dwdp[18];
  udata->dxdotdp[108 + ip*udata->nx] = -udata->p[0]*x_tmp[108]*udata->dwdp[18];
  udata->dxdotdp[109 + ip*udata->nx] = udata->p[0]*x_tmp[108]*udata->dwdp[18];
  udata->dxdotdp[117 + ip*udata->nx] = -udata->p[0]*x_tmp[117]*udata->dwdp[18];
  udata->dxdotdp[118 + ip*udata->nx] = udata->p[0]*x_tmp[117]*udata->dwdp[18];
  udata->dxdotdp[126 + ip*udata->nx] = -udata->p[0]*x_tmp[126]*udata->dwdp[18];
  udata->dxdotdp[127 + ip*udata->nx] = udata->p[0]*x_tmp[126]*udata->dwdp[18];
  udata->dxdotdp[135 + ip*udata->nx] = -udata->p[0]*x_tmp[135]*udata->dwdp[18];
  udata->dxdotdp[136 + ip*udata->nx] = udata->p[0]*x_tmp[135]*udata->dwdp[18];
  udata->dxdotdp[144 + ip*udata->nx] = -udata->p[0]*x_tmp[144]*udata->dwdp[18];
  udata->dxdotdp[145 + ip*udata->nx] = udata->p[0]*x_tmp[144]*udata->dwdp[18];
  udata->dxdotdp[153 + ip*udata->nx] = -udata->p[0]*x_tmp[153]*udata->dwdp[18];
  udata->dxdotdp[154 + ip*udata->nx] = udata->p[0]*x_tmp[153]*udata->dwdp[18];

  } break;

  case 9: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->k[0]*udata->p[0]*x_tmp[0]*udata->w[2]*udata->dwdp[24];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[24];
  udata->dxdotdp[9 + ip*udata->nx] = -udata->dwdp[24]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[10 + ip*udata->nx] = udata->dwdp[24]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  udata->dxdotdp[18 + ip*udata->nx] = -udata->p[0]*x_tmp[18]*udata->dwdp[24];
  udata->dxdotdp[19 + ip*udata->nx] = udata->p[0]*x_tmp[18]*udata->dwdp[24];
  udata->dxdotdp[27 + ip*udata->nx] = -udata->p[0]*x_tmp[27]*udata->dwdp[24];
  udata->dxdotdp[28 + ip*udata->nx] = udata->p[0]*x_tmp[27]*udata->dwdp[24];
  udata->dxdotdp[36 + ip*udata->nx] = -udata->p[0]*x_tmp[36]*udata->dwdp[24];
  udata->dxdotdp[37 + ip*udata->nx] = udata->p[0]*x_tmp[36]*udata->dwdp[24];
  udata->dxdotdp[45 + ip*udata->nx] = -udata->p[0]*x_tmp[45]*udata->dwdp[24];
  udata->dxdotdp[46 + ip*udata->nx] = udata->p[0]*x_tmp[45]*udata->dwdp[24];
  udata->dxdotdp[54 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[25]-udata->p[0]*x_tmp[54]*udata->dwdp[24];
  udata->dxdotdp[55 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[25]+udata->p[0]*x_tmp[54]*udata->dwdp[24];
  udata->dxdotdp[63 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[26]-udata->p[0]*x_tmp[63]*udata->dwdp[24];
  udata->dxdotdp[64 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[26]+udata->p[0]*x_tmp[63]*udata->dwdp[24];
  udata->dxdotdp[72 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[27]-udata->p[0]*x_tmp[72]*udata->dwdp[24];
  udata->dxdotdp[73 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[27]+udata->p[0]*x_tmp[72]*udata->dwdp[24];
  udata->dxdotdp[81 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[28]-udata->p[0]*x_tmp[81]*udata->dwdp[24];
  udata->dxdotdp[82 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[28]+udata->p[0]*x_tmp[81]*udata->dwdp[24];
  udata->dxdotdp[90 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[29]-udata->p[0]*x_tmp[90]*udata->dwdp[24];
  udata->dxdotdp[91 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[29]+udata->p[0]*x_tmp[90]*udata->dwdp[24];
  udata->dxdotdp[99 + ip*udata->nx] = -udata->p[0]*x_tmp[99]*udata->dwdp[24];
  udata->dxdotdp[100 + ip*udata->nx] = udata->p[0]*x_tmp[99]*udata->dwdp[24];
  udata->dxdotdp[108 + ip*udata->nx] = -udata->p[0]*x_tmp[108]*udata->dwdp[24];
  udata->dxdotdp[109 + ip*udata->nx] = udata->p[0]*x_tmp[108]*udata->dwdp[24];
  udata->dxdotdp[117 + ip*udata->nx] = -udata->p[0]*x_tmp[117]*udata->dwdp[24];
  udata->dxdotdp[118 + ip*udata->nx] = udata->p[0]*x_tmp[117]*udata->dwdp[24];
  udata->dxdotdp[126 + ip*udata->nx] = -udata->p[0]*x_tmp[126]*udata->dwdp[24];
  udata->dxdotdp[127 + ip*udata->nx] = udata->p[0]*x_tmp[126]*udata->dwdp[24];
  udata->dxdotdp[135 + ip*udata->nx] = -udata->p[0]*x_tmp[135]*udata->dwdp[24];
  udata->dxdotdp[136 + ip*udata->nx] = udata->p[0]*x_tmp[135]*udata->dwdp[24];
  udata->dxdotdp[144 + ip*udata->nx] = -udata->p[0]*x_tmp[144]*udata->dwdp[24];
  udata->dxdotdp[145 + ip*udata->nx] = udata->p[0]*x_tmp[144]*udata->dwdp[24];
  udata->dxdotdp[153 + ip*udata->nx] = -udata->p[0]*x_tmp[153]*udata->dwdp[24];
  udata->dxdotdp[154 + ip*udata->nx] = udata->p[0]*x_tmp[153]*udata->dwdp[24];

  } break;

}
}
for(ip = 0; ip<udata->nplist; ip++) {
   for(ix = 0; ix<udata->nx; ix++) {
       if(amiIsNaN(udata->dxdotdp[ix+ip*udata->nx])) {
           udata->dxdotdp[ix+ip*udata->nx] = 0;
           if(!udata->nan_dxdotdp) {
               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");
               udata->nan_dxdotdp = TRUE;
           }
       }
       if(amiIsInf(udata->dxdotdp[ix+ip*udata->nx])) {
           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");
           return(-1);
       }
   }
}
return(status);

}


