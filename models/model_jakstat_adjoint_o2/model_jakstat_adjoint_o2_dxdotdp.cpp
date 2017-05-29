
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_dwdp.h"
#include "model_jakstat_adjoint_o2_w.h"

int dxdotdp_model_jakstat_adjoint_o2(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(dxdotdp,0,sizeof(realtype)*162*udata->nplist);
status = dwdp_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  dxdotdp[(0+ip * udata->nx)] = -udata->k[0]*udata->w[0]*udata->w[2]*x_tmp[0];
  dxdotdp[(1+ip * udata->nx)] = udata->w[0]*x_tmp[0];
  dxdotdp[(9+ip * udata->nx)] = -udata->w[0]*x_tmp[9];
  dxdotdp[(10+ip * udata->nx)] = udata->w[0]*x_tmp[9];
  dxdotdp[(18+ip * udata->nx)] = -udata->w[0]*x_tmp[18];
  dxdotdp[(19+ip * udata->nx)] = udata->w[0]*x_tmp[18];
  dxdotdp[(27+ip * udata->nx)] = -udata->w[0]*x_tmp[27];
  dxdotdp[(28+ip * udata->nx)] = udata->w[0]*x_tmp[27];
  dxdotdp[(36+ip * udata->nx)] = -udata->w[0]*x_tmp[36];
  dxdotdp[(37+ip * udata->nx)] = udata->w[0]*x_tmp[36];
  dxdotdp[(45+ip * udata->nx)] = -udata->w[0]*x_tmp[45];
  dxdotdp[(46+ip * udata->nx)] = udata->w[0]*x_tmp[45];
  dxdotdp[(54+ip * udata->nx)] = -udata->w[5]*x_tmp[0]-udata->w[0]*x_tmp[54];
  dxdotdp[(55+ip * udata->nx)] = udata->w[5]*x_tmp[0]+udata->w[0]*x_tmp[54];
  dxdotdp[(63+ip * udata->nx)] = -udata->w[6]*x_tmp[0]-udata->w[0]*x_tmp[63];
  dxdotdp[(64+ip * udata->nx)] = udata->w[6]*x_tmp[0]+udata->w[0]*x_tmp[63];
  dxdotdp[(72+ip * udata->nx)] = -udata->w[7]*x_tmp[0]-udata->w[0]*x_tmp[72];
  dxdotdp[(73+ip * udata->nx)] = udata->w[7]*x_tmp[0]+udata->w[0]*x_tmp[72];
  dxdotdp[(81+ip * udata->nx)] = -udata->w[8]*x_tmp[0]-udata->w[0]*x_tmp[81];
  dxdotdp[(82+ip * udata->nx)] = udata->w[8]*x_tmp[0]+udata->w[0]*x_tmp[81];
  dxdotdp[(90+ip * udata->nx)] = -udata->w[9]*x_tmp[0]-udata->w[0]*x_tmp[90];
  dxdotdp[(91+ip * udata->nx)] = udata->w[9]*x_tmp[0]+udata->w[0]*x_tmp[90];
  dxdotdp[(99+ip * udata->nx)] = -udata->w[0]*x_tmp[99];
  dxdotdp[(100+ip * udata->nx)] = udata->w[0]*x_tmp[99];
  dxdotdp[(108+ip * udata->nx)] = -udata->w[0]*x_tmp[108];
  dxdotdp[(109+ip * udata->nx)] = udata->w[0]*x_tmp[108];
  dxdotdp[(117+ip * udata->nx)] = -udata->w[0]*x_tmp[117];
  dxdotdp[(118+ip * udata->nx)] = udata->w[0]*x_tmp[117];
  dxdotdp[(126+ip * udata->nx)] = -udata->w[0]*x_tmp[126];
  dxdotdp[(127+ip * udata->nx)] = udata->w[0]*x_tmp[126];
  dxdotdp[(135+ip * udata->nx)] = -udata->w[0]*x_tmp[135];
  dxdotdp[(136+ip * udata->nx)] = udata->w[0]*x_tmp[135];
  dxdotdp[(144+ip * udata->nx)] = -udata->w[0]*x_tmp[144];
  dxdotdp[(145+ip * udata->nx)] = udata->w[0]*x_tmp[144];
  dxdotdp[(153+ip * udata->nx)] = -udata->w[0]*x_tmp[153];
  dxdotdp[(154+ip * udata->nx)] = udata->w[0]*x_tmp[153];

  } break;

  case 1: {
  dxdotdp[(1+ip * udata->nx)] = udata->w[1]*-2.0;
  dxdotdp[(2+ip * udata->nx)] = udata->w[1];
  dxdotdp[(10+ip * udata->nx)] = x_tmp[1]*x_tmp[10]*-4.0;
  dxdotdp[(11+ip * udata->nx)] = x_tmp[1]*x_tmp[10]*2.0;
  dxdotdp[(19+ip * udata->nx)] = x_tmp[1]*x_tmp[19]*-4.0;
  dxdotdp[(20+ip * udata->nx)] = x_tmp[1]*x_tmp[19]*2.0;
  dxdotdp[(28+ip * udata->nx)] = x_tmp[1]*x_tmp[28]*-4.0;
  dxdotdp[(29+ip * udata->nx)] = x_tmp[1]*x_tmp[28]*2.0;
  dxdotdp[(37+ip * udata->nx)] = x_tmp[1]*x_tmp[37]*-4.0;
  dxdotdp[(38+ip * udata->nx)] = x_tmp[1]*x_tmp[37]*2.0;
  dxdotdp[(46+ip * udata->nx)] = x_tmp[1]*x_tmp[46]*-4.0;
  dxdotdp[(47+ip * udata->nx)] = x_tmp[1]*x_tmp[46]*2.0;
  dxdotdp[(55+ip * udata->nx)] = x_tmp[1]*x_tmp[55]*-4.0;
  dxdotdp[(56+ip * udata->nx)] = x_tmp[1]*x_tmp[55]*2.0;
  dxdotdp[(64+ip * udata->nx)] = x_tmp[1]*x_tmp[64]*-4.0;
  dxdotdp[(65+ip * udata->nx)] = x_tmp[1]*x_tmp[64]*2.0;
  dxdotdp[(73+ip * udata->nx)] = x_tmp[1]*x_tmp[73]*-4.0;
  dxdotdp[(74+ip * udata->nx)] = x_tmp[1]*x_tmp[73]*2.0;
  dxdotdp[(82+ip * udata->nx)] = x_tmp[1]*x_tmp[82]*-4.0;
  dxdotdp[(83+ip * udata->nx)] = x_tmp[1]*x_tmp[82]*2.0;
  dxdotdp[(91+ip * udata->nx)] = x_tmp[1]*x_tmp[91]*-4.0;
  dxdotdp[(92+ip * udata->nx)] = x_tmp[1]*x_tmp[91]*2.0;
  dxdotdp[(100+ip * udata->nx)] = x_tmp[1]*x_tmp[100]*-4.0;
  dxdotdp[(101+ip * udata->nx)] = x_tmp[1]*x_tmp[100]*2.0;
  dxdotdp[(109+ip * udata->nx)] = x_tmp[1]*x_tmp[109]*-4.0;
  dxdotdp[(110+ip * udata->nx)] = x_tmp[1]*x_tmp[109]*2.0;
  dxdotdp[(118+ip * udata->nx)] = x_tmp[1]*x_tmp[118]*-4.0;
  dxdotdp[(119+ip * udata->nx)] = x_tmp[1]*x_tmp[118]*2.0;
  dxdotdp[(127+ip * udata->nx)] = x_tmp[1]*x_tmp[127]*-4.0;
  dxdotdp[(128+ip * udata->nx)] = x_tmp[1]*x_tmp[127]*2.0;
  dxdotdp[(136+ip * udata->nx)] = x_tmp[1]*x_tmp[136]*-4.0;
  dxdotdp[(137+ip * udata->nx)] = x_tmp[1]*x_tmp[136]*2.0;
  dxdotdp[(145+ip * udata->nx)] = x_tmp[1]*x_tmp[145]*-4.0;
  dxdotdp[(146+ip * udata->nx)] = x_tmp[1]*x_tmp[145]*2.0;
  dxdotdp[(154+ip * udata->nx)] = x_tmp[1]*x_tmp[154]*-4.0;
  dxdotdp[(155+ip * udata->nx)] = x_tmp[1]*x_tmp[154]*2.0;

  } break;

  case 2: {
  dxdotdp[(2+ip * udata->nx)] = -x_tmp[2];
  dxdotdp[(3+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[2];
  dxdotdp[(11+ip * udata->nx)] = -x_tmp[11];
  dxdotdp[(12+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[11];
  dxdotdp[(20+ip * udata->nx)] = -x_tmp[20];
  dxdotdp[(21+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[20];
  dxdotdp[(29+ip * udata->nx)] = -x_tmp[29];
  dxdotdp[(30+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[29];
  dxdotdp[(38+ip * udata->nx)] = -x_tmp[38];
  dxdotdp[(39+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[38];
  dxdotdp[(47+ip * udata->nx)] = -x_tmp[47];
  dxdotdp[(48+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[47];
  dxdotdp[(56+ip * udata->nx)] = -x_tmp[56];
  dxdotdp[(57+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[56];
  dxdotdp[(65+ip * udata->nx)] = -x_tmp[65];
  dxdotdp[(66+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[65];
  dxdotdp[(74+ip * udata->nx)] = -x_tmp[74];
  dxdotdp[(75+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[74];
  dxdotdp[(83+ip * udata->nx)] = -x_tmp[83];
  dxdotdp[(84+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[83];
  dxdotdp[(92+ip * udata->nx)] = -x_tmp[92];
  dxdotdp[(93+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[92];
  dxdotdp[(101+ip * udata->nx)] = -x_tmp[101];
  dxdotdp[(102+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[101];
  dxdotdp[(110+ip * udata->nx)] = -x_tmp[110];
  dxdotdp[(111+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[110];
  dxdotdp[(119+ip * udata->nx)] = -x_tmp[119];
  dxdotdp[(120+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[119];
  dxdotdp[(128+ip * udata->nx)] = -x_tmp[128];
  dxdotdp[(129+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[128];
  dxdotdp[(137+ip * udata->nx)] = -x_tmp[137];
  dxdotdp[(138+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[137];
  dxdotdp[(146+ip * udata->nx)] = -x_tmp[146];
  dxdotdp[(147+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[146];
  dxdotdp[(155+ip * udata->nx)] = -x_tmp[155];
  dxdotdp[(156+ip * udata->nx)] = udata->k[0]*udata->w[3]*x_tmp[155];

  } break;

  case 3: {
  dxdotdp[(0+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[8];
  dxdotdp[(3+ip * udata->nx)] = -udata->k[1]*udata->w[3]*x_tmp[3];
  dxdotdp[(4+ip * udata->nx)] = udata->w[4]-x_tmp[4];
  dxdotdp[(5+ip * udata->nx)] = x_tmp[4]-x_tmp[5];
  dxdotdp[(6+ip * udata->nx)] = x_tmp[5]-x_tmp[6];
  dxdotdp[(7+ip * udata->nx)] = x_tmp[6]-x_tmp[7];
  dxdotdp[(8+ip * udata->nx)] = x_tmp[7]-x_tmp[8];
  dxdotdp[(9+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[17];
  dxdotdp[(12+ip * udata->nx)] = -x_tmp[12];
  dxdotdp[(13+ip * udata->nx)] = x_tmp[12]*2.0-x_tmp[13];
  dxdotdp[(14+ip * udata->nx)] = x_tmp[13]-x_tmp[14];
  dxdotdp[(15+ip * udata->nx)] = x_tmp[14]-x_tmp[15];
  dxdotdp[(16+ip * udata->nx)] = x_tmp[15]-x_tmp[16];
  dxdotdp[(17+ip * udata->nx)] = x_tmp[16]-x_tmp[17];
  dxdotdp[(18+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[26];
  dxdotdp[(21+ip * udata->nx)] = -x_tmp[21];
  dxdotdp[(22+ip * udata->nx)] = x_tmp[21]*2.0-x_tmp[22];
  dxdotdp[(23+ip * udata->nx)] = x_tmp[22]-x_tmp[23];
  dxdotdp[(24+ip * udata->nx)] = x_tmp[23]-x_tmp[24];
  dxdotdp[(25+ip * udata->nx)] = x_tmp[24]-x_tmp[25];
  dxdotdp[(26+ip * udata->nx)] = x_tmp[25]-x_tmp[26];
  dxdotdp[(27+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[35];
  dxdotdp[(30+ip * udata->nx)] = -x_tmp[30];
  dxdotdp[(31+ip * udata->nx)] = x_tmp[30]*2.0-x_tmp[31];
  dxdotdp[(32+ip * udata->nx)] = x_tmp[31]-x_tmp[32];
  dxdotdp[(33+ip * udata->nx)] = x_tmp[32]-x_tmp[33];
  dxdotdp[(34+ip * udata->nx)] = x_tmp[33]-x_tmp[34];
  dxdotdp[(35+ip * udata->nx)] = x_tmp[34]-x_tmp[35];
  dxdotdp[(36+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[44];
  dxdotdp[(39+ip * udata->nx)] = -x_tmp[39];
  dxdotdp[(40+ip * udata->nx)] = x_tmp[39]*2.0-x_tmp[40];
  dxdotdp[(41+ip * udata->nx)] = x_tmp[40]-x_tmp[41];
  dxdotdp[(42+ip * udata->nx)] = x_tmp[41]-x_tmp[42];
  dxdotdp[(43+ip * udata->nx)] = x_tmp[42]-x_tmp[43];
  dxdotdp[(44+ip * udata->nx)] = x_tmp[43]-x_tmp[44];
  dxdotdp[(45+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[53];
  dxdotdp[(48+ip * udata->nx)] = -x_tmp[48];
  dxdotdp[(49+ip * udata->nx)] = x_tmp[48]*2.0-x_tmp[49];
  dxdotdp[(50+ip * udata->nx)] = x_tmp[49]-x_tmp[50];
  dxdotdp[(51+ip * udata->nx)] = x_tmp[50]-x_tmp[51];
  dxdotdp[(52+ip * udata->nx)] = x_tmp[51]-x_tmp[52];
  dxdotdp[(53+ip * udata->nx)] = x_tmp[52]-x_tmp[53];
  dxdotdp[(54+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[62];
  dxdotdp[(57+ip * udata->nx)] = -x_tmp[57];
  dxdotdp[(58+ip * udata->nx)] = x_tmp[57]*2.0-x_tmp[58];
  dxdotdp[(59+ip * udata->nx)] = x_tmp[58]-x_tmp[59];
  dxdotdp[(60+ip * udata->nx)] = x_tmp[59]-x_tmp[60];
  dxdotdp[(61+ip * udata->nx)] = x_tmp[60]-x_tmp[61];
  dxdotdp[(62+ip * udata->nx)] = x_tmp[61]-x_tmp[62];
  dxdotdp[(63+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[71];
  dxdotdp[(66+ip * udata->nx)] = -x_tmp[66];
  dxdotdp[(67+ip * udata->nx)] = x_tmp[66]*2.0-x_tmp[67];
  dxdotdp[(68+ip * udata->nx)] = x_tmp[67]-x_tmp[68];
  dxdotdp[(69+ip * udata->nx)] = x_tmp[68]-x_tmp[69];
  dxdotdp[(70+ip * udata->nx)] = x_tmp[69]-x_tmp[70];
  dxdotdp[(71+ip * udata->nx)] = x_tmp[70]-x_tmp[71];
  dxdotdp[(72+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[80];
  dxdotdp[(75+ip * udata->nx)] = -x_tmp[75];
  dxdotdp[(76+ip * udata->nx)] = x_tmp[75]*2.0-x_tmp[76];
  dxdotdp[(77+ip * udata->nx)] = x_tmp[76]-x_tmp[77];
  dxdotdp[(78+ip * udata->nx)] = x_tmp[77]-x_tmp[78];
  dxdotdp[(79+ip * udata->nx)] = x_tmp[78]-x_tmp[79];
  dxdotdp[(80+ip * udata->nx)] = x_tmp[79]-x_tmp[80];
  dxdotdp[(81+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[89];
  dxdotdp[(84+ip * udata->nx)] = -x_tmp[84];
  dxdotdp[(85+ip * udata->nx)] = x_tmp[84]*2.0-x_tmp[85];
  dxdotdp[(86+ip * udata->nx)] = x_tmp[85]-x_tmp[86];
  dxdotdp[(87+ip * udata->nx)] = x_tmp[86]-x_tmp[87];
  dxdotdp[(88+ip * udata->nx)] = x_tmp[87]-x_tmp[88];
  dxdotdp[(89+ip * udata->nx)] = x_tmp[88]-x_tmp[89];
  dxdotdp[(90+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[98];
  dxdotdp[(93+ip * udata->nx)] = -x_tmp[93];
  dxdotdp[(94+ip * udata->nx)] = x_tmp[93]*2.0-x_tmp[94];
  dxdotdp[(95+ip * udata->nx)] = x_tmp[94]-x_tmp[95];
  dxdotdp[(96+ip * udata->nx)] = x_tmp[95]-x_tmp[96];
  dxdotdp[(97+ip * udata->nx)] = x_tmp[96]-x_tmp[97];
  dxdotdp[(98+ip * udata->nx)] = x_tmp[97]-x_tmp[98];
  dxdotdp[(99+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[107];
  dxdotdp[(102+ip * udata->nx)] = -x_tmp[102];
  dxdotdp[(103+ip * udata->nx)] = x_tmp[102]*2.0-x_tmp[103];
  dxdotdp[(104+ip * udata->nx)] = x_tmp[103]-x_tmp[104];
  dxdotdp[(105+ip * udata->nx)] = x_tmp[104]-x_tmp[105];
  dxdotdp[(106+ip * udata->nx)] = x_tmp[105]-x_tmp[106];
  dxdotdp[(107+ip * udata->nx)] = x_tmp[106]-x_tmp[107];
  dxdotdp[(108+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[116];
  dxdotdp[(111+ip * udata->nx)] = -x_tmp[111];
  dxdotdp[(112+ip * udata->nx)] = x_tmp[111]*2.0-x_tmp[112];
  dxdotdp[(113+ip * udata->nx)] = x_tmp[112]-x_tmp[113];
  dxdotdp[(114+ip * udata->nx)] = x_tmp[113]-x_tmp[114];
  dxdotdp[(115+ip * udata->nx)] = x_tmp[114]-x_tmp[115];
  dxdotdp[(116+ip * udata->nx)] = x_tmp[115]-x_tmp[116];
  dxdotdp[(117+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[125];
  dxdotdp[(120+ip * udata->nx)] = -x_tmp[120];
  dxdotdp[(121+ip * udata->nx)] = x_tmp[120]*2.0-x_tmp[121];
  dxdotdp[(122+ip * udata->nx)] = x_tmp[121]-x_tmp[122];
  dxdotdp[(123+ip * udata->nx)] = x_tmp[122]-x_tmp[123];
  dxdotdp[(124+ip * udata->nx)] = x_tmp[123]-x_tmp[124];
  dxdotdp[(125+ip * udata->nx)] = x_tmp[124]-x_tmp[125];
  dxdotdp[(126+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[134];
  dxdotdp[(129+ip * udata->nx)] = -x_tmp[129];
  dxdotdp[(130+ip * udata->nx)] = x_tmp[129]*2.0-x_tmp[130];
  dxdotdp[(131+ip * udata->nx)] = x_tmp[130]-x_tmp[131];
  dxdotdp[(132+ip * udata->nx)] = x_tmp[131]-x_tmp[132];
  dxdotdp[(133+ip * udata->nx)] = x_tmp[132]-x_tmp[133];
  dxdotdp[(134+ip * udata->nx)] = x_tmp[133]-x_tmp[134];
  dxdotdp[(135+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[143];
  dxdotdp[(138+ip * udata->nx)] = -x_tmp[138];
  dxdotdp[(139+ip * udata->nx)] = x_tmp[138]*2.0-x_tmp[139];
  dxdotdp[(140+ip * udata->nx)] = x_tmp[139]-x_tmp[140];
  dxdotdp[(141+ip * udata->nx)] = x_tmp[140]-x_tmp[141];
  dxdotdp[(142+ip * udata->nx)] = x_tmp[141]-x_tmp[142];
  dxdotdp[(143+ip * udata->nx)] = x_tmp[142]-x_tmp[143];
  dxdotdp[(144+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[152];
  dxdotdp[(147+ip * udata->nx)] = -x_tmp[147];
  dxdotdp[(148+ip * udata->nx)] = x_tmp[147]*2.0-x_tmp[148];
  dxdotdp[(149+ip * udata->nx)] = x_tmp[148]-x_tmp[149];
  dxdotdp[(150+ip * udata->nx)] = x_tmp[149]-x_tmp[150];
  dxdotdp[(151+ip * udata->nx)] = x_tmp[150]-x_tmp[151];
  dxdotdp[(152+ip * udata->nx)] = x_tmp[151]-x_tmp[152];
  dxdotdp[(153+ip * udata->nx)] = udata->k[1]*udata->w[2]*x_tmp[161];
  dxdotdp[(156+ip * udata->nx)] = -x_tmp[156];
  dxdotdp[(157+ip * udata->nx)] = x_tmp[156]*2.0-x_tmp[157];
  dxdotdp[(158+ip * udata->nx)] = x_tmp[157]-x_tmp[158];
  dxdotdp[(159+ip * udata->nx)] = x_tmp[158]-x_tmp[159];
  dxdotdp[(160+ip * udata->nx)] = x_tmp[159]-x_tmp[160];
  dxdotdp[(161+ip * udata->nx)] = x_tmp[160]-x_tmp[161];

  } break;

  case 5: {
  dxdotdp[(0+ip * udata->nx)] = -udata->dwdp[0]*udata->k[0]*udata->p[0]*udata->w[2]*x_tmp[0];
  dxdotdp[(1+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[0];
  dxdotdp[(9+ip * udata->nx)] = -udata->dwdp[0]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(10+ip * udata->nx)] = udata->dwdp[0]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(18+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[18];
  dxdotdp[(19+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[18];
  dxdotdp[(27+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[27];
  dxdotdp[(28+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[27];
  dxdotdp[(36+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[36];
  dxdotdp[(37+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[36];
  dxdotdp[(45+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[45];
  dxdotdp[(46+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[45];
  dxdotdp[(54+ip * udata->nx)] = -udata->dwdp[1]*udata->p[0]*x_tmp[0]-udata->dwdp[0]*udata->p[0]*x_tmp[54];
  dxdotdp[(55+ip * udata->nx)] = udata->dwdp[1]*udata->p[0]*x_tmp[0]+udata->dwdp[0]*udata->p[0]*x_tmp[54];
  dxdotdp[(63+ip * udata->nx)] = -udata->dwdp[2]*udata->p[0]*x_tmp[0]-udata->dwdp[0]*udata->p[0]*x_tmp[63];
  dxdotdp[(64+ip * udata->nx)] = udata->dwdp[2]*udata->p[0]*x_tmp[0]+udata->dwdp[0]*udata->p[0]*x_tmp[63];
  dxdotdp[(72+ip * udata->nx)] = -udata->dwdp[3]*udata->p[0]*x_tmp[0]-udata->dwdp[0]*udata->p[0]*x_tmp[72];
  dxdotdp[(73+ip * udata->nx)] = udata->dwdp[3]*udata->p[0]*x_tmp[0]+udata->dwdp[0]*udata->p[0]*x_tmp[72];
  dxdotdp[(81+ip * udata->nx)] = -udata->dwdp[4]*udata->p[0]*x_tmp[0]-udata->dwdp[0]*udata->p[0]*x_tmp[81];
  dxdotdp[(82+ip * udata->nx)] = udata->dwdp[4]*udata->p[0]*x_tmp[0]+udata->dwdp[0]*udata->p[0]*x_tmp[81];
  dxdotdp[(90+ip * udata->nx)] = -udata->dwdp[5]*udata->p[0]*x_tmp[0]-udata->dwdp[0]*udata->p[0]*x_tmp[90];
  dxdotdp[(91+ip * udata->nx)] = udata->dwdp[5]*udata->p[0]*x_tmp[0]+udata->dwdp[0]*udata->p[0]*x_tmp[90];
  dxdotdp[(99+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[99];
  dxdotdp[(100+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[99];
  dxdotdp[(108+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[108];
  dxdotdp[(109+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[108];
  dxdotdp[(117+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[117];
  dxdotdp[(118+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[117];
  dxdotdp[(126+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[126];
  dxdotdp[(127+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[126];
  dxdotdp[(135+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[135];
  dxdotdp[(136+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[135];
  dxdotdp[(144+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[144];
  dxdotdp[(145+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[144];
  dxdotdp[(153+ip * udata->nx)] = -udata->dwdp[0]*udata->p[0]*x_tmp[153];
  dxdotdp[(154+ip * udata->nx)] = udata->dwdp[0]*udata->p[0]*x_tmp[153];

  } break;

  case 6: {
  dxdotdp[(0+ip * udata->nx)] = -udata->dwdp[6]*udata->k[0]*udata->p[0]*udata->w[2]*x_tmp[0];
  dxdotdp[(1+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[0];
  dxdotdp[(9+ip * udata->nx)] = -udata->dwdp[6]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(10+ip * udata->nx)] = udata->dwdp[6]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(18+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[18];
  dxdotdp[(19+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[18];
  dxdotdp[(27+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[27];
  dxdotdp[(28+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[27];
  dxdotdp[(36+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[36];
  dxdotdp[(37+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[36];
  dxdotdp[(45+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[45];
  dxdotdp[(46+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[45];
  dxdotdp[(54+ip * udata->nx)] = -udata->dwdp[7]*udata->p[0]*x_tmp[0]-udata->dwdp[6]*udata->p[0]*x_tmp[54];
  dxdotdp[(55+ip * udata->nx)] = udata->dwdp[7]*udata->p[0]*x_tmp[0]+udata->dwdp[6]*udata->p[0]*x_tmp[54];
  dxdotdp[(63+ip * udata->nx)] = -udata->dwdp[8]*udata->p[0]*x_tmp[0]-udata->dwdp[6]*udata->p[0]*x_tmp[63];
  dxdotdp[(64+ip * udata->nx)] = udata->dwdp[8]*udata->p[0]*x_tmp[0]+udata->dwdp[6]*udata->p[0]*x_tmp[63];
  dxdotdp[(72+ip * udata->nx)] = -udata->dwdp[9]*udata->p[0]*x_tmp[0]-udata->dwdp[6]*udata->p[0]*x_tmp[72];
  dxdotdp[(73+ip * udata->nx)] = udata->dwdp[9]*udata->p[0]*x_tmp[0]+udata->dwdp[6]*udata->p[0]*x_tmp[72];
  dxdotdp[(81+ip * udata->nx)] = -udata->dwdp[10]*udata->p[0]*x_tmp[0]-udata->dwdp[6]*udata->p[0]*x_tmp[81];
  dxdotdp[(82+ip * udata->nx)] = udata->dwdp[10]*udata->p[0]*x_tmp[0]+udata->dwdp[6]*udata->p[0]*x_tmp[81];
  dxdotdp[(90+ip * udata->nx)] = -udata->dwdp[11]*udata->p[0]*x_tmp[0]-udata->dwdp[6]*udata->p[0]*x_tmp[90];
  dxdotdp[(91+ip * udata->nx)] = udata->dwdp[11]*udata->p[0]*x_tmp[0]+udata->dwdp[6]*udata->p[0]*x_tmp[90];
  dxdotdp[(99+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[99];
  dxdotdp[(100+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[99];
  dxdotdp[(108+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[108];
  dxdotdp[(109+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[108];
  dxdotdp[(117+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[117];
  dxdotdp[(118+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[117];
  dxdotdp[(126+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[126];
  dxdotdp[(127+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[126];
  dxdotdp[(135+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[135];
  dxdotdp[(136+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[135];
  dxdotdp[(144+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[144];
  dxdotdp[(145+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[144];
  dxdotdp[(153+ip * udata->nx)] = -udata->dwdp[6]*udata->p[0]*x_tmp[153];
  dxdotdp[(154+ip * udata->nx)] = udata->dwdp[6]*udata->p[0]*x_tmp[153];

  } break;

  case 7: {
  dxdotdp[(0+ip * udata->nx)] = -udata->dwdp[12]*udata->k[0]*udata->p[0]*udata->w[2]*x_tmp[0];
  dxdotdp[(1+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[0];
  dxdotdp[(9+ip * udata->nx)] = -udata->dwdp[12]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(10+ip * udata->nx)] = udata->dwdp[12]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(18+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[18];
  dxdotdp[(19+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[18];
  dxdotdp[(27+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[27];
  dxdotdp[(28+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[27];
  dxdotdp[(36+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[36];
  dxdotdp[(37+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[36];
  dxdotdp[(45+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[45];
  dxdotdp[(46+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[45];
  dxdotdp[(54+ip * udata->nx)] = -udata->dwdp[13]*udata->p[0]*x_tmp[0]-udata->dwdp[12]*udata->p[0]*x_tmp[54];
  dxdotdp[(55+ip * udata->nx)] = udata->dwdp[13]*udata->p[0]*x_tmp[0]+udata->dwdp[12]*udata->p[0]*x_tmp[54];
  dxdotdp[(63+ip * udata->nx)] = -udata->dwdp[14]*udata->p[0]*x_tmp[0]-udata->dwdp[12]*udata->p[0]*x_tmp[63];
  dxdotdp[(64+ip * udata->nx)] = udata->dwdp[14]*udata->p[0]*x_tmp[0]+udata->dwdp[12]*udata->p[0]*x_tmp[63];
  dxdotdp[(72+ip * udata->nx)] = -udata->dwdp[15]*udata->p[0]*x_tmp[0]-udata->dwdp[12]*udata->p[0]*x_tmp[72];
  dxdotdp[(73+ip * udata->nx)] = udata->dwdp[15]*udata->p[0]*x_tmp[0]+udata->dwdp[12]*udata->p[0]*x_tmp[72];
  dxdotdp[(81+ip * udata->nx)] = -udata->dwdp[16]*udata->p[0]*x_tmp[0]-udata->dwdp[12]*udata->p[0]*x_tmp[81];
  dxdotdp[(82+ip * udata->nx)] = udata->dwdp[16]*udata->p[0]*x_tmp[0]+udata->dwdp[12]*udata->p[0]*x_tmp[81];
  dxdotdp[(90+ip * udata->nx)] = -udata->dwdp[17]*udata->p[0]*x_tmp[0]-udata->dwdp[12]*udata->p[0]*x_tmp[90];
  dxdotdp[(91+ip * udata->nx)] = udata->dwdp[17]*udata->p[0]*x_tmp[0]+udata->dwdp[12]*udata->p[0]*x_tmp[90];
  dxdotdp[(99+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[99];
  dxdotdp[(100+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[99];
  dxdotdp[(108+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[108];
  dxdotdp[(109+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[108];
  dxdotdp[(117+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[117];
  dxdotdp[(118+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[117];
  dxdotdp[(126+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[126];
  dxdotdp[(127+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[126];
  dxdotdp[(135+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[135];
  dxdotdp[(136+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[135];
  dxdotdp[(144+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[144];
  dxdotdp[(145+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[144];
  dxdotdp[(153+ip * udata->nx)] = -udata->dwdp[12]*udata->p[0]*x_tmp[153];
  dxdotdp[(154+ip * udata->nx)] = udata->dwdp[12]*udata->p[0]*x_tmp[153];

  } break;

  case 8: {
  dxdotdp[(0+ip * udata->nx)] = -udata->dwdp[18]*udata->k[0]*udata->p[0]*udata->w[2]*x_tmp[0];
  dxdotdp[(1+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[0];
  dxdotdp[(9+ip * udata->nx)] = -udata->dwdp[18]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(10+ip * udata->nx)] = udata->dwdp[18]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(18+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[18];
  dxdotdp[(19+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[18];
  dxdotdp[(27+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[27];
  dxdotdp[(28+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[27];
  dxdotdp[(36+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[36];
  dxdotdp[(37+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[36];
  dxdotdp[(45+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[45];
  dxdotdp[(46+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[45];
  dxdotdp[(54+ip * udata->nx)] = -udata->dwdp[19]*udata->p[0]*x_tmp[0]-udata->dwdp[18]*udata->p[0]*x_tmp[54];
  dxdotdp[(55+ip * udata->nx)] = udata->dwdp[19]*udata->p[0]*x_tmp[0]+udata->dwdp[18]*udata->p[0]*x_tmp[54];
  dxdotdp[(63+ip * udata->nx)] = -udata->dwdp[20]*udata->p[0]*x_tmp[0]-udata->dwdp[18]*udata->p[0]*x_tmp[63];
  dxdotdp[(64+ip * udata->nx)] = udata->dwdp[20]*udata->p[0]*x_tmp[0]+udata->dwdp[18]*udata->p[0]*x_tmp[63];
  dxdotdp[(72+ip * udata->nx)] = -udata->dwdp[21]*udata->p[0]*x_tmp[0]-udata->dwdp[18]*udata->p[0]*x_tmp[72];
  dxdotdp[(73+ip * udata->nx)] = udata->dwdp[21]*udata->p[0]*x_tmp[0]+udata->dwdp[18]*udata->p[0]*x_tmp[72];
  dxdotdp[(81+ip * udata->nx)] = -udata->dwdp[22]*udata->p[0]*x_tmp[0]-udata->dwdp[18]*udata->p[0]*x_tmp[81];
  dxdotdp[(82+ip * udata->nx)] = udata->dwdp[22]*udata->p[0]*x_tmp[0]+udata->dwdp[18]*udata->p[0]*x_tmp[81];
  dxdotdp[(90+ip * udata->nx)] = -udata->dwdp[23]*udata->p[0]*x_tmp[0]-udata->dwdp[18]*udata->p[0]*x_tmp[90];
  dxdotdp[(91+ip * udata->nx)] = udata->dwdp[23]*udata->p[0]*x_tmp[0]+udata->dwdp[18]*udata->p[0]*x_tmp[90];
  dxdotdp[(99+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[99];
  dxdotdp[(100+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[99];
  dxdotdp[(108+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[108];
  dxdotdp[(109+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[108];
  dxdotdp[(117+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[117];
  dxdotdp[(118+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[117];
  dxdotdp[(126+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[126];
  dxdotdp[(127+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[126];
  dxdotdp[(135+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[135];
  dxdotdp[(136+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[135];
  dxdotdp[(144+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[144];
  dxdotdp[(145+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[144];
  dxdotdp[(153+ip * udata->nx)] = -udata->dwdp[18]*udata->p[0]*x_tmp[153];
  dxdotdp[(154+ip * udata->nx)] = udata->dwdp[18]*udata->p[0]*x_tmp[153];

  } break;

  case 9: {
  dxdotdp[(0+ip * udata->nx)] = -udata->dwdp[24]*udata->k[0]*udata->p[0]*udata->w[2]*x_tmp[0];
  dxdotdp[(1+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[0];
  dxdotdp[(9+ip * udata->nx)] = -udata->dwdp[24]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(10+ip * udata->nx)] = udata->dwdp[24]*(x_tmp[0]+udata->p[0]*x_tmp[9]);
  dxdotdp[(18+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[18];
  dxdotdp[(19+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[18];
  dxdotdp[(27+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[27];
  dxdotdp[(28+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[27];
  dxdotdp[(36+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[36];
  dxdotdp[(37+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[36];
  dxdotdp[(45+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[45];
  dxdotdp[(46+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[45];
  dxdotdp[(54+ip * udata->nx)] = -udata->dwdp[25]*udata->p[0]*x_tmp[0]-udata->dwdp[24]*udata->p[0]*x_tmp[54];
  dxdotdp[(55+ip * udata->nx)] = udata->dwdp[25]*udata->p[0]*x_tmp[0]+udata->dwdp[24]*udata->p[0]*x_tmp[54];
  dxdotdp[(63+ip * udata->nx)] = -udata->dwdp[26]*udata->p[0]*x_tmp[0]-udata->dwdp[24]*udata->p[0]*x_tmp[63];
  dxdotdp[(64+ip * udata->nx)] = udata->dwdp[26]*udata->p[0]*x_tmp[0]+udata->dwdp[24]*udata->p[0]*x_tmp[63];
  dxdotdp[(72+ip * udata->nx)] = -udata->dwdp[27]*udata->p[0]*x_tmp[0]-udata->dwdp[24]*udata->p[0]*x_tmp[72];
  dxdotdp[(73+ip * udata->nx)] = udata->dwdp[27]*udata->p[0]*x_tmp[0]+udata->dwdp[24]*udata->p[0]*x_tmp[72];
  dxdotdp[(81+ip * udata->nx)] = -udata->dwdp[28]*udata->p[0]*x_tmp[0]-udata->dwdp[24]*udata->p[0]*x_tmp[81];
  dxdotdp[(82+ip * udata->nx)] = udata->dwdp[28]*udata->p[0]*x_tmp[0]+udata->dwdp[24]*udata->p[0]*x_tmp[81];
  dxdotdp[(90+ip * udata->nx)] = -udata->dwdp[29]*udata->p[0]*x_tmp[0]-udata->dwdp[24]*udata->p[0]*x_tmp[90];
  dxdotdp[(91+ip * udata->nx)] = udata->dwdp[29]*udata->p[0]*x_tmp[0]+udata->dwdp[24]*udata->p[0]*x_tmp[90];
  dxdotdp[(99+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[99];
  dxdotdp[(100+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[99];
  dxdotdp[(108+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[108];
  dxdotdp[(109+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[108];
  dxdotdp[(117+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[117];
  dxdotdp[(118+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[117];
  dxdotdp[(126+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[126];
  dxdotdp[(127+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[126];
  dxdotdp[(135+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[135];
  dxdotdp[(136+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[135];
  dxdotdp[(144+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[144];
  dxdotdp[(145+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[144];
  dxdotdp[(153+ip * udata->nx)] = -udata->dwdp[24]*udata->p[0]*x_tmp[153];
  dxdotdp[(154+ip * udata->nx)] = udata->dwdp[24]*udata->p[0]*x_tmp[153];

  } break;

}
}
for(ip = 0; ip<udata->nplist; ip++) {
   for(ix = 0; ix<162; ix++) {
       if(amiIsNaN(dxdotdp[ix+ip*162])) {
           dxdotdp[ix+ip*162] = 0;
           if(!udata->nan_dxdotdp) {
               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");
               udata->nan_dxdotdp = TRUE;
           }
       }
       if(amiIsInf(dxdotdp[ix+ip*162])) {
           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");
           return(-1);
       }
   }
}
return(status);

}


