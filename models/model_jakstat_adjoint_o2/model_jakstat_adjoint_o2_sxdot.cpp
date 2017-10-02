
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_JSparse.h"
#include "model_jakstat_adjoint_o2_dxdotdp.h"
#include "model_jakstat_adjoint_o2_w.h"

int sxdot_model_jakstat_adjoint_o2(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
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
realtype *sx_tmp = nullptr;
if(sx)
    sx_tmp = N_VGetArrayPointer(sx);
realtype *sdx_tmp = nullptr;
if(sdx)
    sdx_tmp = N_VGetArrayPointer(sdx);
realtype *sxdot_tmp = nullptr;
if(sxdot)
    sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*162);
if(ip == 0) {
    status = JSparse_model_jakstat_adjoint_o2(t,0.0,x,NULL,xdot,tdata->J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_jakstat_adjoint_o2(t,x,NULL,user_data);
}
  sxdot_tmp[0] = tdata->dxdotdp[0 + ip*model->nx]+tdata->J->data[0]*sx_tmp[0]+tdata->J->data[74]*sx_tmp[8];
  sxdot_tmp[1] = tdata->dxdotdp[1 + ip*model->nx]+tdata->J->data[1]*sx_tmp[0]+tdata->J->data[14]*sx_tmp[1];
  sxdot_tmp[2] = tdata->dxdotdp[2 + ip*model->nx]+tdata->J->data[15]*sx_tmp[1]+tdata->J->data[50]*sx_tmp[2];
  sxdot_tmp[3] = tdata->dxdotdp[3 + ip*model->nx]+tdata->J->data[51]*sx_tmp[2]+tdata->J->data[54]*sx_tmp[3];
  sxdot_tmp[4] = tdata->dxdotdp[4 + ip*model->nx]+tdata->J->data[55]*sx_tmp[3]+tdata->J->data[58]*sx_tmp[4];
  sxdot_tmp[5] = tdata->dxdotdp[5 + ip*model->nx]+tdata->J->data[59]*sx_tmp[4]+tdata->J->data[62]*sx_tmp[5];
  sxdot_tmp[6] = tdata->dxdotdp[6 + ip*model->nx]+tdata->J->data[63]*sx_tmp[5]+tdata->J->data[66]*sx_tmp[6];
  sxdot_tmp[7] = tdata->dxdotdp[7 + ip*model->nx]+tdata->J->data[67]*sx_tmp[6]+tdata->J->data[70]*sx_tmp[7];
  sxdot_tmp[8] = tdata->dxdotdp[8 + ip*model->nx]+tdata->J->data[71]*sx_tmp[7]+tdata->J->data[75]*sx_tmp[8];
  sxdot_tmp[9] = tdata->dxdotdp[9 + ip*model->nx]+tdata->J->data[2]*sx_tmp[0]+tdata->J->data[78]*sx_tmp[9]+tdata->J->data[94]*sx_tmp[17];
  sxdot_tmp[10] = tdata->dxdotdp[10 + ip*model->nx]+tdata->J->data[3]*sx_tmp[0]+tdata->J->data[16]*sx_tmp[1]+tdata->J->data[79]*sx_tmp[9]+tdata->J->data[80]*sx_tmp[10];
  sxdot_tmp[11] = tdata->dxdotdp[11 + ip*model->nx]+tdata->J->data[17]*sx_tmp[1]+tdata->J->data[81]*sx_tmp[10]+tdata->J->data[82]*sx_tmp[11];
  sxdot_tmp[12] = tdata->dxdotdp[12 + ip*model->nx]+tdata->J->data[83]*sx_tmp[11]+tdata->J->data[84]*sx_tmp[12];
  sxdot_tmp[13] = tdata->dxdotdp[13 + ip*model->nx]+tdata->J->data[85]*sx_tmp[12]+tdata->J->data[86]*sx_tmp[13];
  sxdot_tmp[14] = tdata->dxdotdp[14 + ip*model->nx]+tdata->J->data[87]*sx_tmp[13]+tdata->J->data[88]*sx_tmp[14];
  sxdot_tmp[15] = tdata->dxdotdp[15 + ip*model->nx]+tdata->J->data[89]*sx_tmp[14]+tdata->J->data[90]*sx_tmp[15];
  sxdot_tmp[16] = tdata->dxdotdp[16 + ip*model->nx]+tdata->J->data[91]*sx_tmp[15]+tdata->J->data[92]*sx_tmp[16];
  sxdot_tmp[17] = tdata->dxdotdp[17 + ip*model->nx]+tdata->J->data[93]*sx_tmp[16]+tdata->J->data[95]*sx_tmp[17];
  sxdot_tmp[18] = tdata->dxdotdp[18 + ip*model->nx]+tdata->J->data[96]*sx_tmp[18]+tdata->J->data[112]*sx_tmp[26];
  sxdot_tmp[19] = tdata->dxdotdp[19 + ip*model->nx]+tdata->J->data[18]*sx_tmp[1]+tdata->J->data[97]*sx_tmp[18]+tdata->J->data[98]*sx_tmp[19];
  sxdot_tmp[20] = tdata->dxdotdp[20 + ip*model->nx]+tdata->J->data[19]*sx_tmp[1]+tdata->J->data[99]*sx_tmp[19]+tdata->J->data[100]*sx_tmp[20];
  sxdot_tmp[21] = tdata->dxdotdp[21 + ip*model->nx]+tdata->J->data[101]*sx_tmp[20]+tdata->J->data[102]*sx_tmp[21];
  sxdot_tmp[22] = tdata->dxdotdp[22 + ip*model->nx]+tdata->J->data[103]*sx_tmp[21]+tdata->J->data[104]*sx_tmp[22];
  sxdot_tmp[23] = tdata->dxdotdp[23 + ip*model->nx]+tdata->J->data[105]*sx_tmp[22]+tdata->J->data[106]*sx_tmp[23];
  sxdot_tmp[24] = tdata->dxdotdp[24 + ip*model->nx]+tdata->J->data[107]*sx_tmp[23]+tdata->J->data[108]*sx_tmp[24];
  sxdot_tmp[25] = tdata->dxdotdp[25 + ip*model->nx]+tdata->J->data[109]*sx_tmp[24]+tdata->J->data[110]*sx_tmp[25];
  sxdot_tmp[26] = tdata->dxdotdp[26 + ip*model->nx]+tdata->J->data[111]*sx_tmp[25]+tdata->J->data[113]*sx_tmp[26];
  sxdot_tmp[27] = tdata->dxdotdp[27 + ip*model->nx]+tdata->J->data[114]*sx_tmp[27]+tdata->J->data[130]*sx_tmp[35];
  sxdot_tmp[28] = tdata->dxdotdp[28 + ip*model->nx]+tdata->J->data[20]*sx_tmp[1]+tdata->J->data[115]*sx_tmp[27]+tdata->J->data[116]*sx_tmp[28];
  sxdot_tmp[29] = tdata->dxdotdp[29 + ip*model->nx]+tdata->J->data[21]*sx_tmp[1]+tdata->J->data[52]*sx_tmp[2]+tdata->J->data[117]*sx_tmp[28]+tdata->J->data[118]*sx_tmp[29];
  sxdot_tmp[30] = tdata->dxdotdp[30 + ip*model->nx]+tdata->J->data[53]*sx_tmp[2]+tdata->J->data[119]*sx_tmp[29]+tdata->J->data[120]*sx_tmp[30];
  sxdot_tmp[31] = tdata->dxdotdp[31 + ip*model->nx]+tdata->J->data[121]*sx_tmp[30]+tdata->J->data[122]*sx_tmp[31];
  sxdot_tmp[32] = tdata->dxdotdp[32 + ip*model->nx]+tdata->J->data[123]*sx_tmp[31]+tdata->J->data[124]*sx_tmp[32];
  sxdot_tmp[33] = tdata->dxdotdp[33 + ip*model->nx]+tdata->J->data[125]*sx_tmp[32]+tdata->J->data[126]*sx_tmp[33];
  sxdot_tmp[34] = tdata->dxdotdp[34 + ip*model->nx]+tdata->J->data[127]*sx_tmp[33]+tdata->J->data[128]*sx_tmp[34];
  sxdot_tmp[35] = tdata->dxdotdp[35 + ip*model->nx]+tdata->J->data[129]*sx_tmp[34]+tdata->J->data[131]*sx_tmp[35];
  sxdot_tmp[36] = tdata->dxdotdp[36 + ip*model->nx]+tdata->J->data[76]*sx_tmp[8]+tdata->J->data[132]*sx_tmp[36]+tdata->J->data[148]*sx_tmp[44];
  sxdot_tmp[37] = tdata->dxdotdp[37 + ip*model->nx]+tdata->J->data[22]*sx_tmp[1]+tdata->J->data[133]*sx_tmp[36]+tdata->J->data[134]*sx_tmp[37];
  sxdot_tmp[38] = tdata->dxdotdp[38 + ip*model->nx]+tdata->J->data[23]*sx_tmp[1]+tdata->J->data[135]*sx_tmp[37]+tdata->J->data[136]*sx_tmp[38];
  sxdot_tmp[39] = tdata->dxdotdp[39 + ip*model->nx]+tdata->J->data[56]*sx_tmp[3]+tdata->J->data[137]*sx_tmp[38]+tdata->J->data[138]*sx_tmp[39];
  sxdot_tmp[40] = tdata->dxdotdp[40 + ip*model->nx]+tdata->J->data[57]*sx_tmp[3]+tdata->J->data[60]*sx_tmp[4]+tdata->J->data[139]*sx_tmp[39]+tdata->J->data[140]*sx_tmp[40];
  sxdot_tmp[41] = tdata->dxdotdp[41 + ip*model->nx]+tdata->J->data[61]*sx_tmp[4]+tdata->J->data[64]*sx_tmp[5]+tdata->J->data[141]*sx_tmp[40]+tdata->J->data[142]*sx_tmp[41];
  sxdot_tmp[42] = tdata->dxdotdp[42 + ip*model->nx]+tdata->J->data[65]*sx_tmp[5]+tdata->J->data[68]*sx_tmp[6]+tdata->J->data[143]*sx_tmp[41]+tdata->J->data[144]*sx_tmp[42];
  sxdot_tmp[43] = tdata->dxdotdp[43 + ip*model->nx]+tdata->J->data[69]*sx_tmp[6]+tdata->J->data[72]*sx_tmp[7]+tdata->J->data[145]*sx_tmp[42]+tdata->J->data[146]*sx_tmp[43];
  sxdot_tmp[44] = tdata->dxdotdp[44 + ip*model->nx]+tdata->J->data[73]*sx_tmp[7]+tdata->J->data[77]*sx_tmp[8]+tdata->J->data[147]*sx_tmp[43]+tdata->J->data[149]*sx_tmp[44];
  sxdot_tmp[45] = tdata->dxdotdp[45 + ip*model->nx]+tdata->J->data[150]*sx_tmp[45]+tdata->J->data[166]*sx_tmp[53];
  sxdot_tmp[46] = tdata->dxdotdp[46 + ip*model->nx]+tdata->J->data[24]*sx_tmp[1]+tdata->J->data[151]*sx_tmp[45]+tdata->J->data[152]*sx_tmp[46];
  sxdot_tmp[47] = tdata->dxdotdp[47 + ip*model->nx]+tdata->J->data[25]*sx_tmp[1]+tdata->J->data[153]*sx_tmp[46]+tdata->J->data[154]*sx_tmp[47];
  sxdot_tmp[48] = tdata->dxdotdp[48 + ip*model->nx]+tdata->J->data[155]*sx_tmp[47]+tdata->J->data[156]*sx_tmp[48];
  sxdot_tmp[49] = tdata->dxdotdp[49 + ip*model->nx]+tdata->J->data[157]*sx_tmp[48]+tdata->J->data[158]*sx_tmp[49];
  sxdot_tmp[50] = tdata->dxdotdp[50 + ip*model->nx]+tdata->J->data[159]*sx_tmp[49]+tdata->J->data[160]*sx_tmp[50];
  sxdot_tmp[51] = tdata->dxdotdp[51 + ip*model->nx]+tdata->J->data[161]*sx_tmp[50]+tdata->J->data[162]*sx_tmp[51];
  sxdot_tmp[52] = tdata->dxdotdp[52 + ip*model->nx]+tdata->J->data[163]*sx_tmp[51]+tdata->J->data[164]*sx_tmp[52];
  sxdot_tmp[53] = tdata->dxdotdp[53 + ip*model->nx]+tdata->J->data[165]*sx_tmp[52]+tdata->J->data[167]*sx_tmp[53];
  sxdot_tmp[54] = tdata->dxdotdp[54 + ip*model->nx]+tdata->J->data[4]*sx_tmp[0]+tdata->J->data[168]*sx_tmp[54]+tdata->J->data[184]*sx_tmp[62];
  sxdot_tmp[55] = tdata->dxdotdp[55 + ip*model->nx]+tdata->J->data[5]*sx_tmp[0]+tdata->J->data[26]*sx_tmp[1]+tdata->J->data[169]*sx_tmp[54]+tdata->J->data[170]*sx_tmp[55];
  sxdot_tmp[56] = tdata->dxdotdp[56 + ip*model->nx]+tdata->J->data[27]*sx_tmp[1]+tdata->J->data[171]*sx_tmp[55]+tdata->J->data[172]*sx_tmp[56];
  sxdot_tmp[57] = tdata->dxdotdp[57 + ip*model->nx]+tdata->J->data[173]*sx_tmp[56]+tdata->J->data[174]*sx_tmp[57];
  sxdot_tmp[58] = tdata->dxdotdp[58 + ip*model->nx]+tdata->J->data[175]*sx_tmp[57]+tdata->J->data[176]*sx_tmp[58];
  sxdot_tmp[59] = tdata->dxdotdp[59 + ip*model->nx]+tdata->J->data[177]*sx_tmp[58]+tdata->J->data[178]*sx_tmp[59];
  sxdot_tmp[60] = tdata->dxdotdp[60 + ip*model->nx]+tdata->J->data[179]*sx_tmp[59]+tdata->J->data[180]*sx_tmp[60];
  sxdot_tmp[61] = tdata->dxdotdp[61 + ip*model->nx]+tdata->J->data[181]*sx_tmp[60]+tdata->J->data[182]*sx_tmp[61];
  sxdot_tmp[62] = tdata->dxdotdp[62 + ip*model->nx]+tdata->J->data[183]*sx_tmp[61]+tdata->J->data[185]*sx_tmp[62];
  sxdot_tmp[63] = tdata->dxdotdp[63 + ip*model->nx]+tdata->J->data[6]*sx_tmp[0]+tdata->J->data[186]*sx_tmp[63]+tdata->J->data[202]*sx_tmp[71];
  sxdot_tmp[64] = tdata->dxdotdp[64 + ip*model->nx]+tdata->J->data[7]*sx_tmp[0]+tdata->J->data[28]*sx_tmp[1]+tdata->J->data[187]*sx_tmp[63]+tdata->J->data[188]*sx_tmp[64];
  sxdot_tmp[65] = tdata->dxdotdp[65 + ip*model->nx]+tdata->J->data[29]*sx_tmp[1]+tdata->J->data[189]*sx_tmp[64]+tdata->J->data[190]*sx_tmp[65];
  sxdot_tmp[66] = tdata->dxdotdp[66 + ip*model->nx]+tdata->J->data[191]*sx_tmp[65]+tdata->J->data[192]*sx_tmp[66];
  sxdot_tmp[67] = tdata->dxdotdp[67 + ip*model->nx]+tdata->J->data[193]*sx_tmp[66]+tdata->J->data[194]*sx_tmp[67];
  sxdot_tmp[68] = tdata->dxdotdp[68 + ip*model->nx]+tdata->J->data[195]*sx_tmp[67]+tdata->J->data[196]*sx_tmp[68];
  sxdot_tmp[69] = tdata->dxdotdp[69 + ip*model->nx]+tdata->J->data[197]*sx_tmp[68]+tdata->J->data[198]*sx_tmp[69];
  sxdot_tmp[70] = tdata->dxdotdp[70 + ip*model->nx]+tdata->J->data[199]*sx_tmp[69]+tdata->J->data[200]*sx_tmp[70];
  sxdot_tmp[71] = tdata->dxdotdp[71 + ip*model->nx]+tdata->J->data[201]*sx_tmp[70]+tdata->J->data[203]*sx_tmp[71];
  sxdot_tmp[72] = tdata->dxdotdp[72 + ip*model->nx]+tdata->J->data[8]*sx_tmp[0]+tdata->J->data[204]*sx_tmp[72]+tdata->J->data[220]*sx_tmp[80];
  sxdot_tmp[73] = tdata->dxdotdp[73 + ip*model->nx]+tdata->J->data[9]*sx_tmp[0]+tdata->J->data[30]*sx_tmp[1]+tdata->J->data[205]*sx_tmp[72]+tdata->J->data[206]*sx_tmp[73];
  sxdot_tmp[74] = tdata->dxdotdp[74 + ip*model->nx]+tdata->J->data[31]*sx_tmp[1]+tdata->J->data[207]*sx_tmp[73]+tdata->J->data[208]*sx_tmp[74];
  sxdot_tmp[75] = tdata->dxdotdp[75 + ip*model->nx]+tdata->J->data[209]*sx_tmp[74]+tdata->J->data[210]*sx_tmp[75];
  sxdot_tmp[76] = tdata->dxdotdp[76 + ip*model->nx]+tdata->J->data[211]*sx_tmp[75]+tdata->J->data[212]*sx_tmp[76];
  sxdot_tmp[77] = tdata->dxdotdp[77 + ip*model->nx]+tdata->J->data[213]*sx_tmp[76]+tdata->J->data[214]*sx_tmp[77];
  sxdot_tmp[78] = tdata->dxdotdp[78 + ip*model->nx]+tdata->J->data[215]*sx_tmp[77]+tdata->J->data[216]*sx_tmp[78];
  sxdot_tmp[79] = tdata->dxdotdp[79 + ip*model->nx]+tdata->J->data[217]*sx_tmp[78]+tdata->J->data[218]*sx_tmp[79];
  sxdot_tmp[80] = tdata->dxdotdp[80 + ip*model->nx]+tdata->J->data[219]*sx_tmp[79]+tdata->J->data[221]*sx_tmp[80];
  sxdot_tmp[81] = tdata->dxdotdp[81 + ip*model->nx]+tdata->J->data[10]*sx_tmp[0]+tdata->J->data[222]*sx_tmp[81]+tdata->J->data[238]*sx_tmp[89];
  sxdot_tmp[82] = tdata->dxdotdp[82 + ip*model->nx]+tdata->J->data[11]*sx_tmp[0]+tdata->J->data[32]*sx_tmp[1]+tdata->J->data[223]*sx_tmp[81]+tdata->J->data[224]*sx_tmp[82];
  sxdot_tmp[83] = tdata->dxdotdp[83 + ip*model->nx]+tdata->J->data[33]*sx_tmp[1]+tdata->J->data[225]*sx_tmp[82]+tdata->J->data[226]*sx_tmp[83];
  sxdot_tmp[84] = tdata->dxdotdp[84 + ip*model->nx]+tdata->J->data[227]*sx_tmp[83]+tdata->J->data[228]*sx_tmp[84];
  sxdot_tmp[85] = tdata->dxdotdp[85 + ip*model->nx]+tdata->J->data[229]*sx_tmp[84]+tdata->J->data[230]*sx_tmp[85];
  sxdot_tmp[86] = tdata->dxdotdp[86 + ip*model->nx]+tdata->J->data[231]*sx_tmp[85]+tdata->J->data[232]*sx_tmp[86];
  sxdot_tmp[87] = tdata->dxdotdp[87 + ip*model->nx]+tdata->J->data[233]*sx_tmp[86]+tdata->J->data[234]*sx_tmp[87];
  sxdot_tmp[88] = tdata->dxdotdp[88 + ip*model->nx]+tdata->J->data[235]*sx_tmp[87]+tdata->J->data[236]*sx_tmp[88];
  sxdot_tmp[89] = tdata->dxdotdp[89 + ip*model->nx]+tdata->J->data[237]*sx_tmp[88]+tdata->J->data[239]*sx_tmp[89];
  sxdot_tmp[90] = tdata->dxdotdp[90 + ip*model->nx]+tdata->J->data[12]*sx_tmp[0]+tdata->J->data[240]*sx_tmp[90]+tdata->J->data[256]*sx_tmp[98];
  sxdot_tmp[91] = tdata->dxdotdp[91 + ip*model->nx]+tdata->J->data[13]*sx_tmp[0]+tdata->J->data[34]*sx_tmp[1]+tdata->J->data[241]*sx_tmp[90]+tdata->J->data[242]*sx_tmp[91];
  sxdot_tmp[92] = tdata->dxdotdp[92 + ip*model->nx]+tdata->J->data[35]*sx_tmp[1]+tdata->J->data[243]*sx_tmp[91]+tdata->J->data[244]*sx_tmp[92];
  sxdot_tmp[93] = tdata->dxdotdp[93 + ip*model->nx]+tdata->J->data[245]*sx_tmp[92]+tdata->J->data[246]*sx_tmp[93];
  sxdot_tmp[94] = tdata->dxdotdp[94 + ip*model->nx]+tdata->J->data[247]*sx_tmp[93]+tdata->J->data[248]*sx_tmp[94];
  sxdot_tmp[95] = tdata->dxdotdp[95 + ip*model->nx]+tdata->J->data[249]*sx_tmp[94]+tdata->J->data[250]*sx_tmp[95];
  sxdot_tmp[96] = tdata->dxdotdp[96 + ip*model->nx]+tdata->J->data[251]*sx_tmp[95]+tdata->J->data[252]*sx_tmp[96];
  sxdot_tmp[97] = tdata->dxdotdp[97 + ip*model->nx]+tdata->J->data[253]*sx_tmp[96]+tdata->J->data[254]*sx_tmp[97];
  sxdot_tmp[98] = tdata->dxdotdp[98 + ip*model->nx]+tdata->J->data[255]*sx_tmp[97]+tdata->J->data[257]*sx_tmp[98];
  sxdot_tmp[99] = tdata->dxdotdp[99 + ip*model->nx]+tdata->J->data[258]*sx_tmp[99]+tdata->J->data[274]*sx_tmp[107];
  sxdot_tmp[100] = tdata->dxdotdp[100 + ip*model->nx]+tdata->J->data[36]*sx_tmp[1]+tdata->J->data[259]*sx_tmp[99]+tdata->J->data[260]*sx_tmp[100];
  sxdot_tmp[101] = tdata->dxdotdp[101 + ip*model->nx]+tdata->J->data[37]*sx_tmp[1]+tdata->J->data[261]*sx_tmp[100]+tdata->J->data[262]*sx_tmp[101];
  sxdot_tmp[102] = tdata->dxdotdp[102 + ip*model->nx]+tdata->J->data[263]*sx_tmp[101]+tdata->J->data[264]*sx_tmp[102];
  sxdot_tmp[103] = tdata->dxdotdp[103 + ip*model->nx]+tdata->J->data[265]*sx_tmp[102]+tdata->J->data[266]*sx_tmp[103];
  sxdot_tmp[104] = tdata->dxdotdp[104 + ip*model->nx]+tdata->J->data[267]*sx_tmp[103]+tdata->J->data[268]*sx_tmp[104];
  sxdot_tmp[105] = tdata->dxdotdp[105 + ip*model->nx]+tdata->J->data[269]*sx_tmp[104]+tdata->J->data[270]*sx_tmp[105];
  sxdot_tmp[106] = tdata->dxdotdp[106 + ip*model->nx]+tdata->J->data[271]*sx_tmp[105]+tdata->J->data[272]*sx_tmp[106];
  sxdot_tmp[107] = tdata->dxdotdp[107 + ip*model->nx]+tdata->J->data[273]*sx_tmp[106]+tdata->J->data[275]*sx_tmp[107];
  sxdot_tmp[108] = tdata->dxdotdp[108 + ip*model->nx]+tdata->J->data[276]*sx_tmp[108]+tdata->J->data[292]*sx_tmp[116];
  sxdot_tmp[109] = tdata->dxdotdp[109 + ip*model->nx]+tdata->J->data[38]*sx_tmp[1]+tdata->J->data[277]*sx_tmp[108]+tdata->J->data[278]*sx_tmp[109];
  sxdot_tmp[110] = tdata->dxdotdp[110 + ip*model->nx]+tdata->J->data[39]*sx_tmp[1]+tdata->J->data[279]*sx_tmp[109]+tdata->J->data[280]*sx_tmp[110];
  sxdot_tmp[111] = tdata->dxdotdp[111 + ip*model->nx]+tdata->J->data[281]*sx_tmp[110]+tdata->J->data[282]*sx_tmp[111];
  sxdot_tmp[112] = tdata->dxdotdp[112 + ip*model->nx]+tdata->J->data[283]*sx_tmp[111]+tdata->J->data[284]*sx_tmp[112];
  sxdot_tmp[113] = tdata->dxdotdp[113 + ip*model->nx]+tdata->J->data[285]*sx_tmp[112]+tdata->J->data[286]*sx_tmp[113];
  sxdot_tmp[114] = tdata->dxdotdp[114 + ip*model->nx]+tdata->J->data[287]*sx_tmp[113]+tdata->J->data[288]*sx_tmp[114];
  sxdot_tmp[115] = tdata->dxdotdp[115 + ip*model->nx]+tdata->J->data[289]*sx_tmp[114]+tdata->J->data[290]*sx_tmp[115];
  sxdot_tmp[116] = tdata->dxdotdp[116 + ip*model->nx]+tdata->J->data[291]*sx_tmp[115]+tdata->J->data[293]*sx_tmp[116];
  sxdot_tmp[117] = tdata->dxdotdp[117 + ip*model->nx]+tdata->J->data[294]*sx_tmp[117]+tdata->J->data[310]*sx_tmp[125];
  sxdot_tmp[118] = tdata->dxdotdp[118 + ip*model->nx]+tdata->J->data[40]*sx_tmp[1]+tdata->J->data[295]*sx_tmp[117]+tdata->J->data[296]*sx_tmp[118];
  sxdot_tmp[119] = tdata->dxdotdp[119 + ip*model->nx]+tdata->J->data[41]*sx_tmp[1]+tdata->J->data[297]*sx_tmp[118]+tdata->J->data[298]*sx_tmp[119];
  sxdot_tmp[120] = tdata->dxdotdp[120 + ip*model->nx]+tdata->J->data[299]*sx_tmp[119]+tdata->J->data[300]*sx_tmp[120];
  sxdot_tmp[121] = tdata->dxdotdp[121 + ip*model->nx]+tdata->J->data[301]*sx_tmp[120]+tdata->J->data[302]*sx_tmp[121];
  sxdot_tmp[122] = tdata->dxdotdp[122 + ip*model->nx]+tdata->J->data[303]*sx_tmp[121]+tdata->J->data[304]*sx_tmp[122];
  sxdot_tmp[123] = tdata->dxdotdp[123 + ip*model->nx]+tdata->J->data[305]*sx_tmp[122]+tdata->J->data[306]*sx_tmp[123];
  sxdot_tmp[124] = tdata->dxdotdp[124 + ip*model->nx]+tdata->J->data[307]*sx_tmp[123]+tdata->J->data[308]*sx_tmp[124];
  sxdot_tmp[125] = tdata->dxdotdp[125 + ip*model->nx]+tdata->J->data[309]*sx_tmp[124]+tdata->J->data[311]*sx_tmp[125];
  sxdot_tmp[126] = tdata->dxdotdp[126 + ip*model->nx]+tdata->J->data[312]*sx_tmp[126]+tdata->J->data[328]*sx_tmp[134];
  sxdot_tmp[127] = tdata->dxdotdp[127 + ip*model->nx]+tdata->J->data[42]*sx_tmp[1]+tdata->J->data[313]*sx_tmp[126]+tdata->J->data[314]*sx_tmp[127];
  sxdot_tmp[128] = tdata->dxdotdp[128 + ip*model->nx]+tdata->J->data[43]*sx_tmp[1]+tdata->J->data[315]*sx_tmp[127]+tdata->J->data[316]*sx_tmp[128];
  sxdot_tmp[129] = tdata->dxdotdp[129 + ip*model->nx]+tdata->J->data[317]*sx_tmp[128]+tdata->J->data[318]*sx_tmp[129];
  sxdot_tmp[130] = tdata->dxdotdp[130 + ip*model->nx]+tdata->J->data[319]*sx_tmp[129]+tdata->J->data[320]*sx_tmp[130];
  sxdot_tmp[131] = tdata->dxdotdp[131 + ip*model->nx]+tdata->J->data[321]*sx_tmp[130]+tdata->J->data[322]*sx_tmp[131];
  sxdot_tmp[132] = tdata->dxdotdp[132 + ip*model->nx]+tdata->J->data[323]*sx_tmp[131]+tdata->J->data[324]*sx_tmp[132];
  sxdot_tmp[133] = tdata->dxdotdp[133 + ip*model->nx]+tdata->J->data[325]*sx_tmp[132]+tdata->J->data[326]*sx_tmp[133];
  sxdot_tmp[134] = tdata->dxdotdp[134 + ip*model->nx]+tdata->J->data[327]*sx_tmp[133]+tdata->J->data[329]*sx_tmp[134];
  sxdot_tmp[135] = tdata->dxdotdp[135 + ip*model->nx]+tdata->J->data[330]*sx_tmp[135]+tdata->J->data[346]*sx_tmp[143];
  sxdot_tmp[136] = tdata->dxdotdp[136 + ip*model->nx]+tdata->J->data[44]*sx_tmp[1]+tdata->J->data[331]*sx_tmp[135]+tdata->J->data[332]*sx_tmp[136];
  sxdot_tmp[137] = tdata->dxdotdp[137 + ip*model->nx]+tdata->J->data[45]*sx_tmp[1]+tdata->J->data[333]*sx_tmp[136]+tdata->J->data[334]*sx_tmp[137];
  sxdot_tmp[138] = tdata->dxdotdp[138 + ip*model->nx]+tdata->J->data[335]*sx_tmp[137]+tdata->J->data[336]*sx_tmp[138];
  sxdot_tmp[139] = tdata->dxdotdp[139 + ip*model->nx]+tdata->J->data[337]*sx_tmp[138]+tdata->J->data[338]*sx_tmp[139];
  sxdot_tmp[140] = tdata->dxdotdp[140 + ip*model->nx]+tdata->J->data[339]*sx_tmp[139]+tdata->J->data[340]*sx_tmp[140];
  sxdot_tmp[141] = tdata->dxdotdp[141 + ip*model->nx]+tdata->J->data[341]*sx_tmp[140]+tdata->J->data[342]*sx_tmp[141];
  sxdot_tmp[142] = tdata->dxdotdp[142 + ip*model->nx]+tdata->J->data[343]*sx_tmp[141]+tdata->J->data[344]*sx_tmp[142];
  sxdot_tmp[143] = tdata->dxdotdp[143 + ip*model->nx]+tdata->J->data[345]*sx_tmp[142]+tdata->J->data[347]*sx_tmp[143];
  sxdot_tmp[144] = tdata->dxdotdp[144 + ip*model->nx]+tdata->J->data[348]*sx_tmp[144]+tdata->J->data[364]*sx_tmp[152];
  sxdot_tmp[145] = tdata->dxdotdp[145 + ip*model->nx]+tdata->J->data[46]*sx_tmp[1]+tdata->J->data[349]*sx_tmp[144]+tdata->J->data[350]*sx_tmp[145];
  sxdot_tmp[146] = tdata->dxdotdp[146 + ip*model->nx]+tdata->J->data[47]*sx_tmp[1]+tdata->J->data[351]*sx_tmp[145]+tdata->J->data[352]*sx_tmp[146];
  sxdot_tmp[147] = tdata->dxdotdp[147 + ip*model->nx]+tdata->J->data[353]*sx_tmp[146]+tdata->J->data[354]*sx_tmp[147];
  sxdot_tmp[148] = tdata->dxdotdp[148 + ip*model->nx]+tdata->J->data[355]*sx_tmp[147]+tdata->J->data[356]*sx_tmp[148];
  sxdot_tmp[149] = tdata->dxdotdp[149 + ip*model->nx]+tdata->J->data[357]*sx_tmp[148]+tdata->J->data[358]*sx_tmp[149];
  sxdot_tmp[150] = tdata->dxdotdp[150 + ip*model->nx]+tdata->J->data[359]*sx_tmp[149]+tdata->J->data[360]*sx_tmp[150];
  sxdot_tmp[151] = tdata->dxdotdp[151 + ip*model->nx]+tdata->J->data[361]*sx_tmp[150]+tdata->J->data[362]*sx_tmp[151];
  sxdot_tmp[152] = tdata->dxdotdp[152 + ip*model->nx]+tdata->J->data[363]*sx_tmp[151]+tdata->J->data[365]*sx_tmp[152];
  sxdot_tmp[153] = tdata->dxdotdp[153 + ip*model->nx]+tdata->J->data[366]*sx_tmp[153]+tdata->J->data[382]*sx_tmp[161];
  sxdot_tmp[154] = tdata->dxdotdp[154 + ip*model->nx]+tdata->J->data[48]*sx_tmp[1]+tdata->J->data[367]*sx_tmp[153]+tdata->J->data[368]*sx_tmp[154];
  sxdot_tmp[155] = tdata->dxdotdp[155 + ip*model->nx]+tdata->J->data[49]*sx_tmp[1]+tdata->J->data[369]*sx_tmp[154]+tdata->J->data[370]*sx_tmp[155];
  sxdot_tmp[156] = tdata->dxdotdp[156 + ip*model->nx]+tdata->J->data[371]*sx_tmp[155]+tdata->J->data[372]*sx_tmp[156];
  sxdot_tmp[157] = tdata->dxdotdp[157 + ip*model->nx]+tdata->J->data[373]*sx_tmp[156]+tdata->J->data[374]*sx_tmp[157];
  sxdot_tmp[158] = tdata->dxdotdp[158 + ip*model->nx]+tdata->J->data[375]*sx_tmp[157]+tdata->J->data[376]*sx_tmp[158];
  sxdot_tmp[159] = tdata->dxdotdp[159 + ip*model->nx]+tdata->J->data[377]*sx_tmp[158]+tdata->J->data[378]*sx_tmp[159];
  sxdot_tmp[160] = tdata->dxdotdp[160 + ip*model->nx]+tdata->J->data[379]*sx_tmp[159]+tdata->J->data[380]*sx_tmp[160];
  sxdot_tmp[161] = tdata->dxdotdp[161 + ip*model->nx]+tdata->J->data[381]*sx_tmp[160]+tdata->J->data[383]*sx_tmp[161];
return(status);

}


