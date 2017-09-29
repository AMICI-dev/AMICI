
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_dwdx.h"
#include "model_jakstat_adjoint_o2_w.h"

int JDiag_model_jakstat_adjoint_o2(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *JDiag_tmp = N_VGetArrayPointer(JDiag);
int ix;
memset(JDiag_tmp,0,sizeof(realtype)*162);
status = w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
status = dwdx_model_jakstat_adjoint_o2(t,x,NULL,user_data);
  JDiag_tmp[0+0*162] = -udata->k[0]*tdata->p[0]*tdata->w[0]*tdata->w[2];
  JDiag_tmp[1+0*162] = tdata->p[1]*tdata->dwdx[0]*-2.0;
  JDiag_tmp[2+0*162] = -tdata->p[2];
  JDiag_tmp[3+0*162] = -udata->k[1]*tdata->p[3]*tdata->w[3];
  JDiag_tmp[4+0*162] = -tdata->p[3];
  JDiag_tmp[5+0*162] = -tdata->p[3];
  JDiag_tmp[6+0*162] = -tdata->p[3];
  JDiag_tmp[7+0*162] = -tdata->p[3];
  JDiag_tmp[8+0*162] = -tdata->p[3];
  JDiag_tmp[9+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[10+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[11+0*162] = -tdata->p[2];
  JDiag_tmp[12+0*162] = -tdata->p[3];
  JDiag_tmp[13+0*162] = -tdata->p[3];
  JDiag_tmp[14+0*162] = -tdata->p[3];
  JDiag_tmp[15+0*162] = -tdata->p[3];
  JDiag_tmp[16+0*162] = -tdata->p[3];
  JDiag_tmp[17+0*162] = -tdata->p[3];
  JDiag_tmp[18+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[19+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[20+0*162] = -tdata->p[2];
  JDiag_tmp[21+0*162] = -tdata->p[3];
  JDiag_tmp[22+0*162] = -tdata->p[3];
  JDiag_tmp[23+0*162] = -tdata->p[3];
  JDiag_tmp[24+0*162] = -tdata->p[3];
  JDiag_tmp[25+0*162] = -tdata->p[3];
  JDiag_tmp[26+0*162] = -tdata->p[3];
  JDiag_tmp[27+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[28+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[29+0*162] = -tdata->p[2];
  JDiag_tmp[30+0*162] = -tdata->p[3];
  JDiag_tmp[31+0*162] = -tdata->p[3];
  JDiag_tmp[32+0*162] = -tdata->p[3];
  JDiag_tmp[33+0*162] = -tdata->p[3];
  JDiag_tmp[34+0*162] = -tdata->p[3];
  JDiag_tmp[35+0*162] = -tdata->p[3];
  JDiag_tmp[36+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[37+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[38+0*162] = -tdata->p[2];
  JDiag_tmp[39+0*162] = -tdata->p[3];
  JDiag_tmp[40+0*162] = -tdata->p[3];
  JDiag_tmp[41+0*162] = -tdata->p[3];
  JDiag_tmp[42+0*162] = -tdata->p[3];
  JDiag_tmp[43+0*162] = -tdata->p[3];
  JDiag_tmp[44+0*162] = -tdata->p[3];
  JDiag_tmp[45+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[46+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[47+0*162] = -tdata->p[2];
  JDiag_tmp[48+0*162] = -tdata->p[3];
  JDiag_tmp[49+0*162] = -tdata->p[3];
  JDiag_tmp[50+0*162] = -tdata->p[3];
  JDiag_tmp[51+0*162] = -tdata->p[3];
  JDiag_tmp[52+0*162] = -tdata->p[3];
  JDiag_tmp[53+0*162] = -tdata->p[3];
  JDiag_tmp[54+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[55+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[56+0*162] = -tdata->p[2];
  JDiag_tmp[57+0*162] = -tdata->p[3];
  JDiag_tmp[58+0*162] = -tdata->p[3];
  JDiag_tmp[59+0*162] = -tdata->p[3];
  JDiag_tmp[60+0*162] = -tdata->p[3];
  JDiag_tmp[61+0*162] = -tdata->p[3];
  JDiag_tmp[62+0*162] = -tdata->p[3];
  JDiag_tmp[63+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[64+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[65+0*162] = -tdata->p[2];
  JDiag_tmp[66+0*162] = -tdata->p[3];
  JDiag_tmp[67+0*162] = -tdata->p[3];
  JDiag_tmp[68+0*162] = -tdata->p[3];
  JDiag_tmp[69+0*162] = -tdata->p[3];
  JDiag_tmp[70+0*162] = -tdata->p[3];
  JDiag_tmp[71+0*162] = -tdata->p[3];
  JDiag_tmp[72+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[73+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[74+0*162] = -tdata->p[2];
  JDiag_tmp[75+0*162] = -tdata->p[3];
  JDiag_tmp[76+0*162] = -tdata->p[3];
  JDiag_tmp[77+0*162] = -tdata->p[3];
  JDiag_tmp[78+0*162] = -tdata->p[3];
  JDiag_tmp[79+0*162] = -tdata->p[3];
  JDiag_tmp[80+0*162] = -tdata->p[3];
  JDiag_tmp[81+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[82+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[83+0*162] = -tdata->p[2];
  JDiag_tmp[84+0*162] = -tdata->p[3];
  JDiag_tmp[85+0*162] = -tdata->p[3];
  JDiag_tmp[86+0*162] = -tdata->p[3];
  JDiag_tmp[87+0*162] = -tdata->p[3];
  JDiag_tmp[88+0*162] = -tdata->p[3];
  JDiag_tmp[89+0*162] = -tdata->p[3];
  JDiag_tmp[90+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[91+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[92+0*162] = -tdata->p[2];
  JDiag_tmp[93+0*162] = -tdata->p[3];
  JDiag_tmp[94+0*162] = -tdata->p[3];
  JDiag_tmp[95+0*162] = -tdata->p[3];
  JDiag_tmp[96+0*162] = -tdata->p[3];
  JDiag_tmp[97+0*162] = -tdata->p[3];
  JDiag_tmp[98+0*162] = -tdata->p[3];
  JDiag_tmp[99+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[100+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[101+0*162] = -tdata->p[2];
  JDiag_tmp[102+0*162] = -tdata->p[3];
  JDiag_tmp[103+0*162] = -tdata->p[3];
  JDiag_tmp[104+0*162] = -tdata->p[3];
  JDiag_tmp[105+0*162] = -tdata->p[3];
  JDiag_tmp[106+0*162] = -tdata->p[3];
  JDiag_tmp[107+0*162] = -tdata->p[3];
  JDiag_tmp[108+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[109+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[110+0*162] = -tdata->p[2];
  JDiag_tmp[111+0*162] = -tdata->p[3];
  JDiag_tmp[112+0*162] = -tdata->p[3];
  JDiag_tmp[113+0*162] = -tdata->p[3];
  JDiag_tmp[114+0*162] = -tdata->p[3];
  JDiag_tmp[115+0*162] = -tdata->p[3];
  JDiag_tmp[116+0*162] = -tdata->p[3];
  JDiag_tmp[117+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[118+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[119+0*162] = -tdata->p[2];
  JDiag_tmp[120+0*162] = -tdata->p[3];
  JDiag_tmp[121+0*162] = -tdata->p[3];
  JDiag_tmp[122+0*162] = -tdata->p[3];
  JDiag_tmp[123+0*162] = -tdata->p[3];
  JDiag_tmp[124+0*162] = -tdata->p[3];
  JDiag_tmp[125+0*162] = -tdata->p[3];
  JDiag_tmp[126+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[127+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[128+0*162] = -tdata->p[2];
  JDiag_tmp[129+0*162] = -tdata->p[3];
  JDiag_tmp[130+0*162] = -tdata->p[3];
  JDiag_tmp[131+0*162] = -tdata->p[3];
  JDiag_tmp[132+0*162] = -tdata->p[3];
  JDiag_tmp[133+0*162] = -tdata->p[3];
  JDiag_tmp[134+0*162] = -tdata->p[3];
  JDiag_tmp[135+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[136+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[137+0*162] = -tdata->p[2];
  JDiag_tmp[138+0*162] = -tdata->p[3];
  JDiag_tmp[139+0*162] = -tdata->p[3];
  JDiag_tmp[140+0*162] = -tdata->p[3];
  JDiag_tmp[141+0*162] = -tdata->p[3];
  JDiag_tmp[142+0*162] = -tdata->p[3];
  JDiag_tmp[143+0*162] = -tdata->p[3];
  JDiag_tmp[144+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[145+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[146+0*162] = -tdata->p[2];
  JDiag_tmp[147+0*162] = -tdata->p[3];
  JDiag_tmp[148+0*162] = -tdata->p[3];
  JDiag_tmp[149+0*162] = -tdata->p[3];
  JDiag_tmp[150+0*162] = -tdata->p[3];
  JDiag_tmp[151+0*162] = -tdata->p[3];
  JDiag_tmp[152+0*162] = -tdata->p[3];
  JDiag_tmp[153+0*162] = -tdata->p[0]*tdata->w[0];
  JDiag_tmp[154+0*162] = tdata->p[1]*x_tmp[1]*-4.0;
  JDiag_tmp[155+0*162] = -tdata->p[2];
  JDiag_tmp[156+0*162] = -tdata->p[3];
  JDiag_tmp[157+0*162] = -tdata->p[3];
  JDiag_tmp[158+0*162] = -tdata->p[3];
  JDiag_tmp[159+0*162] = -tdata->p[3];
  JDiag_tmp[160+0*162] = -tdata->p[3];
  JDiag_tmp[161+0*162] = -tdata->p[3];
for(ix = 0; ix<162; ix++) {
   if(amiIsNaN(JDiag_tmp[ix])) {
       JDiag_tmp[ix] = 0;
       if(!tdata->nan_JDiag) {
           warnMsgIdAndTxt("AMICI:mex:fJDiag:NaN","AMICI replaced a NaN value on Jacobian diagonal and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_JDiag = TRUE;
       }
   }
   if(amiIsInf(JDiag_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJDiag:Inf","AMICI encountered an Inf value on Jacobian diagonal! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


