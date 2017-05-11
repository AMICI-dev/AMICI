
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_w.h"

int dxdotdp_model_jakstat_adjoint(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(dxdotdp,0,sizeof(realtype)*9*nplist);
status = dwdp_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
  case 0: {
  dxdotdp[(0+ip*9)] = -w_tmp[0]*x_tmp[0];
  dxdotdp[(1+ip*9)] = w_tmp[0]*x_tmp[0];

  } break;

  case 1: {
  dxdotdp[(1+ip*9)] = w_tmp[1]*-2.0;
  dxdotdp[(2+ip*9)] = w_tmp[1];

  } break;

  case 2: {
  dxdotdp[(2+ip*9)] = -x_tmp[2];
  dxdotdp[(3+ip*9)] = (k[0]*x_tmp[2])/k[1];

  } break;

  case 3: {
  dxdotdp[(0+ip*9)] = (k[1]*x_tmp[8])/k[0];
  dxdotdp[(3+ip*9)] = -x_tmp[3];
  dxdotdp[(4+ip*9)] = x_tmp[3]*2.0-x_tmp[4];
  dxdotdp[(5+ip*9)] = x_tmp[4]-x_tmp[5];
  dxdotdp[(6+ip*9)] = x_tmp[5]-x_tmp[6];
  dxdotdp[(7+ip*9)] = x_tmp[6]-x_tmp[7];
  dxdotdp[(8+ip*9)] = x_tmp[7]-x_tmp[8];

  } break;

  case 5: {
  dxdotdp[(0+ip*9)] = -dwdp_tmp[0]*p[0]*x_tmp[0];
  dxdotdp[(1+ip*9)] = dwdp_tmp[0]*p[0]*x_tmp[0];

  } break;

  case 6: {
  dxdotdp[(0+ip*9)] = -dwdp_tmp[1]*p[0]*x_tmp[0];
  dxdotdp[(1+ip*9)] = dwdp_tmp[1]*p[0]*x_tmp[0];

  } break;

  case 7: {
  dxdotdp[(0+ip*9)] = -dwdp_tmp[2]*p[0]*x_tmp[0];
  dxdotdp[(1+ip*9)] = dwdp_tmp[2]*p[0]*x_tmp[0];

  } break;

  case 8: {
  dxdotdp[(0+ip*9)] = -dwdp_tmp[3]*p[0]*x_tmp[0];
  dxdotdp[(1+ip*9)] = dwdp_tmp[3]*p[0]*x_tmp[0];

  } break;

  case 9: {
  dxdotdp[(0+ip*9)] = -dwdp_tmp[4]*p[0]*x_tmp[0];
  dxdotdp[(1+ip*9)] = dwdp_tmp[4]*p[0]*x_tmp[0];

  } break;

}
}
for(ip = 0; ip<nplist; ip++) {
   for(ix = 0; ix<9; ix++) {
       if(amiIsNaN(dxdotdp[ix+ip*9])) {
           dxdotdp[ix+ip*9] = 0;
           if(!udata->am_nan_dxdotdp) {
               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");
               udata->am_nan_dxdotdp = TRUE;
           }
       }
       if(amiIsInf(dxdotdp[ix+ip*9])) {
           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");
           return(-1);
       }
   }
}
return(status);

}


