
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_dwdp.h"
#include "model_steadystate_w.h"

int dxdotdp_model_steadystate(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(dxdotdp,0,sizeof(realtype)*3*nplist);
status = dwdp_model_steadystate(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
  case 0: {
  dxdotdp[(0+ip*3)] = w_tmp[1]*-2.0;
  dxdotdp[(1+ip*3)] = w_tmp[1];

  } break;

  case 1: {
  dxdotdp[(0+ip*3)] = -x_tmp[0]*x_tmp[1];
  dxdotdp[(1+ip*3)] = -x_tmp[0]*x_tmp[1];
  dxdotdp[(2+ip*3)] = x_tmp[0]*x_tmp[1];

  } break;

  case 2: {
  dxdotdp[(0+ip*3)] = x_tmp[1]*2.0;
  dxdotdp[(1+ip*3)] = -x_tmp[1];

  } break;

  case 3: {
  dxdotdp[(0+ip*3)] = dwdp_tmp[0];
  dxdotdp[(1+ip*3)] = dwdp_tmp[0];
  dxdotdp[(2+ip*3)] = -dwdp_tmp[0];

  } break;

  case 4: {
  dxdotdp[(0+ip*3)] = 1.0;

  } break;

}
}
for(ip = 0; ip<nplist; ip++) {
   for(ix = 0; ix<3; ix++) {
       if(amiIsNaN(dxdotdp[ix+ip*3])) {
           dxdotdp[ix+ip*3] = 0;
           if(!udata->am_nan_dxdotdp) {
               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");
               udata->am_nan_dxdotdp = TRUE;
           }
       }
       if(amiIsInf(dxdotdp[ix+ip*3])) {
           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");
           return(-1);
       }
   }
}
return(status);

}


