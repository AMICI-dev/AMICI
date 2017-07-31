
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdp.h"
#include "model_neuron_o2_w.h"

int dxdotdp_model_neuron_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(udata->dxdotdp,0,sizeof(realtype)*10*udata->nplist);
status = dwdp_model_neuron_o2(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  udata->dxdotdp[1 + ip*udata->nx] = -x_tmp[1]+udata->p[1]*x_tmp[0];
  udata->dxdotdp[3 + ip*udata->nx] = -x_tmp[3]+udata->p[1]*x_tmp[2];
  udata->dxdotdp[5 + ip*udata->nx] = x_tmp[0]-x_tmp[5]+udata->p[1]*x_tmp[4];
  udata->dxdotdp[7 + ip*udata->nx] = -x_tmp[7]+udata->p[1]*x_tmp[6];
  udata->dxdotdp[9 + ip*udata->nx] = -x_tmp[9]+udata->p[1]*x_tmp[8];

  } break;

  case 1: {
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0];
  udata->dxdotdp[3 + ip*udata->nx] = x_tmp[0]+udata->p[0]*x_tmp[2];
  udata->dxdotdp[5 + ip*udata->nx] = udata->p[0]*x_tmp[4];
  udata->dxdotdp[7 + ip*udata->nx] = udata->p[0]*x_tmp[6];
  udata->dxdotdp[9 + ip*udata->nx] = udata->p[0]*x_tmp[8];

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


