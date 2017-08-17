
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_w.h"

int dxdotdp_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(udata->dxdotdp,0,sizeof(realtype)*9*udata->nplist);
status = dwdp_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->w[0]*x_tmp[0];
  udata->dxdotdp[1 + ip*udata->nx] = udata->w[0]*x_tmp[0];

  } break;

  case 1: {
  udata->dxdotdp[1 + ip*udata->nx] = udata->w[1]*-2.0;
  udata->dxdotdp[2 + ip*udata->nx] = udata->w[1];

  } break;

  case 2: {
  udata->dxdotdp[2 + ip*udata->nx] = -x_tmp[2];
  udata->dxdotdp[3 + ip*udata->nx] = (udata->k[0]*x_tmp[2])/udata->k[1];

  } break;

  case 3: {
  udata->dxdotdp[0 + ip*udata->nx] = (udata->k[1]*x_tmp[8])/udata->k[0];
  udata->dxdotdp[3 + ip*udata->nx] = -x_tmp[3];
  udata->dxdotdp[4 + ip*udata->nx] = x_tmp[3]*2.0-x_tmp[4];
  udata->dxdotdp[5 + ip*udata->nx] = x_tmp[4]-x_tmp[5];
  udata->dxdotdp[6 + ip*udata->nx] = x_tmp[5]-x_tmp[6];
  udata->dxdotdp[7 + ip*udata->nx] = x_tmp[6]-x_tmp[7];
  udata->dxdotdp[8 + ip*udata->nx] = x_tmp[7]-x_tmp[8];

  } break;

  case 5: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[0];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[0];

  } break;

  case 6: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[1];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[1];

  } break;

  case 7: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[2];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[2];

  } break;

  case 8: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[3];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[3];

  } break;

  case 9: {
  udata->dxdotdp[0 + ip*udata->nx] = -udata->p[0]*x_tmp[0]*udata->dwdp[4];
  udata->dxdotdp[1 + ip*udata->nx] = udata->p[0]*x_tmp[0]*udata->dwdp[4];

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


