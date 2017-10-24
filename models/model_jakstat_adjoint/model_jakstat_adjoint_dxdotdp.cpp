
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_w.h"

using namespace amici;

int dxdotdp_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, void *user_data) {
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
memset(tdata->dxdotdp,0,sizeof(realtype)*9*udata->nplist);
status = dwdp_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->w[0]*x_tmp[0];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->w[0]*x_tmp[0];

  } break;

  case 1: {
  tdata->dxdotdp[1 + ip*model->nx] = tdata->w[1]*-2.0;
  tdata->dxdotdp[2 + ip*model->nx] = tdata->w[1];

  } break;

  case 2: {
  tdata->dxdotdp[2 + ip*model->nx] = -x_tmp[2];
  tdata->dxdotdp[3 + ip*model->nx] = (udata->k[0]*x_tmp[2])/udata->k[1];

  } break;

  case 3: {
  tdata->dxdotdp[0 + ip*model->nx] = (udata->k[1]*x_tmp[8])/udata->k[0];
  tdata->dxdotdp[3 + ip*model->nx] = -x_tmp[3];
  tdata->dxdotdp[4 + ip*model->nx] = x_tmp[3]*2.0-x_tmp[4];
  tdata->dxdotdp[5 + ip*model->nx] = x_tmp[4]-x_tmp[5];
  tdata->dxdotdp[6 + ip*model->nx] = x_tmp[5]-x_tmp[6];
  tdata->dxdotdp[7 + ip*model->nx] = x_tmp[6]-x_tmp[7];
  tdata->dxdotdp[8 + ip*model->nx] = x_tmp[7]-x_tmp[8];

  } break;

  case 5: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[0];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[0];

  } break;

  case 6: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[1];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[1];

  } break;

  case 7: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[2];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[2];

  } break;

  case 8: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[3];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[3];

  } break;

  case 9: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->p[0]*x_tmp[0]*tdata->dwdp[4];
  tdata->dxdotdp[1 + ip*model->nx] = tdata->p[0]*x_tmp[0]*tdata->dwdp[4];

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


