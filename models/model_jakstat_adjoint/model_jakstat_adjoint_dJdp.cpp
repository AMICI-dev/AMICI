
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int dJdp_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, void *user_data) {
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
memset(tdata->dJdp,0,sizeof(realtype)*81*udata->nplist);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  case 0: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 0] = -tdata->w[0];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 1] = tdata->w[0];

  } break;

  case 1: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 1*model->nx + 1] = tdata->dwdx[0]*-2.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 1*model->nx + 2] = tdata->dwdx[0];

  } break;

  case 2: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 2*model->nx + 2] = -1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 2*model->nx + 3] = udata->k[0]/udata->k[1];

  } break;

  case 3: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 8*model->nx + 0] = udata->k[1]/udata->k[0];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 3*model->nx + 3] = -1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 3*model->nx + 4] = 2.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 4*model->nx + 4] = -1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 4*model->nx + 5] = 1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 5*model->nx + 5] = -1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 5*model->nx + 6] = 1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 6*model->nx + 6] = -1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 6*model->nx + 7] = 1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 7*model->nx + 7] = -1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 7*model->nx + 8] = 1.0;
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 8*model->nx + 8] = -1.0;

  } break;

  case 5: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 0] = -tdata->p[0]*tdata->dwdp[0];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 1] = tdata->p[0]*tdata->dwdp[0];

  } break;

  case 6: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 0] = -tdata->p[0]*tdata->dwdp[1];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 1] = tdata->p[0]*tdata->dwdp[1];

  } break;

  case 7: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 0] = -tdata->p[0]*tdata->dwdp[2];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 1] = tdata->p[0]*tdata->dwdp[2];

  } break;

  case 8: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 0] = -tdata->p[0]*tdata->dwdp[3];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 1] = tdata->p[0]*tdata->dwdp[3];

  } break;

  case 9: {
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 0] = -tdata->p[0]*tdata->dwdp[4];
  tdata->dJdp[udata->plist[ip]*model->nx*model->nx + 0*model->nx + 1] = tdata->p[0]*tdata->dwdp[4];

  } break;

}
}
return(status);

}


