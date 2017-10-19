
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_steadystate_dwdp.h"
#include "model_steadystate_w.h"

using namespace amici;

int qBdot_model_steadystate(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data) {
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
realtype *xB_tmp = nullptr;
if(xB)
    xB_tmp = N_VGetArrayPointer(xB);
realtype *dxB_tmp = nullptr;
if(dxB)
    dxB_tmp = N_VGetArrayPointer(dxB);
realtype *qBdot_tmp = nullptr;
if(qBdot)
    qBdot_tmp = N_VGetArrayPointer(qBdot);
int ip;
memset(qBdot_tmp,0,sizeof(realtype)*udata->nplist*model->nJ);
status = dwdp_model_steadystate(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  qBdot_tmp[ip + udata->nplist*0] = tdata->w[1]*xB_tmp[0]*2.0-tdata->w[1]*xB_tmp[1];

  } break;

  case 1: {
  qBdot_tmp[ip + udata->nplist*0] = x_tmp[0]*x_tmp[1]*xB_tmp[0]+x_tmp[0]*x_tmp[1]*xB_tmp[1]-x_tmp[0]*x_tmp[1]*xB_tmp[2];

  } break;

  case 2: {
  qBdot_tmp[ip + udata->nplist*0] = x_tmp[1]*xB_tmp[0]*-2.0+x_tmp[1]*xB_tmp[1];

  } break;

  case 3: {
  qBdot_tmp[ip + udata->nplist*0] = -xB_tmp[0]*tdata->dwdp[0]-xB_tmp[1]*tdata->dwdp[0]+xB_tmp[2]*tdata->dwdp[0];

  } break;

  case 4: {
  qBdot_tmp[ip + udata->nplist*0] = -xB_tmp[0];

  } break;

}
}
for(ip = 0; ip<udata->nplist*model->nJ; ip++) {
   if(amiIsNaN(qBdot_tmp[ip])) {
       qBdot_tmp[ip] = 0;       if(!tdata->nan_qBdot) {
           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_qBdot = TRUE;
       }
   }   if(amiIsInf(qBdot_tmp[ip])) {
       warnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


