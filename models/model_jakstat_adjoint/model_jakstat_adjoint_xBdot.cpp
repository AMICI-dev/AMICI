
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int xBdot_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, void *user_data) {
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
realtype *xBdot_tmp = nullptr;
if(xBdot)
    xBdot_tmp = N_VGetArrayPointer(xBdot);
int ix;
memset(xBdot_tmp,0,sizeof(realtype)*9);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  xBdot_tmp[0] = tdata->p[0]*tdata->w[0]*xB_tmp[0]-tdata->p[0]*tdata->w[0]*xB_tmp[1];
  xBdot_tmp[1] = tdata->p[1]*xB_tmp[1]*tdata->dwdx[0]*2.0-tdata->p[1]*xB_tmp[2]*tdata->dwdx[0];
  xBdot_tmp[2] = tdata->p[2]*xB_tmp[2]-(udata->k[0]*tdata->p[2]*xB_tmp[3])/udata->k[1];
  xBdot_tmp[3] = tdata->p[3]*xB_tmp[3]-tdata->p[3]*xB_tmp[4]*2.0;
  xBdot_tmp[4] = tdata->p[3]*xB_tmp[4]-tdata->p[3]*xB_tmp[5];
  xBdot_tmp[5] = tdata->p[3]*xB_tmp[5]-tdata->p[3]*xB_tmp[6];
  xBdot_tmp[6] = tdata->p[3]*xB_tmp[6]-tdata->p[3]*xB_tmp[7];
  xBdot_tmp[7] = tdata->p[3]*xB_tmp[7]-tdata->p[3]*xB_tmp[8];
  xBdot_tmp[8] = tdata->p[3]*xB_tmp[8]-(udata->k[1]*tdata->p[3]*xB_tmp[0])/udata->k[0];
for(ix = 0; ix<9; ix++) {
   if(amiIsNaN(xBdot_tmp[ix])) {
       xBdot_tmp[ix] = 0;       if(!tdata->nan_xBdot) {
           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_xBdot = TRUE;
       }
   }   if(amiIsInf(xBdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


