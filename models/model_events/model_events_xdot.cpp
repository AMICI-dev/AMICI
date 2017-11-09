
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

using namespace amici;

void xdot_model_events(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(xdot_tmp,0,sizeof(realtype)*3);
w_model_events(t,x,NULL,tdata);
  xdot_tmp[0] = -tdata->h[3]*tdata->p[0]*x_tmp[0];
  xdot_tmp[1] = -tdata->p[2]*x_tmp[1]+tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1));
  xdot_tmp[2] = tdata->h[2]-x_tmp[2];
for(ix = 0; ix<3; ix++) {
   if(amiIsNaN(xdot_tmp[ix])) {
       xdot_tmp[ix] = 0;
       if(!tdata->nan_xdot) {
           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_xdot = TRUE;
       }
   }
   if(amiIsInf(xdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");
       return;
   }   if(udata->qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {
       xdot_tmp[ix] = -xdot_tmp[ix];
   }
}
return;

}


