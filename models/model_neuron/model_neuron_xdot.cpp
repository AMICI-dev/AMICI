
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_w.h"

int xdot_model_neuron(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(xdot_tmp,0,sizeof(realtype)*2);
status = w_model_neuron(t,x,NULL,tdata);
  xdot_tmp[0] = udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2;
  xdot_tmp[1] = -udata->p[0]*(x_tmp[1]-udata->p[1]*x_tmp[0]);
for(ix = 0; ix<2; ix++) {
   if(amiIsNaN(xdot_tmp[ix])) {
       xdot_tmp[ix] = 0;
       if(!udata->nan_xdot) {
           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->nan_xdot = TRUE;
       }
   }
   if(amiIsInf(xdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");
       return(-1);
   }   if(udata->qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {
       xdot_tmp[ix] = -xdot_tmp[ix];
   }
}
return(status);

}


