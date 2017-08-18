
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

int xdot_model_neuron_o2(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(xdot_tmp,0,sizeof(realtype)*10);
status = w_model_neuron_o2(t,x,NULL,tdata);
  xdot_tmp[0] = udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2;
  xdot_tmp[1] = -udata->p[0]*(x_tmp[1]-udata->p[1]*x_tmp[0]);
  xdot_tmp[2] = -x_tmp[3]+tdata->w[1]*x_tmp[2];
  xdot_tmp[3] = -x_tmp[1]+udata->p[1]*x_tmp[0]-udata->p[0]*x_tmp[3]+udata->p[0]*udata->p[1]*x_tmp[2];
  xdot_tmp[4] = -x_tmp[5]+tdata->w[1]*x_tmp[4];
  xdot_tmp[5] = udata->p[0]*x_tmp[0]-udata->p[0]*x_tmp[5]+udata->p[0]*udata->p[1]*x_tmp[4];
  xdot_tmp[6] = -x_tmp[7]+tdata->w[1]*x_tmp[6];
  xdot_tmp[7] = -udata->p[0]*x_tmp[7]+udata->p[0]*udata->p[1]*x_tmp[6];
  xdot_tmp[8] = -x_tmp[9]+tdata->w[1]*x_tmp[8];
  xdot_tmp[9] = -udata->p[0]*x_tmp[9]+udata->p[0]*udata->p[1]*x_tmp[8];
for(ix = 0; ix<10; ix++) {
   if(amiIsNaN(xdot_tmp[ix])) {
       xdot_tmp[ix] = 0;
       if(!tdata->nan_xdot) {
           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_xdot = TRUE;
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


