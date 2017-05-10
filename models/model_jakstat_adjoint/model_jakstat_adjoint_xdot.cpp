
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int xdot_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(xdot_tmp,0,sizeof(realtype)*9);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  xdot_tmp[0] = (k[1]*p[3]*x_tmp[8]-k[0]*p[0]*w_tmp[0]*x_tmp[0])/k[0];
  xdot_tmp[1] = p[1]*w_tmp[1]*-2.0+p[0]*w_tmp[0]*x_tmp[0];
  xdot_tmp[2] = p[1]*w_tmp[1]-p[2]*x_tmp[2];
  xdot_tmp[3] = (k[0]*p[2]*x_tmp[2]-k[1]*p[3]*x_tmp[3])/k[1];
  xdot_tmp[4] = p[3]*(x_tmp[3]*2.0-x_tmp[4]);
  xdot_tmp[5] = p[3]*(x_tmp[4]-x_tmp[5]);
  xdot_tmp[6] = p[3]*(x_tmp[5]-x_tmp[6]);
  xdot_tmp[7] = p[3]*(x_tmp[6]-x_tmp[7]);
  xdot_tmp[8] = p[3]*(x_tmp[7]-x_tmp[8]);
for(ix = 0; ix<9; ix++) {
   if(amiIsNaN(xdot_tmp[ix])) {
       xdot_tmp[ix] = 0;
       if(!udata->am_nan_xdot) {
           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_xdot = TRUE;
       }
   }
   if(amiIsInf(xdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");
       return(-1);
   }   if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {
       xdot_tmp[ix] = -xdot_tmp[ix];
   }
}
return(status);

}


