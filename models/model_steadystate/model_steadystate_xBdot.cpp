
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_dwdx.h"
#include "model_steadystate_w.h"

int xBdot_model_steadystate(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
int ix;
memset(xBdot_tmp,0,sizeof(realtype)*3);
status = w_model_steadystate(t,x,NULL,user_data);
status = dwdx_model_steadystate(t,x,NULL,user_data);
  xBdot_tmp[0] = xB_tmp[0]*(dwdx_tmp[0]*p[0]*2.0+p[1]*x_tmp[1])-xB_tmp[1]*(dwdx_tmp[0]*p[0]-p[1]*x_tmp[1])-p[1]*xB_tmp[2]*x_tmp[1];
  xBdot_tmp[1] = -xB_tmp[0]*(p[2]*2.0-p[1]*x_tmp[0])+xB_tmp[1]*(p[2]+p[1]*x_tmp[0])-p[1]*xB_tmp[2]*x_tmp[0];
  xBdot_tmp[2] = xB_tmp[2]*(dwdx_tmp[1]+k[3])-dwdx_tmp[1]*xB_tmp[0]-dwdx_tmp[1]*xB_tmp[1];
for(ix = 0; ix<3; ix++) {
   if(amiIsNaN(xBdot_tmp[ix])) {
       xBdot_tmp[ix] = 0;       if(!udata->am_nan_xBdot) {
           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_xBdot = TRUE;
       }
   }   if(amiIsInf(xBdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


