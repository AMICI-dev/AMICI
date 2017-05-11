
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int xBdot_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
int ix;
memset(xBdot_tmp,0,sizeof(realtype)*9);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  xBdot_tmp[0] = p[0]*w_tmp[0]*xB_tmp[0]-p[0]*w_tmp[0]*xB_tmp[1];
  xBdot_tmp[1] = dwdx_tmp[0]*p[1]*xB_tmp[1]*2.0-dwdx_tmp[0]*p[1]*xB_tmp[2];
  xBdot_tmp[2] = p[2]*xB_tmp[2]-(k[0]*p[2]*xB_tmp[3])/k[1];
  xBdot_tmp[3] = p[3]*xB_tmp[3]-p[3]*xB_tmp[4]*2.0;
  xBdot_tmp[4] = p[3]*xB_tmp[4]-p[3]*xB_tmp[5];
  xBdot_tmp[5] = p[3]*xB_tmp[5]-p[3]*xB_tmp[6];
  xBdot_tmp[6] = p[3]*xB_tmp[6]-p[3]*xB_tmp[7];
  xBdot_tmp[7] = p[3]*xB_tmp[7]-p[3]*xB_tmp[8];
  xBdot_tmp[8] = p[3]*xB_tmp[8]-(k[1]*p[3]*xB_tmp[0])/k[0];
for(ix = 0; ix<9; ix++) {
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


