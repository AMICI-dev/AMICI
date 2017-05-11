
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int J_model_jakstat_adjoint(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(J->data,0,sizeof(realtype)*81);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  J->data[0+0*9] = -p[0]*w_tmp[0];
  J->data[0+8*9] = (k[1]*p[3])/k[0];
  J->data[1+0*9] = p[0]*w_tmp[0];
  J->data[1+1*9] = dwdx_tmp[0]*p[1]*-2.0;
  J->data[2+1*9] = dwdx_tmp[0]*p[1];
  J->data[2+2*9] = -p[2];
  J->data[3+2*9] = (k[0]*p[2])/k[1];
  J->data[3+3*9] = -p[3];
  J->data[4+3*9] = p[3]*2.0;
  J->data[4+4*9] = -p[3];
  J->data[5+4*9] = p[3];
  J->data[5+5*9] = -p[3];
  J->data[6+5*9] = p[3];
  J->data[6+6*9] = -p[3];
  J->data[7+6*9] = p[3];
  J->data[7+7*9] = -p[3];
  J->data[8+7*9] = p[3];
  J->data[8+8*9] = -p[3];
for(ix = 0; ix<81; ix++) {
   if(amiIsNaN(J->data[ix])) {
       J->data[ix] = 0;
       if(!udata->am_nan_J) {
           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_J = TRUE;
       }
   }
   if(amiIsInf(J->data[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


