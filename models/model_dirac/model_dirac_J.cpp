
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_dwdx.h"
#include "model_dirac_w.h"

int J_model_dirac(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(J->data,0,sizeof(realtype)*4);
status = w_model_dirac(t,x,NULL,user_data);
status = dwdx_model_dirac(t,x,NULL,user_data);
  J->data[0+0*2] = -p[0];
  J->data[1+0*2] = p[2];
  J->data[1+1*2] = -p[3];
for(ix = 0; ix<4; ix++) {
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


