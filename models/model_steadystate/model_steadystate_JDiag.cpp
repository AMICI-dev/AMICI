
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_dwdx.h"
#include "model_steadystate_w.h"

int JDiag_model_steadystate(realtype t, N_Vector JDiag, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *JDiag_tmp = N_VGetArrayPointer(JDiag);
int ix;
memset(JDiag_tmp,0,sizeof(realtype)*3);
status = w_model_steadystate(t,x,NULL,user_data);
status = dwdx_model_steadystate(t,x,NULL,user_data);
  JDiag_tmp[0+0*3] = -udata->p[1]*x_tmp[1]-udata->p[0]*udata->dwdx[0]*2.0;
  JDiag_tmp[1+0*3] = -udata->p[2]-udata->p[1]*x_tmp[0];
  JDiag_tmp[2+0*3] = -udata->k[3]-udata->dwdx[1];
for(ix = 0; ix<3; ix++) {
   if(amiIsNaN(JDiag_tmp[ix])) {
       JDiag_tmp[ix] = 0;
       if(!udata->nan_JDiag) {
           warnMsgIdAndTxt("AMICI:mex:fJDiag:NaN","AMICI replaced a NaN value on Jacobian diagonal and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->nan_JDiag = TRUE;
       }
   }
   if(amiIsInf(JDiag_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJDiag:Inf","AMICI encountered an Inf value on Jacobian diagonal! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


