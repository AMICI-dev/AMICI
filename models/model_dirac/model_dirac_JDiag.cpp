
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_JDiag.h"
#include "model_dirac_w.h"

int JDiag_model_dirac(realtype t, N_Vector JDiag, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *JDiag_tmp = N_VGetArrayPointer(JDiag);
int ix;
memset(JDiag_tmp,0,sizeof(realtype)*2);
status = w_model_dirac(t,x,NULL,user_data);
  JDiag_tmp[0+0*2] = -p[0];
  JDiag_tmp[1+0*2] = -p[3];
for(ix = 0; ix<2; ix++) {
   if(amiIsNaN(JDiag_tmp[ix])) {
       JDiag_tmp[ix] = 0;
       if(!udata->am_nan_JDiag) {
           warnMsgIdAndTxt("AMICI:mex:fJDiag:NaN","AMICI replaced a NaN value on Jacobian diagonal and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_JDiag = TRUE;
       }
   }
   if(amiIsInf(JDiag_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJDiag:Inf","AMICI encountered an Inf value on Jacobian diagonal! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


