
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

int JDiag_model_neuron_o2(realtype t, N_Vector JDiag, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *JDiag_tmp = N_VGetArrayPointer(JDiag);
int ix;
memset(JDiag_tmp,0,sizeof(realtype)*10);
status = w_model_neuron_o2(t,x,NULL,user_data);
status = dwdx_model_neuron_o2(t,x,NULL,user_data);
  JDiag_tmp[0+0*10] = x_tmp[0]*(2.0/2.5E1)+5.0;
  JDiag_tmp[1+0*10] = -udata->p[0];
  JDiag_tmp[2+0*10] = udata->w[1];
  JDiag_tmp[3+0*10] = -udata->p[0];
  JDiag_tmp[4+0*10] = udata->w[1];
  JDiag_tmp[5+0*10] = -udata->p[0];
  JDiag_tmp[6+0*10] = udata->w[1];
  JDiag_tmp[7+0*10] = -udata->p[0];
  JDiag_tmp[8+0*10] = udata->w[1];
  JDiag_tmp[9+0*10] = -udata->p[0];
for(ix = 0; ix<10; ix++) {
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


