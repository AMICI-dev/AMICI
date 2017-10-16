
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_nested_events_dwdx.h"
#include "model_nested_events_w.h"

int JDiag_model_nested_events(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *JDiag_tmp = nullptr;
if(JDiag)
    JDiag_tmp = N_VGetArrayPointer(JDiag);
int ix;
memset(JDiag_tmp,0,sizeof(realtype)*1);
status = w_model_nested_events(t,x,NULL,tdata);
status = dwdx_model_nested_events(t,x,NULL,user_data);
  JDiag_tmp[0] = -tdata->p[4]+tdata->h[1]*tdata->p[3];
for(ix = 0; ix<1; ix++) {
   if(amiIsNaN(JDiag_tmp[ix])) {
       JDiag_tmp[ix] = 0;
       if(!tdata->nan_JDiag) {
           warnMsgIdAndTxt("AMICI:mex:fJDiag:NaN","AMICI replaced a NaN value on Jacobian diagonal and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_JDiag = TRUE;
       }
   }
   if(amiIsInf(JDiag_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJDiag:Inf","AMICI encountered an Inf value on Jacobian diagonal! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


