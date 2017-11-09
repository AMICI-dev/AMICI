
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_dwdx.h"
#include "model_events_w.h"

using namespace amici;

void J_model_events(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(J->data,0,sizeof(realtype)*9);
w_model_events(t,x,NULL,tdata);
dwdx_model_events(t,x,NULL,user_data);
  J->data[0+0*3] = -tdata->h[3]*tdata->p[0];
  J->data[1+0*3] = tdata->p[1]*exp(t*(-1.0/1.0E1));
  J->data[1+1*3] = -tdata->p[2];
  J->data[2+2*3] = -1.0;
for(ix = 0; ix<9; ix++) {
   if(amiIsNaN(J->data[ix])) {
       J->data[ix] = 0;
       if(!tdata->nan_J) {
           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_J = TRUE;
       }
   }
   if(amiIsInf(J->data[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");
       return;
   }
}
return;

}


