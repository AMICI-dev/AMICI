
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int J_model_jakstat_adjoint(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
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
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(J->data,0,sizeof(realtype)*81);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  J->data[0+0*9] = -tdata->p[0]*tdata->w[0];
  J->data[0+8*9] = (udata->k[1]*tdata->p[3])/udata->k[0];
  J->data[1+0*9] = tdata->p[0]*tdata->w[0];
  J->data[1+1*9] = tdata->p[1]*tdata->dwdx[0]*-2.0;
  J->data[2+1*9] = tdata->p[1]*tdata->dwdx[0];
  J->data[2+2*9] = -tdata->p[2];
  J->data[3+2*9] = (udata->k[0]*tdata->p[2])/udata->k[1];
  J->data[3+3*9] = -tdata->p[3];
  J->data[4+3*9] = tdata->p[3]*2.0;
  J->data[4+4*9] = -tdata->p[3];
  J->data[5+4*9] = tdata->p[3];
  J->data[5+5*9] = -tdata->p[3];
  J->data[6+5*9] = tdata->p[3];
  J->data[6+6*9] = -tdata->p[3];
  J->data[7+6*9] = tdata->p[3];
  J->data[7+7*9] = -tdata->p[3];
  J->data[8+7*9] = tdata->p[3];
  J->data[8+8*9] = -tdata->p[3];
for(ix = 0; ix<81; ix++) {
   if(amiIsNaN(J->data[ix])) {
       J->data[ix] = 0;
       if(!tdata->nan_J) {
           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_J = TRUE;
       }
   }
   if(amiIsInf(J->data[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


