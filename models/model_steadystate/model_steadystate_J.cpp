
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_steadystate_dwdx.h"
#include "model_steadystate_w.h"

int J_model_steadystate(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(J->data,0,sizeof(realtype)*9);
status = w_model_steadystate(t,x,NULL,tdata);
status = dwdx_model_steadystate(t,x,NULL,user_data);
  J->data[0+0*3] = -tdata->p[1]*x_tmp[1]-tdata->p[0]*tdata->dwdx[0]*2.0;
  J->data[0+1*3] = tdata->p[2]*2.0-tdata->p[1]*x_tmp[0];
  J->data[0+2*3] = tdata->dwdx[1];
  J->data[1+0*3] = -tdata->p[1]*x_tmp[1]+tdata->p[0]*tdata->dwdx[0];
  J->data[1+1*3] = -tdata->p[2]-tdata->p[1]*x_tmp[0];
  J->data[1+2*3] = tdata->dwdx[1];
  J->data[2+0*3] = tdata->p[1]*x_tmp[1];
  J->data[2+1*3] = tdata->p[1]*x_tmp[0];
  J->data[2+2*3] = -udata->k[3]-tdata->dwdx[1];
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
       return(-1);
   }
}
return(status);

}


