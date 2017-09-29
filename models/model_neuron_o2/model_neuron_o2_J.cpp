
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

int J_model_neuron_o2(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int ix;
memset(J->data,0,sizeof(realtype)*100);
status = w_model_neuron_o2(t,x,NULL,tdata);
status = dwdx_model_neuron_o2(t,x,NULL,user_data);
  J->data[0+0*10] = x_tmp[0]*(2.0/2.5E1)+5.0;
  J->data[0+1*10] = -1.0;
  J->data[1+0*10] = tdata->p[0]*tdata->p[1];
  J->data[1+1*10] = -tdata->p[0];
  J->data[2+0*10] = x_tmp[2]*tdata->dwdx[1];
  J->data[2+2*10] = tdata->w[1];
  J->data[2+3*10] = -1.0;
  J->data[3+0*10] = tdata->p[1];
  J->data[3+1*10] = -1.0;
  J->data[3+2*10] = tdata->p[0]*tdata->p[1];
  J->data[3+3*10] = -tdata->p[0];
  J->data[4+0*10] = x_tmp[4]*tdata->dwdx[1];
  J->data[4+4*10] = tdata->w[1];
  J->data[4+5*10] = -1.0;
  J->data[5+0*10] = tdata->p[0];
  J->data[5+4*10] = tdata->p[0]*tdata->p[1];
  J->data[5+5*10] = -tdata->p[0];
  J->data[6+0*10] = x_tmp[6]*tdata->dwdx[1];
  J->data[6+6*10] = tdata->w[1];
  J->data[6+7*10] = -1.0;
  J->data[7+6*10] = tdata->p[0]*tdata->p[1];
  J->data[7+7*10] = -tdata->p[0];
  J->data[8+0*10] = x_tmp[8]*tdata->dwdx[1];
  J->data[8+8*10] = tdata->w[1];
  J->data[8+9*10] = -1.0;
  J->data[9+8*10] = tdata->p[0]*tdata->p[1];
  J->data[9+9*10] = -tdata->p[0];
for(ix = 0; ix<100; ix++) {
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


