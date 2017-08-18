
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

int xBdot_model_neuron_o2(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
int ix;
memset(xBdot_tmp,0,sizeof(realtype)*10);
status = w_model_neuron_o2(t,x,NULL,tdata);
status = dwdx_model_neuron_o2(t,x,NULL,user_data);
  xBdot_tmp[0] = -xB_tmp[0]*(x_tmp[0]*(2.0/2.5E1)+5.0)-udata->p[0]*udata->p[1]*xB_tmp[1];
  xBdot_tmp[1] = xB_tmp[0]+udata->p[0]*xB_tmp[1];
  xBdot_tmp[2] = -udata->p[1]*xB_tmp[1]-tdata->w[1]*xB_tmp[2]-udata->p[0]*udata->p[1]*xB_tmp[3]-x_tmp[2]*xB_tmp[0]*tdata->dwdx[1];
  xBdot_tmp[3] = xB_tmp[1]+xB_tmp[2]+udata->p[0]*xB_tmp[3];
  xBdot_tmp[4] = -udata->p[0]*xB_tmp[1]-tdata->w[1]*xB_tmp[4]-udata->p[0]*udata->p[1]*xB_tmp[5]-x_tmp[4]*xB_tmp[0]*tdata->dwdx[1];
  xBdot_tmp[5] = xB_tmp[4]+udata->p[0]*xB_tmp[5];
  xBdot_tmp[6] = -tdata->w[1]*xB_tmp[6]-udata->p[0]*udata->p[1]*xB_tmp[7]-x_tmp[6]*xB_tmp[0]*tdata->dwdx[1];
  xBdot_tmp[7] = xB_tmp[6]+udata->p[0]*xB_tmp[7];
  xBdot_tmp[8] = -tdata->w[1]*xB_tmp[8]-udata->p[0]*udata->p[1]*xB_tmp[9]-x_tmp[8]*xB_tmp[0]*tdata->dwdx[1];
  xBdot_tmp[9] = xB_tmp[8]+udata->p[0]*xB_tmp[9];
for(ix = 0; ix<10; ix++) {
   if(amiIsNaN(xBdot_tmp[ix])) {
       xBdot_tmp[ix] = 0;       if(!tdata->nan_xBdot) {
           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_xBdot = TRUE;
       }
   }   if(amiIsInf(xBdot_tmp[ix])) {
       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


