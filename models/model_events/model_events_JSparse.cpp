
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_dwdx.h"
#include "model_events_w.h"

int JSparse_model_events(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
int inz;
SparseSetMatToZero(J);
J->indexvals[0] = 0;
J->indexvals[1] = 1;
J->indexvals[2] = 1;
J->indexvals[3] = 2;
J->indexptrs[0] = 0;
J->indexptrs[1] = 2;
J->indexptrs[2] = 3;
J->indexptrs[3] = 4;
status = w_model_events(t,x,NULL,tdata);
status = dwdx_model_events(t,x,NULL,user_data);
  J->data[0] = -tdata->h[3]*tdata->p[0];
  J->data[1] = tdata->p[1]*exp(t*(-1.0/1.0E1));
  J->data[2] = -tdata->p[2];
  J->data[3] = -1.0;
for(inz = 0; inz<4; inz++) {
   if(amiIsNaN(J->data[inz])) {
       J->data[inz] = 0;
       if(!tdata->nan_JSparse) {
           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");
           tdata->nan_JSparse = TRUE;
       }
   }
   if(amiIsInf(J->data[inz])) {
       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


