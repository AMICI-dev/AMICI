
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_dwdx.h"
#include "model_neuron_w.h"

using namespace amici;

void JSparse_model_neuron(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
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
int inz;
SparseSetMatToZero(J);
J->indexvals[0] = 0;
J->indexvals[1] = 1;
J->indexvals[2] = 0;
J->indexvals[3] = 1;
J->indexptrs[0] = 0;
J->indexptrs[1] = 2;
J->indexptrs[2] = 4;
w_model_neuron(t,x,NULL,tdata);
dwdx_model_neuron(t,x,NULL,user_data);
  J->data[0] = x_tmp[0]*(2.0/2.5E1)+5.0;
  J->data[1] = tdata->p[0]*tdata->p[1];
  J->data[2] = -1.0;
  J->data[3] = -tdata->p[0];
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
       return;
   }
}
return;

}


