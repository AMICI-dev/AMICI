
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_dwdx.h"
#include "model_neuron_o2_w.h"

using namespace amici;

int JSparse_model_neuron_o2(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
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
int inz;
SparseSetMatToZero(J);
J->indexvals[0] = 0;
J->indexvals[1] = 1;
J->indexvals[2] = 2;
J->indexvals[3] = 3;
J->indexvals[4] = 4;
J->indexvals[5] = 5;
J->indexvals[6] = 6;
J->indexvals[7] = 8;
J->indexvals[8] = 0;
J->indexvals[9] = 1;
J->indexvals[10] = 3;
J->indexvals[11] = 2;
J->indexvals[12] = 3;
J->indexvals[13] = 2;
J->indexvals[14] = 3;
J->indexvals[15] = 4;
J->indexvals[16] = 5;
J->indexvals[17] = 4;
J->indexvals[18] = 5;
J->indexvals[19] = 6;
J->indexvals[20] = 7;
J->indexvals[21] = 6;
J->indexvals[22] = 7;
J->indexvals[23] = 8;
J->indexvals[24] = 9;
J->indexvals[25] = 8;
J->indexvals[26] = 9;
J->indexptrs[0] = 0;
J->indexptrs[1] = 8;
J->indexptrs[2] = 11;
J->indexptrs[3] = 13;
J->indexptrs[4] = 15;
J->indexptrs[5] = 17;
J->indexptrs[6] = 19;
J->indexptrs[7] = 21;
J->indexptrs[8] = 23;
J->indexptrs[9] = 25;
J->indexptrs[10] = 27;
status = w_model_neuron_o2(t,x,NULL,tdata);
status = dwdx_model_neuron_o2(t,x,NULL,user_data);
  J->data[0] = x_tmp[0]*(2.0/2.5E1)+5.0;
  J->data[1] = tdata->p[0]*tdata->p[1];
  J->data[2] = x_tmp[2]*tdata->dwdx[1];
  J->data[3] = tdata->p[1];
  J->data[4] = x_tmp[4]*tdata->dwdx[1];
  J->data[5] = tdata->p[0];
  J->data[6] = x_tmp[6]*tdata->dwdx[1];
  J->data[7] = x_tmp[8]*tdata->dwdx[1];
  J->data[8] = -1.0;
  J->data[9] = -tdata->p[0];
  J->data[10] = -1.0;
  J->data[11] = tdata->w[1];
  J->data[12] = tdata->p[0]*tdata->p[1];
  J->data[13] = -1.0;
  J->data[14] = -tdata->p[0];
  J->data[15] = tdata->w[1];
  J->data[16] = tdata->p[0]*tdata->p[1];
  J->data[17] = -1.0;
  J->data[18] = -tdata->p[0];
  J->data[19] = tdata->w[1];
  J->data[20] = tdata->p[0]*tdata->p[1];
  J->data[21] = -1.0;
  J->data[22] = -tdata->p[0];
  J->data[23] = tdata->w[1];
  J->data[24] = tdata->p[0]*tdata->p[1];
  J->data[25] = -1.0;
  J->data[26] = -tdata->p[0];
for(inz = 0; inz<27; inz++) {
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


