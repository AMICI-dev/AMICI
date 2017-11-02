
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_robertson_dwdx.h"
#include "model_robertson_w.h"

using namespace amici;

int JSparse_model_robertson(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
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
J->indexvals[3] = 0;
J->indexvals[4] = 1;
J->indexvals[5] = 2;
J->indexvals[6] = 0;
J->indexvals[7] = 1;
J->indexvals[8] = 2;
J->indexptrs[0] = 0;
J->indexptrs[1] = 3;
J->indexptrs[2] = 6;
J->indexptrs[3] = 9;
status = w_model_robertson(t,x,dx,tdata);
status = dwdx_model_robertson(t,x,dx,user_data);
  J->data[0] = -cj-tdata->p[0];
  J->data[1] = tdata->p[0];
  J->data[2] = 1.0;
  J->data[3] = tdata->dwdx[0];
  J->data[4] = -cj-tdata->dwdx[0]-tdata->p[2]*x_tmp[1]*2.0;
  J->data[5] = 1.0;
  J->data[6] = tdata->dwdx[1];
  J->data[7] = -tdata->dwdx[1];
  J->data[8] = 1.0;
for(inz = 0; inz<9; inz++) {
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


