
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_dwdx.h"
#include "model_steadystate_w.h"

int JSparse_model_steadystate(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
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
status = w_model_steadystate(t,x,NULL,user_data);
status = dwdx_model_steadystate(t,x,NULL,user_data);
  J->data[0] = dwdx_tmp[0]*p[0]*-2.0-p[1]*x_tmp[1];
  J->data[1] = dwdx_tmp[0]*p[0]-p[1]*x_tmp[1];
  J->data[2] = p[1]*x_tmp[1];
  J->data[3] = p[2]*2.0-p[1]*x_tmp[0];
  J->data[4] = -p[2]-p[1]*x_tmp[0];
  J->data[5] = p[1]*x_tmp[0];
  J->data[6] = dwdx_tmp[1];
  J->data[7] = dwdx_tmp[1];
  J->data[8] = -dwdx_tmp[1]-k[3];
for(inz = 0; inz<9; inz++) {
   if(amiIsNaN(J->data[inz])) {
       J->data[inz] = 0;
       if(!udata->am_nan_JSparse) {
           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_JSparse = TRUE;
       }
   }
   if(amiIsInf(J->data[inz])) {
       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");
       return(-1);
   }
}
return(status);

}


