
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_w.h"

int qBdot_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
int ip;
memset(qBdot_tmp,0,sizeof(realtype)*nplist*ng);
status = dwdp_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
  case 0: {
  qBdot_tmp[ip+nplist*0] = w_tmp[0]*xB_tmp[0]*x_tmp[0]-w_tmp[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 1: {
  qBdot_tmp[ip+nplist*0] = w_tmp[1]*xB_tmp[1]*2.0-w_tmp[1]*xB_tmp[2];

  } break;

  case 2: {
  qBdot_tmp[ip+nplist*0] = xB_tmp[2]*x_tmp[2]-(k[0]*xB_tmp[3]*x_tmp[2])/k[1];

  } break;

  case 3: {
  qBdot_tmp[ip+nplist*0] = xB_tmp[3]*x_tmp[3]-xB_tmp[5]*(x_tmp[4]-x_tmp[5])-xB_tmp[6]*(x_tmp[5]-x_tmp[6])-xB_tmp[7]*(x_tmp[6]-x_tmp[7])-xB_tmp[8]*(x_tmp[7]-x_tmp[8])-xB_tmp[4]*(x_tmp[3]*2.0-x_tmp[4])-(k[1]*xB_tmp[0]*x_tmp[8])/k[0];

  } break;

  case 5: {
  qBdot_tmp[ip+nplist*0] = dwdp_tmp[0]*p[0]*xB_tmp[0]*x_tmp[0]-dwdp_tmp[0]*p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 6: {
  qBdot_tmp[ip+nplist*0] = dwdp_tmp[1]*p[0]*xB_tmp[0]*x_tmp[0]-dwdp_tmp[1]*p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 7: {
  qBdot_tmp[ip+nplist*0] = dwdp_tmp[2]*p[0]*xB_tmp[0]*x_tmp[0]-dwdp_tmp[2]*p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 8: {
  qBdot_tmp[ip+nplist*0] = dwdp_tmp[3]*p[0]*xB_tmp[0]*x_tmp[0]-dwdp_tmp[3]*p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 9: {
  qBdot_tmp[ip+nplist*0] = dwdp_tmp[4]*p[0]*xB_tmp[0]*x_tmp[0]-dwdp_tmp[4]*p[0]*xB_tmp[1]*x_tmp[0];

  } break;

}
}
for(ip = 0; ip<nplist*ng; ip++) {
   if(amiIsNaN(qBdot_tmp[ip])) {
       qBdot_tmp[ip] = 0;       if(!udata->am_nan_qBdot) {
           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->am_nan_qBdot = TRUE;
       }
   }   if(amiIsInf(qBdot_tmp[ip])) {
       warnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


