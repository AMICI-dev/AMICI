
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_w.h"

int qBdot_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
int ip;
memset(qBdot_tmp,0,sizeof(realtype)*udata->nplist*udata->ng);
status = dwdp_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
  qBdot_tmp[ip + udata->nplist*0] = udata->w[0]*xB_tmp[0]*x_tmp[0]-udata->w[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 1: {
  qBdot_tmp[ip + udata->nplist*0] = udata->w[1]*xB_tmp[1]*2.0-udata->w[1]*xB_tmp[2];

  } break;

  case 2: {
  qBdot_tmp[ip + udata->nplist*0] = xB_tmp[2]*x_tmp[2]-(udata->k[0]*xB_tmp[3]*x_tmp[2])/udata->k[1];

  } break;

  case 3: {
  qBdot_tmp[ip + udata->nplist*0] = xB_tmp[3]*x_tmp[3]-xB_tmp[5]*(x_tmp[4]-x_tmp[5])-xB_tmp[6]*(x_tmp[5]-x_tmp[6])-xB_tmp[7]*(x_tmp[6]-x_tmp[7])-xB_tmp[8]*(x_tmp[7]-x_tmp[8])-xB_tmp[4]*(x_tmp[3]*2.0-x_tmp[4])-(udata->k[1]*xB_tmp[0]*x_tmp[8])/udata->k[0];

  } break;

  case 5: {
  qBdot_tmp[ip + udata->nplist*0] = udata->dwdp[0]*udata->p[0]*xB_tmp[0]*x_tmp[0]-udata->dwdp[0]*udata->p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 6: {
  qBdot_tmp[ip + udata->nplist*0] = udata->dwdp[1]*udata->p[0]*xB_tmp[0]*x_tmp[0]-udata->dwdp[1]*udata->p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 7: {
  qBdot_tmp[ip + udata->nplist*0] = udata->dwdp[2]*udata->p[0]*xB_tmp[0]*x_tmp[0]-udata->dwdp[2]*udata->p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 8: {
  qBdot_tmp[ip + udata->nplist*0] = udata->dwdp[3]*udata->p[0]*xB_tmp[0]*x_tmp[0]-udata->dwdp[3]*udata->p[0]*xB_tmp[1]*x_tmp[0];

  } break;

  case 9: {
  qBdot_tmp[ip + udata->nplist*0] = udata->dwdp[4]*udata->p[0]*xB_tmp[0]*x_tmp[0]-udata->dwdp[4]*udata->p[0]*xB_tmp[1]*x_tmp[0];

  } break;

}
}
for(ip = 0; ip<udata->nplist*udata->ng; ip++) {
   if(amiIsNaN(qBdot_tmp[ip])) {
       qBdot_tmp[ip] = 0;       if(!udata->nan_qBdot) {
           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");
           udata->nan_qBdot = TRUE;
       }
   }   if(amiIsInf(qBdot_tmp[ip])) {
       warnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");
       return(-1);
   }}
return(status);

}


