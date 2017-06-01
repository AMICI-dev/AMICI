
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int dydp_model_jakstat_adjoint(realtype t, int it, realtype *dydp, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 4: {
  dydp[ip*3 + 0] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  dydp[ip*3 + 1] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);

  } break;

  case 5: {
  dydp[ip*3 + 2] = am_Dspline_pos(4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 6: {
  dydp[ip*3 + 2] = am_Dspline_pos(6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 7: {
  dydp[ip*3 + 2] = am_Dspline_pos(8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 8: {
  dydp[ip*3 + 2] = am_Dspline_pos(10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 9: {
  dydp[ip*3 + 2] = am_Dspline_pos(12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 10: {
  dydp[ip*3 + 1] = 1.0;

  } break;

  case 11: {
  dydp[ip*3 + 0] = 1.0;

  } break;

  case 12: {
  dydp[ip*3 + 1] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/udata->p[4];

  } break;

  case 13: {
  dydp[ip*3 + 0] = (x_tmp[1]+x_tmp[2]*2.0)/udata->p[4];

  } break;

}
}
return(status);

}


