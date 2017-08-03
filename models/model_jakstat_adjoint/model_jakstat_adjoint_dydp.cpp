
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_jakstat_adjoint_w.h"

int dydp_model_jakstat_adjoint(realtype t, int it, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 4: {
  tdata->dydp[ip*udata->ny + 0] = -1.0/(udata->p[4]*udata->p[4])*udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  tdata->dydp[ip*udata->ny + 1] = -1.0/(udata->p[4]*udata->p[4])*udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);

  } break;

  case 5: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 6: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 7: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 8: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 9: {
  tdata->dydp[ip*udata->ny + 2] = am_Dspline_pos(12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);

  } break;

  case 10: {
  tdata->dydp[ip*udata->ny + 1] = 1.0;

  } break;

  case 11: {
  tdata->dydp[ip*udata->ny + 0] = 1.0;

  } break;

  case 12: {
  tdata->dydp[ip*udata->ny + 1] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/udata->p[4];

  } break;

  case 13: {
  tdata->dydp[ip*udata->ny + 0] = (x_tmp[1]+x_tmp[2]*2.0)/udata->p[4];

  } break;

}
}
return(status);

}


