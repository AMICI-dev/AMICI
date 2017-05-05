
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int sy_model_steadystate(realtype t, int it, realtype *sy, realtype *dydx, realtype *dydp, N_Vector *sx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *sx_tmp;
int ip;
for(ip = 0; ip<nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
switch (plist[ip]) {
  case 0: {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
  sy[it+nt*((0)+ip*3)] = dydx[0]*sx_tmp[0];
  sy[it+nt*((1)+ip*3)] = dydx[4]*sx_tmp[1];
  sy[it+nt*((2)+ip*3)] = dydx[8]*sx_tmp[2];

  } break;

  case 1: {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
  sy[it+nt*((0)+ip*3)] = dydx[0]*sx_tmp[0];
  sy[it+nt*((1)+ip*3)] = dydx[4]*sx_tmp[1];
  sy[it+nt*((2)+ip*3)] = dydx[8]*sx_tmp[2];

  } break;

  case 2: {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
  sy[it+nt*((0)+ip*3)] = dydx[0]*sx_tmp[0];
  sy[it+nt*((1)+ip*3)] = dydx[4]*sx_tmp[1];
  sy[it+nt*((2)+ip*3)] = dydx[8]*sx_tmp[2];

  } break;

  case 3: {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
  sy[it+nt*((0)+ip*3)] = dydx[0]*sx_tmp[0];
  sy[it+nt*((1)+ip*3)] = dydx[4]*sx_tmp[1];
  sy[it+nt*((2)+ip*3)] = dydx[8]*sx_tmp[2];

  } break;

  case 4: {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
  sy[it+nt*((0)+ip*3)] = dydx[0]*sx_tmp[0];
  sy[it+nt*((1)+ip*3)] = dydx[4]*sx_tmp[1];
  sy[it+nt*((2)+ip*3)] = dydx[8]*sx_tmp[2];

  } break;

}
}
return(status);

}


