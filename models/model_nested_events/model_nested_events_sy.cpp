
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_nested_events_w.h"

int sy_model_nested_events(realtype t, int it, realtype *sy, realtype *dydx, realtype *dydp, N_Vector *sx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *sx_tmp;
int ip;
for(ip = 0; ip<nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (plist[ip]) {
  case 0: {
  sy[it+nt*((0)+ip*1)] = dydx[0]*sx_tmp[0];

  } break;

  case 1: {
  sy[it+nt*((0)+ip*1)] = dydx[0]*sx_tmp[0];

  } break;

  case 2: {
  sy[it+nt*((0)+ip*1)] = dydx[0]*sx_tmp[0];

  } break;

  case 3: {
  sy[it+nt*((0)+ip*1)] = dydx[0]*sx_tmp[0];

  } break;

  case 4: {
  sy[it+nt*((0)+ip*1)] = dydx[0]*sx_tmp[0];

  } break;

}
}
return(status);

}


