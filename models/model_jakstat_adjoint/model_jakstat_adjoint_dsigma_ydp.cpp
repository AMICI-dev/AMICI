
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int dsigma_ydp_model_jakstat_adjoint(realtype t, realtype *dsigma_ydp, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
int ip;
memset(dsigma_ydp,0,sizeof(realtype)*3*udata->nplist);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 14: {
  dsigma_ydp[ip*3 + 0] = 1.0;

  } break;

  case 15: {
  dsigma_ydp[ip*3 + 1] = 1.0;

  } break;

  case 16: {
  dsigma_ydp[ip*3 + 2] = 1.0;

  } break;

}
}
return(status);

}


