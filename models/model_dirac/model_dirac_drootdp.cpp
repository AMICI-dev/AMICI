
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_dirac_w.h"

int drootdp_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_dirac(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 1: {
  tdata->drzdp[ip*udata->nz + 0] = 1.0;
  tdata->drzdp[ip*udata->nz + 1] = -1.0;

  } break;

}
}
return(status);

}


