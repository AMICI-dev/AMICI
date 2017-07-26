
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_steadystate_w.h"

int dydx_model_steadystate(realtype t, int it, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,user_data);
  tdata->dydx[0+0*3] = 1.0;
  tdata->dydx[1+1*3] = 1.0;
  tdata->dydx[2+2*3] = 1.0;
return(status);

}


