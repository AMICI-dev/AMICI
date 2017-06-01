
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

int dydx_model_steadystate(realtype t, int it, realtype *dydx, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,user_data);
  dydx[0+0*3] = 1.0;
  dydx[1+1*3] = 1.0;
  dydx[2+2*3] = 1.0;
return(status);

}


