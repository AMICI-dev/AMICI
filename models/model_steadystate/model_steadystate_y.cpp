
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int y_model_steadystate(realtype t, int it, realtype *y, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,user_data);
  y[it+nt*(0)] = x_tmp[0];
  y[it+nt*(1)] = x_tmp[1];
  y[it+nt*(2)] = x_tmp[2];
return(status);

}


