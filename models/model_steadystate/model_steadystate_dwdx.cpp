
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

int dwdx_model_steadystate(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(udata->dwdx,0,sizeof(realtype)*2);
status = w_model_steadystate(t,x,NULL,user_data);
  udata->dwdx[0] = x_tmp[0]*2.0;
  udata->dwdx[1] = udata->p[3];
return(status);

}


