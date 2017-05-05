
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int w_model_steadystate(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(w_tmp,0,sizeof(realtype)*2);
  w_tmp[0] = p[3]*x_tmp[2];
  w_tmp[1] = x_tmp[0]*x_tmp[0];
return(status);

}


