
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int x0_model_steadystate(N_Vector x0, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x0_tmp = N_VGetArrayPointer(x0);
memset(x0_tmp,0,sizeof(realtype)*3);
realtype t = tstart;
  x0_tmp[0] = k[0];
  x0_tmp[1] = k[1];
  x0_tmp[2] = k[2];
return(status);

}


