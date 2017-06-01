
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

int dzdx_model_steadystate(realtype t, int ie, realtype *dzdx, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,user_data);
return(status);

}


