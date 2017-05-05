
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int dJydp_model_steadystate(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dJydp,0,sizeof(realtype)*nytrue*nplist*ng);
status = w_model_steadystate(t,x,NULL,user_data);
return(status);

}


