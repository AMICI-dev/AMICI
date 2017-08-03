
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int root_model_jakstat_adjoint(realtype t, N_Vector x, realtype *root, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
return(status);

}


