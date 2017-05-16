
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int root_model_jakstat_adjoint_o2(realtype t, N_Vector x, realtype *root, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
return(status);

}


