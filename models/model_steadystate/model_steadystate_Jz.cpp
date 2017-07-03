
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_steadystate_w.h"

int Jz_model_steadystate(realtype t, int ie, realtype *Jz, realtype *z, N_Vector x, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,user_data);
return(status);

}


