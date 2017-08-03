
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_dirac_w.h"

int y_model_dirac(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_dirac(t,x,NULL,user_data);
  rdata->y[it + udata->nt*0] = x_tmp[1];
return(status);

}


