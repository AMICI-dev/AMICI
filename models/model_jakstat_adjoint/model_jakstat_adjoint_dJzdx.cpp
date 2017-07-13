
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_jakstat_adjoint_w.h"

int dJzdx_model_jakstat_adjoint(realtype t, int ie, N_Vector x, realtype *z, realtype *mz, realtype *dzdx,  void *user_data, TempData *tdata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
return(status);

}


