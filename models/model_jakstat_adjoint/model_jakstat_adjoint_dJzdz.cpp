
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_jakstat_adjoint_w.h"

int dJzdz_model_jakstat_adjoint(realtype t, int ie, N_Vector x, realtype *z, realtype *mz, void *user_data, TempData *tdata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJzdz,0,sizeof(realtype)*udata->nz*udata->nztrue*udata->nJ);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
return(status);

}


