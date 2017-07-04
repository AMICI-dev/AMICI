
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_steadystate_w.h"

int dJydy_model_steadystate(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydy,0,sizeof(realtype)*udata->nytrue*udata->nytrue*udata->nJ);
status = w_model_steadystate(t,x,NULL,user_data);
return(status);

}


