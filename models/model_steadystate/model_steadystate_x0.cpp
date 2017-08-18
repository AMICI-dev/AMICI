
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

int x0_model_steadystate(N_Vector x0, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x0_tmp = N_VGetArrayPointer(x0);
memset(x0_tmp,0,sizeof(realtype)*3);
realtype t = udata->tstart;
  x0_tmp[0] = udata->k[0];
  x0_tmp[1] = udata->k[1];
  x0_tmp[2] = udata->k[2];
return(status);

}


