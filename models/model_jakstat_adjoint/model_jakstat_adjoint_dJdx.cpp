
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdx.h"
#include "model_jakstat_adjoint_w.h"

int dJdx_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
int ix;
memset(tdata->dJdx,0,sizeof(realtype)*729);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
status = dwdx_model_jakstat_adjoint(t,x,NULL,user_data);
  tdata->dJdx[1*model->nx*model->nx + 1*model->nx + 1] = tdata->p[1]*-4.0;
  tdata->dJdx[1*model->nx*model->nx + 1*model->nx + 2] = tdata->p[1]*2.0;
return(status);

}


