
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int sigma_y_model_jakstat_adjoint(realtype t, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
memset(tdata->sigmay,0,sizeof(realtype)*3);
  tdata->sigmay[0] = tdata->p[14];
  tdata->sigmay[1] = tdata->p[15];
  tdata->sigmay[2] = tdata->p[16];
return(status);

}


