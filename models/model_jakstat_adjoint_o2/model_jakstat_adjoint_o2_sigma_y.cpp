
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

void sigma_y_model_jakstat_adjoint_o2(realtype t, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
memset(tdata->sigmay,0,sizeof(realtype)*54);
  tdata->sigmay[0] = tdata->p[14];
  tdata->sigmay[1] = tdata->p[15];
  tdata->sigmay[2] = tdata->p[16];
  tdata->sigmay[17] = 1.0;
  tdata->sigmay[35] = 1.0;
  tdata->sigmay[53] = 1.0;
return;

}


