
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

using namespace amici;

void sigma_z_model_jakstat_adjoint(realtype t, int ie, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
memset(tdata->sigmaz,0,sizeof(realtype)*0);
return;

}


