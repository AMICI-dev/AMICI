
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_events_w.h"

using namespace amici;

int rz_model_events(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,tdata);
  rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0] = x_tmp[0]-x_tmp[2];
return(status);

}


