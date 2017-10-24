
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_events_w.h"

using namespace amici;

int z_model_events(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,tdata);
    switch(ie) { 
        case 0: {
  rdata->z[tdata->nroots[ie]+udata->nmaxevent*0] = t;

        } break;

        case 1: {
  rdata->z[tdata->nroots[ie]+udata->nmaxevent*1] = t;

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 
return(status);

}


