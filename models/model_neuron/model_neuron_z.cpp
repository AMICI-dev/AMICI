
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_neuron_w.h"

using namespace amici;

void z_model_neuron(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
w_model_neuron(t,x,NULL,tdata);
    switch(ie) { 
        case 0: {
  rdata->z[tdata->nroots[ie]+udata->nmaxevent*0] = t;

        } break;

    } 
return;

}


