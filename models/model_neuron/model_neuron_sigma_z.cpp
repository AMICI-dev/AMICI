
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include "model_neuron_w.h"

int sigma_z_model_neuron(realtype t, int ie, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
memset(tdata->sigmaz,0,sizeof(realtype)*1);
  tdata->sigmaz[0] = 1.0;
return(status);

}


