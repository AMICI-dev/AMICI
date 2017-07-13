
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_neuron_w.h"

int sigma_y_model_neuron(realtype t, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
memset(tdata->sigmay,0,sizeof(realtype)*1);
  tdata->sigmay[0] = 1.0;
return(status);

}


