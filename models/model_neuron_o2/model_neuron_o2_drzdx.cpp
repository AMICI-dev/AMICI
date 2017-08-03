
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_neuron_o2_w.h"

int drzdx_model_neuron_o2(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_neuron_o2(t,x,NULL,user_data);
  tdata->drzdx[0+0*5] = 1.0;
  tdata->drzdx[1+2*5] = 1.0;
  tdata->drzdx[2+4*5] = 1.0;
  tdata->drzdx[3+6*5] = 1.0;
  tdata->drzdx[4+8*5] = 1.0;
return(status);

}


