
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_nested_events_w.h"

int drootdx_model_nested_events(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_nested_events(t,x,NULL,user_data);
  tdata->drzdx[0+0*3] = 1.0;
return(status);

}


