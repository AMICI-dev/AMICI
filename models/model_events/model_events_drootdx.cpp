
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_events_w.h"

int drootdx_model_events(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,user_data);
  tdata->drzdx[0+1*4] = 1.0;
  tdata->drzdx[0+2*4] = -1.0;
  tdata->drzdx[1+0*4] = 1.0;
  tdata->drzdx[1+2*4] = -1.0;
return(status);

}


