
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_events_w.h"

int dJzdp_model_events(realtype t, int ie, realtype *dJzdp, realtype *z, N_Vector x, realtype *dzdp, realtype *mz, realtype *sigma_z, realtype *dsigma_zdp, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,user_data);
return(status);

}


