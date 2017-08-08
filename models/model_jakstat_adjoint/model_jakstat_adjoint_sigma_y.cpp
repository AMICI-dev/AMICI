
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_jakstat_adjoint_w.h"

int sigma_y_model_jakstat_adjoint(realtype t, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
memset(tdata->sigmay,0,sizeof(realtype)*3);
  tdata->sigmay[0] = udata->p[14];
  tdata->sigmay[1] = udata->p[15];
  tdata->sigmay[2] = udata->p[16];
return(status);

}


