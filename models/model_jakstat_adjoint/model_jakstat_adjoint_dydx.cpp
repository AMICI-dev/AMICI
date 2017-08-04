
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_jakstat_adjoint_w.h"

int dydx_model_jakstat_adjoint(realtype t, int it, N_Vector x, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  tdata->dydx[0+1*3] = udata->p[13]/udata->p[4];
  tdata->dydx[0+2*3] = (udata->p[13]*2.0)/udata->p[4];
  tdata->dydx[1+0*3] = udata->p[12]/udata->p[4];
  tdata->dydx[1+1*3] = udata->p[12]/udata->p[4];
  tdata->dydx[1+2*3] = (udata->p[12]*2.0)/udata->p[4];
return(status);

}


