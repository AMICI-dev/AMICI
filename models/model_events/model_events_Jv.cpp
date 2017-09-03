
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

int Jv_model_events(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *v_tmp = N_VGetArrayPointer(v);
realtype *Jv_tmp = N_VGetArrayPointer(Jv);
memset(Jv_tmp,0,sizeof(realtype)*3);
status = w_model_events(t,x,NULL,tdata);
  Jv_tmp[0] = -tdata->h[3]*tdata->p[0]*v_tmp[0];
  Jv_tmp[1] = -tdata->p[2]*v_tmp[1]+tdata->p[1]*v_tmp[0]*exp(t*(-1.0/1.0E1));
  Jv_tmp[2] = -v_tmp[2];
return(status);

}


