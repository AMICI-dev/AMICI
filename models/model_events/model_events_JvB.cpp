
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

int JvB_model_events(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
realtype *vB_tmp = N_VGetArrayPointer(vB);
realtype *JvB_tmp = N_VGetArrayPointer(JvB);
memset(JvB_tmp,0,sizeof(realtype)*3);
status = w_model_events(t,x,NULL,tdata);
  JvB_tmp[0] = -tdata->p[1]*vB_tmp[1]*exp(t*(-1.0/1.0E1))+tdata->h[3]*tdata->p[0]*vB_tmp[0];
  JvB_tmp[1] = tdata->p[2]*vB_tmp[1];
  JvB_tmp[2] = vB_tmp[2];
return(status);

}


