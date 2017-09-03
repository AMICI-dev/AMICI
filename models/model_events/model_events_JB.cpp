
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_dwdx.h"
#include "model_events_w.h"

int JB_model_events(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(JB->data,0,sizeof(realtype)*9);
status = w_model_events(t,x,NULL,tdata);
status = dwdx_model_events(t,x,NULL,user_data);
  JB->data[0+0*3] = tdata->h[3]*tdata->p[0];
  JB->data[0+1*3] = -tdata->p[1]*exp(t*(-1.0/1.0E1));
  JB->data[1+1*3] = tdata->p[2];
  JB->data[2+2*3] = 1.0;
return(status);

}


