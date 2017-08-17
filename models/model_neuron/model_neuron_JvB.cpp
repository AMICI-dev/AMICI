
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_w.h"

int JvB_model_neuron(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
realtype *vB_tmp = N_VGetArrayPointer(vB);
realtype *JvB_tmp = N_VGetArrayPointer(JvB);
memset(JvB_tmp,0,sizeof(realtype)*2);
status = w_model_neuron(t,x,NULL,tdata);
  JvB_tmp[0] = -vB_tmp[0]*(x_tmp[0]*(2.0/2.5E1)+5.0)-udata->p[0]*udata->p[1]*vB_tmp[1];
  JvB_tmp[1] = vB_tmp[0]+udata->p[0]*vB_tmp[1];
return(status);

}


