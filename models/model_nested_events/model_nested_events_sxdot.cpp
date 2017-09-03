
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_nested_events_JSparse.h"
#include "model_nested_events_dxdotdp.h"
#include "model_nested_events_w.h"

int sxdot_model_nested_events(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp = N_VGetArrayPointer(sx);
realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*1);
if(ip == 0) {
    status = JSparse_model_nested_events(t,x,xdot,tdata->J,user_data,NULL,NULL,NULL);
    status = dxdotdp_model_nested_events(t,x,NULL,user_data);
}
  sxdot_tmp[0] = tdata->dxdotdp[0 + ip*model->nx]+tdata->J->data[0]*sx_tmp[0];
return(status);

}


