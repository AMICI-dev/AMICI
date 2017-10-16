
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_jakstat_adjoint_w.h"

int y_model_jakstat_adjoint(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
  rdata->y[it + udata->nt*0] = tdata->p[11]+(tdata->p[13]*(x_tmp[1]+x_tmp[2]*2.0))/tdata->p[4];
  rdata->y[it + udata->nt*1] = tdata->p[10]+(tdata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0))/tdata->p[4];
  rdata->y[it + udata->nt*2] = am_spline_pos(t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
return(status);

}


