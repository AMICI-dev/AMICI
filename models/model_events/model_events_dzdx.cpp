
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

using namespace amici;

void dzdx_model_events(realtype t, int ie, N_Vector x, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
w_model_events(t,x,NULL,tdata);
  tdata->dzdx[0+1*2] = 1.0/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));
  tdata->dzdx[0+2*2] = -1.0/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));
  tdata->dzdx[1+0*2] = 1.0/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);
  tdata->dzdx[1+2*2] = -1.0/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);
return;

}


