
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

using namespace amici;

int Jv_model_steadystate(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
realtype *v_tmp = nullptr;
if(v)
    v_tmp = N_VGetArrayPointer(v);
realtype *Jv_tmp = nullptr;
if(Jv)
    Jv_tmp = N_VGetArrayPointer(Jv);
memset(Jv_tmp,0,sizeof(realtype)*3);
status = w_model_steadystate(t,x,NULL,tdata);
  Jv_tmp[0] = v_tmp[2]*tdata->dwdx[1]+v_tmp[1]*(tdata->p[2]*2.0-tdata->p[1]*x_tmp[0])-v_tmp[0]*(tdata->p[1]*x_tmp[1]+tdata->p[0]*tdata->dwdx[0]*2.0);
  Jv_tmp[1] = v_tmp[2]*tdata->dwdx[1]-v_tmp[0]*(tdata->p[1]*x_tmp[1]-tdata->p[0]*tdata->dwdx[0])-v_tmp[1]*(tdata->p[2]+tdata->p[1]*x_tmp[0]);
  Jv_tmp[2] = -v_tmp[2]*(udata->k[3]+tdata->dwdx[1])+tdata->p[1]*v_tmp[0]*x_tmp[1]+tdata->p[1]*v_tmp[1]*x_tmp[0];
return(status);

}


