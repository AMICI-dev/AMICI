
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_robertson_JSparse.h"
#include "model_robertson_dxdotdp.h"
#include "model_robertson_dfdx.h"
#include "model_robertson_M.h"
#include "model_robertson_w.h"

using namespace amici;

void sxdot_model_robertson(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *sx_tmp = nullptr;
if(sx)
    sx_tmp = N_VGetArrayPointer(sx);
realtype *sdx_tmp = nullptr;
if(sdx)
    sdx_tmp = N_VGetArrayPointer(sdx);
realtype *sxdot_tmp = nullptr;
if(sxdot)
    sxdot_tmp = N_VGetArrayPointer(sxdot);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
memset(sxdot_tmp,0,sizeof(realtype)*3);
if(ip == 0) {
    dfdx_model_robertson(t,x,dx,user_data);
    M_model_robertson(t,x,dx,user_data);
    dxdotdp_model_robertson(t,x,dx,user_data);
}
  sxdot_tmp[0] = tdata->dxdotdp[0 + ip*model->nx]-tdata->M[0]*sdx_tmp[0]+tdata->dfdx[0]*sx_tmp[0]+tdata->dfdx[3]*sx_tmp[1]+tdata->dfdx[6]*sx_tmp[2];
  sxdot_tmp[1] = tdata->dxdotdp[1 + ip*model->nx]-tdata->M[4]*sdx_tmp[1]+tdata->dfdx[1]*sx_tmp[0]+tdata->dfdx[4]*sx_tmp[1]+tdata->dfdx[7]*sx_tmp[2];
  sxdot_tmp[2] = tdata->dxdotdp[2 + ip*model->nx]+tdata->dfdx[2]*sx_tmp[0]+tdata->dfdx[5]*sx_tmp[1]+tdata->dfdx[8]*sx_tmp[2];
return;

}


