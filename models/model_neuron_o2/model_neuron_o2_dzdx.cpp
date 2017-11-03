
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

using namespace amici;

void dzdx_model_neuron_o2(realtype t, int ie, N_Vector x, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
w_model_neuron_o2(t,x,NULL,tdata);
  tdata->dzdx[0+0*5] = -1.0/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->dzdx[1+0*5] = x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[1+1*5] = -x_tmp[2]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[1+2*5] = -1.0/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->dzdx[2+0*5] = x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[2+1*5] = -x_tmp[4]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[2+4*5] = -1.0/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->dzdx[3+0*5] = x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[3+1*5] = -x_tmp[6]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[3+6*5] = -1.0/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->dzdx[4+0*5] = x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[4+1*5] = -x_tmp[8]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0);
  tdata->dzdx[4+8*5] = -1.0/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
return;

}


