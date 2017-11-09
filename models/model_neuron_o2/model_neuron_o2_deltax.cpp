
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_o2_w.h"

using namespace amici;

void deltax_model_neuron_o2(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = nullptr;
if(xdot_old)
    xdot_old_tmp = N_VGetArrayPointer(xdot_old);
memset(tdata->deltax,0,sizeof(realtype)*10);
w_model_neuron_o2(t,x,NULL,tdata);
              switch(ie) { 
              case 0: {
  tdata->deltax[0] = -tdata->p[2]-x_tmp[0];
  tdata->deltax[1] = tdata->p[3];
  tdata->deltax[2] = -(x_tmp[2]*(tdata->p[2]*5.0+tdata->p[3]+x_tmp[0]*5.0-(tdata->p[2]*tdata->p[2])*(1.0/2.5E1)+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->deltax[3] = -(x_tmp[2]*(tdata->p[0]*(tdata->p[3]+x_tmp[1]+tdata->p[1]*tdata->p[2])-tdata->p[0]*(x_tmp[1]-tdata->p[1]*x_tmp[0])))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->deltax[4] = -(x_tmp[4]*(tdata->p[2]*5.0+tdata->p[3]+x_tmp[0]*5.0-(tdata->p[2]*tdata->p[2])*(1.0/2.5E1)+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->deltax[5] = -(x_tmp[4]*(tdata->p[0]*(tdata->p[3]+x_tmp[1]+tdata->p[1]*tdata->p[2])-tdata->p[0]*(x_tmp[1]-tdata->p[1]*x_tmp[0])))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->deltax[6] = -(x_tmp[6]*(tdata->p[2]*5.0+tdata->p[3]+x_tmp[0]*5.0-(tdata->p[2]*tdata->p[2])*(1.0/2.5E1)+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-1.0;
  tdata->deltax[7] = -(x_tmp[6]*(tdata->p[0]*(tdata->p[3]+x_tmp[1]+tdata->p[1]*tdata->p[2])-tdata->p[0]*(x_tmp[1]-tdata->p[1]*x_tmp[0])))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->deltax[8] = -(x_tmp[8]*(tdata->p[2]*5.0+tdata->p[3]+x_tmp[0]*5.0-(tdata->p[2]*tdata->p[2])*(1.0/2.5E1)+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  tdata->deltax[9] = -(x_tmp[8]*(tdata->p[0]*(tdata->p[3]+x_tmp[1]+tdata->p[1]*tdata->p[2])-tdata->p[0]*(x_tmp[1]-tdata->p[1]*x_tmp[0])))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+1.0;

              } break;

              } 
return;

}


