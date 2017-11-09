
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_neuron_w.h"

using namespace amici;

void deltaxB_model_neuron(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = nullptr;
if(xB)
    xB_tmp = N_VGetArrayPointer(xB);
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = nullptr;
if(xdot_old)
    xdot_old_tmp = N_VGetArrayPointer(xdot_old);
memset(tdata->deltaxB,0,sizeof(realtype)*2);
w_model_neuron(t,x,NULL,tdata);
              switch(ie) { 
              case 0: {
  tdata->deltaxB[0] = xB_tmp[0];

              } break;

              } 
return;

}


