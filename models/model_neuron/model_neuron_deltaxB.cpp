
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include "model_neuron_w.h"

int deltaxB_model_neuron(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
memset(tdata->deltaxB,0,sizeof(realtype)*2);
status = w_model_neuron(t,x,NULL,tdata);
              switch(ie) { 
              case 0: {
  tdata->deltaxB[0] = xB_tmp[0];

              } break;

              } 
return(status);

}


