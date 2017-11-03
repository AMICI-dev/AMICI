
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_neuron_w.h"

using namespace amici;

void sz_model_neuron(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
w_model_neuron(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

}
}
return;

}


