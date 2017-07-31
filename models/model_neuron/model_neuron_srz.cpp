
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_neuron_w.h"

int srz_model_neuron(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
status = w_model_neuron(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*udata->nz + 0)] = sx_tmp[0];

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*udata->nz + 0)] = sx_tmp[0];

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*udata->nz + 0)] = sx_tmp[0];

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*udata->nz + 0)] = sx_tmp[0];

        } break;

    } 

  } break;

}
}
return(status);

}


