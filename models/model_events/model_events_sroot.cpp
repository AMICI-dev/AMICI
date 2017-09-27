
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_events_w.h"

int sroot_model_events(realtype t, int ie, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
status = w_model_events(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((0) + ip*2)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((1) + ip*2)] = sx_tmp[0]-sx_tmp[2];

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((0) + ip*2)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((1) + ip*2)] = sx_tmp[0]-sx_tmp[2];

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((0) + ip*2)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((1) + ip*2)] = sx_tmp[0]-sx_tmp[2];

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((0) + ip*2)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  sroot[tdata->nroots[ie] + udata->nmaxevent*((1) + ip*2)] = sx_tmp[0]-sx_tmp[2];

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

}
}
return(status);

}


