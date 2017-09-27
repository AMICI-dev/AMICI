
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_events_w.h"

int sJz_model_events(realtype t, int ie, realtype *sJz, realtype *s2Jz, realtype *dJzdz, realtype *dJzdp, realtype *sz, realtype *dzdp, realtype *mz, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
int ip;
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
int iz;
if(!amiIsNaN(mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*0]*sz[tdata->nroots[ie]+udata->nmaxevent*(0+ip*2)];
}
if(!amiIsNaN(mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*3]*sz[tdata->nroots[ie]+udata->nmaxevent*(1+ip*2)];
}

  } break;

  case 1: {
int iz;
if(!amiIsNaN(mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*0]*sz[tdata->nroots[ie]+udata->nmaxevent*(0+ip*2)];
}
if(!amiIsNaN(mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*3]*sz[tdata->nroots[ie]+udata->nmaxevent*(1+ip*2)];
}

  } break;

  case 2: {
int iz;
if(!amiIsNaN(mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*0]*sz[tdata->nroots[ie]+udata->nmaxevent*(0+ip*2)];
}
if(!amiIsNaN(mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*3]*sz[tdata->nroots[ie]+udata->nmaxevent*(1+ip*2)];
}

  } break;

  case 3: {
int iz;
if(!amiIsNaN(mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*0]*sz[tdata->nroots[ie]+udata->nmaxevent*(0+ip*2)];
}
if(!amiIsNaN(mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  sJz[0] += dJzdz[tdata->nroots[ie]+udata->nmaxevent*3]*sz[tdata->nroots[ie]+udata->nmaxevent*(1+ip*2)];
}

  } break;

}
}
return(status);

}


