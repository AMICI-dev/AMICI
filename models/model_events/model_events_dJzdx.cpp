
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_events_w.h"

int dJzdx_model_events(realtype t, int ie, realtype *dJzdx, realtype *z, N_Vector x, realtype *dzdx, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,user_data);
int iz;
if(!amiIsNaN(mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  dJzdx[tdata->nroots[ie]+(0+1*1)*udata->nmaxevent] += -dzdx[2]*1.0/(sigma_z[0]*sigma_z[0])*(mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0);
  dJzdx[tdata->nroots[ie]+(0+2*1)*udata->nmaxevent] += -dzdx[4]*1.0/(sigma_z[0]*sigma_z[0])*(mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0);
}
if(!amiIsNaN(mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  dJzdx[tdata->nroots[ie]+(0+0*1)*udata->nmaxevent] += -dzdx[1]*1.0/(sigma_z[1]*sigma_z[1])*(mz[tdata->nroots[ie]+udata->nmaxevent*1]*2.0-z[tdata->nroots[ie]+udata->nmaxevent*1]*2.0);
  dJzdx[tdata->nroots[ie]+(0+2*1)*udata->nmaxevent] += -dzdx[5]*1.0/(sigma_z[1]*sigma_z[1])*(mz[tdata->nroots[ie]+udata->nmaxevent*1]*2.0-z[tdata->nroots[ie]+udata->nmaxevent*1]*2.0);
}
return(status);

}


