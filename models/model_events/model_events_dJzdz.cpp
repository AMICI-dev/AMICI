
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_events_w.h"

int dJzdz_model_events(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJzdz,0,sizeof(realtype)*model->nz*model->nztrue*model->nJ);
status = w_model_events(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJzdz[iz+(0+0*1)*model->nztrue] = 1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  tdata->dJzdz[iz+(0+1*1)*model->nztrue] = 1.0/(tdata->sigmaz[1]*tdata->sigmaz[1])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*1]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*1]*2.0)*-5.0E-1;
}
return(status);

}


