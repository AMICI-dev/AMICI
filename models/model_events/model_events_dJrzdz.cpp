
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_events_w.h"

int dJrzdz_model_events(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJrzdz,0,sizeof(realtype)*model->nz*model->nztrue*model->nJ);
status = w_model_events(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJrzdz[iz+(0+0*1)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(0+1*1)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
}
return(status);

}


