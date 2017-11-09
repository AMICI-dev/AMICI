
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_events_w.h"

using namespace amici;

void dJrzdsigma_model_events(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJrzdsigma,0,sizeof(realtype)*model->nztrue*model->nz*model->nJ);
w_model_events(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJrzdsigma[iz+(0+0*1)*model->nztrue] = (rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0])*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-1.0+1.0/tdata->sigmaz[0];
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  tdata->dJrzdsigma[iz+(0+1*1)*model->nztrue] = (rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1])*1.0/(tdata->sigmaz[1]*tdata->sigmaz[1]*tdata->sigmaz[1])*-1.0+1.0/tdata->sigmaz[1];
}
return;

}


