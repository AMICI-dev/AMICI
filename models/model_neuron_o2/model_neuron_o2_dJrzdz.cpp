
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_neuron_o2_w.h"

using namespace amici;

void dJrzdz_model_neuron_o2(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJrzdz,0,sizeof(realtype)*model->nz*model->nztrue*model->nJ);
w_model_neuron_o2(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJrzdz[iz+(0+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(1+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(1+1*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(2+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*2]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(2+2*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(3+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*3]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(3+3*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(4+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*4]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(4+4*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
}
return;

}


