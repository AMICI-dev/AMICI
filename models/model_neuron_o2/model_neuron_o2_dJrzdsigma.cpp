
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

int dJrzdsigma_model_neuron_o2(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJrzdsigma,0,sizeof(realtype)*model->nztrue*model->nz*model->nJ);
status = w_model_neuron_o2(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJrzdsigma[iz+(0+0*5)*model->nztrue] = (rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0])*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-1.0+1.0/tdata->sigmaz[0];
  tdata->dJrzdsigma[iz+(1+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-2.0;
  tdata->dJrzdsigma[iz+(2+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*2]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-2.0;
  tdata->dJrzdsigma[iz+(3+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*3]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-2.0;
  tdata->dJrzdsigma[iz+(4+0*5)*model->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*4]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-2.0;
}
return(status);

}


