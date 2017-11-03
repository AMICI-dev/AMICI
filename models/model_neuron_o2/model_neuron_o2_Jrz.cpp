
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

void Jrz_model_neuron_o2(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
w_model_neuron_o2(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->Jz[0] += amilog((tdata->sigmaz[0]*tdata->sigmaz[0])*3.141592653589793*2.0)*5.0E-1+(rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0])*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*5.0E-1;
  tdata->Jz[1] += rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->Jz[2] += rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*2]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->Jz[3] += rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*3]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->Jz[4] += rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*4]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
}
return;

}


