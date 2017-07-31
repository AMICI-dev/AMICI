
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_neuron_o2_w.h"

int dJrzdz_model_neuron_o2(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJrzdz,0,sizeof(realtype)*udata->nz*udata->nztrue*udata->nJ);
status = w_model_neuron_o2(t,x,NULL,user_data);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJrzdz[iz+(0+0*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(1+0*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(1+1*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(2+0*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*2]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(2+2*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(3+0*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*3]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(3+3*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(4+0*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*4]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
  tdata->dJrzdz[iz+(4+4*5)*udata->nztrue] = rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*1.0;
}
return(status);

}


