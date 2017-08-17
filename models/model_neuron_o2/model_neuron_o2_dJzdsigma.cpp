
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_neuron_o2_w.h"

int dJzdsigma_model_neuron_o2(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJzdsigma,0,sizeof(realtype)*udata->nztrue*udata->nz*udata->nJ);
status = w_model_neuron_o2(t,x,NULL,tdata);
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJzdsigma[iz+(0+0*5)*udata->nztrue] = 1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*pow(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0],2.0)*-1.0+1.0/tdata->sigmaz[0];
  tdata->dJzdsigma[iz+(1+0*5)*udata->nztrue] = rdata->z[tdata->nroots[ie]+udata->nmaxevent*1]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0)*1.0;
  tdata->dJzdsigma[iz+(2+0*5)*udata->nztrue] = rdata->z[tdata->nroots[ie]+udata->nmaxevent*2]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0)*1.0;
  tdata->dJzdsigma[iz+(3+0*5)*udata->nztrue] = rdata->z[tdata->nroots[ie]+udata->nmaxevent*3]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0)*1.0;
  tdata->dJzdsigma[iz+(4+0*5)*udata->nztrue] = rdata->z[tdata->nroots[ie]+udata->nmaxevent*4]*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0)*1.0;
}
return(status);

}


