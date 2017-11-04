
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dJrzdsigma_model_events(double *dJrzdsigma, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJrzdsigma[iz+(0+0*1)*model->nztrue] = (rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0])*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*-1.0+1.0/tdata->sigmaz[0];
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  tdata->dJrzdsigma[iz+(0+1*1)*model->nztrue] = (rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1])*1.0/(tdata->sigmaz[1]*tdata->sigmaz[1]*tdata->sigmaz[1])*-1.0+1.0/tdata->sigmaz[1];
}
}

