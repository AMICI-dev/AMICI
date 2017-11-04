
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dJzdsigma_model_events(double *dJzdsigma, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJzdsigma[iz+(0+0*1)*model->nztrue] = 1.0/(tdata->sigmaz[0]*tdata->sigmaz[0]*tdata->sigmaz[0])*pow(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0],2.0)*-1.0+1.0/tdata->sigmaz[0];
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  tdata->dJzdsigma[iz+(0+1*1)*model->nztrue] = 1.0/(tdata->sigmaz[1]*tdata->sigmaz[1]*tdata->sigmaz[1])*pow(edata->mz[tdata->nroots[ie]+udata->nmaxevent*1]-rdata->z[tdata->nroots[ie]+udata->nmaxevent*1],2.0)*-1.0+1.0/tdata->sigmaz[1];
}
}

