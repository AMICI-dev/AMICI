
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dJzdz_model_events(double *dJzdz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->dJzdz[iz+(0+0*1)*model->nztrue] = 1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*0]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  tdata->dJzdz[iz+(0+1*1)*model->nztrue] = 1.0/(tdata->sigmaz[1]*tdata->sigmaz[1])*(edata->mz[tdata->nroots[ie]+udata->nmaxevent*1]*2.0-rdata->z[tdata->nroots[ie]+udata->nmaxevent*1]*2.0)*-5.0E-1;
}
}

