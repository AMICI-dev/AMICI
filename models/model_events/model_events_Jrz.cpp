
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void Jrz_model_events(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz) {
int iz;
if(!amiIsNaN(edata->mz[0*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 0;
  tdata->Jz[0] += amilog((tdata->sigmaz[0]*tdata->sigmaz[0])*3.141592653589793*2.0)*5.0E-1+(rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0])*1.0/(tdata->sigmaz[0]*tdata->sigmaz[0])*5.0E-1;
}
if(!amiIsNaN(edata->mz[1*udata->nmaxevent+tdata->nroots[ie]])){
    iz = 1;
  tdata->Jz[0] += amilog((tdata->sigmaz[1]*tdata->sigmaz[1])*3.141592653589793*2.0)*5.0E-1+(rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1]*rdata->rz[tdata->nroots[ie]+udata->nmaxevent*1])*1.0/(tdata->sigmaz[1]*tdata->sigmaz[1])*5.0E-1;
}
}

