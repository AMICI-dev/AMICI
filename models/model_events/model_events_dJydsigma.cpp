
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dJydsigma_model_events(double *dJydsigma, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydsigma[iy+(0)*model->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0]*1.0,2.0)*-1.0+1.0/tdata->sigmay[0];
}
}

