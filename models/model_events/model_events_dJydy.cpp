
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dJydy_model_events(double *dJydy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydy[iy+(0)*model->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
}
}

