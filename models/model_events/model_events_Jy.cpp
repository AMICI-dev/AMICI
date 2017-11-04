
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void Jy_model_events(double *nllh, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->Jy[0] += amilog((tdata->sigmay[0]*tdata->sigmay[0])*3.141592653589793*2.0)*5.0E-1+1.0/(tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0],2.0)*5.0E-1;
}
}

