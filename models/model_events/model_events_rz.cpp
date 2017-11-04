
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void rz_model_events(double *rz, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
  rdata->rz[tdata->nroots[ie]+udata->nmaxevent*0] = x_tmp[0]-x_tmp[2];
}

